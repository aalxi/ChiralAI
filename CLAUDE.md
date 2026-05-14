# CLAUDE.md — ChiralAI Project Brief

This file is the standing brief for any AI assistant working in this repository. Read it before touching any code.

---

## What This Project Is

ChiralAI is an **AI-guided discovery, validation, and route-prediction engine for biocatalytic chiral molecule targets**. It is not a general chemistry tool.

The pipeline answers: *given a natural-language research goal, which chiral molecules are strong biocatalytic targets, which enzymes produce them with high enantioselectivity, what is the shortest known biosynthetic route from central metabolites, and is the biosynthesis metabolically feasible in E. coli?*

For compounds in KEGG (~12k reactions), it now **does** enumerate the enzymatic steps needed to reach a target from central metabolites via Tier 1 route prediction. Novel-target retrobiosynthesis via SMARTS pattern matching (Tier 2) is planned for the next sprint.

The scientific premise: biomanufacturing is inherently enantioselective because enzymes are chiral catalysts. This system helps researchers find molecules that are (1) chiral with defined stereochemistry, (2) reachable through known metabolic pathways, (3) producible with high enantioselectivity by known or engineerable enzymes, and (4) metabolically feasible in a target host organism.

---

## Scientific Domain Context

**Chirality / Enantioselectivity:** Chiral molecules exist as non-superimposable mirror images (enantiomers). Enantiomers are designated R or S by CIP priority rules. Enzymes produce single enantiomers due to their chiral active sites. Enantiomeric excess (ee%) measures selectivity: ee = |[R] - [S]| / ([R] + [S]) × 100. Pharmaceutical APIs typically require >99% ee; >80% ee is considered synthetically useful.

**Biocatalysis vs. chemical asymmetric synthesis:** Biocatalytic routes are preferred in pharmaceutical and green chemistry contexts. Key enzyme classes for enantioselective synthesis: ketoreductases (KRED), transaminases (TA), lipases, cytochrome P450s, epoxide hydrolases, lyases. These all operate with natural enantioselectivity that can be shifted by directed evolution or rational design.

**KEGG flat-file format:** KEGG REST API returns compound entries as plain text with labeled fields: ENTRY, NAME, FORMULA, EXACT_MASS, MOL_WEIGHT, REACTION, PATHWAY, ENZYME, BRITE, DBLINKS. Parse these by splitting on newline and extracting field labels.

**COBRA / FBA:** Flux Balance Analysis maximizes an objective function (e.g., production of target metabolite) subject to stoichiometric constraints of a metabolic network. `cobra.io.load_model('iJO1366')` loads the standard *E. coli* K-12 genome-scale model. Use `model.optimize()` and check `solution.fluxes` for pathway flux.

**BRENDA:** The gold-standard enzyme database. Accessible via SOAP API (`brenda-enzymes.org/soap/brenda_zeep.wsdl`). Enantioselectivity data is stored in free-text commentary fields (`commentarySubstrates`, `commentaryProducts`) — there is no structured `getEnantioselectivity` method. Use `getSubstratesProducts` with positional param-strings and parse ee% from the commentary via regex.

---

## The Competitive Landscape (Do Not Duplicate These)

- **ASKCOS** — organic retrosynthesis, trained on USPTO. Do not try to compete here.
- **RetroBioCat** — biocatalytic route planning by enzyme class. Our reference for what "route planning" looks like; ChiralAI is complementary (we validate and score targets; we don't enumerate routes).
- **COBRApy / StrainDesign / Cameo** — FBA and strain optimization. Use COBRApy as a library; do not reinvent it.
- **Rosetta / ProteinMPNN** — de novo enzyme design. Do not try to do this computationally — flag when engineering is needed and leave it to the user.
- **ChemCrow** — LLM + chemistry tools for organic synthesis. Does not address biocatalysis or enantioselectivity.

ChiralAI's value is **integration and natural-language accessibility**, not outcompeting any of the above at their core function.

---

## Architecture

All 8 modules are implemented and wired. `main.py` calls them in sequence.

```
main.py                               Orchestrator — runs the full pipeline, thin logic only.

ChiraLLM/query_handler.py             GPT-4.1 interface. Returns 5 ranked candidates as JSON.
                                      System prompt explicitly requests: R/S-specified SMILES,
                                      enzyme class, KEGG ID, known ee. JSON schema:
                                      {"suggestions": [{name, SMILES, KEGG_ID, enzyme_class,
                                       known_ee, rationale, applications}, ...]}

ChiraLLM/chirality_checker.py         RDKit chirality validation. validate_chirality(smiles)
                                      returns {valid, chiral_centers, n_centers, all_assigned}.
                                      chiral_centers is a list of (atom_idx, 'R'|'S'|'?') tuples.

ChiraLLM/database_validator.py        KEGG REST flat-file parser. query_kegg(compound_id)
                                      returns {status, name, formula, pathway_names, enzymes,
                                      reactions}. Parses structured fields; does not return
                                      raw text.

ChiraLLM/brenda_client.py             BRENDA SOAP API via zeep (Settings(strict=False) required).
                                      query_enantioselectivity(ec_number) calls getSubstratesProducts
                                      with positional param-strings. Extracts ee% via regex from
                                      commentary fields. Filters out substrates='more' sentinel.
                                      Returns {status, entries: [{substrates, enantioselectivity,
                                      stereo, organism, commentary}]}.

ChiraLLM/feasibility_checker.py       COBRApy FBA on iJO1366 (E. coli K-12). Matches KEGG ID
                                      via kegg.compound annotation. Adds a demand reaction,
                                      optimizes, reports flux + cofactor requirements.
                                      Returns {status, flux, cofactors_required, metabolite_id}.

ChiraLLM/enantioselectivity_scorer.py Composite scoring. score_suggestion(suggestion) synthesizes
                                      BRENDA ee, Tanimoto substrate similarity (PubChem SMILES +
                                      Morgan fingerprints), stereo confirmation, and FBA flux into
                                      a 0–1 composite score with confidence tier (high/medium/low).
                                      Degrades gracefully to LLM-claim ee when BRENDA is absent.

ChiraLLM/route_predictor.py           Tier 1 biosynthesis route predictor. predict_route(compound_id,
                                      mode='top_n'|'full_tree'|'shortest_plus_diverse', n, budget)
                                      runs weighted A* backward through KEGG reaction graph from
                                      target to curated central metabolites. Edge cost is a
                                      transparent breakdown: thermodynamic (eQuilibrator ΔG),
                                      directionality (KEGG <=> vs =>), industrial-reversibility
                                      override (KREDs/TAs/IREDs/lipases/BVMOs). Heuristic is
                                      Morgan-fingerprint Tanimoto distance to nearest central.
                                      Disk-cached at ~/.cache/chiralai/. Returns RouteResult
                                      with status field; never raises.

utils/file_saver.py                   Timestamped CSV + JSON sidecar. Flattens scoring dict into
                                      flat columns; preserves full nested structure in JSON.
```

`filter.py` at repo root is a legacy prototype (pre-module architecture). It is not wired into `main.py` and should be ignored.

---

## Development Conventions

- **The LLM is the orchestration layer, not the ground truth.** Every molecular claim made by GPT must be verifiable by a database call or RDKit computation. Never present an LLM hallucination as a fact.
- **Database-first for enantioselectivity.** Use BRENDA ee values (`ee_source: "brenda_verified"`) when available. LLM-claimed ee is discounted 0.6× and labeled `"llm_claim"` — it must be visually distinguishable from verified data in any output.
- **Fail informatively.** If KEGG returns no data, or BRENDA has no ee data for a substrate, say so in the output — don't silently skip. Researchers need to know the confidence level of each suggestion.
- **SMILES must be stereospecific.** All SMILES strings in the pipeline should use `@` and `@@` to encode defined stereochemistry. Flag SMILES without stereo specification — they indicate the LLM failed to commit to an enantiomer.
- **No hardcoded fallbacks.** If a module fails, it should raise or return an error dict, not silently return a default molecule.
- **Keep main.py thin.** Orchestration logic only. All science belongs in `ChiraLLM/` modules.
- **BRENDA SOAP calling convention is positional, not keyword.** Arguments must be passed as `"paramName*filterValue"` strings. Keyword arguments silently return 0 results.

---

## Key External APIs

| Resource | Access | Notes |
|----------|--------|-------|
| KEGG REST | `http://rest.kegg.jp/get/{id}` | Free, no auth. Parse flat-file text response. |
| BRENDA SOAP | `brenda-enzymes.org/soap/brenda_zeep.wsdl` | Requires free registration. Use zeep with `Settings(strict=False)`. Positional param-strings only. |
| PubChem | `pubchempy` library | Used for substrate SMILES lookup in scorer (Tanimoto computation). |
| OpenAI | `openai` client, model `gpt-4.1` | Key in `.env` as `OPENAI_API_KEY`. |
| COBRApy models | `cobra.io.load_model('iJO1366')` | Downloads automatically on first run; cached in `_model_cache`. |

---

## Current Module Status (April 2026)

| Module | Status | Notes |
|--------|--------|-------|
| `query_handler.py` | Working | GPT-4.1, returns 5 ranked candidates as JSON |
| `chirality_checker.py` | Working | RDKit, `validate_chirality()` returns chiral_centers with R/S tags |
| `database_validator.py` | Working | KEGG REST flat-file parser; structured output |
| `brenda_client.py` | Working | Full BRENDA SOAP with correct positional calling convention |
| `feasibility_checker.py` | Working | COBRApy FBA on iJO1366; cofactor flagging |
| `enantioselectivity_scorer.py` | Working | Composite scoring with BRENDA-verified ee + Tanimoto similarity |
| `route_predictor.py` | Working | Tier 1 (KEGG traversal); Tier 2 RetroRules SMARTS deferred |
| `file_saver.py` | Working | CSV + JSON sidecar with flat scoring columns |

## Known Gaps

- **BRENDA credentials**: Free-tier API access requires `BRENDA_EMAIL` and `BRENDA_PASSWORD` in `.env`. Without credentials, every suggestion falls back to `ee_source: llm_claim` and `confidence: medium`. Acquiring credentials is the single highest-leverage improvement.
- **Non-E. coli host support**: `feasibility_checker.py` only loads iJO1366 (E. coli K-12). Secondary metabolites and many pharmaceutical targets are `not_in_model`. Non-native pathway hosts (S. cerevisiae iMM904, P. putida) require additional models.
- **`engineered_variants` field**: Placeholder `None` in all scorer output. `getEngineering` BRENDA SOAP method works but is not yet wired.
- **Tier 2 retrobiosynthesis**: Tier 1 covers compounds KEGG already knows. Novel-target retrobiosynthesis (target SMILES → enzymatic disconnection via RetroRules SMARTS + RDKit `RunReactants`) is planned for Sprint 2.
