# CLAUDE.md — ChiralAI Project Brief

This file is the standing brief for any AI assistant working in this repository. Read it before touching any code.

---

## What This Project Is

ChiralAI is an **AI-guided discovery engine for biosynthetically-accessible chiral molecules**. It is not a general chemistry tool. It targets a specific, unoccupied gap in the field: no existing software integrates retrosynthetic biocatalysis planning with enantioselectivity data and metabolic feasibility in a single pipeline.

The scientific premise: biomanufacturing is inherently enantioselective because enzymes are chiral catalysts. This system helps researchers find molecules that are (1) chiral with defined stereochemistry, (2) reachable through known biosynthetic pathways, (3) producible with high enantioselectivity by known or engineerable enzymes, and (4) metabolically feasible in a target host organism.

---

## Scientific Domain Context

**Chirality / Enantioselectivity:** Chiral molecules exist as non-superimposable mirror images (enantiomers). Enantiomers are designated R or S by CIP priority rules. Enzymes produce single enantiomers due to their chiral active sites. Enantiomeric excess (ee%) measures selectivity: ee = |[R] - [S]| / ([R] + [S]) × 100. Pharmaceutical APIs typically require >99% ee; >80% ee is considered synthetically useful.

**Biocatalysis vs. chemical asymmetric synthesis:** Biocatalytic routes are preferred in pharmaceutical and green chemistry contexts. Key enzyme classes for enantioselective synthesis: ketoreductases (KRED), transaminases (TA), lipases, cytochrome P450s, epoxide hydrolases, lyases. These all operate with natural enantioselectivity that can be shifted by directed evolution or rational design.

**KEGG flat-file format:** KEGG REST API returns compound entries as plain text with labeled fields: ENTRY, NAME, FORMULA, EXACT_MASS, MOL_WEIGHT, REACTION, PATHWAY, ENZYME, BRITE, DBLINKS. Parse these by splitting on newline and extracting field labels.

**COBRA / FBA:** Flux Balance Analysis maximizes an objective function (e.g., production of target metabolite) subject to stoichiometric constraints of a metabolic network. `cobra.io.load_model('iJO1366')` loads the standard *E. coli* K-12 genome-scale model. Use `model.optimize()` and check `solution.fluxes` for pathway flux.

**BRENDA:** The gold-standard enzyme database. Accessible via SOAP API (`brenda-enzymes.org/soap/brenda_zeep.wsdl`). Stores chiral SMILES for stereospecific substrates. Key query: `getKmValue(enzyme_name, smiles)` and `getSpecificActivity`. Has ee data for many well-characterized enzymes.

---

## The Competitive Landscape (Do Not Duplicate These)

- **ASKCOS** — organic retrosynthesis, trained on USPTO. Do not try to compete here.
- **RetroBioCat** — biocatalytic route planning by enzyme class. Our reference for what "route planning" looks like; we extend it by adding enantioselectivity.
- **COBRApy / StrainDesign / Cameo** — FBA and strain optimization. Use COBRApy as a library; do not reinvent it.
- **Rosetta / ProteinMPNN** — de novo enzyme design. Do not try to do this computationally — flag when engineering is needed and leave it to the user.
- **ChemCrow** — LLM + chemistry tools for organic synthesis. Does not address biocatalysis or enantioselectivity.

ChiralAI's value is **integration and natural-language accessibility**, not outcompeting any of the above at their core function.

---

## Architecture

```
main.py                          Orchestrator. Calls modules in sequence, handles JSON parsing,
                                 saves results. Keep this file thin — logic belongs in modules.

ChiraLLM/query_handler.py        GPT-4.1 interface. System prompt must explicitly request:
                                 - Defined stereochemistry (R/S in SMILES and by name)
                                 - Known enantioselective enzyme class
                                 - KEGG ID and BRENDA-referenced ee data if available
                                 - 5 ranked candidates, not 1
                                 Response format: JSON array of molecule objects.

ChiraLLM/chirality_checker.py    RDKit chirality validation. validate_chirality(smiles) is the
                                 only function. Returns: valid (bool), chiral_centers (list of
                                 (atom_idx, R/S) tuples), n_centers (int), all_assigned (bool).
                                 NOTHING ELSE belongs in this file.

ChiraLLM/database_validator.py   KEGG REST API. Parse flat-file format: extract NAME, FORMULA,
                                 PATHWAY, ENZYME, REACTION fields. Return a structured dict,
                                 not raw trimmed text.

utils/file_saver.py              CSV export. Add columns for: chiral_centers, R_S_config,
                                 known_ee, enantioselective_enzyme, pathway_names,
                                 feasibility_score when those fields exist on suggestions.

ChiraLLM/brenda_client.py        [TO BUILD] BRENDA SOAP API client. Query by EC number or
                                 substrate SMILES. Return: enzyme names, ee values, substrate
                                 specificity data, known organisms.

ChiraLLM/feasibility_checker.py  [TO BUILD] COBRApy FBA layer. Load iJO1366 or user-specified
                                 model. Check if target KEGG compound ID maps to a model
                                 metabolite. Run FBA. Flag cofactor dependencies (NADPH, PLP).

ChiraLLM/enantioselectivity_scorer.py  [TO BUILD] Rank enzyme candidates by known ee from
                                       BRENDA, substrate structural similarity (RDKit Tanimoto),
                                       and availability of engineered variants.
```

---

## Development Conventions

- **The LLM is the orchestration layer, not the ground truth.** Every molecular claim made by GPT must be verifiable by a database call or RDKit computation. Never present an LLM hallucination as a fact.
- **Database-first for enantioselectivity.** Use BRENDA ee values when available. ML prediction is a fallback for novel substrates, not the primary path.
- **Fail informatively.** If KEGG returns no data, or BRENDA has no ee data for a substrate, say so in the output — don't silently skip. A researcher needs to know the confidence level of each suggestion.
- **SMILES must be stereospecific.** All SMILES strings in the pipeline should use `@` and `@@` to encode defined stereochemistry. Flag SMILES without stereo specification — they indicate the LLM failed to commit to an enantiomer.
- **No hardcoded fallbacks.** If a module fails, it should raise or return an error dict, not silently return a default molecule.
- **Keep main.py thin.** Orchestration logic only. All science belongs in `ChiraLLM/` modules.

---

## Key External APIs

| Resource | Access | Notes |
|----------|--------|-------|
| KEGG REST | `http://rest.kegg.jp/get/{id}` | Free, no auth. Parse flat-file text response. |
| BRENDA SOAP | `brenda-enzymes.org/soap/brenda_zeep.wsdl` | Requires free registration for API key. |
| PubChem | `pubchempy` library | Already in requirements; use for SMILES validation / CID lookup. |
| OpenAI | `openai` client, model `gpt-4.1` | Key in `.env` as `OPENAI_API_KEY`. |
| COBRApy models | `cobra.io.load_model('iJO1366')` | Downloads automatically; cache locally. |

---

## Current Status (April 2026)

| Module | Status | Blocker |
|--------|--------|---------|
| `query_handler.py` | Working | Prompt too weak on stereochemistry; returns 1 molecule |
| `chirality_checker.py` | Broken | Lines 13–142 are garbage scratch-pad code — delete them |
| `database_validator.py` | Working but shallow | Returns raw text; needs flat-file parser |
| `file_saver.py` | Working | Missing stereochem columns |
| BRENDA client | Not started | — |
| COBRApy layer | Not started | — |
| Enantioselectivity scorer | Not started | — |

**Immediate unblock:** Delete `chirality_checker.py` lines 13–142, uncomment `main.py:39-40`.
