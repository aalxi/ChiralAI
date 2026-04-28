# ChiralAI

**AI-guided discovery and validation engine for biocatalytic chiral molecule targets.**

---

## The Problem

Enantiopure compounds are essential in pharmaceutical synthesis, agrochemistry, and advanced materials. ~80% of chiral active pharmaceuticals now have to be single enantiomers. Biocatalysis is the preferred route: enzymes are inherently chiral, operate under mild conditions, and can achieve >99% ee (enantiomeric excess). But identifying the right enzyme, for the right substrate, with the right stereochemical outcome, in a host that can actually produce the compound, that requires manually stitching together five different tools and databases that were never designed work compatibly.

No existing software does this end-to-end. 

---

## What ChiralAI Does

Given a natural-language query, ChiralAI runs a grounded discovery and validation pipeline:

```
User query (e.g., "enantiopure amine building block for beta-lactam synthesis")
    ↓
GPT-4.1 — suggests 5 candidate chiral molecules with defined R/S stereochemistry
    ↓
RDKit — validates chirality; identifies and assigns R/S stereocenters
    ↓
KEGG — maps compounds to known metabolic pathways and enzyme classes
    ↓
BRENDA — retrieves known ee values from enzyme substrate data (requires API credentials)
    ↓
COBRApy FBA — checks metabolic feasibility in E. coli iJO1366; flags cofactor requirements
    ↓
Composite scorer — ranks candidates by ee source, Tanimoto substrate similarity, and feasibility
    ↓
Timestamped CSV + JSON output with full provenance and confidence tiers
```

The LLM is the **orchestration and reasoning layer**, not the scientific ground truth. Every suggestion is grounded in a database call or computational result. LLM-claimed ee values are labeled `"llm_claim"` and discounted; BRENDA-verified ee is labeled `"brenda_verified"` and weighted higher.

---

## Why This Is a Real Gap

| Tool | What It Does | What It Misses |
|------|-------------|----------------|
| RetroBioCat | Biocatalytic route planning | No stereochemistry or enantioselectivity awareness |
| ASKCOS | Organic retrosynthesis | Not built for enzymatic pathways |
| ChemCrow | GPT-4 + chemistry tools | Organic synthesis only; no biocatalysis |
| COBRApy | Genome-scale metabolic FBA | Ignores stereochemistry entirely |
| BRENDA | Gold-standard enzyme database with ee values | A database, not a discovery tool |

ChiralAI's contribution is integration — connecting chiral validation, pathway context, enantioselectivity data, and metabolic feasibility in a single workflow accessible via natural language.

---

## Quickstart

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Set credentials
cp env.example .env
# Add: OPENAI_API_KEY, BRENDA_EMAIL, BRENDA_PASSWORD

# 3. Run
python3 main.py
```

Enter a query like:
- `"suggest a chiral amino acid precursor for asymmetric synthesis"`
- `"enantiopure lactone building blocks for biodegradable polymers"`
- `"(R)-selective secondary alcohol for pharmaceutical synthesis via E. coli fermentation"`

Results are saved to a timestamped CSV and JSON sidecar in the project directory.

---

## Architecture

```
ChiralAI/
├── main.py                            # Orchestrator — runs the full pipeline
├── ChiraLLM/
│   ├── query_handler.py               # GPT-4.1 — molecule suggestion (5 ranked candidates)
│   ├── chirality_checker.py           # RDKit — stereocenter detection and R/S assignment
│   ├── database_validator.py          # KEGG REST API — pathway and enzyme lookup
│   ├── brenda_client.py               # BRENDA SOAP — ee values from enzyme substrate data
│   ├── feasibility_checker.py         # COBRApy FBA — metabolic feasibility in iJO1366
│   └── enantioselectivity_scorer.py   # Composite scorer — ranks by ee, Tanimoto, feasibility
└── utils/
    └── file_saver.py                  # Timestamped CSV + JSON export
```

---

## Output Columns (CSV)

| Column | Source | Notes |
|--------|--------|-------|
| `scoring_composite_score` | Scorer | 0–1; weighted ee + Tanimoto + feasibility |
| `scoring_confidence` | Scorer | `high` / `medium` / `low` |
| `scoring_top_enzyme_ec` | BRENDA / KEGG | Best-ranked EC number |
| `scoring_top_enzyme_ee` | BRENDA / LLM | ee% value |
| `scoring_top_enzyme_source` | Scorer | `brenda_verified` or `llm_claim` |
| `scoring_stereo_confirmed` | RDKit | True only if all stereocenters are R/S assigned |
| `scoring_feasibility_flux` | COBRApy | mmol/gDW/h in iJO1366; None if not in model |
| `scoring_notes` | Scorer | Human-readable provenance and caveats |

---

## Scientific Grounding

- **RDKit** — stereocenter detection, CIP R/S assignment, Morgan fingerprints for Tanimoto
- **KEGG** — compound, pathway, and enzyme commission data
- **BRENDA** — 112k enzymes, 5.8M data points; ee% extracted from substrate commentary fields
- **PubChem** — substrate SMILES resolution for Tanimoto computation
- **COBRApy + iJO1366** — genome-scale *E. coli* K-12 metabolic model for flux analysis

---

## Known Limitations

- **BRENDA credentials required for verified ee.** Without `BRENDA_EMAIL` / `BRENDA_PASSWORD` in `.env`, all ee values fall back to LLM claims labeled `llm_claim`. BRENDA registration is free.
- **E. coli only.** The FBA layer uses iJO1366 (E. coli K-12). Secondary metabolites and many pharmaceutical targets return `not_in_model`. Other host models (S. cerevisiae, P. putida) are not yet supported.
- **No retrosynthetic route planning.** ChiralAI validates and scores named targets; it does not enumerate the enzymatic steps needed to build a molecule from simpler precursors.

---

## Roadmap

- [ ] Retrosynthetic route prediction — enumerate enzymatic steps from target to precursors
- [ ] Non-E. coli host models — S. cerevisiae (iMM904), P. putida support in FBA layer
- [ ] Engineered variant data — wire BRENDA `getEngineering` for directed evolution candidates
- [ ] Name↔SMILES stereo consistency check — programmatic CIP verification against molecule name

---

Questions or ideas? Connect on LinkedIn: https://www.linkedin.com/in/alexeimanuel/
