# ChiralAI

**AI-guided discovery of biosynthetically-accessible chiral molecules.**

---

## The Problem

Enantiopure compounds are essential in pharmaceutical synthesis, agrochemistry, and advanced materials — ~80% of chiral APIs are now required as single enantiomers. Biocatalysis is the preferred route: enzymes are inherently chiral, operate under mild conditions, and can achieve >99% ee. But identifying the right enzyme, for the right substrate, with the right stereochemical outcome, in a host that can actually produce the compound — that requires manually stitching together five different tools and databases that were never designed to talk to each other.

No existing software does this end-to-end. **ChiralAI does.**

---

## What ChiralAI Does

Given a natural-language query, ChiralAI runs a grounded discovery pipeline:

```
User query (e.g., "enantiopure amine building block for beta-lactam synthesis")
    ↓
GPT-4.1 — suggests candidate chiral molecules with defined stereochemistry
    ↓
RDKit — validates chirality; identifies and assigns R/S stereocenters
    ↓
KEGG — maps compounds to known metabolic pathways and enzyme classes
    ↓
[BRENDA — retrieves known ee values and enantioselective enzyme data]     ← in progress
    ↓
[COBRApy FBA — checks metabolic feasibility in a target host organism]    ← in progress
    ↓
Ranked CSV output with stereochemistry, pathway, enzyme, and feasibility data
```

The LLM is the **orchestration and reasoning layer** — not the scientific ground truth. Every suggestion is grounded in a database call or computational result.

---

## Why This Is a Real Gap

| Tool | What It Does | What It Misses |
|------|-------------|----------------|
| RetroBioCat | Biocatalytic route planning | No stereochemistry or enantioselectivity awareness |
| ASKCOS | Organic retrosynthesis | Not built for enzymatic pathways |
| ChemCrow | GPT-4 + chemistry tools | Organic synthesis only; no biocatalysis |
| COBRApy | Genome-scale metabolic FBA | Ignores stereochemistry entirely |
| BRENDA | Gold-standard enzyme database with ee values | A database, not a discovery tool |

ChiralAI's contribution is integration — connecting retrosynthetic reasoning, chiral validation, pathway context, and metabolic feasibility in a single workflow accessible via natural language.

---

## Current State

**Working:**
- GPT-4.1 → structured JSON molecule suggestions
- RDKit chirality validation (stereocenters, R/S assignments)
- KEGG compound and pathway lookup
- Timestamped CSV export

**In Progress:**
- BRENDA integration (enantioselective enzyme data, known ee values)
- Structured KEGG data parsing (pathways, enzyme classes, reactions)
- COBRApy metabolic feasibility layer (*E. coli* iJO1366 model)
- Enantioselectivity scoring from database-retrieved ee data

---

## Quickstart

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Set your OpenAI API key
cp env.example .env
# Edit .env and add: OPENAI_API_KEY=your_key_here

# 3. Run
python3 main.py
```

Enter a query like:
- `"suggest a chiral amino acid precursor for asymmetric synthesis"`
- `"enantiopure lactone building blocks for biodegradable polymers"`
- `"(R)-selective secondary alcohol for pharmaceutical synthesis via E. coli fermentation"`

Results are saved to a timestamped CSV in the project directory.

---

## Architecture

```
ChiralAI/
├── main.py                      # Orchestrator — runs the full pipeline
├── ChiraLLM/
│   ├── query_handler.py         # GPT-4.1 interface — molecule suggestion
│   ├── chirality_checker.py     # RDKit — stereocenter detection and R/S assignment
│   └── database_validator.py   # KEGG REST API — pathway and enzyme lookup
└── utils/
    └── file_saver.py            # Timestamped CSV export
```

**Planned additions:**
- `ChiraLLM/brenda_client.py` — BRENDA SOAP API for ee values and substrate specificity
- `ChiraLLM/feasibility_checker.py` — COBRApy FBA for host organism metabolic feasibility
- `ChiraLLM/enantioselectivity_scorer.py` — enzyme ranking by predicted/known ee

---

## Scientific Grounding

- **KEGG** — compound, pathway, and enzyme commission data
- **BRENDA** (planned) — 112k enzymes, 5.8M data points, chiral SMILES, ee values
- **MetaCyc** (planned) — 3,284 curated biosynthetic pathways
- **COBRApy + iJO1366** (planned) — genome-scale *E. coli* metabolic model for flux analysis
- **RDKit** — open-source cheminformatics; stereocenters, SMILES validation, R/S assignment

---

## Roadmap

- [ ] BRENDA integration — enantioselective enzyme lookup by substrate
- [ ] KEGG flat-file parsing — structured pathway/enzyme/reaction extraction
- [ ] COBRApy feasibility layer — FBA with target metabolite, cofactor flagging
- [ ] Enantioselectivity scoring — rank enzymes by known ee, substrate similarity
- [ ] Multi-candidate output — ranked list of 5–10 molecules per query
- [ ] Enzyme engineering flags — identify substrates requiring directed evolution

---

Questions or ideas? Connect on LinkedIn: https://www.linkedin.com/in/alexeimanuel/
