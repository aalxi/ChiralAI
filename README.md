# ChiralAI

**A natural-language discovery engine for biocatalytic chiral targets — from query to ranked, route-resolved, feasibility-checked candidates in one pipeline.**

Ask in English. Get back a ranked list of chiral molecules with: defined R/S stereochemistry, the enzyme classes that make them, a predicted biosynthetic route from central metabolism, BRENDA-verified enantioselectivity at each step, and metabolic feasibility in *E. coli* — with every claim traceable to a database call.

---

## The Problem

~80% of new chiral pharmaceuticals must be sold as single enantiomers. Biocatalysis is the preferred route — enzymes are inherently chiral and routinely hit >99% ee — but answering "which chiral target is actually buildable, by which enzyme, in which host?" today means manually stitching KEGG, BRENDA, PubChem, COBRApy, and RDKit, none of which were designed to talk to each other. ChiralAI is the missing connective tissue.

---

## Pipeline

```
Natural-language query
    │
    ▼
GPT-4.1                  → 5 candidate chiral molecules with R/S-specified SMILES
    │
    ▼
RDKit                    → CIP stereocenter detection; flags unassigned centers
    │
    ▼
KEGG                     → pathway + EC enzyme lookup
    │
    ▼
Route predictor (Tier 1) → weighted A* over KEGG reactions; backward search to
                           central metabolites; transparent edge costs
    │
    ▼
BRENDA                   → ee% per EC, deduplicated across all route steps
    │
    ▼
COBRApy / iJO1366        → FBA on terminal precursor in E. coli K-12
    │
    ▼
Composite scorer         → weighted ee × Tanimoto × stereo × feasibility,
                           with confidence tier (high / medium / low)
    │
    ▼
Timestamped CSV + JSON   → flat scoring columns + nested provenance
```

---

## Design choices that matter

These are the decisions that make the pipeline useful rather than impressive-looking:

- **Database-first, with provenance labels.** Every ee value is tagged `brenda_verified` or `llm_claim`. LLM claims are discounted 0.6× and visually distinguishable in output. The LLM proposes; the databases dispose.

- **Weighted A* with a chemical-distance heuristic.** The route predictor isn't blind BFS — it uses Morgan-fingerprint Tanimoto distance to the nearest curated central metabolite as the heuristic, dramatically pruning the KEGG search space while keeping it provably better than greedy.

- **Transparent edge costs, not a black box.** Each route step's cost decomposes into `base + thermodynamic (eQuilibrator ΔG) + directionality (KEGG ⇌ vs →) + industrial reversibility override`. A user can read the breakdown and disagree with any term.

- **Industrial reversibility override.** Pure ΔG would penalize KREDs, transaminases, IREDs, BVMOs, and lipases for running "uphill" — but in practice these enzyme classes are routinely reversed in industrial biocatalysis via cofactor recycling and substrate engineering. The override encodes that domain knowledge as an EC-prefix lookup.

- **Multiplicative ee composition.** Whole-route ee is `∏(ee_i / 100) × 100` across all steps. Surfaces the compounding stereochemistry loss that single-step ee numbers hide.

- **Curated central metabolite stop set.** 41 nodes — full TCA, glycolysis, PPP, the 20 amino acids, KIV. Backward search terminates at biological reality, not arbitrary depth.

- **Three output modes.** `top_n` for ranked alternatives, `full_tree` for the complete branching DAG, `shortest_plus_diverse` for a designable portfolio. Researchers ask different questions; the predictor answers all three.

- **Disk caching at `~/.cache/chiralai/`.** KEGG reactions, MOL files, and eQuilibrator ΔG values cache with a 30-day TTL. Repeat runs are essentially free.

---

## What it is not (yet)

- **Tier 2 retrobiosynthesis (novel targets).** Route prediction currently covers compounds KEGG already knows (~12k reactions). Targets outside KEGG need RetroRules SMARTS + RDKit `RunReactants` — planned for Sprint 2.
- **Non-*E. coli* hosts.** Feasibility uses iJO1366 only. iMM904 (yeast), P. putida, etc., are not yet wired.
- **Enzyme engineering.** ChiralAI flags when a wild-type enzyme's ee is insufficient and leaves variant design to the user (Rosetta, ProteinMPNN). It does not attempt this computationally.

---

## The wedge

ChiralAI does not compete with the tools below at their core functions. It runs the chain none of them span — chiral validation, route prediction, enantioselectivity, host feasibility — with explicit provenance throughout.

| Tool | Core strength | What it leaves out |
|------|---------------|--------------------|
| RetroBioCat | Biocatalytic route planning | Stereochemistry / ee |
| ASKCOS | Organic retrosynthesis | Enzymes |
| ChemCrow | LLM + chemistry tools | Biocatalysis |
| COBRApy | Genome-scale FBA | Stereochemistry |
| BRENDA | Enzyme database | Discovery |

The discipline none of them enforce: labeling LLM speculation distinctly from database-verified fact.

---

## Quickstart

```bash
pip install -r requirements.txt
cp env.example .env
# Add OPENAI_API_KEY (required), BRENDA_EMAIL + BRENDA_PASSWORD (recommended)
python3 main.py
```

Try:
- `"enantiopure amine for a beta-lactam side chain"`
- `"chiral lactone monomer for biodegradable polymers"`
- `"(R)-secondary alcohol producible in E. coli fermentation"`

Outputs land as `suggestions_<timestamp>.csv` plus a JSON sidecar with the full route trees and scoring breakdown.

---

## Architecture

```
ChiralAI/
├── main.py                            Orchestrator — thin
└── ChiraLLM/
    ├── query_handler.py               GPT-4.1 → 5 candidates as JSON
    ├── chirality_checker.py           RDKit stereocenter detection
    ├── database_validator.py          KEGG REST flat-file parser
    ├── route_predictor.py             Tier 1 A* over KEGG; transparent costs
    ├── brenda_client.py               BRENDA SOAP; ee% from commentary
    ├── feasibility_checker.py         COBRApy FBA on iJO1366
    └── enantioselectivity_scorer.py   Composite score + confidence tier
└── utils/file_saver.py                Timestamped CSV + JSON
```

---

## CSV output

| Column | Source | Notes |
|--------|--------|-------|
| `scoring_composite_score` | Scorer | 0–1; weighted ee + Tanimoto + feasibility |
| `scoring_confidence` | Scorer | `high` / `medium` / `low` |
| `scoring_top_enzyme_ec` | BRENDA / KEGG | Best-ranked EC |
| `scoring_top_enzyme_ee` | BRENDA / LLM | ee% value |
| `scoring_top_enzyme_source` | Scorer | `brenda_verified` or `llm_claim` |
| `scoring_stereo_confirmed` | RDKit | True only if all centers R/S assigned |
| `scoring_feasibility_flux` | COBRApy | mmol/gDW/h; None if not in model |
| `route_top1_step_count` | Route predictor | Steps in shortest predicted route |
| `route_top1_terminal_precursor` | Route predictor | KEGG ID of the central metabolite reached |
| `route_top1_total_cost` | Route predictor | Summed transparent edge cost |
| `route_top1_composed_ee` | Route predictor + BRENDA | Multiplicative ee across the route |
| `scoring_notes` | Scorer | Human-readable provenance and caveats |

The JSON sidecar preserves the full nested structure: every route step, every BRENDA hit, every cost component.

---

## Known limitations

- **BRENDA credentials matter.** Without them, every suggestion falls back to `llm_claim` ee. Registration is free; this is the single highest-leverage configuration step.
- **iJO1366 only.** Secondary metabolites and many pharmaceutical targets return `not_in_model`.
- **KEGG coverage ceiling.** Tier 1 cannot route to compounds KEGG doesn't know about. Tier 2 will close this.

---

## Roadmap

- [ ] **Tier 2** — novel-target retrobiosynthesis via RetroRules SMARTS + RDKit `RunReactants`
- [ ] Non-*E. coli* host models — iMM904 (S. cerevisiae), P. putida
- [ ] BRENDA `getEngineering` — surface known directed-evolution variants per EC
- [ ] Programmatic name ↔ SMILES CIP consistency check

---

Questions or ideas? [LinkedIn](https://www.linkedin.com/in/alexeimanuel/)
