# Retrospective benchmark: design

## Why this test exists

ChiralAI claims to tell a metabolic engineer whether a chiral target can be made
biologically. The only honest way to support that claim is to run it on
compounds where the answer is already known from the literature, and to check
whether the system's answer tracks the truth *better than a trivial rule would*.

Every metric below is therefore paired with a baseline. A number with no
baseline is not evidence.

## Three separable claims

The pipeline makes three claims that fail independently, so they are scored
independently. Reporting one aggregate "accuracy" would hide which part is broken.

### Claim 1 — Identifier grounding
"Given a compound name, the system retrieves the right database entry."

Under test: `find_compound_kegg()`, which takes the **first line** of
KEGG `/find/compound/{name}`.

- **Coverage**: fraction of names returning any hit.
- **Top-1 accuracy**: fraction where the hit is the correct compound.
- **Stereo accuracy**: fraction where the hit is the correct *enantiomer*,
  not merely the right skeleton. This is the metric the project's name implies.
- Baseline: the same query with URL encoding plus exact-synonym re-ranking.
  The gap between the two separates a fixable bug from an intrinsic limit.

### Claim 2 — Feasibility discrimination
"The feasibility verdict separates makeable from not-makeable."

Under test: `combined_feasibility_analysis()` → `biosynthetic_feasibility`, `confidence`.

Positives are compounds with documented production in an engineered or
industrial microbial host. Negatives are **hard** negatives: chiral,
drug-sized molecules made industrially by chemical synthesis, with no reported
de novo microbial route. Hard negatives matter because a set of obviously
non-biological negatives (polymers, perfluorinated solvents) lets any
classifier score well and teaches nothing.

- **AUROC** on `confidence`, with bootstrap confidence interval.
- **Balanced accuracy** at the system's own operating point.
- **Precision** at that point, because a false "feasible" costs a wet-lab
  campaign while a false "infeasible" costs only an idea.
- Baselines: constant predictor; molecular weight under 500 alone; "present in
  KEGG at all"; "has at least one KEGG pathway".

### Claim 3 — Stereochemical handling
"ChiralAI is about chirality."

- Fraction of verdicts where the requested enantiomer was checked at all.
- Fraction where retrieved structure's stereocentres match the request.

## Stratification, and what it controls for

- **Stratum A1, classic positives** (pre-2015, heavily cited). A language model
  has read these.
- **Stratum A2, recent positives** (2015+). Performance that holds on A1 but
  collapses on A2 indicates recall of literature, not inference.
- **Stratum B, hard negatives.** Chiral synthetic drugs, mostly under 500 Da,
  so they pass the size rule and probe real biosynthetic knowledge.
- **Stratum C, near-miss.** Made by plants or mammals but not by a microbial
  host. These test the distinction a metabolic engineer actually cares about:
  "biology makes this somewhere" versus "a microbe has been engineered to make it."

## Power, stated before running

At roughly 20 positives and 16 negatives, the 95% interval on an AUROC near 0.8
spans about ±0.15. That resolves "informative" from "coin flip". It cannot rank
two similar methods, and no claim of that kind will be made from this run.

## Ground truth provenance

Structures and formulas are pulled from KEGG's own MOL records rather than
written from memory, so the reference set is auditable and not a second source
of error. Labels come from the published production literature and are cited
per row in `dataset.csv`.

## Pre-registered predictions

Recorded before running, so the results cannot be read as post-hoc narrative.

1. Retrieval coverage collapses for multi-word names. `/find/compound/` is
   called with raw spaces in the URL.
2. Where a name is enantiomer-ambiguous, first-hit retrieval returns an
   arbitrary enantiomer.
3. `confidence` carries at most two bits. BRENDA needs a licensed multi-GB dump
   that is absent, and the thermodynamic call is handed pathway map identifiers
   rather than reactions, so both terms are constant. The score should reduce
   exactly to size-rule AND database-presence.
4. Hard negatives will be called feasible, because many drugs appear in KEGG via
   xenobiotic degradation and drug-metabolism maps. Presence in KEGG means the
   cell *encounters* the compound, not that it *makes* it.
