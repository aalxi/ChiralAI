# Retrospective benchmark: results

Run 2026-10-07. Design and pre-registered predictions: [DESIGN.md](DESIGN.md).
Reference set: [dataset.csv](dataset.csv). Raw metric dump: [METRICS.md](METRICS.md).

Reproduce with `run_benchmark.py`, `run_stereo_naming.py`, `pathway_classes.py`,
`structure_first.py`. All four call the repository's own functions. Structures and
pathway annotations come from KEGG, so the reference set is auditable rather than
asserted.

## What was tested

36 compounds: 20 with documented production in a microbial host, 16 without.
35 of the 36 have at least one stereocentre. Negatives are deliberately hard,
mostly chiral drugs between 200 and 410 Da, so they pass the pipeline's size
filter and probe whether it knows any biology. One stratum of four plant natural
products above 780 Da is reported separately because it is the only group the
size filter can act on.

At this sample size the 95% interval on an AUROC near 0.86 spans roughly ±0.12.
That separates "informative" from "coin flip". It cannot rank two close methods,
and no such claim is made below.

## Headline

The feasibility module reaches 86% balanced accuracy, and that number is real.
But its verdict is **exactly** reproduced, on all 36 compounds with no
exceptions, by two lines of database lookup:

```
(molecular_weight <= 500) AND (KEGG lists at least one pathway)
```

Everything layered above those two lookups is provably inert on this set. The
enzyme lookup returned nothing 36 times out of 36. The thermodynamic function
was called 299 times and returned a number zero times. The `confidence` output
took two distinct values. None of those components could have influenced a
single verdict.

So the honest reading is not that the pipeline fails. It is that the pipeline's
accuracy belongs to KEGG, and the parts that would make it more than a KEGG
query are not running. A differently chosen pair of lookups from the same free
database does better on the metric that matters most here, cutting false
"feasible" calls from three to one.

## Claim 1 — identifier grounding

Scored on the 32 compounds that have a KEGG COMPOUND entry.

| retrieval method | top-1 exact | coverage |
|---|---|---|
| name-first, canonical name (as shipped) | 30/32 = 94% | 32/32 |
| name-first, loose name | 28/32 = 88% | 32/32 |
| name-first, stereochemical name | — | **7/32 = 22%** |
| **structure-first (formula + InChIKey)** | **32/32 = 100%** | 32/32 |

Retrieval works well when the user happens to type the name KEGG chose. It
fails when they specify which enantiomer they want. Asking for "(S)-tryptophan"
returns nothing, because KEGG indexes that molecule as "L-Tryptophan" and the
lookup compares strings, with no notion that the two names denote the same
compound. 25 of 32 stereochemical queries returned nothing at all.

The two outright errors come from taking the first line of the result. The query
"farnesene" returns alpha-farnesene (C09665) when the compound produced
industrially in yeast is the beta isomer (C09666). The query "strictosidine"
returns C03309 rather than C03470.

Structure-first retrieval resolved all 32 compounds, and all 32 to the correct
stereoisomer, because an InChIKey distinguishes enantiomers by construction.
This is the same free KEGG API, queried by molecule instead of by name.

## Claim 2 — feasibility discrimination

Two comparisons are needed, and they say different things.

**At matched binary operating points**, which is the apples-to-apples view:

| decision rule | sens | spec | precision | balanced acc | false "feasible" |
|---|---|---|---|---|---|
| ChiralAI as shipped | 90% | 81% | 86% | 86% | 3 |
| has ≥1 KEGG pathway (no size gate) | 95% | 62% | 76% | 79% | 6 |
| MW ≤ 500 only | 95% | 19% | 59% | 57% | 13 |
| has ≥1 biosynthesis map only | 85% | 81% | 85% | 83% | 3 |
| **biosynthesis map AND MW ≤ 500** | 80% | **94%** | **94%** | **87%** | **1** |

The shipped combination beats either of its two halves alone, so the size gate
is doing real work: it removes the four large plant natural products that a raw
pathway lookup waves through. Filtering pathways by what the map describes then
trades four points of sensitivity for thirteen of specificity, which is the
right direction when a false "feasible" costs a wet-lab campaign.

**As a ranking score**, measured by AUROC:

| model | full set | size-matched |
|---|---|---|
| ChiralAI `confidence` | 0.856 (95% CI 0.73–0.97) | 0.867 (0.72–0.98) |
| KEGG pathway count (continuous) | 0.939 | 0.971 |
| biosynthesis-map count (continuous) | 0.887 | — |
| molecular weight alone | 0.850 | 0.808 |
| the code's own MW ≤ 500 rule | 0.569 | 0.475 |
| constant predictor | 0.500 | 0.500 |

The raw pathway count ranks better than `confidence`. That gap is not evidence
of a cleverer model, it is mostly an artifact: a score with two distinct values
cannot rank 36 compounds, so it forfeits AUROC by construction. The finding is
about `confidence` being uninformative as a score, not about pathway counting
being a good predictor. Read the operating-point table for how the rules
actually compare.

### The verdict is two database lookups

For all 36 compounds, with no exceptions, `biosynthetic_feasibility` equals

```
(molecular_weight <= 500) AND (KEGG lists at least one pathway)
```

That identity is exact, 36/36. The remaining machinery contributes nothing
because each of its terms is constant across the whole set:

| component | measured behaviour |
|---|---|
| BRENDA enzyme lookup | returned a hit for **0 of 36** compounds |
| thermodynamics | called **299 times**, returned a number **0 times** |
| `confidence` | took exactly **two** distinct values, 0.88 and 0.95 |

The enzyme layer is silent because the download URL hardcoded in
`utils/brenda.py` returns HTTP 404. No user of this repository has ever had that
feature work.

The thermodynamic layer is silent because it is handed the wrong kind of object.
It receives pathway map identifiers such as `path:map00350`, where
eQuilibrator expects a reaction equation. Every call fails into a bare
`except` that returns `None`. The dependency is also absent from the
environment, so this code path has never executed successfully.

### `confidence` is one bit, and it is ambiguous

Two values over 36 compounds is not a probability. Worse, 0.88 is emitted for
two opposite diagnoses: a molecule too large for the size rule that *does* have
a known pathway, and a molecule small enough that has *no* known pathway. A user
reading 0.88 cannot tell whether the system means "too complex" or "no known
route". Those call for different next experiments.

### Why the errors happen

The three false positives are not near-misses. They pass because KEGG links
them to a pathway map that has nothing to do with making them:

| compound | the pathway that passed it |
|---|---|
| levofloxacin | Bile secretion |
| citalopram | Drug metabolism — cytochrome P450 |
| camptothecin | alkaloid biosynthesis (a plant route, not microbial) |

This is a conceptual error, not a tuning problem. Presence in KEGG means a cell
*encounters* a compound. The pipeline reads it as evidence the cell can *make*
it. A drug-degradation map is evidence of the opposite.

The two false negatives are more damaging than the false positives. Artemisinic
acid is rejected because KEGG happens to list no pathway link for C20309,
although semisynthetic artemisinin from engineered yeast is the most cited
success in the field. Strictosidine is rejected for weighing 530 Da, having been
made in yeast in 2015. Meanwhile paclitaxel is rejected for the right answer but
the wrong reason: it is a genuine natural product that the size rule happens to
exclude at 854 Da.

### A better signal exists in the same database

Classifying every pathway link by what the map actually describes, then
requiring a biosynthesis map rather than any map:

| signal | sens | specificity | precision | false positives |
|---|---|---|---|---|
| as shipped | 90% | 81% | 86% | 3 |
| biosynthesis maps only, plus size rule | 80% | **94%** | **94%** | **1** |

Same database, same cost, one filter on the pathway class. The fix is to ask
KEGG whether anything *makes* the compound, instead of whether KEGG has heard
of it. It is a genuine trade: two more compounds are wrongly rejected in
exchange for two fewer wrongly accepted. For a tool that sends people into a
six-month strain-engineering campaign, that trade is worth taking, but it is a
choice and should be exposed as a threshold rather than hardcoded.

## Is any output quantitative?

Asked because the project describes itself as checking "whether the pathway
works in E. coli, how selectively each step makes the correct form". Measured
against the code, neither clause holds.

| output field | distinct values over 36 compounds | is it a prediction? |
|---|---|---|
| `biosynthetic_feasibility` | 2 (true/false) | a verdict, not a quantity |
| `confidence` | **2** (0.88, 0.95) | a 4-item checklist through a sigmoid |
| `molecular_weight` | 36 | a property of the input, computed by RDKit |
| `num_chiral_centers` | 10 | a count of the input's centres, not of selectivity |
| `delta_g` | **1** (always `None`) | never produced a number in 299 calls |

There is no predicted yield, titre, enantiomeric excess, selectivity, flux or
conversion anywhere in the pipeline. A search of the source for those terms
returns no matches outside Python's own `yield` keyword in the BRENDA iterator.
So the output is qualitative: one boolean plus a two-valued score.

### "works in E. coli" is not host-specific

Every compound was run twice, once with host `E. coli` and once with `yeast`.

| measurement | result |
|---|---|
| compounds where any substantive output field changed | **0 / 36** |

The host string is consumed in exactly one place,
`ChiraLLM/feasibility.py:91`, as `host in ORGANISMS`. That is a membership test
on a two-key dictionary, so it checks that the user typed a supported host name.
The KEGG organism codes the dictionary maps to, `eco` and `sce`, are defined and
then never referenced anywhere in the codebase. No organism-specific pathway,
genome or enzyme query is ever issued. KEGG is queried with the species-neutral
`map` identifiers throughout.

### "how selectively each step makes the correct form" is not computed

Selectivity per step would require knowing the enzymes in the route and their
stereochemical preference. The enzyme lookup returned nothing 36 times out of
36, and the `pathway` field is not a route: it renumbers the KEGG pathway maps a
compound appears in as "step 1" through "step N". For L-DOPA that yields 13
steps, one of which is the global metabolism map. No step is an enzymatic
reaction, so no step has a selectivity.

The stereochemical measurement in Claim 3 below closes this off: the pipeline
cannot distinguish two enantiomers from their structures at all, so it has no
basis on which to score which form is produced.

## Claim 3 — stereochemistry

Nine enantiomer pairs, each validated programmatically as a true pair: identical
flat skeleton, inverted stereocentre assignments.

| measurement | result |
|---|---|
| pairs whose structural analysis was byte-identical | **9 / 9** |
| pairs where a correct name with the wrong enantiomer's structure was flagged | **0 / 9** |

The structural layer reports molecular weight and a count of stereocentres. It
never reports which configuration. So (R)- and (S)-limonene receive the same
analysis, as do L- and D-tryptophan. The only field that differs between
enantiomers is the KEGG identifier, and that is derived from the name string,
never from the structure.

The consequence is the practical one. Nothing cross-checks the name against the
structure. Feeding the correct name together with the mirror-image structure
produced output identical to the correct input in all nine pairs. This is the
mechanism behind the repository's own saved limonene error, where the name and
the structure disagreed and nothing noticed. `validate_chirality()` exists in
`ChiraLLM/chirality_checker.py` and is imported by nothing.

For a tool whose premise is that the arrangement of atoms matters, the
stereochemistry is the one thing it does not check.

## Pre-registered predictions, scored honestly

| prediction | outcome |
|---|---|
| 1. Retrieval collapses for multi-word names | **Refuted.** `requests` percent-encodes spaces and KEGG handles them. The real collapse is for stereochemical names, 22% coverage, which prediction 1 missed. |
| 2. First-hit retrieval returns arbitrary enantiomers | **Supported**, and sharper than stated: 22% coverage for stereochemical queries, plus the farnesene isomer error. |
| 3. `confidence` carries at most two bits | **Confirmed exactly.** Two distinct values; BRENDA 0/36; thermodynamics 0 numbers from 299 calls. |
| 4. Hard negatives pass via xenobiotic and drug-metabolism maps | **Confirmed in mechanism, overstated in scale.** Only 2 of 12 hard negatives passed, but both by exactly the predicted route. |

## Defects that block running the pipeline at all

Found while wiring the harness, independent of the benchmark.

1. `main.py:4` imports `save_suggestions_to_csv` from `utils/file_saver.py`,
   which defines only `save_data_to_csv`. The command-line entry point raises
   ImportError before doing anything.
2. The prompt instructs the model to emit a JSON key `SMILES`, and
   `ChiraLLM/feasibility.py:174` reads `r["smiles"]`. When the model honours the
   prompt, every row is discarded as "Missing name or SMILES". This is the
   recorded failure in `suggestions_20250702_155147.csv`. The June run succeeded
   only because that run happened to emit lowercase keys.
3. `ChiraLLM/feasibility.py` lines 1–2 are the literal text `+42` and `-15`,
   leftover diff markers. They parse as valid Python expressions, so they are
   silent.

## What this benchmark does not establish

- Labels come from literature a language model has read. Stratum A1, the classic
  positives, scored 12/12 and stratum A2, the post-2015 ones, scored 6/8. The
  gap is in the right direction for memorisation but far too small a sample to
  claim it.
- The discovery stage was not run. It needs a paid API call per query and is
  stochastic at temperature 0.7, so a fair measurement needs repeated sampling.
  The structural errors in the repository's own saved output, a several-hundred
  character degenerate string offered as limonene, and a tricarboxylic structure
  offered as lactic acid, suggest a validity gate is needed before that stage is
  worth benchmarking.
- Negatives are labelled "no reported microbial route", which is an absence of
  evidence. A compound could move to the positive class next year.

## Where this leaves the project

The integration runs end to end and the database plumbing is sound. The 86%
balanced accuracy is a real result and better than either half of its own rule.
What is missing is that no layer *above* the two database lookups currently
earns its place, and the one thing the project is named for is unchecked.

Three changes are each measured above, not speculative:

1. Retrieve by structure, not by name. 94% to 100% top-1, and stereochemically
   correct by construction.
2. Require a biosynthesis pathway rather than any pathway. Specificity 81% to
   94%, false "feasible" calls 3 to 1, at a cost of two sensitivity points.
3. Cross-check the name against the structure, and report which enantiomer was
   assessed. Currently 0 of 9 mismatches are caught.

Until the enzyme and thermodynamic layers either work or are removed, reporting
a `confidence` number implies evidence that was never gathered.
