# Post-fix measurement

Re-run of the fixed pipeline on the **same 36 compounds, same ground truth** as
the original audit ([DESIGN.md](DESIGN.md), [RESULTS.md](RESULTS.md)). The
baseline is not re-derived: it is read from `stage2.json`, recorded during the
original run, so it cannot drift to flatter the new code.

Reproduce with `python benchmarks/run_fixed.py` (writes `fixed.json`).

## Headline

| metric | shipped | fixed | note |
|---|---|---|---|
| balanced accuracy | 86% | **92%** | |
| sensitivity | 90% | 90% | unchanged by design |
| specificity | 81% | **94%** | false positives 3 → 1 |
| retrieval top-1 | 94% | **100%** | and stereochemically exact, 32/32 |
| distinct output values | 2 | 5 verdict classes + 4 stereo outcomes | |
| compounds given a stereochemical answer | 0 | 35/35 chiral | |
| named selectivity risks | 0 | 7 | |

The sensitivity figure is deliberately flat. The gain is all in specificity:
the original passed three compounds it should not have, and two of those were
passed on evidence that does not mean what the code assumed.

## What fixed each false positive

| compound | why the old code passed it | why the new code does not |
|---|---|---|
| levofloxacin | KEGG lists it on a bile-secretion map, and "any pathway" was the test | that map classifies as `physiology_or_disease`, not biosynthesis |
| citalopram | listed via "Drug metabolism — cytochrome P450" | classifies as `xenobiotic_degradation`: the cell *breaks it down* |
| camptothecin | still a false positive | a plant alkaloid with a curated biosynthesis map, 348 Da, so the size discount does not catch it. Honest miss. |

## The trap in removing the molecular-weight gate

The original verdict was exactly `MW <= 500 AND has-any-KEGG-pathway`. MW alone
scores 57% balanced accuracy with 13 false positives, so it is worthless as a
feasibility signal — but deleting it outright made things worse, not better:

| rule | sens | spec | balanced | FP | FN |
|---|---|---|---|---|---|
| any pathway AND MW ≤ 500 (shipped) | 90% | 81% | 86% | 3 | 2 |
| biosynthesis pathway documented anywhere | 95% | 81% | 88% | 3 | 1 |
| biosynthesis pathway AND present in host | 60% | 100% | 80% | 0 | 8 |
| **biosynthesis AND (in host OR within size gate)** | 90% | 94% | **92%** | 1 | 2 |

Dropping the size gate let paclitaxel (854 Da) and vinblastine (811 Da) through:
both have curated plant biosynthesis and no microbial route. Requiring the route
to be in the host's own KEGG map set fixes specificity completely and destroys
sensitivity, because **most successes in this field are heterologous expression**
and are therefore absent from the host's native maps — the strict rule penalises
exactly the engineering achievements the benchmark counts as positives.

So molecular weight is kept, but demoted: it never asserts feasibility, and is
consulted only to discount a route that exists in nature but not in this host,
as a proxy for how many dedicated enzymes would have to be transplanted.

## Remaining errors, named

Three compounds are still wrong, and all three are explicable rather than
mysterious:

- **artemisinic acid** (false negative). KEGG's record C20309 carries no
  biosynthesis-class pathway, though the yeast route is famous. A gap in KEGG's
  curation, not in the logic.
- **strictosidine** (false negative). 530 Da, route documented but not in the
  host, so the size discount excludes it. This is the known cost of that rule.
- **camptothecin** (false positive). As above.

## Stereochemistry, which the old pipeline had none of

On the 35 chiral compounds: 17 `EXACT`, 1 `WRONG-FORM`, 17 `NONE`, and 7 carry a
named enantiocomplementary selectivity risk. The original gave byte-identical
output to 9/9 enantiomer pairs.

The `WRONG-FORM` class is the one that changes a decision. KEGG pathway links
attach to a **connectivity skeleton**, so L-glucose inherits D-glucose's maps
and the pathway evidence reads as a documented route in host. It is not one.
A stereochemistry-aware verdict is what stops that, and without it the headline
verdict would cover two opposite situations — the same non-identifiability the
original `confidence` field had.

## What this still does not establish

- The LLM discovery stage is **not** benchmarked. It needs paid API calls and
  repeated sampling to characterise. The structure gate now rejects its known
  bad output, which is a different and weaker claim than measuring it.
- n = 36, so each proportion carries roughly ±0.12. A 6-point balanced-accuracy
  gain is inside that band for a single compound's worth of movement; the
  specificity gain rests on two compounds whose mechanism is identified above,
  which is why the mechanism is reported alongside the number.
- Nothing here predicts yield, titre or enantiomeric excess. No such code
  exists, and none is claimed.
