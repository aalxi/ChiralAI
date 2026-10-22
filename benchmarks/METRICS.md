## Claim 1 - identifier grounding

Scored on the 32 of 36 compounds that have a KEGG COMPOUND entry.

| metric | value |
|---|---|
| coverage (any hit) | 32/32 = 100% |
| top-1 exact | 30/32 = 94% |
| right skeleton, wrong stereoisomer | 0 |
| coverage, multi-word names (n=1) | 1/1 |
| coverage, single-word names (n=31) | 31/31 |

Wrong retrievals:

- strictosidine: returned C03309, correct is C03470
- farnesene: returned C09665, correct is C09666

## Claim 2 - feasibility discrimination

### Full set  (n=36, 20 positive / 16 negative)

| model | AUROC |
|---|---|
| **ChiralAI confidence** | **0.856**  (95% CI 0.73-0.97) |
| baseline: present in KEGG | 0.625 |
| baseline: KEGG pathway count | 0.939 |
| baseline: molecular weight alone | 0.850 |
| baseline: the code's own MW<=500 rule | 0.569 |
| baseline: constant | 0.500 |

At the system's own operating point (`biosynthetic_feasibility`):

| | value |
|---|---|
| sensitivity | 90% (18/20) |
| specificity | 81% (13/16) |
| precision | 86% |
| balanced accuracy | 86% |
| false 'feasible' calls | 3 |

### Size-matched subset (drops stratum C, the only group above 500 Da)  (n=32, 20 positive / 12 negative)

| model | AUROC |
|---|---|
| **ChiralAI confidence** | **0.867**  (95% CI 0.72-0.98) |
| baseline: present in KEGG | 0.667 |
| baseline: KEGG pathway count | 0.971 |
| baseline: molecular weight alone | 0.808 |
| baseline: the code's own MW<=500 rule | 0.475 |
| baseline: constant | 0.500 |

At the system's own operating point (`biosynthetic_feasibility`):

| | value |
|---|---|
| sensitivity | 90% (18/20) |
| specificity | 83% (10/12) |
| precision | 90% |
| balanced accuracy | 87% |
| false 'feasible' calls | 2 |

Distinct `confidence` values across all 36 compounds: 2 -> [0.88, 0.95]

Mean confidence, positives 0.943 / negatives 0.893

BRENDA returned a hit for 0/36 compounds.
Thermodynamic function called 299 times; returned a number 0 times.
Example arguments it received: ['path:map00350', 'path:map00950', 'path:map00965', 'path:map01060']

False 'feasible' calls, by stratum:

- stratum B: 2/12 called feasible -> levofloxacin, citalopram
- stratum C: 1/4 called feasible -> camptothecin

- stratum A1 (true positives): 12/12 called feasible
- stratum A2 (true positives): 6/8 called feasible

## Claim 3 - stereochemistry

Enantiomer pairs tested: 6. Each pair is two different KEGG compounds with mirror-image structures.

| pair | inputs differ | identical verdict | fields that differed |
|---|---|---|---|
| limonene | yes | no | kegg_compound_id |
| DOPA | yes | no | biosynthetic_feasibility, bottlenecks, confidence, kegg_compound_id |
| lactate | yes | no | kegg_compound_id |
| alanine | yes | no | kegg_compound_id |
| malate | yes | no | kegg_compound_id |
| tryptophan | yes | no | biosynthetic_feasibility, bottlenecks, confidence, kegg_compound_id |

Pairs given an identical assessment despite being different molecules: 0/6.
Pairs whose structural analysis block was byte-identical: 5/6.

