# ChiralAI — operating notes

Natural-language query → candidate chiral molecules → biosynthetic feasibility
in a microbial host. Stage A is an LLM call; stage B is database-grounded.

## Locked findings — audit of the ORIGINAL pipeline (2026-10-07)

All of these were fixed in the 2026-10-21 rewrite below. Kept because they are
what the defects looked like, and the regression tests are named after them.

Full write-up: [benchmarks/RESULTS.md](benchmarks/RESULTS.md).
Design and pre-registered predictions: [benchmarks/DESIGN.md](benchmarks/DESIGN.md).

| finding | evidence |
|---|---|
| Feasibility verdict is exactly `MW<=500 AND has-any-KEGG-pathway` | identity holds 36/36 compounds |
| Its 86% balanced accuracy belongs to KEGG, not to the pipeline | every layer above the two lookups is inert on this set |
| `confidence` cannot rank, so it trails a raw pathway count | AUROC 0.86 vs 0.94; a 2-valued score forfeits AUROC by construction |
| `confidence` is one ambiguous bit, not a probability | 2 distinct values over 36; 0.88 means both "too large" and "no route" |
| BRENDA layer has never worked for anyone | hardcoded download URL returns 404 |
| Thermodynamics has never produced a number | 299 calls, 0 results; receives `path:mapXXXXX` where a reaction is expected |
| System is stereo-blind | 9/9 enantiomer pairs get byte-identical analysis; 0/9 name/structure mismatches caught |
| No quantitative output exists | no yield/titre/ee/selectivity/flux code anywhere; output is 1 boolean + a 2-valued score |
| Host choice has no effect | 0/36 compounds changed any field between `E. coli` and `yeast`; `eco`/`sce` codes defined but never used |
| Name-based retrieval fails exactly on chirality | stereochemical queries return nothing 25/32 times |

## Measured improvements identified by that audit (all now shipped)

1. **Retrieve by structure, not name.** Exact formula search then InChIKey match:
   94% → 100% top-1, and stereochemically correct by construction.
2. **Require a biosynthesis pathway, not any pathway.** KEGG links drugs via
   degradation and physiology maps; presence means the cell *encounters* the
   compound, not that it *makes* it. Specificity 81% → 94%, false "feasible"
   3 → 1, costing 2 sensitivity points. Expose it as a threshold, not a hardcode.
3. **Cross-check name against structure** and report which enantiomer was assessed.

## Why we don't do X

- **Don't trust first-hit KEGG name search.** It returned alpha-farnesene for
  "farnesene" when the yeast product is beta. Take all candidates, pick by structure.
- **Don't treat KEGG presence as biosynthesability.** That is what passes
  levofloxacin (via a bile-secretion map) and citalopram (via a P450
  drug-degradation map).
- **Don't report `confidence` until a layer under it works.** It currently
  implies evidence that was never gathered.
- **Don't describe the tool as host-specific or as scoring selectivity.**
  Measured: host changes 0/36 outputs, and selectivity is not computed at all.
  Both claims are falsifiable in one command.
- **Don't rely on MW<=500 as a biosynthesis signal.** Alone it is worthless
  (balanced accuracy 57%, 13 false positives) and it wrongly rejects
  strictosidine (530 Da, made in yeast 2015). It does, however, carry the size
  gate that stops large plant natural products slipping through a raw pathway
  lookup, so replace it deliberately rather than just deleting it.

## The wedge: stereochemistry-resolved enzyme selection

Built after the benchmark showed the feasibility module adds nothing over two
database lookups and is stereo-blind. Full write-up: [stereo/README.md](stereo/README.md).

The useful question is not "is this feasible" but **"which enzyme makes THIS
enantiomer, and which makes the one I don't want?"** Nobody had published that
join: Rhea curates stereochemically explicit reactions with EC/UniProt links but
is organised by reaction, KEGG cannot retrieve a specified enantiomer by name,
BRENDA's stereoselectivity is free text behind a licence, and retrosynthesis
rule sets are usually stereo-stripped.

| artifact | contents |
|---|---|
| `stereo/data/stereo_enzyme_index.tsv` | 3,357 stereo-defined compounds with an EC acting on them |
| `stereo/data/enantiocomplementary_pairs.tsv` | 174 pairs, of which 57 are true enantiomer pairs |
| `stereo/query.py` | SMILES in, enzymes out, wrong-enantiomer enzymes flagged |

Validation: 13/15 recall on literature pairs held out before inspecting output
(one of the two misses was an error in my own held-out list, so 13/14 on
correctly specified pairs); 57/57 enantiomer pairs independently confirmed as
mirror images by RDKit inversion; both negative controls absent.

The four query outcomes matter more than the counts. `WRONG-FORM` means the
scaffold is known but only the other stereoisomer, which is precisely what every
name-matching tool reports as success. That is the limonene error from this
repo's own saved output, now caught.

### Two traps worth remembering

- **Match InChIKey blocks 1+2, never the full key.** Block 3 is protonation.
  ChEBI stores physiological ions, chemists draw neutral molecules, so
  (S)-lactic acid and ChEBI's (S)-lactate differ in the last block. Full-key
  matching returns a false `WRONG-FORM` on a correct query.
- **Chiral cofactors break naive stereo-pair finding.** Requiring one chiral
  participant per reaction excluded every dehydrogenase, because NAD/NADH have
  their own stereocentres. Lactate dehydrogenase went missing. Compare which
  single participant was *swapped* instead: 108 pairs became 175. Always assert
  the canonical case explicitly in validation.

## The rewrite (2026-10-21) — all blockers cleared

Post-fix measurement: [benchmarks/FIXED.md](benchmarks/FIXED.md). Same 36
compounds, same ground truth; baseline read from the recorded `stage2.json` so
it cannot drift.

| | shipped | now |
|---|---|---|
| balanced accuracy | 86% | 92% |
| specificity | 81% | 94% (FP 3 → 1) |
| retrieval top-1 | 94% | 100%, stereo-exact 32/32 |
| stereochemical answers | 0 | 35/35 chiral compounds, 7 named selectivity risks |
| thermodynamic values | 0 of 299 calls | works, verified against literature sign |

What changed, by defect:

| was | now |
|---|---|
| CLI could not start (bad import, `+42`/`-15` diff markers, `openai.ChatCompletion`) | runs; `tests/test_pipeline.py` pins each |
| LLM output trusted | `ChiraLLM/structure_gate.py` checks formula vs KEGG *before* assessment |
| retrieval by name | `ChiraLLM/retrieval.py`: exact formula → InChIKey |
| any KEGG pathway = feasible | pathways classified; biosynthesis-only by default, `--any-pathway` restores the old behaviour |
| `confidence` (2 values, non-identifiable) | **deleted**; explicit `evidence` record + named `reasons` |
| host argument inert | host's own KEGG map set intersected with the compound's biosynthesis maps |
| stereo-blind | `stereo/api.py` wired in; verdict is stereochemistry-aware |
| BRENDA/thermo silently empty | report `unavailable` with the reason in `layers_unavailable` |
| 3 modules duplicating KEGG logic, imported by nothing | deleted; one `utils/kegg.py` |
| README claimed CRISPR/plasmid/ML output | rewritten: measured vs roadmap |

### Why we don't do X (new entries)

- **Don't delete the MW gate, demote it.** Removing it outright let paclitaxel
  (854 Da) and vinblastine (811 Da) through on curated *plant* biosynthesis.
  Requiring the route to be in the host's native maps gives 100% specificity and
  60% sensitivity, because most wins in this field *are* heterologous expression
  and so are absent from native maps. Final rule: biosynthesis route AND
  (in host OR within size gate) → 92% balanced. Table in FIXED.md.
- **Don't judge a structure against the FIRST KEGG name hit**, and don't judge it
  against a substring match either. "glucose" returns UDP-glucose first (rejects
  correct glucose); "L-glucose" only matches GDP-L-glucose, so a substring
  comparison rejects correct L-glucose against a nucleotide sugar's formula. Only
  compare against candidates where the queried name is a **full synonym**; with
  none, report UNVERIFIED rather than reject.
- **Don't let a pathway verdict outrank stereochemistry.** KEGG pathway links
  attach to a connectivity skeleton, so L-glucose inherits D-glucose's maps and
  reads as "route documented in host". `WRONG-FORM` must override the verdict, or
  one output covers two opposite situations — the exact fault of the old
  `confidence`.
- **Don't wrap a whole function in `except Exception: return None`.** `thermo.py`
  hid three separate bugs that way (wrong input type, unit-less `p_h`, a
  module-level `ComponentContribution()` that downloaded a huge cache on import).
  A function that cannot fail visibly cannot be debugged.
- **Don't trust a mock over a contract.** The old `test_thermo.py` defined a
  `DummyCC` accepting bare floats for `p_h` and returning `.to(unit).m`. Both are
  wrong for the real API, so the suite stayed green through all three bugs. The
  mock encoded the bug.
- **Don't report a central value with a huge uncertainty.** eQuilibrator returns
  e.g. −87 kJ/mol ±100000 when a participant is outside its coverage. Reporting
  the number alone manufactures precision; `status: not-estimable` instead.
- **`*.csv` in .gitignore was excluding `benchmarks/dataset.csv`** — the ground
  truth the whole benchmark rests on was never committable. Negated explicitly.

## Running the benchmark

```bash
python -m venv venv && ./venv/bin/pip install -r requirements.txt
./venv/bin/python -m pytest tests/ -q                  # 36 regression tests
./venv/bin/python benchmarks/run_fixed.py              # fixed vs shipped
./venv/bin/python benchmarks/run_benchmark.py          # claims 1–3
./venv/bin/python benchmarks/run_stereo_naming.py      # validated enantiomers + naming
./venv/bin/python benchmarks/pathway_classes.py        # pathway-class arm
./venv/bin/python benchmarks/structure_first.py        # structure-first retrieval
./venv/bin/python benchmarks/analyze.py                # writes METRICS.md
```

KEGG responses cache under `benchmarks/.kegg_cache` and `benchmarks/.http_cache`,
so reruns are offline and free. Delete to refresh.
