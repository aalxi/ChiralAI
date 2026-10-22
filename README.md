# ChiralAI

**For a chiral target: which enzyme makes the form you want, and which one makes
the mirror image you don't?**

---

## The problem

Most drugs, agrochemicals and advanced materials are chiral: the molecule exists
in two mirror-image forms, and usually only one of them works. The other is
inactive, hard to separate, and occasionally toxic. Enantiopurity is the whole
commercial point.

Biology is extremely good at this — enzymes are chiral catalysts and routinely
deliver one enantiomer cleanly. But that only helps if you pick the right
enzyme. Many scaffolds have **enantiocomplementary** enzyme families: two
different EC numbers acting on the same skeleton, each producing the opposite
form. Pick the wrong one and you don't get less product, you get the wrong
product.

So the question that decides a project is narrow and specific: *for this exact
stereoisomer, which enzyme families are documented to make it, and which
documented families would make its mirror image instead?*

## What you get

```
$ python -m stereo.query --smiles 'CC1=CC[C@@H](CC1)C(=C)C'

target            : C=C(C)[C@H]1CC=C(C)CC1
InChIKey          : XMGQYMWWDOXHJM-JTQLQIEISA-N
stereocentres     : 1 total, 1 defined

RESULT: EXACT      — a curated reaction acts on this exact stereoisomer.

  compound      : (4R)-limonene (CHEBI:15382)
  makes it      : EC:1.14.13.107, EC:1.14.14.53, EC:1.17.99.8, EC:4.2.3.20
  Rhea          : RHEA:10940;RHEA:18957;RHEA:26093;RHEA:26097;RHEA:62560

  SELECTIVITY RISK — the enantiomer is also enzymatically accessible:
    (4S)-limonene (CHEBI:15383)
    made by     : EC:4.2.3.16
    Choosing between these two EC families IS the stereochemical decision.
```

Four outcomes, and the distinction between them is the product:

| outcome | meaning | what to do |
|---|---|---|
| `EXACT` | a curated reaction acts on your exact stereoisomer | pick among the EC families listed; check the selectivity risk |
| `WRONG-FORM` | the scaffold is curated, but **only as a different stereoisomer** | there is no documented natural route to your target |
| `SKELETON` | you didn't specify stereochemistry | specify it; no enantiomer-specific route can be chosen without it |
| `NONE` | no curated reaction for this connectivity | absence of evidence, not evidence of absence |

`WRONG-FORM` is the one that changes a decision, because it is exactly what a
name-matching tool reports as success. Ask for L-glucose and a name-based
lookup finds D-glucose's enzymes and pathways and tells you it's feasible. It
isn't.

## The gap this fills

The data exists. The *join* didn't. Each resource below is excellent and none of
them answers "which enzyme makes this enantiomer":

| resource | what it does well | where it falls short for this question |
|---|---|---|
| **Rhea** | reactions with stereochemically explicit participants, linked to EC and UniProt | organised **by reaction**. You can ask what a reaction does; you cannot ask which form of a compound you can make |
| **KEGG** | pathway context, host organism maps, free REST API | can't reliably retrieve a *named* enantiomer — only 22% of stereochemical names resolve (measured). Worse, pathway links attach to a **connectivity skeleton**, so enantiomers share them and the two forms look identical |
| **BRENDA** | the deepest enzyme-property resource there is | stereoselectivity is largely free text, and the bulk download sits behind a registered licence |
| **MetaCyc** | high-quality curated pathways | licensed for bulk use, and not organised by product stereochemistry |
| **RetroRules / RetroPath** | retrosynthesis rules for pathway design | reaction rules are generally stereochemistry-stripped or generalised, which is the information you need |

ChiralAI joins Rhea's stereochemically explicit reactions to ChEBI structures,
then **inverts the index**: from "what does this reaction do" to "which enzymes
make this exact stereoisomer, and which make its mirror image". That inversion
is the contribution. The underlying curation is EBI's, and every row traces back
to a Rhea reaction and ChEBI identifier a reader can check.

**3,357** stereo-defined compounds with an enzyme acting on them · **174**
stereocomplementary pairs, of which **57** are structurally verified enantiomer
pairs across 50 scaffolds.

Validation: 13/15 recall on literature enantiocomplementary pairs written down
*before* inspecting the output; 57/57 enantiomer labels independently re-derived
by inverting every stereocentre in RDKit, so nothing is labelled an enantiomer
pair on the strength of its name; both negative controls absent. Method and
limits: [stereo/README.md](stereo/README.md).

## Install

```bash
git clone https://github.com/aalxi/ChiralAI.git && cd ChiralAI
python -m venv venv && ./venv/bin/pip install -r requirements.txt
```

Nothing above needs an API key or a licence. Only the optional
natural-language stage calls OpenAI; only `--thermo` needs `equilibrator-api`.

## Use

```bash
# 1. which enzyme makes this enantiomer?
./venv/bin/python -m stereo.query --smiles 'CC1=CC[C@@H](CC1)C(=C)C'

# ...and the UniProt accessions to order
./venv/bin/python -m stereo.query --smiles 'CC1=CC[C@@H](CC1)C(=C)C' --sequences

# browse what's indexed under a name
./venv/bin/python -m stereo.query --name limonene

# 2. full assessment: structure check, pathway evidence, host, stereochemistry
./venv/bin/python main.py --smiles 'C[C@H](O)C(=O)O' --name 'L-lactic acid' --host 'E. coli'

# with reaction thermodynamics
./venv/bin/python main.py --smiles 'C[C@@H](O)C(=O)O' --name 'D-lactic acid' --thermo

# 3. optional: natural-language discovery (needs OPENAI_API_KEY in .env)
./venv/bin/python main.py --query 'chiral monomer for a biodegradable polyester'
```

The full assessment adds, on top of the stereochemical answer:

- **A structure check.** Does the molecule match the name it was given? Verified
  against KEGG by molecular formula. The language-model stage is never trusted —
  it has produced an 896-character degenerate string labelled "limonene" and a
  tricarboxylic acid labelled "lactic acid", both of which are now rejected.
- **Pathway evidence that distinguishes making from encountering.** A compound
  on a drug-degradation map is *broken down* by the cell, not made by it.
- **Host specificity that is real.** The compound's biosynthesis maps are
  intersected with the maps the host organism actually encodes. E. coli encodes
  138, yeast 142, and only 95 are shared.
- **An explicit list of what was not consulted, and why.**

No field is a prediction. Each is either a curated database record or a stated
absence of one.

## Measured

Retrospective benchmark, 36 compounds, pre-registered predictions and hard
negatives. Baseline figures are read from the original recorded run, so they
can't drift to flatter the current code.

| | shipped (July 2026) | now |
|---|---|---|
| balanced accuracy | 86% | **92%** |
| specificity | 81% | **94%** (false positives 3 → 1) |
| KEGG retrieval, top-1 correct | 94% | **100%**, stereochemically exact |
| enantiomer pairs given different answers | 0 of 9 | all |
| thermodynamic values ever returned | 0 of 299 calls | works |
| distinct output values | 2 | 5 verdict classes × 4 stereo outcomes |

Audit of the original, including a refuted prediction of mine:
[benchmarks/RESULTS.md](benchmarks/RESULTS.md) ·
design and power: [benchmarks/DESIGN.md](benchmarks/DESIGN.md) ·
post-fix measurement and the three errors that remain:
[benchmarks/FIXED.md](benchmarks/FIXED.md)

## What it does not do

Stated plainly, because an earlier version of this README implied otherwise:

- **No yield, titre, flux or enantiomeric-excess prediction.** No such code
  exists and no number is claimed.
- **No plasmid, CRISPR or codon-optimisation output.**
- **No trained ML model.** The optional LLM suggests candidate molecules;
  everything after that is database lookup and cheminformatics.
- **EC-family resolution, not variant resolution.** Kinetics and measured ee for
  a specific enzyme variant still need BRENDA or the primary literature.
- **No benchmark of the discovery stage.** Its known bad output is rejected by
  the structure gate, which is a weaker claim than measuring it.
- **The pair table is small** — 57 enantiomer pairs over 50 scaffolds. That
  ceiling is set by Rhea's curation, not by the method, and it grows as Rhea does.

## Layout

```
main.py                      CLI
ChiraLLM/query_handler.py    stage A: LLM suggestions (untrusted)
ChiraLLM/structure_gate.py   does the structure match the name?
ChiraLLM/retrieval.py        KEGG lookup by structure, not by name
ChiraLLM/feasibility.py      the assessment and its evidence record
stereo/                      Rhea x ChEBI stereochemistry index + query tool
utils/kegg.py                cached KEGG client
utils/thermo.py              eQuilibrator wrapper
utils/brenda.py              BRENDA (needs a local licensed dump)
benchmarks/                  dataset, harness, results
tests/                       one regression test per defect found
```

```bash
./venv/bin/python -m pytest tests/ -q                   # 36 tests
./venv/bin/python -m pytest tests/ -q -m "not network"  # offline subset
```

Engineering decisions, failure modes already ruled out and the rationale for
each are logged in [CLAUDE.md](CLAUDE.md).

## Roadmap

Not built yet, and listed here rather than above:

- Multi-step retrosynthesis to a host's native precursors with stereochemistry
  retained at every step.
- Generic Rhea participants (`an L-alpha-amino acid`) are invisible to
  structure-exact matching, and shouldn't be.
- Expression-level host matching: cofactor availability, toxicity, codon usage.
- Enzyme-variant resolution — kinetics and measured ee, not just EC family.
- Target classes: arylpropionates, beta-lactams, chiral amines and alcohols.

## Developed with BettermindLabs

Questions or collaborations: https://www.linkedin.com/in/alexeimanuel/

> *Biology doesn't ask how to make more of X — it asks what X's it could make best.*
