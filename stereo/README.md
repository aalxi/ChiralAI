# Stereochemistry-resolved enzyme index

Answers one question, from curated data only:

> I want **this enantiomer**. Which enzymes make it, and which enzyme makes the
> one I do not want?

## Why this is the useful question

A metabolic engineer choosing a route to a chiral target does not need a
feasibility score. They need the EC family and the sequences, and they need to
know the selectivity risk, because an enzyme that makes the mirror image of
your target is worse than no enzyme at all. Enantiopurity is usually the whole
commercial point: the wrong enantiomer of a drug is an impurity to be purified
away at best.

Existing resources each stop short of this. KEGG's compound pages do not
organise enzymes by product stereochemistry, and its name search cannot even
retrieve a specified enantiomer reliably (measured at 22% coverage in
[../benchmarks/RESULTS.md](../benchmarks/RESULTS.md)). BRENDA holds
stereoselectivity but largely as free text behind a licence. Retrobiosynthesis
tools such as RetroPath and RetroRules expand reaction rules that are often
stereochemistry-stripped. Rhea curates reactions with stereochemically explicit
ChEBI participants and EC/UniProt links, but is organised by reaction, not by
"which form can I make".

This index is the join nobody had published: Rhea reactions, resolved to exact
stereoisomers through ChEBI structures, inverted into a per-enantiomer enzyme
lookup, with the enantiocomplementary partner attached.

## What is in it

| artifact | contents |
|---|---|
| `data/stereo_enzyme_index.tsv` | 3,357 stereochemically-defined compounds with at least one EC acting on them |
| `data/enantiocomplementary_pairs.tsv` | 174 stereocomplementary enzyme pairs, of which **57 are true enantiomer pairs** spanning 50 distinct scaffolds |

Built from 18,611 Rhea reactions and 190,078 ChEBI structures. Of 14,261
distinct Rhea participants, 11,687 resolve to a structure and 7,332 carry at
least one defined stereocentre, falling into 724 stereoisomer families.

## Use

```bash
python query.py --smiles 'C[C@H](O)C(=O)O'      # (S)-lactic acid
python query.py --name limonene
python query.py --smiles '...' --sequences      # add UniProt accessions
```

Four outcomes, and the distinction between them is the product:

| outcome | meaning |
|---|---|
| `EXACT` | an enzyme acting on your exact stereoisomer is documented |
| `WRONG-FORM` | the scaffold is known but **only the other stereoisomer**, so the enzymes listed make the form you do not want |
| `SKELETON` | scaffold known, your query left stereochemistry unspecified |
| `NONE` | no curated reaction matches this connectivity |

`WRONG-FORM` is the outcome that matters, because it is the one every
name-matching tool silently reports as success.

Worked example, the case the original pipeline got wrong:

```
target     : (R)-limonene          InChIKey XMGQYMWWDOXHJM-JTQLQIEISA-N
RESULT     : EXACT  — (4R)-limonene (CHEBI:15382)
made by    : EC 4.2.3.20, plus 1.14.13.107 / 1.14.14.53 / 1.17.99.8
RISK       : (4S)-limonene is made by EC 4.2.3.16
```

The repository's earlier KEGG name lookup returned *(S)*-limonene for the query
"limonene", which is the wrong enantiomer for the industrial target, and nothing
flagged it.

## Validation

`python validate.py`

Recall is measured against 15 enantiocomplementary enzyme pairs written from the
biocatalysis literature before inspecting the output.

| check | result |
|---|---|
| recall on held-out literature pairs | **13 / 15 = 87%** |
| negative controls (same EC both sides; unrelated ECs) | both correctly absent |
| structural audit of every pair labelled `enantiomers` | **57 / 57** confirmed mirror images, 0 false |

The structural audit is the important one: each pair is independently re-derived
by inverting every stereocentre of one member with RDKit and requiring the
canonical SMILES to equal the other member. No pair is labelled an enantiomer
pair on the strength of its name.

Both misses are accounted for, and only one is a real limitation:

- **Amino-acid oxidase, EC 1.4.3.2 vs 1.4.3.3.** Rhea states these over the
  generic classes "an L-alpha-amino acid" and "a D-alpha-amino acid". Generic
  classes have no concrete structure and therefore no InChIKey, so
  structure-exact matching cannot see them. A genuine coverage gap.
- **Mandelate dehydrogenase, EC 1.1.99.31 vs 1.1.3.46.** My held-out entry was
  wrong. EC 1.1.99.31 acts on *(S)*-mandelate and EC 1.1.3.46 on
  *(S)*-4-hydroxymandelate, a different scaffold rather than the enantiomer.
  Excluding it was correct, so true recall on correctly specified pairs is
  13 / 14 = 93%.

## Two design details that decide whether this works

**Match on InChIKey blocks 1 and 2, never the whole key.** The third block
encodes protonation. ChEBI stores physiological ions and chemists draw neutral
molecules, so *(S)*-lactic acid keys as `...-N` while ChEBI's *(S)*-lactate keys
as `...-M`. Matching the full key reports a false `WRONG-FORM` on a correct
query, which is worse than no answer. Blocks 1 and 2 are
stereochemistry-exact and ionisation-insensitive.

**Chiral cofactors break naive pair-finding.** The first version required
exactly one chiral participant per reaction, which silently excluded every
dehydrogenase, because NAD and NADH carry their own stereocentres. Lactate
dehydrogenase, the most obvious pair in all of biochemistry, was missing.
Comparing *which single participant was swapped* between two reactions raised
the yield from 108 pairs to 175. If the canonical validation case had not been
checked explicitly, that bug would have shipped looking like a complete result.

## What this does not claim

- **Curated coverage, not prediction.** `NONE` means no curated evidence, not
  that no enzyme exists. Rhea is incomplete and biased toward well-studied
  metabolism. Nothing here extrapolates to novel substrates.
- **EC-level, not variant-level.** An EC number names a family. Catalytic
  efficiency, expression, and enantiomeric excess for a given variant on a given
  substrate are not in Rhea. Use this to pick the family, then BRENDA or the
  primary literature for kinetics.
- **Enzyme-level, not pathway-level.** This is one step, not a route. Pathway
  assembly, cofactor balance and thermodynamics are separate problems.
- **Generic reaction entries are invisible**, as the amino-acid oxidase miss shows.
- **Diastereomer pairs are reported separately** from enantiomer pairs and are
  not interchangeable. 117 of the 174 pairs are diastereomeric.
- **The table is small.** 57 enantiomer pairs over 50 scaffolds, and 106 of the
  3,357 indexed compounds have a recorded counterpart. That is a ceiling set by
  Rhea's curation, not by the method. It will grow as Rhea does.

## Reproduce

```bash
python -m venv venv && ./venv/bin/pip install rdkit requests
./venv/bin/python build_index.py    # Rhea + ChEBI -> stereoisomer families
./venv/bin/python find_pairs.py     # stereocomplementary enzyme pairs
./venv/bin/python enzyme_index.py   # final TSVs
./venv/bin/python validate.py       # recall + structural audit
```

Source data is downloaded in `data/`. Rhea and ChEBI are both free and
redistributable under CC-BY; the enzyme assignments are EBI's curation, not ours.
