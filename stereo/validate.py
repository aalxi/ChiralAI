"""Validation: recall against enantiocomplementary enzyme pairs specified
independently of the build, plus negative controls.

The held-out list is written from the biocatalysis literature BEFORE looking at
the output table. Each entry names the two EC numbers that act on opposite
stereoisomers of one scaffold. Recall is then whether the pipeline recovered
that EC pair.
"""
import csv, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "data")

# (label, EC on one form, EC on the other form)
KNOWN = [
 ("lactate dehydrogenase, L vs D",            "1.1.1.27", "1.1.1.28"),
 ("limonene synthase, (4S) vs (4R)",          "4.2.3.16", "4.2.3.20"),
 ("serine ammonia-lyase, L vs D",             "4.3.1.17", "4.3.1.18"),
 ("amino-acid oxidase, L vs D",               "1.4.3.2",  "1.4.3.3"),
 ("glutamate oxidase, L vs D",                "1.4.3.11", "1.4.3.7"),
 ("threonine aldolase, L vs D",               "4.1.2.48", "4.1.2.42"),
 ("copalyl diphosphate synthase, ent vs normal","5.5.1.13","5.5.1.12"),
 ("glycerol-3-P vs glycerol-1-P dehydrogenase","1.1.1.8", "1.1.1.261"),
 ("cysteine desulfhydrase, L vs D",           "4.4.1.1",  "4.4.1.15"),
 ("tagatose/galactitol dehydrogenase, L vs D","1.1.1.406","1.1.1.16"),
 ("3-hydroxyaspartate dehydratase, L vs D",   "4.3.1.16", "4.3.1.27"),
 ("mandelate dehydrogenase, S vs R",          "1.1.99.31","1.1.3.46"),
 ("lysine oxidase, L vs D",                   "1.4.3.14", "1.4.3.3"),
 ("germacrene D synthase, (-) vs (+)",        "4.2.3.75", "4.2.3.77"),
 ("6-hydroxynicotine oxidase, L vs D",        "1.5.3.5",  "1.5.3.6"),
]
# pairs that must NOT appear: same enzyme on both, or unrelated ECs
NEG = [("same EC both sides", "1.1.1.27", "1.1.1.27"),
       ("unrelated ECs",      "1.1.1.27", "4.2.3.16")]

rows = list(csv.DictReader(open(f"{D}/enantiocomplementary_pairs.tsv"), delimiter="\t"))
def ecset(s): return {e.replace("EC:", "").strip() for e in s.split(";") if e.strip()}
found = []
for r in rows:
    found.append((ecset(r["a_ec"]), ecset(r["b_ec"]), r["kind"], r["a_name"], r["b_name"]))

def recovered(e1, e2):
    for A, B, kind, an, bn in found:
        if (e1 in A and e2 in B) or (e1 in B and e2 in A):
            return kind, an, bn
    return None

print(f"pairs in output table: {len(rows)}  "
      f"(enantiomers {sum(1 for r in rows if r['kind']=='enantiomers')})\n")
print("=== RECALL on held-out literature pairs ===\n")
hits = 0
for lab, e1, e2 in KNOWN:
    g = recovered(e1, e2)
    if g:
        hits += 1
        print(f"  FOUND   {lab:44s} EC {e1} / {e2}   [{g[0]}]")
    else:
        print(f"  MISSED  {lab:44s} EC {e1} / {e2}")
print(f"\nrecall: {hits}/{len(KNOWN)} = {hits/len(KNOWN):.0%}\n")
print("=== NEGATIVE CONTROLS (must not appear) ===\n")
for lab, e1, e2 in NEG:
    g = recovered(e1, e2)
    print(f"  {'FAIL - appeared' if g else 'ok - absent  '}  {lab:26s} EC {e1} / {e2}")
print()
# independent consistency check: every 'enantiomers' pair must have mirror SMILES
sys.path.insert(0, HERE)
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
idx = {r["chebi"]: r for r in csv.DictReader(open(f"{D}/stereo_enzyme_index.tsv"), delimiter="\t")}
def mir(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    for a in m.GetAtoms():
        t = a.GetChiralTag()
        if t == Chem.ChiralType.CHI_TETRAHEDRAL_CW: a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW: a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return Chem.MolToSmiles(m)
ok = bad = skip = 0
for r in rows:
    if r["kind"] != "enantiomers": continue
    a, b = idx.get(r["a_chebi"]), idx.get(r["b_chebi"])
    if not a or not b or not a["smiles"] or not b["smiles"]: skip += 1; continue
    ma, mb = mir(a["smiles"]), Chem.MolToSmiles(Chem.MolFromSmiles(b["smiles"]))
    if ma == mb: ok += 1
    else: bad += 1
print("=== STRUCTURAL AUDIT of every pair labelled 'enantiomers' ===\n")
print(f"  mirror image confirmed : {ok}")
print(f"  NOT mirror images      : {bad}")
print(f"  unresolvable structure : {skip}")
json.dump(dict(recall=f"{hits}/{len(KNOWN)}", audit_ok=ok, audit_bad=bad, audit_skip=skip),
          open(f"{D}/validation.json", "w"), indent=1)
