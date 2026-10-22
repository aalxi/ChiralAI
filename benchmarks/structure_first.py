"""Measured-improvement arm for Claim 1.

Name-first retrieval (what the repo does) matches a user's string against KEGG's
chosen synonyms. Structure-first retrieval matches the molecule itself:
exact formula search to narrow, then full InChIKey to pick, which is
stereochemistry-sensitive by construction.

Both are scored on the same 32 compounds against the same ground truth.
"""
import sys, os, json, csv
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import keggc_shim as K
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors as rmd, inchi
RDLogger.DisableLog("rdApp.*")

def ikey(smi):
    m = Chem.MolFromSmiles(smi)
    return inchi.MolToInchiKey(m) if m else None

def formula(smi):
    m = Chem.MolFromSmiles(smi)
    return rmd.CalcMolFormula(m).replace("+","").replace("-","") if m else None

def kegg_formula_candidates(f):
    txt = K._get(f"find/compound/{f}/formula")
    out = []
    for line in (txt.strip().split("\n") if txt.strip() else []):
        if "\t" not in line: continue
        cid, fo = line.split("\t")[0], line.split("\t")[1].strip()
        if fo == f:                     # exact formula only
            out.append(cid.replace("cpd:",""))
    return out

def structure_first(smi, cap=25):
    f = formula(smi)
    if not f: return None, 0, "no formula"
    cands = kegg_formula_candidates(f)
    want = ikey(smi)
    want_skel = want.split("-")[0] if want else None
    skel_hit = None
    for cid in cands[:cap]:
        s = K.smiles_for(f"cpd:{cid}")
        if not s: continue
        got = ikey(s)
        if not got: continue
        if got == want:
            return cid, len(cands), "exact inchikey"
        if want_skel and got.split("-")[0] == want_skel and skel_hit is None:
            skel_hit = cid
    if skel_hit: return skel_hit, len(cands), "skeleton only (stereo undetermined)"
    return None, len(cands), "no structural match"

rows = [r for r in csv.DictReader(open(os.path.join(HERE,"dataset.csv"))) if r["truth_kegg_cpd"]]
s3b = json.load(open(os.path.join(HERE,"stage3b.json")))
namearm = {r["query_name"]: r["results"] for r in s3b["naming"]}

res = []
for r in rows:
    cid, n, how = structure_first(r["smiles"])
    nf = namearm.get(r["query_name"], {})
    res.append(dict(query_name=r["query_name"], truth=r["truth_kegg_cpd"],
                    struct_got=cid or "", struct_exact=(cid == r["truth_kegg_cpd"]),
                    struct_how=how, n_formula_candidates=n,
                    name_canonical_exact=nf.get("canonical",{}).get("exact"),
                    name_plain_exact=nf.get("plain",{}).get("exact"),
                    name_stereo_hit=nf.get("stereo_prefixed",{}).get("hit")))
    print(f"{r['query_name'][:20]:20s} truth={r['truth_kegg_cpd']} struct={cid or '---':8s} "
          f"{'OK ' if cid==r['truth_kegg_cpd'] else 'NO '} cands={n:3d} {how}", flush=True)

n = len(res)
print()
print("| retrieval method | exact top-1 | coverage |")
print("|---|---|---|")
for nm, key, hk in [("name-first, canonical name (as shipped)","name_canonical_exact",None),
                    ("name-first, loose name","name_plain_exact",None),
                    ("structure-first (formula + InChIKey)","struct_exact",None)]:
    ok = sum(1 for r in res if r[key]); print(f"| {nm} | {ok}/{n} = {ok/n:.0%} | |")
sh = sum(1 for r in res if r["name_stereo_hit"])
print(f"| name-first, stereo-descriptor name | - | {sh}/{n} = {sh/n:.0%} |")
print()
print("structure-first resolutions that were stereo-exact:",
      sum(1 for r in res if r["struct_how"]=="exact inchikey"), "/", n)
print("structure-first that could only match the flat skeleton:",
      sum(1 for r in res if r["struct_how"].startswith("skeleton")), "/", n)
json.dump(res, open(os.path.join(HERE,"structure_first.json"),"w"), indent=1)
print("\nwrote structure_first.json")
