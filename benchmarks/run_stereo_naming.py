"""Stage 3b: validated enantiomer pairs + name/structure consistency.
   Stage 1b: retrieval robustness to the naming styles an LLM actually emits."""
import sys, os, json, csv, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import run_benchmark as B          # installs shims, imports repo code
F = B.F
import keggc_shim as K
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

def blind(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(m, isomericSmiles=False)))

def tags(smi):
    m = Chem.MolFromSmiles(smi)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return tuple(t for _, t in Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                                         useLegacyImplementation=False))

CAND = [("limonene","cpd:C06099","cpd:C00521"), ("lactate","cpd:C00186","cpd:C00256"),
        ("alanine","cpd:C00041","cpd:C00133"), ("malate","cpd:C00149","cpd:C00497"),
        ("tryptophan","cpd:C00078","cpd:C00525"), ("DOPA","cpd:C00355","cpd:C00822"),
        ("serine","cpd:C00065","cpd:C00740"), ("threonine","cpd:C00188","cpd:C00820"),
        ("glutamate","cpd:C00025","cpd:C00217"), ("valine","cpd:C00183","cpd:C06417")]

out = {"pairs": [], "mismatch": [], "naming": []}
valid = []
for lab, ea, eb in CAND:
    sa, sb = K.smiles_for(ea), K.smiles_for(eb)
    na, nb = K.name_for(ea), K.name_for(eb)
    if not sa or not sb: continue
    same_skel = blind(sa) == blind(sb)
    ta, tb = tags(sa), tags(sb)
    inverted = (len(ta) == len(tb) and len(ta) > 0
                and all({x, y} == {"R", "S"} for x, y in zip(ta, tb)))
    rec = dict(pair=lab, a=na, b=nb, a_id=ea, b_id=eb, same_skeleton=same_skel,
               a_tags="".join(ta), b_tags="".join(tb), true_enantiomers=bool(same_skel and inverted))
    out["pairs"].append(rec)
    if rec["true_enantiomers"]: valid.append((lab, na, sa, nb, sb))
    print(f"{'VALID ' if rec['true_enantiomers'] else 'REJECT'} {lab:11s} {na[:26]:26s}/{nb[:26]:26s} skel={same_skel} {''.join(ta)}|{''.join(tb)}")

print(f"\n{len(valid)} validated enantiomer pairs\n")
HOST = "E. coli"
for lab, na, sa, nb, sb in valid:
    A = F.combined_feasibility_analysis(na, sa, HOST)
    Bv = F.combined_feasibility_analysis(nb, sb, HOST)
    rd_same = A.get("rdkit_analysis") == Bv.get("rdkit_analysis")
    diff = [k for k in ("biosynthetic_feasibility","confidence","kegg_compound_id","bottlenecks")
            if A.get(k) != Bv.get(k)]
    # name/structure mismatch: correct name, WRONG enantiomer's structure
    M = F.combined_feasibility_analysis(na, sb, HOST)
    flagged = (M.get("rdkit_analysis") != A.get("rdkit_analysis")
               or M.get("biosynthetic_feasibility") != A.get("biosynthetic_feasibility")
               or M.get("confidence") != A.get("confidence")
               or M.get("bottlenecks") != A.get("bottlenecks"))
    out["mismatch"].append(dict(pair=lab, name_used=na, structure_used=nb,
                                identical_to_correct=not flagged,
                                mismatch_flagged=bool(flagged),
                                rdkit_identical_between_enantiomers=rd_same,
                                differing_fields=diff))
    print(f"{lab:11s} rdkit identical R vs S: {str(rd_same):5s} | "
          f"name+wrong-structure flagged: {str(flagged):5s} | fields differing: {diff}")

# ---- Stage 1b: naming robustness -------------------------------------------
rows = list(csv.DictReader(open(os.path.join(HERE, "dataset.csv"))))
def variants(name, smi):
    n = name
    v = {"canonical": n}
    base = n.replace("L-", "").replace("(R)-", "").replace("(S)-", "")
    v["plain"] = base
    v["acid_style"] = (base + " acid") if base.endswith(("ate",)) else base
    v["stereo_prefixed"] = "(S)-" + base if not n.startswith("(") else n
    v["title_case"] = " ".join(w.capitalize() for w in base.split())
    return v

print("\n--- naming robustness ---")
styles = ["canonical","plain","acid_style","stereo_prefixed","title_case"]
agg = {s: {"hit":0,"exact":0,"n":0} for s in styles}
for r in rows:
    if not r["truth_kegg_cpd"]: continue
    vs = variants(r["query_name"], r["smiles"])
    rec = {"query_name": r["query_name"], "truth": r["truth_kegg_cpd"], "results": {}}
    for s in styles:
        got = (F.find_compound_kegg(vs[s]) or "").replace("cpd:","").strip()
        rec["results"][s] = {"query": vs[s], "got": got,
                             "hit": bool(got), "exact": got == r["truth_kegg_cpd"]}
        agg[s]["n"] += 1; agg[s]["hit"] += bool(got); agg[s]["exact"] += got == r["truth_kegg_cpd"]
    out["naming"].append(rec)
for s in styles:
    a = agg[s]
    print(f"{s:16s} n={a['n']:2d} coverage={a['hit']/a['n']:.0%} exact={a['exact']/a['n']:.0%}")
out["naming_summary"] = agg
json.dump(out, open(os.path.join(HERE,"stage3b.json"),"w"), indent=1, default=str)
print("\nwrote stage3b.json")
