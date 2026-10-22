"""Emit the final stereochemistry-resolved enzyme index.

For every stereochemically-defined compound in Rhea: which EC numbers act on
*that exact stereoisomer*, and which enzyme (if any) acts on its enantiomer.
"""
import csv, json, os, re, collections
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "data")
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

TAG = re.compile(r"<[^>]+>")
def clean(n): return TAG.sub("", n or "").strip()

cm = json.load(open(f"{D}/chiral_map.json"))
skel_of, chiral, nstereo = cm["skel_of"], cm["chiral"], cm["nstereo"]
rxns = json.load(open(f"{D}/rxns.json"))
names = json.load(open(f"{D}/names_used.json"))
pairs = json.load(open(f"{D}/pairs.json"))
struct = {r["chebi_id"]: (r["smiles"], r["inchikey"])
          for r in csv.DictReader(open(f"{D}/chebi_struct.tsv"), delimiter="\t")}

# compound -> reactions acting on it
acts = collections.defaultdict(list)
for rx in rxns:
    for p in rx["parts"]:
        if p in chiral: acts[p].append(rx)

# enantiomer counterpart, from the validated pairs
counterpart = {}
for p in pairs:
    if p["kind"] != "enantiomers": continue
    counterpart.setdefault(p["a_chebi"], set()).update(p["b_ec"])
    counterpart.setdefault(p["b_chebi"], set()).update(p["a_ec"])
enant_of = {}
for p in pairs:
    if p["kind"] == "enantiomers":
        enant_of[p["a_chebi"]] = p["b_chebi"]; enant_of[p["b_chebi"]] = p["a_chebi"]

rows = []
for cid, rs in sorted(acts.items(), key=lambda kv: -len(kv[1])):
    smi, ik = struct.get(cid, ("", ""))
    ecs = sorted({e for rx in rs for e in rx["ec"]})
    if not ecs: continue
    other = enant_of.get(cid, "")
    rows.append(dict(
        chebi=f"CHEBI:{cid}", name=clean(names.get(cid, cid)), smiles=smi, inchikey=ik,
        skeleton=ik.split("-")[0] if ik else "",
        # blocks 1+2: connectivity + stereochemistry, dropping the final
        # protonation block so a neutral drawing matches ChEBI's physiological ion
        stereo_key="-".join(ik.split("-")[:2]) if ik else "", n_stereocentres=nstereo.get(cid, 0),
        n_reactions=len(rs), ec_acting_on_this_isomer=";".join(ecs),
        rhea_ids=";".join(sorted({rx["rhea"] for rx in rs}))[:300],
        enantiomer_chebi=f"CHEBI:{other}" if other else "",
        enantiomer_name=clean(names.get(other, "")) if other else "",
        ec_acting_on_enantiomer=";".join(sorted(counterpart.get(cid, []))) ,
        enantiocomplementary_pair_known=bool(other)))

with open(f"{D}/stereo_enzyme_index.tsv", "w", newline="") as w:
    wr = csv.DictWriter(w, fieldnames=list(rows[0].keys()), delimiter="\t")
    wr.writeheader(); wr.writerows(rows)

# clean pair table
ptab = []
for p in pairs:
    ptab.append(dict(kind=p["kind"], skeleton=p["skeleton"],
        a_chebi=f"CHEBI:{p['a_chebi']}", a_name=clean(p["a_name"]), a_ec=";".join(p["a_ec"]),
        a_rhea=";".join(p["a_rhea"]), a_equation=clean(p["a_eq"]),
        b_chebi=f"CHEBI:{p['b_chebi']}", b_name=clean(p["b_name"]), b_ec=";".join(p["b_ec"]),
        b_rhea=";".join(p["b_rhea"]), b_equation=clean(p["b_eq"]),
        n_stereocentres=p["n_stereocentres"]))
# dedupe identical skeleton+EC-set rows
seen=set(); ded=[]
for r in ptab:
    k=(r["skeleton"],r["a_ec"],r["b_ec"],r["a_chebi"],r["b_chebi"])
    if k in seen: continue
    seen.add(k); ded.append(r)
with open(f"{D}/enantiocomplementary_pairs.tsv","w",newline="") as w:
    wr=csv.DictWriter(w,fieldnames=list(ded[0].keys()),delimiter="\t")
    wr.writeheader(); wr.writerows(ded)

n_en=sum(1 for r in ded if r["kind"]=="enantiomers")
print(f"index rows (stereo-defined compounds with >=1 EC) : {len(rows)}")
print(f"  with a known enantiocomplementary enzyme        : {sum(1 for r in rows if r['enantiocomplementary_pair_known'])}")
print(f"pair table rows after dedupe                      : {len(ded)}")
print(f"  true enantiomer pairs                           : {n_en}")
print(f"  other stereoisomer (diastereomer) pairs         : {len(ded)-n_en}")
print(f"\nwrote data/stereo_enzyme_index.tsv and data/enantiocomplementary_pairs.tsv")
