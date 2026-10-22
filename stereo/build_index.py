"""Build a stereochemistry-resolved enzyme index from Rhea + ChEBI.

The question a biocatalysis or metabolic-engineering group actually asks is
"which enzyme makes THIS enantiomer, and is there one that makes the other?"
Rhea curates reactions with stereochemically explicit ChEBI participants and
EC/UniProt links, so the answer is derivable rather than predicted.

Method
  1. Resolve every Rhea participant to a standard InChIKey via ChEBI.
  2. The InChIKey's first block is the connectivity skeleton; two stereoisomers
     share it and differ in the second block. Group participants by skeleton to
     find stereoisomer families.
  3. Give each reaction a signature in which every participant belonging to a
     stereoisomer family is replaced by its skeleton. Reactions sharing a
     signature are the same chemistry differing only in stereochemistry, so the
     grouping is O(n) rather than pairwise.
  4. Within a signature bucket, two reactions with different stereoisomers and
     different EC numbers are a stereocomplementary enzyme pair.
"""
import csv, gzip, json, os, sys, collections

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

# ---- ChEBI structures ------------------------------------------------------
struct = {}
with open(os.path.join(D, "chebi_struct.tsv")) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        struct[r["chebi_id"]] = (r["smiles"], r["inchikey"])
print("chebi structures:", len(struct))

names = {}
csv.field_size_limit(10**9)
with gzip.open(os.path.join(D, "compounds.tsv.gz"), "rt", encoding="utf-8",
               errors="replace", newline="") as fh:
    rd = csv.reader(fh, delimiter="\t", quotechar='"'); hdr = next(rd)
    ix = {k: i for i, k in enumerate(hdr)}
    for r in rd:
        if len(r) < len(hdr): continue
        names[r[ix["id"]]] = r[ix["name"]]
print("chebi names:", len(names))

# ---- Rhea reactions --------------------------------------------------------
rxns = []
with open(os.path.join(D, "rhea.tsv")) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        parts = [p.replace("CHEBI:", "").strip()
                 for p in (r["ChEBI identifier"] or "").split(";") if p.strip()]
        if not parts: continue
        rxns.append(dict(rhea=r["Reaction identifier"], eq=r["Equation"],
                         ec=[e.strip() for e in (r["EC number"] or "").split(";") if e.strip()],
                         parts=parts, enzymes=(r["Enzymes"] or "").strip()))
print("rhea reactions with participants:", len(rxns))

# ---- which participants are chiral, and what are their stereo families? ----
used = sorted({p for rx in rxns for p in rx["parts"]})
print("distinct participants:", len(used))

chiral, skel_of, nstereo = {}, {}, {}
for cid in used:
    s = struct.get(cid)
    if not s: continue
    smi, ik = s
    if not ik or "-" not in ik: continue
    skel_of[cid] = ik.split("-")[0]
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None: continue
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    cen = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
    d = [c for c in cen if c[1] in ("R", "S")]
    nstereo[cid] = len(d)
    if d: chiral[cid] = ik
print("participants resolved to a structure:", len(skel_of))
print("participants with >=1 DEFINED stereocentre:", len(chiral))

# a stereoisomer family = one skeleton carrying >1 distinct full InChIKey
by_skel = collections.defaultdict(set)
for cid, ik in chiral.items():
    by_skel[skel_of[cid]].add(ik)
families = {s for s, iks in by_skel.items() if len(iks) > 1}
fam_members = collections.defaultdict(set)
for cid in chiral:
    if skel_of[cid] in families: fam_members[skel_of[cid]].add(cid)
print("stereoisomer families among Rhea participants:", len(families))

# ---- signature grouping ----------------------------------------------------
def signature(rx):
    sig = []
    for p in rx["parts"]:
        sk = skel_of.get(p)
        sig.append(("SKEL", sk) if (sk in families and p in chiral) else ("CPD", p))
    return tuple(sorted(sig))

buckets = collections.defaultdict(list)
for rx in rxns:
    if any(skel_of.get(p) in families and p in chiral for p in rx["parts"]):
        buckets[signature(rx)].append(rx)
multi = {k: v for k, v in buckets.items() if len(v) > 1}
print("signature buckets containing a chiral participant:", len(buckets))
print("buckets holding more than one reaction:", len(multi))

json.dump(dict(
    n_reactions=len(rxns), n_participants=len(used),
    n_resolved=len(skel_of), n_chiral=len(chiral), n_families=len(families),
    n_buckets=len(buckets), n_multi=len(multi)),
    open(os.path.join(D, "build_stats.json"), "w"), indent=1)

# persist what later steps need
json.dump(dict(skel_of=skel_of, chiral=chiral, nstereo=nstereo,
               families=sorted(families),
               fam_members={k: sorted(v) for k, v in fam_members.items()}),
          open(os.path.join(D, "chiral_map.json"), "w"))
with open(os.path.join(D, "buckets.json"), "w") as w:
    json.dump([[list(k), [r["rhea"] for r in v]] for k, v in multi.items()], w)
json.dump(rxns, open(os.path.join(D, "rxns.json"), "w"))
json.dump({k: v for k, v in names.items() if k in used},
          open(os.path.join(D, "names_used.json"), "w"))
print("\nwrote data/chiral_map.json, buckets.json, rxns.json, names_used.json")
