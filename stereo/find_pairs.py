"""Extract stereocomplementary enzyme pairs and classify enantiomer vs diastereomer."""
import csv, json, os, collections
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "data")
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

cm = json.load(open(os.path.join(D, "chiral_map.json")))
skel_of, chiral, nstereo = cm["skel_of"], cm["chiral"], cm["nstereo"]
families = set(cm["families"])
rxns = {r["rhea"]: r for r in json.load(open(os.path.join(D, "rxns.json")))}
names = json.load(open(os.path.join(D, "names_used.json")))
buckets = json.load(open(os.path.join(D, "buckets.json")))
struct = {}
with open(os.path.join(D, "chebi_struct.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        struct[r["chebi_id"]] = r["smiles"]

def mirror(smi):
    """Invert every defined stereocentre: the enantiomer."""
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    for a in m.GetAtoms():
        t = a.GetChiralTag()
        if t == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    try: return Chem.MolToSmiles(m)
    except Exception: return None

def canon(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m else None

pairs = []
for sig, rids in buckets:
    rs = [rxns[r] for r in rids if r in rxns]
    # the chiral family member each reaction uses
    def fam_member(rx):
        got = [p for p in rx["parts"] if p in chiral and skel_of.get(p) in families]
        return tuple(sorted(got))
    by_member = collections.defaultdict(list)
    for rx in rs: by_member[fam_member(rx)].append(rx)
    if len(by_member) < 2: continue            # all reactions use the same isomer
    keys = sorted(by_member)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            A, Bk = keys[i], keys[j]
            # Require a single swapped participant. Cofactors such as NAD(+) are
            # themselves chiral, so compare the family-member sets rather than
            # demanding only one chiral participant in the whole reaction.
            onlyA = sorted(set(A) - set(Bk)); onlyB = sorted(set(Bk) - set(A))
            if len(onlyA) != 1 or len(onlyB) != 1: continue
            a, b = onlyA[0], onlyB[0]
            if skel_of.get(a) != skel_of.get(b): continue
            ecA = sorted({e for rx in by_member[A] for e in rx["ec"]})
            ecB = sorted({e for rx in by_member[Bk] for e in rx["ec"]})
            if not ecA or not ecB: continue
            if set(ecA) == set(ecB): continue          # same enzyme, not complementary
            sa, sb = struct.get(a), struct.get(b)
            kind = "stereoisomers"
            if sa and sb:
                ma, mb = canon(sa), canon(sb)
                mir = mirror(sa)
                if mir and mb and mir == mb and ma != mb:
                    kind = "enantiomers"
                elif ma == mb:
                    kind = "identical"
            if kind == "identical": continue
            pairs.append(dict(
                skeleton=skel_of.get(a), kind=kind,
                a_chebi=a, a_name=names.get(a, a), a_smiles=sa,
                a_ec=ecA, a_rhea=[rx["rhea"] for rx in by_member[A]],
                a_eq=by_member[A][0]["eq"],
                a_enzymes=by_member[A][0]["enzymes"],
                b_chebi=b, b_name=names.get(b, b), b_smiles=sb,
                b_ec=ecB, b_rhea=[rx["rhea"] for rx in by_member[Bk]],
                b_eq=by_member[Bk][0]["eq"],
                b_enzymes=by_member[Bk][0]["enzymes"],
                n_stereocentres=nstereo.get(a)))

enant = [p for p in pairs if p["kind"] == "enantiomers"]
print(f"stereocomplementary enzyme pairs : {len(pairs)}")
print(f"  of which true enantiomer pairs : {len(enant)}")
print(f"  distinct skeletons covered     : {len({p['skeleton'] for p in pairs})}")
print()
print("=== VALIDATION: lactate dehydrogenase must appear as EC 1.1.1.27 vs 1.1.1.28 ===")
hit = [p for p in pairs if "lactate" in (p["a_name"] + p["b_name"]).lower()
       and any("1.1.1.27" in e for e in p["a_ec"] + p["b_ec"])]
for p in hit[:4]:
    print(f"  {p['a_name']} {p['a_ec']}  <->  {p['b_name']} {p['b_ec']}   [{p['kind']}]")
print()
print("=== sample of enantiomer pairs ===")
for p in enant[:14]:
    print(f"  {p['a_name'][:34]:34s} {','.join(p['a_ec'])[:22]:22s} <-> "
          f"{p['b_name'][:34]:34s} {','.join(p['b_ec'])[:22]:22s}")
json.dump(pairs, open(os.path.join(D, "pairs.json"), "w"), indent=1)
print(f"\nwrote data/pairs.json ({len(pairs)} pairs)")
