"""Classify every KEGG pathway link, then test whether a biosynthesis-only
filter discriminates better than the repo's 'any pathway' signal.

This is the measured-improvement arm: it asks not just 'is the current signal
weak' but 'is there a better signal available from the same free database'.
"""
import sys, os, json, csv, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import keggc_shim as K
from run_benchmark import auroc, boot_ci, confusion

GLOBAL_MAPS = {"map01100","map01110","map01120","map01200","map01210","map01212",
               "map01230","map01232","map01250","map01240","map01060","map01061",
               "map01062","map01063","map01064","map01065","map01066","map01070",
               "map01002","map01003"}

def classify(mapid, name):
    m = mapid.replace("path:", ""); n = name.lower()
    if m in GLOBAL_MAPS:                       return "global_overview"
    if "drug metabolism" in n or "xenobiotic" in n: return "xenobiotic_degradation"
    if "degradation" in n or "catabolism" in n: return "degradation"
    if re.match(r"map0[45]\d{4}", m) or re.match(r"map0[45]\d{3}", m): return "physiology_or_disease"
    if m.startswith("map07"):                  return "drug_classification"
    if "biosynthesis" in n or "biosynthe" in n: return "biosynthesis"
    if "metabolism" in n:                      return "metabolism"
    return "other"

rows = list(csv.DictReader(open(os.path.join(HERE, "dataset.csv"))))
s2 = {r["query_name"]: r for r in json.load(open(os.path.join(HERE,"stage2.json")))["stage2_feasibility"]}

recs = []
for r in rows:
    cid = r["truth_kegg_cpd"]
    links = []
    if cid:
        txt = K._get(f"link/pathway/cpd:{cid}")
        links = [l.split("\t")[1] for l in txt.strip().split("\n") if "\t" in l] if txt.strip() else []
    cls = {}
    for m in links:
        info = K._get(f"get/{m.replace('path:','')}")
        nm = next((l.replace("NAME","").strip() for l in info.split("\n") if l.startswith("NAME")), "?")
        c = classify(m, nm)
        cls.setdefault(c, []).append(f"{m}:{nm[:40]}")
    recs.append(dict(query_name=r["query_name"], label=int(r["label"]), stratum=r["stratum"],
                     mw=float(r["mw"]), n_links=len(links), classes={k: len(v) for k,v in cls.items()},
                     detail=cls,
                     n_biosynthesis=len(cls.get("biosynthesis", [])),
                     n_specific_metab=len(cls.get("biosynthesis", []))+len(cls.get("metabolism", [])),
                     n_noise=sum(len(cls.get(k,[])) for k in
                                 ("xenobiotic_degradation","physiology_or_disease",
                                  "drug_classification","global_overview","degradation")),
                     ))
    print(f"{r['query_name'][:20]:20s} lab={r['label']} links={len(links):2d} "
          f"{ {k:len(v) for k,v in cls.items()} }", flush=True)

lab = [x["label"] for x in recs]
def ev(score, nm, pred=None):
    a = auroc(score, lab); lo,hi = boot_ci(score, lab)
    line = f"| {nm} | {a:.3f} | {lo:.2f}-{hi:.2f} |"
    if pred is not None:
        c = confusion(pred, lab)
        line += f" sens {c['sensitivity']:.0%} spec {c['specificity']:.0%} bal {c['balanced_accuracy']:.0%} fp {c['fp']}"
    return line, a

print("\n| signal | AUROC | 95% CI |")
print("|---|---|---|")
cur_pred = [(x["mw"]<=500) and (x["n_links"]>0) for x in recs]
out = []
for nm, sc, pr in [
    ("as shipped: any KEGG pathway AND MW<=500", [1.0 if p else 0.0 for p in cur_pred], cur_pred),
    ("any KEGG pathway (count)", [float(x["n_links"]) for x in recs], None),
    ("biosynthesis-map count only", [float(x["n_biosynthesis"]) for x in recs], None),
    ("has >=1 biosynthesis map", [1.0 if x["n_biosynthesis"]>0 else 0.0 for x in recs],
                                 [x["n_biosynthesis"]>0 for x in recs]),
    ("has >=1 biosynthesis map AND MW<=500",
        [1.0 if (x["n_biosynthesis"]>0 and x["mw"]<=500) else 0.0 for x in recs],
        [x["n_biosynthesis"]>0 and x["mw"]<=500 for x in recs]),
    ("biosynthesis or specific metabolism map", [float(x["n_specific_metab"]) for x in recs], None),
]:
    line, a = ev(sc, nm, pr); print(line); out.append((nm, a))

json.dump(dict(records=recs, auroc=out), open(os.path.join(HERE,"pathway_classes.json"),"w"),
          indent=1, default=str)
print("\nwrote pathway_classes.json")
