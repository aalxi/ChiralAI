"""Re-measure the FIXED pipeline on the same 36 compounds as the original audit.

A number with no baseline is not evidence, so every figure here is paired with
the shipped pipeline's figure on the identical dataset and identical ground
truth (benchmarks/dataset.csv, labels and provenance in DESIGN.md).

The old verdict is read from stage2.json, recorded during the original run.
Nothing is re-derived for the baseline, so it cannot drift to flatter the new
code.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from ChiraLLM.feasibility import assess


def metrics(pred, lab):
    tp = sum(1 for p, l in zip(pred, lab) if p and l)
    tn = sum(1 for p, l in zip(pred, lab) if not p and not l)
    fp = sum(1 for p, l in zip(pred, lab) if p and not l)
    fn = sum(1 for p, l in zip(pred, lab) if not p and l)
    sens = tp / (tp + fn) if tp + fn else 0.0
    spec = tn / (tn + fp) if tn + fp else 0.0
    return dict(tp=tp, tn=tn, fp=fp, fn=fn, sensitivity=sens, specificity=spec,
                balanced_accuracy=(sens + spec) / 2,
                accuracy=(tp + tn) / len(lab))


rows = list(csv.DictReader(open(os.path.join(HERE, "dataset.csv"))))
old = {r["query_name"]: r for r in json.load(open(os.path.join(HERE, "stage2.json")))["stage2_feasibility"]}

recs = []
for i, r in enumerate(rows, 1):
    rec = assess(r["query_name"], r["smiles"], "E. coli", cross_check=False)
    o = old.get(r["query_name"], {})
    recs.append(dict(
        query_name=r["query_name"], label=int(r["label"]), stratum=r["stratum"],
        mw=float(r["mw"]), truth=r["truth_kegg_cpd"],
        n_stereocentres=int(r["n_stereocentres"] or 0),
        old_verdict=bool(o.get("verdict")), old_confidence=o.get("confidence"),
        old_kegg=o.get("kegg_id") or "",
        new_verdict=bool(rec["biosynthetic_feasibility"]), new_label=rec["verdict"],
        new_kegg=rec.get("kegg_compound_id") or "",
        retrieved_by=rec["retrieval"].get("method"),
        stereo_exact=bool(rec["retrieval"].get("stereo_exact")),
        stereo=rec["stereochemistry"]["outcome"],
        stereo_risk=bool((rec["stereochemistry"] or {}).get("selectivity_risk")),
        n_reasons=len(rec["reasons"]),
    ))
    print(f"{i:2d}/{len(rows)} {r['query_name'][:24]:24s} lab={r['label']} "
          f"old={int(bool(o.get('verdict')))} new={int(recs[-1]['new_verdict'])} "
          f"kegg={recs[-1]['new_kegg'] or '--':7s} stereo={recs[-1]['stereo']}", flush=True)

lab = [x["label"] for x in recs]
m_old = metrics([x["old_verdict"] for x in recs], lab)
m_new = metrics([x["new_verdict"] for x in recs], lab)

out = {"n": len(recs), "old": m_old, "new": m_new, "records": recs}

print("\n| metric | shipped pipeline | fixed pipeline |")
print("|---|---|---|")
for k in ("sensitivity", "specificity", "balanced_accuracy", "accuracy"):
    print(f"| {k.replace('_',' ')} | {m_old[k]:.0%} | {m_new[k]:.0%} |")
for k in ("tp", "fp", "tn", "fn"):
    print(f"| {k} | {m_old[k]} | {m_new[k]} |")

truth_known = [x for x in recs if x["truth"]]
oldok = sum(1 for x in truth_known if x["old_kegg"].replace("cpd:", "") == x["truth"])
newok = sum(1 for x in truth_known if x["new_kegg"] == x["truth"])
print(f"\nretrieval top-1 correct: shipped {oldok}/{len(truth_known)} "
      f"({oldok/len(truth_known):.0%}), fixed {newok}/{len(truth_known)} "
      f"({newok/len(truth_known):.0%})")
print("retrieved by structure:", sum(1 for x in recs if x["retrieved_by"] == "structure"), "/", len(recs))
print("stereo-exact retrieval :", sum(1 for x in recs if x["stereo_exact"]), "/", len(recs))

chiral = [x for x in recs if x["n_stereocentres"] > 0]
print(f"\nstereochemistry, on the {len(chiral)} chiral compounds:")
for o in ("EXACT", "WRONG-FORM", "SKELETON", "NONE"):
    n = sum(1 for x in chiral if x["stereo"] == o)
    print(f"  {o:11s} {n:2d}")
print("  with a named selectivity risk:", sum(1 for x in chiral if x["stereo_risk"]))
print("\nold distinct confidence values:",
      sorted({x["old_confidence"] for x in recs if x["old_confidence"] is not None}))
print("new verdict classes:", sorted({x["new_label"] for x in recs}))

json.dump(out, open(os.path.join(HERE, "fixed.json"), "w"), indent=1)
print("\nwrote benchmarks/fixed.json")
