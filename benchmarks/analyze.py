"""Compute metrics and baselines from benchmark results. No metric without a baseline."""
import json, csv, sys, os, math, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_benchmark import auroc, boot_ci, confusion

HERE = os.path.dirname(os.path.abspath(__file__))
ds = {r["query_name"]: r for r in csv.DictReader(open(os.path.join(HERE, "dataset.csv")))}

def load(*names):
    out = {}
    for n in names:
        p = os.path.join(HERE, n)
        if os.path.exists(p): out.update(json.load(open(p)))
    return out

R = load("stage1.json", "stage2.json", "stage3.json", "results.json")
L = []
def say(s=""): L.append(s); print(s)

# ------------------------------------------------- Claim 1
if "stage1_retrieval" in R:
    s1 = R["stage1_retrieval"]
    scored = [r for r in s1 if r["truth"]]
    say("## Claim 1 - identifier grounding")
    say()
    say(f"Scored on the {len(scored)} of {len(s1)} compounds that have a KEGG COMPOUND entry.")
    say()
    hit = sum(r["hit"] for r in scored); ex = sum(r["exact"] for r in scored)
    say(f"| metric | value |")
    say(f"|---|---|")
    say(f"| coverage (any hit) | {hit}/{len(scored)} = {hit/len(scored):.0%} |")
    say(f"| top-1 exact | {ex}/{len(scored)} = {ex/len(scored):.0%} |")
    sw = [r for r in scored if r["same_skeleton_wrong_stereo"]]
    say(f"| right skeleton, wrong stereoisomer | {len(sw)} |")
    mw_ = [r for r in scored if r["multiword"]]; sw_ = [r for r in scored if not r["multiword"]]
    if mw_:
        say(f"| coverage, multi-word names (n={len(mw_)}) | {sum(r['hit'] for r in mw_)}/{len(mw_)} |")
    say(f"| coverage, single-word names (n={len(sw_)}) | {sum(r['hit'] for r in sw_)}/{len(sw_)} |")
    say()
    wrong = [r for r in scored if r["hit"] and not r["exact"]]
    if wrong:
        say("Wrong retrievals:")
        say()
        for r in wrong:
            tag = " (same skeleton, different stereoisomer)" if r["same_skeleton_wrong_stereo"] else ""
            say(f"- {r['query_name']}: returned {r['retrieved']}, correct is {r['truth']}{tag}")
        say()
    miss = [r for r in scored if not r["hit"]]
    if miss:
        say("Returned nothing: " + ", ".join(r["query_name"] for r in miss)); say()

# ------------------------------------------------- Claim 2
if "stage2_feasibility" in R:
    s2 = R["stage2_feasibility"]
    say("## Claim 2 - feasibility discrimination")
    say()
    def block(rows, title):
        lab = [r["label"] for r in rows]
        if not (1 in lab and 0 in lab):
            say(f"### {title}: not scoreable"); return
        conf = [r["confidence"] for r in rows]
        pred = [r["verdict"] for r in rows]
        a = auroc(conf, lab); lo, hi = boot_ci(conf, lab)
        c = confusion(pred, lab)
        # baselines
        b_mw   = auroc([-r["mw"] for r in rows], lab)                      # smaller = feasible
        b_kegg = auroc([1.0 if r["kegg_id"] else 0.0 for r in rows], lab)   # in KEGG at all
        b_path = auroc([float(r["n_pathways"]) for r in rows], lab)         # pathway count
        b_rule = auroc([1.0 if r["mw"] <= 500 else 0.0 for r in rows], lab) # the code's own rule
        say(f"### {title}  (n={len(rows)}, {sum(lab)} positive / {len(lab)-sum(lab)} negative)")
        say()
        say("| model | AUROC |")
        say("|---|---|")
        say(f"| **ChiralAI confidence** | **{a:.3f}**  (95% CI {lo:.2f}-{hi:.2f}) |")
        say(f"| baseline: present in KEGG | {b_kegg:.3f} |")
        say(f"| baseline: KEGG pathway count | {b_path:.3f} |")
        say(f"| baseline: molecular weight alone | {b_mw:.3f} |")
        say(f"| baseline: the code's own MW<=500 rule | {b_rule:.3f} |")
        say(f"| baseline: constant | 0.500 |")
        say()
        say("At the system's own operating point (`biosynthetic_feasibility`):")
        say()
        say("| | value |")
        say("|---|---|")
        say(f"| sensitivity | {c['sensitivity']:.0%} ({c['tp']}/{c['tp']+c['fn']}) |")
        say(f"| specificity | {c['specificity']:.0%} ({c['tn']}/{c['tn']+c['fp']}) |")
        say(f"| precision | {c['precision']:.0%} |")
        say(f"| balanced accuracy | {c['balanced_accuracy']:.0%} |")
        say(f"| false 'feasible' calls | {c['fp']} |")
        say()
        return a
    block(s2, "Full set")
    matched = [r for r in s2 if r["stratum"] in ("A1", "A2", "B")]
    block(matched, "Size-matched subset (drops stratum C, the only group above 500 Da)")

    # confidence distribution
    vals = sorted({round(r["confidence"], 4) for r in s2})
    say(f"Distinct `confidence` values across all {len(s2)} compounds: {len(vals)} -> {vals}")
    say()
    pv = [r['confidence'] for r in s2 if r['label']==1]
    nv = [r['confidence'] for r in s2 if r['label']==0]
    say(f"Mean confidence, positives {statistics.mean(pv):.3f} / negatives {statistics.mean(nv):.3f}")
    say()
    # instrumentation
    nb = sum(1 for r in s2 if r["brenda_found"])
    say(f"BRENDA returned a hit for {nb}/{len(s2)} compounds.")
    dgin, dgval = [], []
    for r in s2:
        dgin += r["dg_inputs"]; dgval += [v for v in r["dg_values"]]
    nonnull = [v for v in dgval if v is not None]
    say(f"Thermodynamic function called {len(dgin)} times; returned a number {len(nonnull)} times.")
    if dgin: say(f"Example arguments it received: {dgin[:4]}")
    say()
    say("False 'feasible' calls, by stratum:")
    say()
    for s in ["B", "C"]:
        sub = [r for r in s2 if r["stratum"] == s]
        if not sub: continue
        fp = [r["query_name"] for r in sub if r["verdict"]]
        say(f"- stratum {s}: {len(fp)}/{len(sub)} called feasible" + (f" -> {', '.join(fp)}" if fp else ""))
    say()
    for s in ["A1", "A2"]:
        sub = [r for r in s2 if r["stratum"] == s]
        if not sub: continue
        tp = sum(1 for r in sub if r["verdict"])
        say(f"- stratum {s} (true positives): {tp}/{len(sub)} called feasible")
    say()

# ------------------------------------------------- Claim 3
if "stage3_enantiomers" in R:
    s3 = R["stage3_enantiomers"]
    say("## Claim 3 - stereochemistry")
    say()
    say(f"Enantiomer pairs tested: {len(s3)}. Each pair is two different KEGG compounds "
        f"with mirror-image structures.")
    say()
    say("| pair | inputs differ | identical verdict | fields that differed |")
    say("|---|---|---|---|")
    for r in s3:
        d = ", ".join(r["differing_fields"]) or "none"
        say(f"| {r['pair']} | {'yes' if r['smiles_differ'] else 'no'} | "
            f"{'yes' if not r['differing_fields'] else 'no'} | {d} |")
    say()
    nd = sum(1 for r in s3 if not r["differing_fields"])
    say(f"Pairs given an identical assessment despite being different molecules: {nd}/{len(s3)}.")
    ri = sum(1 for r in s3 if r["rdkit_identical"])
    say(f"Pairs whose structural analysis block was byte-identical: {ri}/{len(s3)}.")
    say()

open(os.path.join(HERE, "METRICS.md"), "w").write("\n".join(L) + "\n")
print("\n[wrote METRICS.md]")
