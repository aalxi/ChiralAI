"""Retrospective benchmark for the ChiralAI feasibility pipeline.

Runs the repository's OWN functions. Two shims are required and both are
faithful to what a user of this repo actually experiences:

  * equilibrator_api is not installed (absent from the environment), so it is
    replaced by a stub whose parse calls raise. utils.thermo.delta_g wraps
    everything in `except Exception: return None`, so the stub reproduces the
    real return value. Separately, every argument delta_g receives is recorded,
    which is the measurement that matters.
  * requests.get is wrapped in a transparent on-disk cache that preserves URL
    and response bodies, so reruns do not hammer KEGG. Every URL the repo
    constructs is logged verbatim for audit.

utils.brenda.query_brenda is NOT stubbed. It runs for real and fails for real.
"""
import sys, os, json, csv, types, math, random, hashlib, argparse, urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
CACHE = os.path.join(HERE, ".http_cache")
os.makedirs(CACHE, exist_ok=True)

# ---------------------------------------------------------------- shims
DG_CALLS = []            # every argument handed to delta_g
URL_LOG = []             # every URL the repo builds

eq = types.ModuleType("equilibrator_api")
class _CC:
    def __init__(self): self.p_h = 7.0; self.ionic_strength = 0.1; self.ccache = self
    def get_compound(self, *a, **k): raise RuntimeError("equilibrator_api not installed")
    def parse_reaction_formula(self, s): raise RuntimeError("equilibrator_api not installed")
    def standard_dg_prime(self, r): raise RuntimeError("equilibrator_api not installed")
class _Rxn:
    @classmethod
    def parse_formula(cls, f, s, rid=None): raise RuntimeError("equilibrator_api not installed")
eq.ComponentContribution = _CC; eq.Reaction = _Rxn
sys.modules["equilibrator_api"] = eq

import requests as _requests
_real_get = _requests.get
def _cached_get(url, *a, **k):
    URL_LOG.append(url)
    key = hashlib.sha256(url.encode()).hexdigest()[:40]
    fp = os.path.join(CACHE, key)
    if os.path.exists(fp):
        body = open(fp, encoding="utf-8").read()
        r = _requests.Response(); r.status_code = 200
        r._content = body.encode(); r.encoding = "utf-8"
        r.url = url
        return r
    r = _real_get(url, *a, **k)
    if r.ok:
        open(fp, "w", encoding="utf-8").write(r.text)
    return r

import utils.thermo as thermo
_real_dg = thermo.delta_g
def _logged_dg(s, *a, **k):
    DG_CALLS.append(s)
    return _real_dg(s, *a, **k)

import ChiraLLM.feasibility as F
F.requests.get = _cached_get
F.delta_g = _logged_dg

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------- metrics
def auroc(scores, labels):
    """Mann-Whitney U with explicit tie handling (0.5 credit for ties)."""
    pos = [s for s, l in zip(scores, labels) if l == 1]
    neg = [s for s, l in zip(scores, labels) if l == 0]
    if not pos or not neg: return float("nan")
    wins = sum((1.0 if p > n else 0.5 if p == n else 0.0) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))

def boot_ci(scores, labels, n=4000, seed=0):
    rng = random.Random(seed)
    idx = list(range(len(scores))); vals = []
    for _ in range(n):
        smp = [rng.choice(idx) for _ in idx]
        s = [scores[i] for i in smp]; l = [labels[i] for i in smp]
        if 1 in l and 0 in l:
            v = auroc(s, l)
            if not math.isnan(v): vals.append(v)
    if not vals: return (float("nan"), float("nan"))
    vals.sort()
    return (vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals))])

def confusion(pred, labels):
    tp = sum(1 for p, l in zip(pred, labels) if p and l == 1)
    fp = sum(1 for p, l in zip(pred, labels) if p and l == 0)
    tn = sum(1 for p, l in zip(pred, labels) if not p and l == 0)
    fn = sum(1 for p, l in zip(pred, labels) if not p and l == 1)
    sens = tp / (tp + fn) if tp + fn else float("nan")
    spec = tn / (tn + fp) if tn + fp else float("nan")
    prec = tp / (tp + fp) if tp + fp else float("nan")
    return dict(tp=tp, fp=fp, tn=tn, fn=fn, sensitivity=sens, specificity=spec,
                precision=prec, balanced_accuracy=(sens + spec) / 2,
                accuracy=(tp + tn) / len(labels))

# ---------------------------------------------------------------- stages
def stereo_blind_smiles(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(m, isomericSmiles=False)))

def stage1_retrieval(rows):
    """Claim 1: does find_compound_kegg return the right KEGG entry?"""
    res = []
    for r in rows:
        name = r["query_name"]
        got = F.find_compound_kegg(name)
        truth = r["truth_kegg_cpd"]
        got_id = (got or "").replace("cpd:", "").strip()
        hit = bool(got_id)
        exact = bool(truth) and got_id == truth
        # right skeleton, wrong stereochemistry?
        skel = None
        if hit and truth and not exact:
            import keggc_shim as K
            a = K.smiles_for(f"cpd:{got_id}"); b = K.smiles_for(f"cpd:{truth}")
            if a and b:
                skel = (stereo_blind_smiles(a) == stereo_blind_smiles(b))
        res.append(dict(query_name=name, stratum=r["stratum"], truth=truth,
                        retrieved=got_id, multiword=(" " in name),
                        hit=hit, exact=exact, same_skeleton_wrong_stereo=bool(skel)))
    return res

def stage2_feasibility(rows, host):
    res = []
    for r in rows:
        DG_CALLS.clear()
        out = F.combined_feasibility_analysis(r["query_name"], r["smiles"], host)
        res.append(dict(
            query_name=r["query_name"], label=r["label"], stratum=r["stratum"],
            mw=r["mw"], n_stereocentres=r["n_stereocentres"],
            verdict=bool(out.get("biosynthetic_feasibility")),
            confidence=float(out.get("confidence", 0.0)),
            kegg_id=out.get("kegg_compound_id") or "",
            n_pathways=len(out.get("kegg_pathways") or []),
            brenda_found=bool((out.get("brenda") or {}).get("found")),
            brenda_ecs=len((out.get("brenda") or {}).get("ec_numbers") or []),
            dg_values=[t["delta_g"] for t in (out.get("thermo") or [])],
            dg_inputs=list(DG_CALLS),
            bottlenecks=json.dumps(out.get("bottlenecks")),
            rdkit=json.dumps(out.get("rdkit_analysis")),
        ))
    return res

def stage3_enantiomers(pairs, host):
    """Claim 3: can the pipeline tell two enantiomers apart at all?"""
    res = []
    for label, (nameA, smiA), (nameB, smiB) in pairs:
        a = F.combined_feasibility_analysis(nameA, smiA, host)
        b = F.combined_feasibility_analysis(nameB, smiB, host)
        keys = ["biosynthetic_feasibility", "confidence", "kegg_compound_id", "bottlenecks"]
        diffs = {k: (a.get(k), b.get(k)) for k in keys if a.get(k) != b.get(k)}
        ra, rb = a.get("rdkit_analysis", {}), b.get("rdkit_analysis", {})
        res.append(dict(pair=label, a=nameA, b=nameB,
                        smiles_differ=(smiA != smiB),
                        rdkit_identical=(ra == rb),
                        differing_fields=sorted(diffs.keys()),
                        detail=json.dumps(diffs, default=str)[:400]))
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="E. coli")
    ap.add_argument("--dataset", default=os.path.join(HERE, "dataset.csv"))
    ap.add_argument("--out", default=os.path.join(HERE, "results.json"))
    ap.add_argument("--stage", default="all")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.dataset)))
    for r in rows:
        r["label"] = int(r["label"]); r["mw"] = float(r["mw"])
        r["n_stereocentres"] = int(r["n_stereocentres"])

    out = {"host": a.host, "n": len(rows)}
    if a.stage in ("all", "1"):
        out["stage1_retrieval"] = stage1_retrieval(rows)
        print("stage 1 done", flush=True)
    if a.stage in ("all", "2"):
        out["stage2_feasibility"] = stage2_feasibility(rows, a.host)
        print("stage 2 done", flush=True)
    if a.stage in ("all", "3"):
        import keggc_shim as K
        pairs = []
        for lab, ea, eb in [("limonene", "cpd:C06099", "cpd:C00521"),
                            ("DOPA", "cpd:C00355", "cpd:C00822"),
                            ("lactate", "cpd:C00186", "cpd:C00256"),
                            ("alanine", "cpd:C00041", "cpd:C00133"),
                            ("malate", "cpd:C00149", "cpd:C00497"),
                            ("tryptophan", "cpd:C00078", "cpd:C00525")]:
            sa, na = K.smiles_for(ea), K.name_for(ea)
            sb, nb = K.smiles_for(eb), K.name_for(eb)
            if sa and sb: pairs.append((lab, (na, sa), (nb, sb)))
        out["stage3_enantiomers"] = stage3_enantiomers(pairs, a.host)
        print("stage 3 done", flush=True)
    out["url_log_sample"] = URL_LOG[:12]
    out["url_log_count"] = len(URL_LOG)
    json.dump(out, open(a.out, "w"), indent=1, default=str)
    print("wrote", a.out)

if __name__ == "__main__":
    main()
