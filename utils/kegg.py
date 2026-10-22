"""One cached KEGG client.

Replaces three partial, divergent copies of this logic that lived in
``ChiraLLM/feasibility.py``, ``utils/utils.py`` and
``ChiraLLM/database_validator.py``.

Two things here are not cosmetic:

* **HTTPS.** The old code called ``http://rest.kegg.jp``, which 301-redirects.
  Anything not following redirects silently received an empty body and read it
  as "compound not found".
* **Caching with backoff.** KEGG throttles a bare loop of requests and starts
  refusing connections part-way through a run, which looks exactly like a run
  of negative results. Responses are cached on disk, so a rerun is offline,
  free, and reproducible.
"""
from __future__ import annotations
import hashlib, os, re, time
from pathlib import Path

import requests

KEGG_REST = "https://rest.kegg.jp"
CACHE_DIR = Path(os.getenv("CHIRALAI_KEGG_CACHE", Path.home() / ".cache" / "chiralai" / "kegg"))
TIMEOUT = 30
RETRIES = 4

# Maps that link a compound to the whole of metabolism, or to the fact that an
# organism merely encounters it. Presence in one of these is not evidence that
# a cell can make the compound. See benchmarks/RESULTS.md.
GLOBAL_MAPS = {
    "map01100", "map01110", "map01120", "map01200", "map01210", "map01212",
    "map01230", "map01232", "map01250", "map01240", "map01060", "map01061",
    "map01062", "map01063", "map01064", "map01065", "map01066", "map01070",
    "map01002", "map01003",
}

BIOSYNTHETIC_CLASSES = {"biosynthesis", "metabolism"}


def _cache_path(path: str) -> Path:
    return CACHE_DIR / (hashlib.sha256(path.encode()).hexdigest() + ".txt")


def _get(path: str, *, use_cache: bool = True) -> str:
    """GET a KEGG REST path, cached on disk. Returns "" for a 404.

    An empty body is cached too: "this compound has no pathway links" is a real
    answer and re-asking it every run is what gets a client throttled.
    """
    cp = _cache_path(path)
    if use_cache and cp.exists():
        return cp.read_text()
    last = None
    for attempt in range(RETRIES):
        try:
            r = requests.get(f"{KEGG_REST}/{path}", timeout=TIMEOUT)
            if r.status_code == 404:
                text = ""
            elif r.ok:
                text = r.text
            else:
                raise requests.HTTPError(f"HTTP {r.status_code}")
            CACHE_DIR.mkdir(parents=True, exist_ok=True)
            cp.write_text(text)
            return text
        except Exception as exc:          # throttling, DNS, reset connections
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"KEGG unreachable for {path!r} after {RETRIES} tries: {last}")


# ─── compound retrieval ──────────────────────────────────────────────────────

def find_by_name(name: str) -> list[tuple[str, str]]:
    """All name-search candidates as (compound_id, synonyms). Never just the first.

    The first hit is not the best hit: searching "farnesene" returns
    alpha-farnesene, while the yeast product is the beta isomer.
    """
    txt = _get(f"find/compound/{requests.utils.quote(name)}")
    out = []
    for line in txt.strip().split("\n"):
        if "\t" in line:
            cid, syn = line.split("\t", 1)
            out.append((cid.replace("cpd:", ""), syn.strip()))
    return out


def find_by_formula(formula: str) -> list[str]:
    """Compound ids whose formula matches exactly."""
    txt = _get(f"find/compound/{formula}/formula")
    out = []
    for line in txt.strip().split("\n"):
        if "\t" not in line:
            continue
        cid, got = line.split("\t")[0], line.split("\t")[1].strip()
        if got == formula:
            out.append(cid.replace("cpd:", ""))
    return out


def mol_block(compound_id: str) -> str | None:
    txt = _get(f"get/{compound_id.replace('cpd:', '')}/mol")
    return txt if txt.strip() else None


# ─── pathways ────────────────────────────────────────────────────────────────

def pathways_for(compound_id: str) -> list[str]:
    txt = _get(f"link/pathway/cpd:{compound_id.replace('cpd:', '')}")
    return [l.split("\t")[1] for l in txt.strip().split("\n") if "\t" in l]


def pathway_name(map_id: str) -> str:
    txt = _get(f"get/{map_id.replace('path:', '')}")
    for line in txt.split("\n"):
        if line.startswith("NAME"):
            return line.replace("NAME", "").strip()
    return "?"


def classify_pathway(map_id: str, name: str | None = None) -> str:
    """Why is this compound on this map? Encountering it is not making it.

    Ordering matters: the degradation and xenobiotic tests must run before the
    generic "metabolism" test, because "Drug metabolism - cytochrome P450"
    contains the word metabolism but is the reason citalopram was passed as
    biosynthesable.
    """
    m = map_id.replace("path:", "")
    n = (name if name is not None else pathway_name(map_id)).lower()
    if m in GLOBAL_MAPS:
        return "global_overview"
    if "drug metabolism" in n or "xenobiotic" in n:
        return "xenobiotic_degradation"
    if "degradation" in n or "catabolism" in n:
        return "degradation"
    if re.match(r"map0[45]\d{3,4}", m):
        return "physiology_or_disease"
    if m.startswith("map07"):
        return "drug_classification"
    if "biosynthe" in n:
        return "biosynthesis"
    if "metabolism" in n:
        return "metabolism"
    return "other"


def classified_pathways(compound_id: str) -> list[dict]:
    out = []
    for m in pathways_for(compound_id):
        nm = pathway_name(m)
        out.append({"map": m, "name": nm, "class": classify_pathway(m, nm)})
    return out


def host_pathway_ids(org_code: str) -> set[str]:
    """The map numbers an organism actually encodes, e.g. eco00260 -> 00260.

    This is the only genuinely host-specific free signal in KEGG, and the thing
    the shipped code promised but never computed: it defined ``eco``/``sce`` and
    then never referenced them.
    """
    txt = _get(f"list/pathway/{org_code}")
    ids = set()
    for line in txt.strip().split("\n"):
        mid = line.split("\t")[0].replace("path:", "")
        n = re.sub(r"^[a-z]{3,4}", "", mid)
        if n:
            ids.add(n)
    return ids


# ─── reactions (for thermodynamics) ──────────────────────────────────────────

def reactions_for(compound_id: str) -> list[str]:
    txt = _get(f"link/reaction/cpd:{compound_id.replace('cpd:', '')}")
    return [l.split("\t")[1].replace("rn:", "") for l in txt.strip().split("\n") if "\t" in l]


def reaction_equation(rn_id: str) -> str | None:
    """The stoichiometric equation, e.g. 'C00022 + C00003 <=> C00186 + C00004'.

    This is what a thermodynamics engine needs. The shipped code passed it
    ``path:map00260`` instead — a pathway map identifier — which is why 299
    calls produced zero numbers.
    """
    txt = _get(f"get/{rn_id.replace('rn:', '')}")
    for line in txt.split("\n"):
        if line.startswith("EQUATION"):
            return line.replace("EQUATION", "").strip()
    return None
