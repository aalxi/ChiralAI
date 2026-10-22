from __future__ import annotations
import os, gzip
from pathlib import Path

# This URL returns HTTP 404 and has done for the whole life of this repo, so
# the BRENDA layer has never produced a result for anyone. BRENDA's bulk
# download now sits behind a registered-user licence, so there is no drop-in
# replacement URL. The layer therefore reports itself UNAVAILABLE rather than
# returning an empty result, which the caller previously could not distinguish
# from "searched BRENDA, found no enzyme" -- a false negative dressed as data.
# To enable it: register at brenda-enzymes.org, download the dump, and point
# BRENDA_JSON_PATH or BRENDA_TEXT_PATH at it.
BRENDA_URL = os.getenv(
    "BRENDA_URL", "https://www.brenda-enzymes.org/brenda_download/brenda.json.gz")

CACHE_DIR = Path(os.getenv("BRENDA_CACHE_DIR", Path.home() / ".cache" / "chiralai"))
JSON_PATH = Path(os.getenv("BRENDA_JSON_PATH", CACHE_DIR / "brenda.json.gz"))
TXT_PATH = Path(os.getenv("BRENDA_TEXT_PATH", CACHE_DIR / "brenda_download.txt"))


def _download_json(path: Path) -> None:
    import requests
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with requests.get(BRENDA_URL, stream=True) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)


def _json_entries(path: Path):
    import ijson
    with gzip.open(path, "rb") as fh:
        for entry in ijson.items(fh, "item"):
            yield entry


def _txt_entries(path: Path, compound: str, host: str | None):
    from brendapyrser import parser as bp
    parser = bp.BrendaParser(path)
    enzymes = parser.search_substrate(compound)
    for e in enzymes:
        yield {
            "ec_number": e.ec,
            "organism": e.organisms,
            "pathway": [e.pathway] if e.pathway else [],
        }


UNAVAILABLE = {
    "status": "unavailable", "found": False,
    "ec_numbers": [], "host_matches": [], "pathways": [],
}


def query_brenda(compound_name: str, host: str | None = None, refresh: bool = False) -> dict:
    """Search a local BRENDA dump for a substrate.

    Returns a dict whose ``status`` is "ok" or "unavailable". Callers must check
    it: an empty result with status "unavailable" means the database was never
    consulted, which is not the same as finding nothing in it.
    """
    path = JSON_PATH
    if not path.exists() and TXT_PATH.exists():
        path = TXT_PATH
    elif refresh or not path.exists():
        try:
            _download_json(path)
        except Exception as exc:
            if TXT_PATH.exists():
                path = TXT_PATH
            else:
                return dict(UNAVAILABLE, reason=(
                    f"no local BRENDA dump ({JSON_PATH}) and the download failed: "
                    f"{type(exc).__name__}: {str(exc)[:120]}"))

    ec_numbers: set[str] = set()
    host_matches: set[str] = set()
    pathways: set[str] = set()

    if path.suffix == ".gz":
        try:
            entries = _json_entries(path)
        except FileNotFoundError:
            entries = []
    else:
        entries = _txt_entries(path, compound_name, host)

    host_l = host.lower() if host else None
    for entry in entries:
        subs = entry.get("substrate") or entry.get("substrates") or []
        if isinstance(subs, str):
            subs = [subs]
        if not any(compound_name.lower() in str(s).lower() for s in subs):
            continue
        ec = entry.get("ec_number")
        if ec:
            ec_numbers.add(str(ec))
            if host_l:
                orgs = entry.get("organism") or []
                for org in orgs:
                    if host_l in str(org).lower():
                        host_matches.add(str(ec))
        pth = entry.get("pathway") or []
        if isinstance(pth, str):
            pathways.add(pth)
        else:
            pathways.update([str(p) for p in pth])
    return {
        "status": "ok",
        "found": bool(ec_numbers),
        "ec_numbers": sorted(ec_numbers),
        "host_matches": sorted(host_matches),
        "pathways": sorted(pathways),
    }