"""Cached KEGG helper used only to build ground truth, never by the system under test."""
import os, hashlib, time, urllib.parse, requests
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".kegg_cache")
os.makedirs(CACHE, exist_ok=True)
_S = requests.Session()

def _get(path):
    key = hashlib.sha256(path.encode()).hexdigest()[:40]
    fp = os.path.join(CACHE, key)
    if os.path.exists(fp): return open(fp, encoding="utf-8").read()
    seg = path.split("/")
    url = "https://rest.kegg.jp/" + "/".join(seg[:2] + [urllib.parse.quote(s, safe=":+") for s in seg[2:]])
    body = ""
    for i in range(5):
        try:
            r = _S.get(url, timeout=40)
            if r.status_code == 404: break
            if r.ok: body = r.text; break
            time.sleep(1.2 * 2 ** i)
        except Exception:
            time.sleep(1.2 * 2 ** i)
    open(fp, "w", encoding="utf-8").write(body); time.sleep(0.8)
    return body

def smiles_for(eid):
    txt = _get(f"get/{eid}/mol")
    if not txt.strip(): return None
    m = Chem.MolFromMolBlock(txt)
    if m is None: return None
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return Chem.MolToSmiles(m)

def name_for(eid):
    txt = _get(f"get/{eid}")
    for line in txt.split("\n"):
        if line.startswith("NAME"):
            return line.replace("NAME", "").strip().rstrip(";")
    return eid
