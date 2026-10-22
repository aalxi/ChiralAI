"""Stereochemistry-resolved enzyme lookup, as a library.

The question a biocatalysis group actually asks is not "is this feasible" but
"which enzyme makes THIS enantiomer, and which one makes the mirror image I
don't want?" This answers that from curated Rhea reactions joined to ChEBI
structures. See stereo/README.md for provenance and limits.

``lookup(smiles)`` returns one of four outcomes, and the distinction is the
whole product:

  EXACT       Rhea documents an enzyme acting on this exact stereoisomer.
  WRONG-FORM  The scaffold is known, but only as a different stereoisomer.
              This is what a name-matching tool reports as success.
  SKELETON    Scaffold known; the query did not specify stereochemistry.
  NONE        No curated reaction matches this connectivity.
"""
from __future__ import annotations
import csv, os, functools

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
INDEX = os.path.join(DATA, "stereo_enzyme_index.tsv")
PAIRS = os.path.join(DATA, "enantiocomplementary_pairs.tsv")


@functools.lru_cache(maxsize=1)
def _index():
    rows = list(csv.DictReader(open(INDEX), delimiter="\t"))
    by_key, by_skel = {}, {}
    for r in rows:
        if r.get("stereo_key"):
            by_key.setdefault(r["stereo_key"], r)
        by_skel.setdefault(r["skeleton"], []).append(r)
    return rows, by_key, by_skel


@functools.lru_cache(maxsize=1)
def pairs():
    return list(csv.DictReader(open(PAIRS), delimiter="\t"))


def stereo_key(mol) -> str:
    """InChIKey blocks 1+2: connectivity plus stereochemistry, dropping block 3.

    Block 3 encodes protonation. ChEBI stores physiological ions and chemists
    draw neutral molecules, so (S)-lactic acid and ChEBI's (S)-lactate differ in
    that block alone. Matching the full key returns a confident, wrong
    WRONG-FORM on a perfectly correct query.
    """
    return "-".join(inchi.MolToInchiKey(mol).split("-")[:2])


def lookup(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"outcome": "ERROR", "detail": "unparseable SMILES"}
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                        useLegacyImplementation=False)
    defined = [c for c in centres if c[1] in ("R", "S")]
    ik = inchi.MolToInchiKey(mol)
    skel, key = ik.split("-")[0], stereo_key(mol)
    _rows, by_key, by_skel = _index()

    res = {
        "inchikey": ik, "skeleton": skel, "stereo_key": key,
        "n_stereocentres": len(centres), "n_defined": len(defined),
        "enzymes_for_requested_form": [], "enzymes_for_other_forms": [],
        "selectivity_risk": None, "rhea_ids": [], "compound": None,
    }

    hit = by_key.get(key)
    family = by_skel.get(skel, [])

    if hit:
        res.update(
            outcome="EXACT",
            compound=f"{hit['name']} ({hit['chebi']})",
            enzymes_for_requested_form=[e for e in hit["ec_acting_on_this_isomer"].split(";") if e],
            rhea_ids=[r for r in hit["rhea_ids"].split(";") if r][:8],
            detail="a curated reaction acts on this exact stereoisomer",
        )
        if hit["enantiomer_chebi"]:
            res["selectivity_risk"] = {
                "enantiomer": f"{hit['enantiomer_name']} ({hit['enantiomer_chebi']})",
                "made_by_ec": [e for e in hit["ec_acting_on_enantiomer"].split(";") if e],
                "note": "the mirror image is also enzymatically accessible; "
                        "choosing between these EC families IS the stereochemical decision",
            }
        return res

    if family and defined:
        res.update(
            outcome="WRONG-FORM",
            enzymes_for_other_forms=[
                {"compound": f"{r['name']} ({r['chebi']})",
                 "ec": [e for e in r["ec_acting_on_this_isomer"].split(";") if e]}
                for r in family[:6]],
            detail="the scaffold is curated, but only as a different stereoisomer: "
                   "the enzymes below would make the wrong form",
        )
        return res

    if family:
        res.update(
            outcome="SKELETON",
            enzymes_for_other_forms=[
                {"compound": f"{r['name']} ({r['chebi']})",
                 "ec": [e for e in r["ec_acting_on_this_isomer"].split(";") if e]}
                for r in family[:6]],
            detail="scaffold curated, but the query did not specify stereochemistry",
        )
        return res

    res.update(outcome="NONE",
               detail="no curated reaction matches this connectivity; this is an "
                      "absence of evidence, not a prediction that no route exists")
    return res


def available() -> bool:
    return os.path.exists(INDEX)
