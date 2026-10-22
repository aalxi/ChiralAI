"""Retrieve a KEGG compound by structure, not by name.

Benchmarked on 36 compounds (benchmarks/RESULTS.md): name-first retrieval got
94% top-1 and failed on stereochemical queries 25 times out of 32;
structure-first got 100% and is stereochemistry-exact by construction.

Why name matching fails here specifically: KEGG stores one chosen set of
synonyms per compound, and only 22% of the stereochemical names a chemist would
write resolve at all. Worse, the first hit is often the wrong molecule --
"farnesene" returns the alpha isomer when the yeast product is beta, and
"glucose" returns UDP-glucose. A structure cannot be ambiguous in that way.

The ladder is: exact molecular formula narrows the candidate set, then InChIKey
decides, strictest match first.
"""
from __future__ import annotations

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from rdkit.Chem import rdMolDescriptors as rmd

from utils import kegg

RDLogger.DisableLog("rdApp.*")


def _keys(mol) -> tuple[str, str, str]:
    """(full InChIKey, connectivity+stereo, connectivity only)."""
    ik = inchi.MolToInchiKey(mol)
    return ik, "-".join(ik.split("-")[:2]), ik.split("-")[0]


def find_by_structure(smiles: str, cap: int = 25) -> dict:
    """Resolve a SMILES to a KEGG compound id.

    ``how`` records the strictness of the match, and callers must look at it:
    a "skeleton" match means KEGG has the scaffold but could not confirm which
    stereoisomer, which is exactly the case where a stereochemical claim must
    not be made.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"kegg_id": None, "how": "unparseable SMILES", "n_candidates": 0}
    formula = rmd.CalcMolFormula(mol).replace("+", "").replace("-", "")
    want, want_stereo, want_skel = _keys(mol)

    candidates = kegg.find_by_formula(formula)
    skeleton_hit = None
    for cid in candidates[:cap]:
        block = kegg.mol_block(cid)
        if not block:
            continue
        ref = Chem.MolFromMolBlock(block)
        if ref is None:
            continue
        got, got_stereo, got_skel = _keys(ref)
        if got == want:
            return {"kegg_id": cid, "how": "exact InChIKey", "formula": formula,
                    "n_candidates": len(candidates), "stereo_exact": True}
        if got_stereo == want_stereo:
            return {"kegg_id": cid, "how": "InChIKey blocks 1+2 (ionisation differs)",
                    "formula": formula, "n_candidates": len(candidates),
                    "stereo_exact": True}
        if got_skel == want_skel and skeleton_hit is None:
            skeleton_hit = cid
    if skeleton_hit:
        return {"kegg_id": skeleton_hit, "formula": formula,
                "how": "skeleton only - KEGG's record does not pin the stereoisomer",
                "n_candidates": len(candidates), "stereo_exact": False}
    return {"kegg_id": None, "how": "no structural match in KEGG", "formula": formula,
            "n_candidates": len(candidates), "stereo_exact": False}


def resolve(name: str, smiles: str) -> dict:
    """Structure first; fall back to name search only if that finds nothing.

    The fallback is recorded, never hidden, because a name-resolved hit carries
    no guarantee about stereochemistry.
    """
    res = find_by_structure(smiles)
    if res["kegg_id"]:
        res["method"] = "structure"
        return res
    cands = kegg.find_by_name(name)
    if cands:
        return {"kegg_id": cands[0][0], "method": "name (fallback)",
                "how": f"name search, {len(cands)} candidates, first taken - "
                       f"stereochemistry NOT verified",
                "n_candidates": len(cands), "stereo_exact": False,
                "name_candidates": [c for c, _ in cands[:8]]}
    return dict(res, method="none")
