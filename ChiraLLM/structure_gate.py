"""Does the LLM's structure actually correspond to the molecule it named?

Stage A is a language model. It produces, with equal fluency, a correct SMILES
and a several-hundred-character degenerate string labelled "limonene" (that is
in this repo's own saved output, ``suggestions_20250702_155147.csv``). Nothing
downstream can tell the difference, because every downstream check keys on the
*name*.

The gate runs before feasibility, not after, and asks three questions:

1. Does RDKit parse it, and does it survive a canonical round trip?
2. Does its molecular formula agree with the KEGG entry the model cited, or
   failing that with the compound KEGG returns for the name?
3. Is the stereochemistry specified at all?

Question 3 never rejects. A molecule with unspecified stereocentres is a
legitimate query; it just cannot have a stereochemical route selected for it,
and the caller needs to be told that rather than silently handed one.
"""
from __future__ import annotations

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem import rdMolDescriptors as rmd

RDLogger.DisableLog("rdApp.*")


def _formula(mol) -> str:
    """Formula with charge markers stripped, so an anion matches its acid."""
    return rmd.CalcMolFormula(mol).replace("+", "").replace("-", "")


def _norm(text: str) -> str:
    """Compare names ignoring case, spacing and hyphenation: KEGG writes
    "L-Lactic acid" where a user types "l lactic acid"."""
    return "".join(ch for ch in text.lower() if ch.isalnum())


def kegg_name_candidates(name: str):
    from utils import kegg
    return kegg.find_by_name(name)


def _kegg_formula(compound_id: str) -> str | None:
    from utils import kegg
    txt = kegg._get(f"get/{compound_id.replace('cpd:', '')}")
    for line in txt.split("\n"):
        if line.startswith("FORMULA"):
            return line.replace("FORMULA", "").strip()
    return None


def check(name: str, smiles: str, cited_kegg_id: str | None = None,
          *, cross_check: bool = True) -> dict:
    """Return a verdict dict. ``ok`` False means do not analyse this row."""
    out = {
        "name": name, "smiles_in": smiles, "cited_kegg_id": cited_kegg_id or "",
        "ok": False, "rejected_because": None, "warnings": [],
    }
    if not smiles or not isinstance(smiles, str):
        out["rejected_because"] = "no SMILES supplied"
        return out

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        out["rejected_because"] = "RDKit cannot parse the SMILES"
        return out

    canonical = Chem.MolToSmiles(mol)
    if Chem.MolFromSmiles(canonical) is None:
        out["rejected_because"] = "SMILES does not survive a canonical round trip"
        return out

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                        useLegacyImplementation=False)
    defined = [c for c in centres if c[1] in ("R", "S")]

    out.update(
        canonical_smiles=canonical,
        inchikey=inchi.MolToInchiKey(mol),
        formula=_formula(mol),
        molecular_weight=round(Descriptors.MolWt(mol), 2),
        heavy_atoms=mol.GetNumHeavyAtoms(),
        n_stereocentres=len(centres),
        n_defined_stereocentres=len(defined),
        stereo_tags=",".join(f"{i}{t}" for i, t in centres),
    )

    if centres and not defined:
        out["warnings"].append(
            "stereocentres present but none specified: no enantiomer-specific "
            "route can be selected for this structure")
    if not centres:
        out["warnings"].append("no stereocentres: this molecule is not chiral")

    # ── formula agreement against an independent database record ──────────────
    # Two false rejections had to be designed out of this check, and both came
    # from treating a loose name hit as authoritative:
    #
    #  * KEGG returns UDP-glucose FIRST for "glucose", so comparing against the
    #    first hit rejects correct glucose.
    #  * KEGG has no entry for free L-glucose, only GDP-L-glucose, so a substring
    #    match rejects correct L-glucose against a nucleotide sugar's formula.
    #
    # So the comparison is only made against candidates for which the queried
    # name is a genuine full synonym, not merely a substring of one. When no such
    # candidate exists there is nothing to verify against, and the structure is
    # reported as unverified rather than condemned.
    if cross_check:
        wanted = _norm(name)
        strict, loose = [], []
        if cited_kegg_id:
            cid = cited_kegg_id.replace("cpd:", "")
            f = _kegg_formula(cid)
            if f:
                strict.append((cid, f, "cited KEGG id"))
        for cid, syns in kegg_name_candidates(name)[:20]:
            if any(_norm(x) == wanted for x in syns.split(";")):
                f = _kegg_formula(cid)
                if f:
                    strict.append((cid, f, "KEGG synonym match"))
                    if f == out["formula"]:
                        break
            else:
                loose.append(cid)

        out["kegg_synonym_matches"] = [f"{c}:{f}" for c, f, _ in strict]
        out["kegg_loose_matches"] = loose[:8]
        matched = next(((c, f, h) for c, f, h in strict if f == out["formula"]), None)
        if matched:
            out["reference_kegg_id"], out["reference_formula"], out["reference_via"] = matched
            out["warnings"].append(
                f"formula {out['formula']} agrees with KEGG {matched[0]} ({matched[2]})")
        elif strict:
            out["rejected_because"] = (
                f"formula mismatch: the structure is {out['formula']}, but every KEGG "
                f"compound whose synonyms include {name!r} has a different formula "
                f"({', '.join(f for _, f, _ in strict[:6])}). The name and the "
                f"structure are different molecules.")
            return out
        else:
            out["warnings"].append(
                f"no KEGG compound has {name!r} as a synonym, so the structure could "
                f"not be cross-checked and is UNVERIFIED"
                + (f" (loose name hits, not used: {', '.join(loose[:4])})" if loose else ""))

    out["ok"] = True
    return out
