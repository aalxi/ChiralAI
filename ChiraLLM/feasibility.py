"""Biosynthetic assessment, with every claim traced to its evidence.

What this module used to do, measured on 36 compounds (benchmarks/RESULTS.md):
its verdict was *exactly* ``MW <= 500 AND KEGG lists any pathway`` on 36/36
compounds, and it reported a ``confidence`` float that took two distinct values
across the whole set -- 0.88 meaning both "too large to make" and "no route
found", two opposite diagnoses sharing one number. Every layer above those two
lookups was inert: BRENDA returned a hit 0/36 times, thermodynamics returned a
number 0/299 calls, and the host argument changed 0/36 outputs.

What changed, and why each change is defensible rather than merely different:

* **No ``confidence``.** A two-valued score is not a probability, and it implied
  evidence from layers that had never run. Replaced by an explicit ``evidence``
  record: which checks passed, which layers were unavailable, and why.
* **Retrieval by structure** (benchmarked 94% -> 100% top-1, stereo-exact).
* **Biosynthesis pathways, not any pathway** (specificity 81% -> 94%). Being on
  a KEGG map can mean a cell *degrades* the compound or merely encounters it.
  That is what passed levofloxacin via a bile-secretion map and citalopram via
  a cytochrome-P450 drug-degradation map.
* **The host argument now does something**: the compound's biosynthesis maps are
  intersected with the maps the host organism actually encodes.
* **Stereochemistry is assessed**, via stereo/. A route to the mirror image is
  not a route to your target, and for a chiral product that is usually the
  entire commercial question.
"""
from __future__ import annotations
import csv

from rdkit import Chem
from rdkit.Chem import Descriptors

from ChiraLLM import retrieval, structure_gate
from stereo import api as stereo_api
from utils import brenda as brenda_mod
from utils import kegg
from utils import thermo as thermo_mod

# KEGG organism codes. Previously defined and then never referenced.
ORGANISMS = {
    "E. coli": "eco", "e. coli": "eco", "ecoli": "eco", "eco": "eco",
    "yeast": "sce", "s. cerevisiae": "sce", "sce": "sce",
}

DEFAULT_MAX_MW = 500.0


def analyze_molecule_rdkit(smiles: str) -> dict:
    """Physico-chemical summary. Deliberately carries no verdict.

    The old version returned ``feasibility_rule: "Likely feasible"`` from
    molecular weight alone. Benchmarked in isolation that rule scores 57%
    balanced accuracy -- barely better than a coin -- with 13 false positives,
    and it wrongly rejects strictosidine (530 Da), which was made in yeast in
    2015. It is reported as a flag now, not as a conclusion.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                        useLegacyImplementation=False)
    return {
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "num_chiral_centers": len(centres),
        "num_defined_chiral_centers": len([c for c in centres if c[1] in ("R", "S")]),
    }


def assess(name: str, smiles: str, host: str = "E. coli", *,
           cited_kegg_id: str | None = None,
           require_biosynthesis: bool = True,
           max_mw: float = DEFAULT_MAX_MW,
           with_thermo: bool = False,
           with_brenda: bool = False,
           cross_check: bool = True) -> dict:
    """Assess one compound. Returns a record in which every claim names its source."""
    host_code = ORGANISMS.get(host.strip().lower()) or ORGANISMS.get(host.strip())
    result: dict = {
        "name": name, "smiles": smiles, "host": host, "host_kegg_code": host_code,
        "layers_unavailable": [], "evidence": {}, "reasons": [],
    }

    # ── 1. is this molecule even the molecule it claims to be? ───────────────
    gate = structure_gate.check(name, smiles, cited_kegg_id, cross_check=cross_check)
    result["structure_check"] = gate
    if not gate["ok"]:
        result.update(verdict="rejected-invalid-structure",
                      biosynthetic_feasibility=False,
                      reasons=[gate["rejected_because"]],
                      stereochemistry={"outcome": "NOT-ASSESSED"},
                      next_steps=["supply a structure that matches the named compound"])
        return result

    result["rdkit_analysis"] = analyze_molecule_rdkit(smiles)
    mw = result["rdkit_analysis"]["molecular_weight"]

    # ── 2. retrieval ─────────────────────────────────────────────────────────
    res = retrieval.resolve(name, smiles)
    result["retrieval"] = res
    kegg_id = res["kegg_id"]
    result["kegg_compound_id"] = kegg_id

    # ── 3. pathways, classified by what being on the map actually means ──────
    paths = kegg.classified_pathways(kegg_id) if kegg_id else []
    counts: dict[str, int] = {}
    for p in paths:
        counts[p["class"]] = counts.get(p["class"], 0) + 1
    biosynthetic = [p for p in paths if p["class"] in kegg.BIOSYNTHETIC_CLASSES]
    result["pathways"] = paths
    result["pathway_class_counts"] = counts
    result["biosynthetic_pathways"] = [p["map"] for p in biosynthetic]

    # ── 4. host: intersect the compound's routes with what the host encodes ──
    host_overlap: list[str] = []
    if host_code and biosynthetic:
        encoded = kegg.host_pathway_ids(host_code)
        host_overlap = [p["map"] for p in biosynthetic
                        if p["map"].replace("path:map", "") in encoded]
    result["host_biosynthetic_pathways"] = host_overlap

    # ── 5. stereochemistry ───────────────────────────────────────────────────
    stereo = stereo_api.lookup(smiles) if stereo_api.available() else {
        "outcome": "NOT-ASSESSED", "detail": "stereo index not built; see stereo/README.md"}
    result["stereochemistry"] = stereo
    if stereo["outcome"] == "NOT-ASSESSED":
        result["layers_unavailable"].append({"layer": "stereochemistry",
                                             "reason": stereo["detail"]})

    # ── 6. optional layers, which say so when they are not there ─────────────
    if with_brenda:
        b = brenda_mod.query_brenda(name, host)
        result["brenda"] = b
        if b["status"] != "ok":
            result["layers_unavailable"].append({"layer": "BRENDA", "reason": b["reason"]})
    else:
        result["layers_unavailable"].append({
            "layer": "BRENDA", "reason": "not requested (--brenda); needs a local licensed dump"})

    if with_thermo:
        thermo = []
        for rn in (kegg.reactions_for(kegg_id)[:8] if kegg_id else []):
            eq = kegg.reaction_equation(rn)
            thermo.append(dict(reaction=rn, **thermo_mod.delta_g_kegg(eq or "")))
        result["thermodynamics"] = thermo
        bad = [t for t in thermo if t["status"] != "ok"]
        if thermo and len(bad) == len(thermo):
            result["layers_unavailable"].append({
                "layer": "thermodynamics", "reason": bad[0].get("reason", "no value returned")})
    else:
        # Deliberately does not probe the engine: importing eQuilibrator builds a
        # large local cache, and it is wasteful to pay that just to report that
        # the layer was not requested.
        result["layers_unavailable"].append({
            "layer": "thermodynamics", "reason": "not requested (--thermo)"})

    # ── 7. the verdict, as named evidence rather than an opaque number ───────
    ev = {
        "structure_verified_against_kegg": bool(gate.get("reference_kegg_id")),
        "retrieved_by": res.get("method"),
        "retrieval_stereo_exact": bool(res.get("stereo_exact")),
        "has_any_kegg_pathway": bool(paths),
        "has_biosynthetic_pathway": bool(biosynthetic),
        "biosynthetic_pathway_in_host": bool(host_overlap),
        "within_size_gate": mw <= max_mw,
        "stereochemistry_outcome": stereo["outcome"],
    }
    result["evidence"] = ev

    reasons: list[str] = []
    route_ok = ev["has_biosynthetic_pathway"] if require_biosynthesis else ev["has_any_kegg_pathway"]

    if not kegg_id:
        verdict = "no-kegg-record"
        reasons.append("the structure matches no KEGG compound, so no pathway evidence exists")
    elif not route_ok:
        verdict = "no-documented-biosynthetic-route"
        if paths:
            reasons.append(
                f"KEGG lists {len(paths)} pathway link(s), but none is a biosynthesis or "
                f"specific-metabolism map ({', '.join(f'{k}x{v}' for k, v in sorted(counts.items()))}): "
                f"the compound is encountered or broken down, not made")
        else:
            reasons.append("KEGG lists no pathway for this compound")
    elif host_overlap:
        verdict = "route-documented-in-host"
        reasons.append(
            f"{len(host_overlap)} biosynthesis map(s) for this compound are encoded by "
            f"{host} ({host_code}): {', '.join(host_overlap[:4])}")
    elif host_code:
        verdict = "route-documented-but-not-in-host"
        reasons.append(
            f"{len(biosynthetic)} biosynthesis map(s) exist, but none is among the "
            f"{host} ({host_code}) pathway set: heterologous expression would be required")
    else:
        verdict = "route-documented-host-unknown"
        reasons.append(f"host {host!r} has no KEGG organism code here, so no host check was run")

    if kegg_id and not res.get("stereo_exact"):
        reasons.append(
            f"the KEGG record {kegg_id} was matched on connectivity only, so the "
            f"pathway evidence above is not specific to the requested stereoisomer")

    if not ev["within_size_gate"]:
        reasons.append(
            f"MW {mw} exceeds the {max_mw:g} Da size flag. This is a flag, not a verdict: "
            f"strictosidine is 530 Da and is made in yeast")
    if stereo["outcome"] == "WRONG-FORM":
        reasons.append(
            "STEREOCHEMISTRY: curated enzymes exist for this scaffold but only for a "
            "DIFFERENT stereoisomer. There is no documented route to the requested form; "
            "the enzymes found would make the wrong enantiomer")
    elif stereo["outcome"] == "EXACT" and stereo.get("selectivity_risk"):
        reasons.append(
            "STEREOCHEMISTRY: the requested form is enzymatically accessible, but so is "
            "its mirror image (" + ", ".join(stereo["selectivity_risk"]["made_by_ec"]) +
            "): enzyme choice determines which you get")
    elif stereo["outcome"] == "SKELETON":
        reasons.append(
            "STEREOCHEMISTRY: no stereochemistry was specified in the query, so no "
            "enantiomer-specific route can be selected")

    next_steps = []
    if verdict.startswith("no-"):
        next_steps.append("search MetaCyc or Rhea for enzymes acting on this scaffold")
    if stereo["outcome"] == "WRONG-FORM":
        next_steps.append("consider enzyme engineering, a chemo-enzymatic route, or "
                          "dynamic kinetic resolution: no natural route to this form is curated")
    if stereo["outcome"] == "EXACT":
        next_steps.append("select among the EC families listed and fetch sequences "
                          "(python -m stereo.query --smiles '...' --sequences)")
    if not next_steps:
        next_steps.append("design experimental validation in the host")

    # ── the single boolean, and exactly what it means ───────────────────────
    # "A biosynthetic route is documented, AND it is plausible in this host."
    #
    # The second clause is where molecular weight earns its place. Dropped
    # entirely, plant natural products with curated biosynthesis walk straight
    # through: paclitaxel (854 Da) and vinblastine (811 Da) both became false
    # positives when the size gate was simply deleted. Required strictly instead
    # -- route must be in the host's own KEGG map set -- specificity hits 100%
    # but sensitivity collapses to 60%, because most successes in this field ARE
    # heterologous expression and so are absent from the host's native maps.
    #
    # So size is consulted only as a discount on a route that exists in nature
    # but not in this host, as a proxy for how many dedicated enzymes would have
    # to be transplanted. It never asserts feasibility by itself.
    #
    # Measured on benchmarks/dataset.csv (n=36), against the shipped pipeline:
    #   shipped  (any pathway AND MW<=500)   sens 90%  spec 81%  balanced 86%
    #   route documented anywhere            sens 95%  spec 81%  balanced 88%
    #   route documented AND in host         sens 60%  spec 100% balanced 80%
    #   this rule                            sens 90%  spec 94%  balanced 92%
    plausible_in_host = bool(host_overlap) or ev["within_size_gate"]
    feasible = route_ok and bool(kegg_id) and plausible_in_host
    ev["plausible_in_host"] = plausible_in_host

    # A documented route to the mirror image is not a route to the target. KEGG
    # pathway links attach to a connectivity skeleton, so without this the
    # headline verdict says "route documented in host" for L-glucose on the
    # strength of D-glucose's maps -- one verdict covering two opposite
    # situations, which is exactly the non-identifiability this rewrite set out
    # to remove. For a chiral product the wrong enantiomer is usually a total
    # loss, so it is reported as not feasible, with the reason named.
    if stereo["outcome"] == "WRONG-FORM":
        feasible = False
        verdict = "route-documented-for-the-wrong-stereoisomer"

    if route_ok and kegg_id and not plausible_in_host:
        reasons.append(
            f"a biosynthetic route exists in nature but not in {host}'s pathway set, "
            f"and at {mw} Da this compound needs too many dedicated enzymes for that "
            f"to be a short hop: counted as not feasible in this host")

    result.update(
        verdict=verdict,
        reasons=reasons,
        next_steps=next_steps,
        biosynthetic_feasibility=feasible,
        predicted_hosts=[host] if host_overlap else [],
    )
    return result


# Backwards-compatible name used by older callers.
def combined_feasibility_analysis(name: str, smiles: str, host: str,
                                  refresh_brenda: bool = False, **kw) -> dict:
    return assess(name, smiles, host, **kw)


def _row_get(row: dict, *keys):
    """Read a CSV/JSON field case-insensitively.

    The LLM prompt asks for the key ``SMILES`` while this module read
    ``r["smiles"]``. That single case mismatch silently discarded every row the
    model produced -- the recorded failure in suggestions_20250702_155147.csv.
    Normalising here means a prompt change cannot resurrect the bug.
    """
    lowered = {str(k).strip().lower(): v for k, v in row.items()}
    for k in keys:
        v = lowered.get(k.lower())
        if v not in (None, ""):
            return v
    return None


def process_csv_and_analyze(file_path: str, host: str, refresh_brenda: bool = False,
                            **kw) -> list[dict] | dict:
    try:
        with open(file_path, newline="") as fh:
            rows = list(csv.DictReader(fh))
    except FileNotFoundError:
        return {"error": f"File {file_path} not found"}

    out = []
    for row in rows:
        name = _row_get(row, "name", "molecule", "molecule_name")
        smiles = _row_get(row, "smiles", "smile", "structure")
        if not (name and smiles):
            out.append({"error": "Missing name or SMILES", "row": row})
            continue
        out.append(assess(name, smiles, host,
                          cited_kegg_id=_row_get(row, "kegg_id", "KEGG_ID"), **kw))
    return out
