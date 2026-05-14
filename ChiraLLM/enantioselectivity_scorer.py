import math
import re
import logging
from typing import Optional

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, DataStructs

logger = logging.getLogger(__name__)

_smiles_cache: dict[str, str | None] = {}
# includeChirality=False: BRENDA substrate names resolve to connectivity SMILES via PubChem,
# so chirality bits would deflate every Tanimoto score against the stereo-specific query molecule.
_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=False)

W_EE_TAN = 0.70  # combined weight for best enzyme score (ee=0.50 + tanimoto=0.20)
W_FEAS   = 0.20
W_STEREO = 0.10


def _parse_ee_percentage(ee_str: object) -> Optional[float]:
    if ee_str is None:
        return None
    s = str(ee_str).strip().lower()
    if not s or s == "unknown":
        return None
    m = re.search(r"(\d+(?:\.\d+)?)", s)
    if not m:
        return None
    return float(m.group(1))


# Cofactors that appear in BRENDA reaction strings but carry no substrate-identity information.
# Filtering these from "substrate + cofactor + H+" gives the primary substrate compound name.
_COFACTORS = frozenset({
    "nadph", "nadp+", "nadp", "nadh", "nad+", "nad",
    "h+", "h2o", "o2", "co2", "atp", "adp", "amp",
    "fad", "fadh2", "fmn", "fmnh2", "coenzyme a", "coa",
    "pyridoxal 5'-phosphate", "plp", "thiamine diphosphate", "tpp",
})


def _primary_substrate(reaction_string: str) -> str:
    """Extract the primary (non-cofactor) substrate from a BRENDA reaction string.

    BRENDA's substrates field contains the full reaction equation, e.g.:
    "acetophenone + NADPH + H+" → returns "acetophenone"
    Returns the original string unchanged if no cofactor splitting applies.
    """
    parts = [p.strip() for p in reaction_string.split(" + ")]
    non_cofactor = [p for p in parts if p.lower() not in _COFACTORS and p]
    return non_cofactor[0] if len(non_cofactor) == 1 else reaction_string


def _smiles_for_substrate(substrate_name: str) -> Optional[str]:
    # Strip reaction-equation cofactors before PubChem lookup
    clean_name = _primary_substrate(substrate_name)
    key = clean_name.strip().lower()
    if key in _smiles_cache:
        return _smiles_cache[key]
    smiles = None
    try:
        results = pcp.get_properties(["IsomericSMILES", "CanonicalSMILES"], clean_name, "name")
        if results:
            # Key names vary by pubchempy version: try all known aliases
            r = results[0]
            candidate = r.get("SMILES") or r.get("IsomericSMILES") or r.get("CanonicalSMILES")
            if candidate and Chem.MolFromSmiles(candidate) is not None:
                smiles = candidate
    except Exception:
        pass
    _smiles_cache[key] = smiles
    return smiles


def _compute_tanimoto(query_smiles: str, substrate_smiles: str) -> Optional[float]:
    mol_q = Chem.MolFromSmiles(query_smiles)
    mol_s = Chem.MolFromSmiles(substrate_smiles)
    if mol_q is None or mol_s is None:
        return None
    fp_q = _morgan_gen.GetFingerprint(mol_q)
    fp_s = _morgan_gen.GetFingerprint(mol_s)
    return DataStructs.TanimotoSimilarity(fp_q, fp_s)


def _score_enzyme_entry(entry: dict, query_smiles: Optional[str]) -> dict:
    ec = entry.get("ec_number", "unknown")
    ee_raw = entry.get("enantioselectivity")
    # New brenda_client parses ee to float; old path may pass a string — handle both.
    ee_value = ee_raw if isinstance(ee_raw, (int, float)) else _parse_ee_percentage(ee_raw)
    ee_norm = min(ee_value / 100.0, 1.0) if ee_value is not None else 0.0

    tanimoto = None
    # brenda_client uses "substrates" (plural, matching BRENDA field name)
    substrate_name = entry.get("substrates") or entry.get("substrate")
    if query_smiles and substrate_name:
        sub_smiles = _smiles_for_substrate(substrate_name)
        if sub_smiles:
            tanimoto = _compute_tanimoto(query_smiles, sub_smiles)

    if tanimoto is not None:
        score = (ee_norm * 0.50 + tanimoto * 0.20) / 0.70
    else:
        score = ee_norm

    return {
        "ec_number": ec,
        "score": round(score, 4),
        "ee_value": ee_value,
        "ee_source": "brenda_verified",
        "tanimoto_similarity": round(tanimoto, 4) if tanimoto is not None else None,
        "organism": entry.get("organism"),
        "engineered_variants": None,  # requires BRENDA getEngineeringInformation; deferred to v2
    }


def _compose_route_ee(route: dict, brenda_data: dict) -> Optional[float]:
    """Computes whole-route ee composition for a multi-step biosynthetic route.

    Formula: ee_overall = ∏(ee_step_i / 100) * 100, taking the best per-step ee
    across all ECs assigned to that step.

    Returns None if any step has no BRENDA-verified ee, or if the route has no steps.

    Known limitation: this assumes step independence. Routes with dynamic kinetic
    resolution (DKR) — where an upstream racemization step combined with a downstream
    enantioselective step produces an artificially-high terminal ee — are underestimated.
    See the route_predictor design spec §7.1.
    """
    steps = route.get("steps") or []
    if not steps:
        return None

    step_ees: list[float] = []
    for step in steps:
        ec_numbers = step.get("ec_numbers", [])
        best_ee = None
        for ec in ec_numbers:
            ec_data = brenda_data.get(ec)
            if not isinstance(ec_data, dict) or ec_data.get("status") != "success":
                continue
            for entry in ec_data.get("entries", []):
                ee_val = entry.get("enantioselectivity")
                if isinstance(ee_val, (int, float)):
                    if best_ee is None or ee_val > best_ee:
                        best_ee = float(ee_val)
        if best_ee is None:
            return None
        step_ees.append(best_ee)

    return math.prod(e / 100.0 for e in step_ees) * 100.0


def score_suggestion(suggestion: dict) -> dict:
    """
    Synthesize all enrichment data for one molecule suggestion into a ranked,
    interpretable enantioselectivity score.

    Degrades gracefully when BRENDA credentials are absent, KEGG data is missing,
    or feasibility could not be computed. Never raises — returns low-confidence
    score with explanatory notes instead.

    Returns a 'scoring' dict with composite_score (0–1), confidence tier,
    top_enzyme breakdown, enzyme_rankings list, stereo_confirmed flag,
    feasibility_flux, and human-readable scoring_notes for researchers.
    """
    notes: list[str] = []
    query_smiles = suggestion.get("SMILES")

    # --- Stereo confirmed (stricter than chirality_validation.valid alone) ---
    chir = suggestion.get("chirality_validation", {})
    centers = chir.get("chiral_centers", [])
    stereo_confirmed = bool(centers) and all(tag in ("R", "S") for _, tag in centers)
    if centers and not stereo_confirmed:
        notes.append(
            "Chiral centers found but stereochemistry unassigned in SMILES — LLM may have omitted @ notation"
        )
    stereo_component = 1.0 if stereo_confirmed else 0.0

    # --- Feasibility component ---
    feas = suggestion.get("feasibility", {})
    feas_status = feas.get("status") if isinstance(feas, dict) else None
    if feas_status == "feasible":
        feas_component: Optional[float] = 1.0
        feasibility_flux: Optional[float] = feas.get("flux")
    elif feas_status == "infeasible":
        feas_component = 0.0
        feasibility_flux = 0.0
    else:
        feas_component = None  # not_in_model / error / missing — excluded from denominator
        feasibility_flux = None
        if feas_status:
            notes.append(f"Feasibility: {feas_status} — excluded from composite score")

    # --- Build enzyme_rankings ---
    brenda_data = suggestion.get("brenda_data", {})
    kegg_data = suggestion.get("kegg_data", {})

    # Collect BRENDA entries that actually have data
    verified_entries: list[dict] = []
    absent_ecs: list[str] = []

    if isinstance(brenda_data, dict) and "status" not in brenda_data:
        # Normal case: keyed by EC number
        for ec, ec_result in brenda_data.items():
            if not isinstance(ec_result, dict):
                continue
            if ec_result.get("status") == "success" and ec_result.get("entries"):
                for entry in ec_result["entries"]:
                    verified_entries.append({**entry, "ec_number": ec})
            else:
                absent_ecs.append(f"EC {ec} (status: {ec_result.get('status', 'unknown')})")
    else:
        # Top-level status dict like {"status": "no_ec_numbers"} or {"status": "no_credentials"}
        status_val = brenda_data.get("status", "unknown") if isinstance(brenda_data, dict) else "unknown"
        absent_ecs.append(f"all ECs (brenda_data status: {status_val})")

    if verified_entries:
        # Case A: at least one BRENDA-verified entry
        enzyme_rankings = sorted(
            [_score_enzyme_entry(e, query_smiles) for e in verified_entries],
            key=lambda r: r["score"],
            reverse=True,
        )
        for note in absent_ecs:
            notes.append(f"No BRENDA data for {note}")
        has_brenda = True
    else:
        # Case B: no BRENDA data — synthetic entries from LLM known_ee
        has_brenda = False
        notes.append("No BRENDA credentials or data available; ee sourced from LLM claim (unverified)")

        llm_ee_value = _parse_ee_percentage(suggestion.get("known_ee"))
        # Discount by 0.6: an unverified LLM-claimed 98% ee must score below a
        # BRENDA-verified 98% ee — otherwise the composite is indistinguishable
        # and researchers sorting by score will treat hallucinations as facts.
        llm_ee_norm = min(llm_ee_value / 100.0, 1.0) * 0.6 if llm_ee_value is not None else 0.0

        # EC numbers from KEGG; fall back to the LLM enzyme_class string
        ec_list: list[str] = []
        if isinstance(kegg_data, dict) and kegg_data.get("status") == "success":
            ec_list = kegg_data.get("enzymes", [])
        if not ec_list:
            fallback_class = suggestion.get("enzyme_class", "unknown")
            ec_list = [fallback_class]

        enzyme_rankings = [
            {
                "ec_number": ec,
                "score": round(llm_ee_norm, 4),
                "ee_value": llm_ee_value,
                "ee_source": "llm_claim",
                "tanimoto_similarity": None,
                "organism": None,
                "engineered_variants": None,
            }
            for ec in ec_list
        ]

    # --- Composite score ---
    best_enzyme_score = max((r["score"] for r in enzyme_rankings), default=0.0)

    numerator   = W_EE_TAN * best_enzyme_score + W_STEREO * stereo_component
    denominator = W_EE_TAN + W_STEREO

    if feas_component is not None:
        numerator   += W_FEAS * feas_component
        denominator += W_FEAS

    composite_score = round(numerator / denominator, 4)

    # --- Confidence tier ---
    has_tanimoto = any(r.get("tanimoto_similarity") is not None for r in enzyme_rankings)
    if has_brenda and has_tanimoto:
        confidence = "high"
    elif has_brenda or (not has_brenda and feas_component == 1.0):
        confidence = "medium"
    else:
        confidence = "low"

    # --- Top enzyme ---
    if enzyme_rankings:
        best = enzyme_rankings[0]
        top_enzyme = {
            "ec_number": best["ec_number"],
            "ee_value": best["ee_value"],
            "ee_source": best["ee_source"],
            "tanimoto_similarity": best["tanimoto_similarity"],
            "organism": best["organism"],
        }
    else:
        top_enzyme = {
            "ec_number": None,
            "ee_value": None,
            "ee_source": "unknown",
            "tanimoto_similarity": None,
            "organism": None,
        }
        notes.append("No enzyme candidates found in KEGG or BRENDA data")

    # Whole-route ee composition (per-route, attached back into the route_prediction structure).
    # See route_predictor spec §7.1 for the multiplicative formula and DKR limitation.
    route_pred = suggestion.get("route_prediction") or {}
    routes = route_pred.get("routes") or []
    for route in routes:
        composed = _compose_route_ee(route, brenda_data)
        route["composed_ee"] = composed

    return {
        "composite_score": composite_score,
        "confidence": confidence,
        "top_enzyme": top_enzyme,
        "enzyme_rankings": enzyme_rankings,
        "stereo_confirmed": stereo_confirmed,
        "feasibility_flux": feasibility_flux,
        "scoring_notes": notes,
    }
