import cobra
import cobra.io

# BiGG metabolite IDs for cofactors worth flagging to the researcher
COFACTORS = {
    "nadph_c":   "NADPH",
    "nadp_c":    "NADP+",
    "nadh_c":    "NADH",
    "nad_c":     "NAD+",
    "pydx5p_c":  "PLP (pyridoxal-5'-phosphate)",
    "fad_c":     "FAD",
    "fadh2_c":   "FADH2",
    "coa_c":     "CoA",
    "atp_c":     "ATP",
    "fmn_c":     "FMN",
    "thf_c":     "THF",
}

_model_cache: dict = {}


def _load_model(model_id: str):
    if model_id not in _model_cache:
        _model_cache[model_id] = cobra.io.load_model(model_id)
    return _model_cache[model_id]


def _find_metabolite(model, kegg_id: str):
    """
    Search iJO1366 metabolites for a matching KEGG compound annotation.
    iJO1366 stores KEGG IDs under the 'kegg.compound' annotation key,
    either as a string or list of strings.
    """
    matches = []
    for met in model.metabolites:
        refs = met.annotation.get("kegg.compound", [])
        if isinstance(refs, str):
            refs = [refs]
        if kegg_id in refs:
            matches.append(met)
    return matches


def check_feasibility(kegg_compound_id: str, model_id: str = "iJO1366") -> dict:
    """
    Run FBA on iJO1366 (or user-specified model) to assess whether a target
    compound is producible in E. coli and which cofactors its biosynthesis requires.

    Returns a dict with:
      status         — 'feasible' | 'infeasible' | 'not_in_model' | 'error'
      flux           — optimal production flux in mmol/gDW/h
      cofactors      — list of cofactors involved in producing reactions
      metabolite_id  — BiGG metabolite ID matched from KEGG annotation
    """
    try:
        model = _load_model(model_id)
    except Exception as e:
        return {"status": "error", "message": f"Failed to load {model_id}: {e}"}

    matches = _find_metabolite(model, kegg_compound_id)
    if not matches:
        return {
            "status": "not_in_model",
            "kegg_id": kegg_compound_id,
            "model": model_id,
            "message": f"{kegg_compound_id} has no matching metabolite in {model_id} — "
                       "compound may be outside E. coli metabolism or use a non-standard ID",
        }

    # Prefer cytoplasmic form; fall back to first match
    cytoplasmic = [m for m in matches if m.id.endswith("_c")]
    target = cytoplasmic[0] if cytoplasmic else matches[0]

    # Use context manager so demand reaction + objective change are reverted on exit
    with model:
        demand = cobra.Reaction("DM_CHIRALAI_TARGET")
        demand.lower_bound = 0
        demand.upper_bound = 1000.0
        demand.add_metabolites({target: -1.0})
        model.add_reactions([demand])
        model.objective = demand

        solution = model.optimize()

    # Cofactors: scan all reactions that stoichiometrically produce the target
    producing = [
        rxn for rxn in target.reactions
        if rxn.metabolites.get(target, 0) > 0
    ]
    cofactors_found = sorted({
        COFACTORS[met.id]
        for rxn in producing
        for met in rxn.metabolites
        if met.id in COFACTORS
    })

    base = {
        "kegg_id": kegg_compound_id,
        "metabolite_id": target.id,
        "metabolite_name": target.name,
        "model": model_id,
        "cofactors_required": cofactors_found,
        "n_producing_reactions": len(producing),
    }

    if solution.status == "optimal" and solution.objective_value > 1e-6:
        return {**base, "status": "feasible", "flux": round(solution.objective_value, 4)}
    else:
        return {**base, "status": "infeasible", "flux": 0.0, "solver_status": solution.status}
