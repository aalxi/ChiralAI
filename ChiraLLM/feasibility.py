import csv
# Import our new centralized utility functions from the utils folder
from utils.utils import get_compound_properties, find_kegg_compound_id, get_kegg_pathways_for_compound

# Supported hosts for the MVP
SUPPORTED_HOSTS = {
    "e. coli": "eco",
    "yeast": "sce"
}

def analyze_feasibility(name: str, smiles: str, host: str) -> dict:
    """
    Performs a lightweight feasibility analysis for a single molecule.
    This is the core of Module 2.
    """
    # 1. Get basic properties from our utils
    properties = get_compound_properties(smiles)
    if "error" in properties:
        # If SMILES is invalid, we can't proceed.
        return {
            "name": name,
            "smiles": smiles,
            "biosynthetic_feasibility": False,
            "confidence": 0.0,
            "bottlenecks": [properties["error"]],
            "next_steps": ["Fix the provided SMILES string."],
            "pathway": [],
            "references": []
        }

    # 2. Check KEGG for pathways
    kegg_id = find_kegg_compound_id(name)
    pathways = get_kegg_pathways_for_compound(kegg_id)

    # 3. Apply MVP Heuristics (Simple Rules)
    is_host_supported = host.lower() in SUPPORTED_HOSTS
    is_pathway_known = bool(pathways)
    is_mw_reasonable = properties.get("molecular_weight", 999) < 500

    # The final feasibility decision for the MVP
    is_feasible = is_host_supported and is_pathway_known and is_mw_reasonable

    # 4. Identify Bottlenecks
    bottlenecks = []
    if not is_host_supported:
        bottlenecks.append(f"Host '{host}' is not a supported chassis in this MVP.")
    if not is_pathway_known:
        bottlenecks.append("No known biosynthetic pathway found in KEGG.")
    if not is_mw_reasonable:
        bottlenecks.append("High molecular weight (>500 Da) may challenge expression/transport.")
    if properties.get("num_chiral_centers", 0) > 4:
        bottlenecks.append("High stereochemical complexity.")

    # 5. Calculate a simple, rule-based confidence score
    confidence = 0.1
    if is_host_supported: confidence += 0.2
    if is_pathway_known:  confidence += 0.4
    if is_mw_reasonable:  confidence += 0.2

    # 6. Suggest Next Steps
    next_steps = []
    if not is_feasible:
        next_steps.append("Re-evaluate with different hosts or consider pathway engineering.")
    else:
        next_steps.append("Proceed to thermodynamic analysis and experimental design.")

    # 7. Assemble the final, structured output
    return {
        "name": name,
        "smiles": smiles,
        "biosynthetic_feasibility": is_feasible,
        "confidence": round(confidence, 2),
        "predicted_host": host if is_host_supported else "N/A",
        "kegg_compound_id": kegg_id or "N/A",
        "bottlenecks": bottlenecks or ["None identified"],
        "next_steps": next_steps,
        "references": pathways,
        "pathway": [
            {"step": i+1, "enzyme_status": "known_pathway", "enzyme_id": path}
            for i, path in enumerate(pathways)
        ]
    }