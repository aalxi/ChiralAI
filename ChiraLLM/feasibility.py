import requests, csv
from rdkit import Chem
from rdkit.Chem import Descriptors

# ――― KEGG constants ――― #
KEGG_REST = "http://rest.kegg.jp"
ORGANISMS = {
    "E. coli": "eco",
    "yeast":   "sce",
}

# ――― RDKit helper ――― #
def analyze_molecule_rdkit(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}
    mw  = Descriptors.MolWt(mol)
    n_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    feasibility = "Hard to make" if mw > 500 else "Likely feasible"
    return {
        "molecular_weight": round(mw, 2),
        "num_chiral_centers": n_chiral,
        "feasibility_rule": feasibility,
    }

# ――― KEGG helpers ――― #
def find_compound_kegg(name: str) -> str | None:
    r = requests.get(f"{KEGG_REST}/find/compound/{name}")
    if r.ok and r.text.strip():
        return r.text.strip().split('\n')[0].split('\t')[0]
    return None

def get_kegg_pathways(compound_id: str) -> list[str]:
    r = requests.get(f"{KEGG_REST}/link/pathway/{compound_id}")
    if r.ok and r.text.strip():
        return [ln.split('\t')[1] for ln in r.text.strip().split('\n')]
    return []

# ――― public API ――― #
def combined_feasibility_analysis(name: str, smiles: str, host: str) -> dict:
    """
    Returns an enriched feasibility assessment that follows the ChiralAI
    Module‑2 schema.  Uses lightweight heuristics until deeper ML / thermodynamic
    integrations (eQuilibrator, ProGen, etc.) are wired in.

    Schema keys (all present):
        biosynthetic_feasibility : bool
        predicted_hosts          : list[str]
        pathway                  : list[dict]  # at most a KEGG‑ID placeholder per step
        bottlenecks              : list[str]
        confidence               : float 0‑1
        references               : list[str]
        next_steps               : list[str]

    The function still returns a superset with the legacy keys so nothing else
    downstream breaks.
    """
    # 1) RDKit phys‑chem quick check
    rdkit_out = analyze_molecule_rdkit(smiles)
    if "error" in rdkit_out:
        base_result = {"error": f"RDKit failed for {name}: {rdkit_out['error']}"}
        # keep legacy fields so callers relying on them do not explode
        base_result.update(
            biosynthetic_feasibility=False,
            predicted_hosts=[],
            pathway=[],
            bottlenecks=["Invalid SMILES"],
            confidence=0.0,
            references=[],
            next_steps=["provide valid SMILES"],
        )
        return base_result

    # 2) KEGG lookup
    kegg_id = find_compound_kegg(name)
    pathways = []
    if kegg_id:
        pathways = get_kegg_pathways(kegg_id)

    # 3) Heuristic feasibility decision
    rule_ok      = rdkit_out["feasibility_rule"] == "Likely feasible"
    pathway_ok   = bool(pathways)
    organism_ok  = host in ORGANISMS

    biosyn_feasible = rule_ok and pathway_ok and organism_ok

    # 4) Bottleneck reasoning
    bottlenecks = []
    if not pathway_ok:
        bottlenecks.append("no documented pathway")
    if rdkit_out["num_chiral_centers"] > 5:
        bottlenecks.append("high stereochemical complexity")
    if not organism_ok:
        bottlenecks.append(f"host '{host}' unsupported")

    # 5) Confidence – toy weighting until ML model replaces it
    confidence = 0.2
    if rule_ok:        confidence += 0.3
    if pathway_ok:     confidence += 0.4
    if organism_ok:    confidence += 0.1
    confidence = round(min(confidence, 0.99), 2)

    # 6) Assemble pathway detail placeholder
    pathway_detail = [
        {"step": idx + 1,
         "enzyme": "known",
         "notes": pid}
        for idx, pid in enumerate(pathways)
    ]

    # 7) High‑level next steps
    next_steps = []
    if not pathway_ok:
        next_steps.append("search MetaCyc/BRENDA for missing enzymes")
    if "high stereochemical complexity" in bottlenecks:
        next_steps.append("evaluate chemo‑enzymatic or dynamic‑kinetic routes")
    if not next_steps:
        next_steps.append("design experimental validation in host")

    # 8) Final result (includes legacy keys for backward compatibility)
    return {
        "name":               name,
        "smiles":             smiles,
        "rdkit_analysis":     rdkit_out,
        "kegg_compound_id":   kegg_id,
        "kegg_pathways":      pathways,
        # --- new schema ---
        "biosynthetic_feasibility": biosyn_feasible,
        "predicted_hosts":    [host] if organism_ok else [],
        "pathway":            pathway_detail,
        "bottlenecks":        bottlenecks,
        "confidence":         confidence,
        "references":         pathways if pathways else [],
        "next_steps":         next_steps,
    }

# ――― CSV helper (still optional) ――― #
def process_csv_and_analyze(file_path: str, host: str) -> list[dict] | dict:
    try:
        with open(file_path, newline='') as f:
            rdr = csv.DictReader(f)
            return [
                combined_feasibility_analysis(r["name"], r["smiles"], host)
                if r.get("name") and r.get("smiles")
                else {"error": "Missing name or SMILES"}
                for r in rdr
            ]
    except FileNotFoundError:
        return {"error": f"File {file_path} not found"}
    
