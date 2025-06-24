import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
# KEGG constants
KEGG_REST = "http://rest.kegg.jp"
ORGANISMS = {
    "E. coli": "eco",
    "yeast": "sce"
}
def analyze_molecule_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}
    mw = Descriptors.MolWt(mol)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral_centers = len(chiral_centers)
    feasibility = "Hard to make" if mw > 500 else "Likely feasible"
    return {
        "molecular_weight": round(mw, 2),
        "num_chiral_centers": num_chiral_centers,
        "feasibility_rule": feasibility
    }
def find_compound_kegg(name):
    url = f"{KEGG_REST}/find/compound/{name}"
    r = requests.get(url)
    if r.ok and r.text.strip():
        lines = r.text.strip().split('\n')
        entry = lines[0].split('\t')[0]
        return entry
    return None
def get_kegg_pathways(compound_id, host_code):
    url = f"{KEGG_REST}/link/pathway/{compound_id}"
    r = requests.get(url)
    if r.ok and r.text.strip():
        lines = r.text.strip().split('\n')
        pathways = [line.split('\t')[1] for line in lines]
        return [p for p in pathways if p.startswith(f"path:{host_code}")]
    return []
def combined_feasibility_analysis(name, smiles, host):
    host_code = ORGANISMS.get(host)
    if not host_code:
        return {"error": f"Unknown host: {host}"}
    rdkit_result = analyze_molecule_rdkit(smiles)
    kegg_id = find_compound_kegg(name)
    if not kegg_id:
        kegg_status = "Not found in KEGG"
        pathways = []
    else:
        kegg_status = f"Found as {kegg_id}"
        pathways = get_kegg_pathways(kegg_id, host_code)
    feasibility_kegg = "YES" if pathways else "NO"
    return {
        "input": {"name": name, "host": host},
        "rdkit_analysis": rdkit_result,
        "kegg": {
            "compound_status": kegg_status,
            "feasible_in_host": feasibility_kegg,
            "pathways": pathways
        }
    }
# Example usage:
if __name__ == "__main__":
    name = input("Molecule name (e.g., glucose): ")
    smiles = input("SMILES string: ")
    host = input("Host organism (E. coli or yeast): ")
    result = combined_feasibility_analysis(name, smiles, host)
    from pprint import pprint
    pprint(result)