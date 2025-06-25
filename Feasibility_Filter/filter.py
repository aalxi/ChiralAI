import requests
import csv
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
        pathways = []
        for line in r.text.strip().split('\n'):
            pathway_id = line.split('\t')[1]
            pathways.append(pathway_id)
        return pathways
    return []

def combined_feasibility_analysis(name, smiles, host):
    # Analyze molecule with RDKit
    rdkit_analysis = analyze_molecule_rdkit(smiles)
    if "error" in rdkit_analysis:
        return {"error": f"RDKit analysis failed for {name}: {rdkit_analysis['error']}"}
    
    # Find compound in KEGG
    compound_id = find_compound_kegg(name)
    if not compound_id:
        return {"error": f"Compound {name} not found in KEGG"}
    
    # Get KEGG pathways
    host_code = ORGANISMS.get(host)
    if not host_code:
        return {"error": f"Host organism {host} not supported"}
    pathways = get_kegg_pathways(compound_id, host_code)
    
    # Combine results
    return {
        "name": name,
        "smiles": smiles,
        "rdkit_analysis": rdkit_analysis,
        "kegg_compound_id": compound_id,
        "kegg_pathways": pathways
    }

def process_csv_and_analyze(file_path, host):
    results = []
    try:
        with open(file_path, mode='r') as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                name = row.get("name")
                smiles = row.get("smiles")
                if name and smiles:
                    result = combined_feasibility_analysis(name, smiles, host)
                    results.append(result)
                else:
                    results.append({"error": "Missing name or SMILES in CSV row"})
    except FileNotFoundError:
        return {"error": f"File {file_path} not found"}
    return results