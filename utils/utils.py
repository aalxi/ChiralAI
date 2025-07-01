import requests
from rdkit import Chem
from rdkit.Chem import Descriptors

# --- KEGG API Constants ---
KEGG_REST_URL = "http://rest.kegg.jp"

def get_compound_properties(smiles: str) -> dict:
    """
    Analyzes a SMILES string to get basic molecular properties using RDKit.
    Returns a dictionary with properties or an error.
    """
    if not smiles or not isinstance(smiles, str):
        return {"error": "Invalid or empty SMILES string"}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "RDKit could not parse SMILES string"}
    
    # Calculate properties
    molecular_weight = Descriptors.MolWt(mol)
    num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    
    return {
        "molecular_weight": round(molecular_weight, 2),
        "num_chiral_centers": num_chiral_centers,
        "is_chiral": num_chiral_centers > 0
    }

def find_kegg_compound_id(compound_name: str) -> str | None:
    """Finds the KEGG ID for a given compound name."""
    try:
        response = requests.get(f"{KEGG_REST_URL}/find/compound/{compound_name}")
        response.raise_for_status()  # Raises an exception for bad status codes
        if response.text.strip():
            # The first line, first entry is usually the best match
            return response.text.strip().split('\n')[0].split('\t')[0]
    except requests.RequestException:
        return None
    return None

def get_kegg_pathways_for_compound(compound_id: str) -> list[str]:
    """Gets all KEGG pathway maps associated with a compound ID."""
    if not compound_id:
        return []
    try:
        response = requests.get(f"{KEGG_REST_URL}/link/pathway/{compound_id}")
        response.raise_for_status()
        if response.text.strip():
            # Return a list of pathway IDs (e.g., 'path:eco00260')
            return [line.split('\t')[1] for line in response.text.strip().split('\n')]
    except requests.RequestException:
        return []
    return []