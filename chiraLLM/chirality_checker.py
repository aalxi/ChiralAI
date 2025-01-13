from rdkit import Chem

def validate_chirality(smiles):
    """
    Checks if the molecule represented by the given SMILES string is chiral.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return {"valid": False, "error": "Invalid SMILES"}
    
    chiral_centers = Chem.FindMolChiralCenters(molecule, includeUnassigned=True)
    return {"valid": len(chiral_centers) > 0, "chiral_centers": chiral_centers}
