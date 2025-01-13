import requests

def query_kegg(compound_id):
    """
    Queries KEGG for compound data using the compound ID.
    """
    url = f"http://rest.kegg.jp/get/{compound_id}"
    response = requests.get(url)
    if response.status_code == 200:
        # Parse KEGG data to extract useful information
        data = response.text
        return {"status": "success", "data": data[:500]}  # Trimmed for brevity
    else:
        return {"status": "error", "message": f"Failed to retrieve data for {compound_id}"}
