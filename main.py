from ChiraLLM.query_handler import ask_gpt_chirality
from ChiraLLM.database_validator import query_kegg
from ChiraLLM.chirality_checker import validate_chirality
from utils.file_saver import save_suggestions_to_csv

def main():
    print("Welcome to ChiraLLM (Discovery Engine of ChiralAI)!")
    query = input("Enter your query (e.g., 'suggest a biodegradable polymer precursor'): ")
    print(f"Processing query: {query}")

    # Step 1: Query GPT
    response = ask_gpt_chirality(query)
    print(f"Response from GPT: {response}")

    # Step 2: Parse and validate suggestions
    suggestions = [
        {
            "name": "L-glutamate",
            "SMILES": "C(CC(=O)O)C(N)=O",
            "KEGG_ID": "C00025",
            "applications": "Used in agriculture as a nitrogen source."
        }  # Example suggestion; replace with parsed response
    ]

    for suggestion in suggestions:
        smiles = suggestion.get("SMILES")
        if smiles:
            suggestion["chirality_validation"] = validate_chirality(smiles)

        compound_id = suggestion.get("KEGG_ID")
        if compound_id:
            suggestion["kegg_data"] = query_kegg(compound_id)

    # Step 3: Save and display results
    print(f"Processed suggestions: {suggestions}")
    filename = save_suggestions_to_csv(suggestions)
    print(f"Results saved to {filename}")

if __name__ == "__main__":
    main()
