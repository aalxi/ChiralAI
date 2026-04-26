import json
import sys
from ChiraLLM.query_handler import ask_gpt_chirality
from ChiraLLM.database_validator import query_kegg
from ChiraLLM.chirality_checker import validate_chirality
from ChiraLLM.brenda_client import query_enantioselectivity_batch
from ChiraLLM.feasibility_checker import check_feasibility
from ChiraLLM.enantioselectivity_scorer import score_suggestion
from utils.file_saver import save_suggestions_to_csv

def main():
    print("Welcome to ChiraLLM (Discovery Engine of ChiralAI)!")
    query = input("Enter your query (e.g., 'suggest a biodegradable polymer precursor'): ")
    print(f"Processing query: {query}")

    # Step 1: Query GPT
    response = ask_gpt_chirality(query)
    print("Raw GPT output:", response)
    # Optionally, remove one of the duplicate prints if not needed:
    # print(f"Response from GPT: {response}")

    # Step 2: Parse and validate suggestions
    try:
        parsed_response = json.loads(response)
        if isinstance(parsed_response, dict):
            # New prompt returns {"suggestions": [...]}; fall back to wrapping bare dict
            suggestions = parsed_response.get("suggestions", [parsed_response])
        elif isinstance(parsed_response, list):
            suggestions = parsed_response
        else:
            suggestions = []
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)
        sys.exit(-1)
        suggestions = []

    # Process each suggestion: validate chirality and fetch KEGG data if applicable
    for suggestion in suggestions:
        smiles = suggestion.get("SMILES")
        if smiles:
            suggestion["chirality_validation"] = validate_chirality(smiles)

        compound_id = suggestion.get("KEGG_ID")
        if compound_id:
            kegg = query_kegg(compound_id)
            suggestion["kegg_data"] = kegg
            # Use EC numbers from KEGG to pull enantioselectivity from BRENDA
            ec_numbers = kegg.get("enzymes", []) if kegg.get("status") == "success" else []
            if ec_numbers:
                suggestion["brenda_data"] = query_enantioselectivity_batch(ec_numbers[:5])
            else:
                suggestion["brenda_data"] = {"status": "no_ec_numbers"}
            suggestion["feasibility"] = check_feasibility(compound_id)

        suggestion["scoring"] = score_suggestion(suggestion)

    # Step 3: Save and display results
    # print(f"Processed suggestions: {suggestions}")
    filename = save_suggestions_to_csv(suggestions)
    print(f"Results saved to {filename}")

if __name__ == "__main__":
    main()
