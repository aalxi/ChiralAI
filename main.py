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
    print("Raw GPT output:", response)
    print(f"Response from GPT: {response}")

    # Step 2: Parse and validate suggestions
    # Don't offer a fallback suggestion - the ChatGPT call is failing, and we didn't realize that at first.
import json

# Step 2: Parse and validate suggestions
try:
    parsed_response = json.loads(response)
    # If GPT returns a single dictionary, put it in a list so we can iterate.
    if isinstance(parsed_response, dict):
        suggestions = [parsed_response]
    # If GPT returns a list of dicts, just use it directly.
    elif isinstance(parsed_response, list):
        suggestions = parsed_response
    else:
        # If it's neither a dict nor a list, fallback to empty list.
        suggestions = []
except Exception as e:
    print(f"Error parsing GPT response as JSON: {e}")
    suggestions = []


    for suggestion in suggestions:
        smiles = suggestion.get("SMILES")
        if smiles:
            suggestion["chirality_validation"] = validate_chirality(smiles)

        compound_id = suggestion.get("KEGG_ID")
        if compound_id:
            suggestion["kegg_data"] = query_kegg(compound_id) # Formatting for this is coming back as a KEY (spaces) VALUE\n KEY (spaces) VALUE \n thing, and we're looking for CSVs.
            print("kegg_data = " + str(suggestion["kegg_data"]))

    # Step 3: Save and display results
    print(f"Processed suggestions: {suggestions}")
    filename = save_suggestions_to_csv(suggestions)
    print(f"Results saved to {filename}")

if __name__ == "__main__":
    main()
