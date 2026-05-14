import argparse
import json
import sys
from dataclasses import asdict
from ChiraLLM.query_handler import ask_gpt_chirality
from ChiraLLM.database_validator import query_kegg
from ChiraLLM.chirality_checker import validate_chirality
from ChiraLLM.brenda_client import query_enantioselectivity_batch
from ChiraLLM.feasibility_checker import check_feasibility
from ChiraLLM.enantioselectivity_scorer import score_suggestion
from ChiraLLM.route_predictor import predict_route
from utils.file_saver import save_suggestions_to_csv

def main():
    parser = argparse.ArgumentParser(description="ChiraLLM — AI-guided chiral molecule discovery")
    parser.add_argument("--query", "-q", type=str, default=None,
                        help="Discovery query (if omitted, prompts interactively)")
    parser.add_argument("--out-dir", type=str, default=".",
                        help="Directory for output CSV/JSON files (default: current dir)")
    args = parser.parse_args()

    print("Welcome to ChiraLLM (Discovery Engine of ChiralAI)!")
    if args.query:
        query = args.query
    else:
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

            # NEW: route prediction runs after KEGG validation, before BRENDA enrichment.
            # The full set of ECs across all routes' steps becomes the BRENDA query input.
            route_result = predict_route(compound_id, mode="top_n", n=3)
            suggestion["route_prediction"] = asdict(route_result)

            if route_result.status == "success":
                # Deduplicate ECs across all routes' steps before the BRENDA batch call.
                all_ecs = sorted({
                    ec
                    for route in route_result.routes
                    for step in route.steps
                    for ec in step.ec_numbers
                })
                if all_ecs:
                    suggestion["brenda_data"] = query_enantioselectivity_batch(all_ecs[:20])
                else:
                    suggestion["brenda_data"] = {"status": "no_ec_numbers"}

                # Per-route terminal-precursor feasibility (call site change per spec §2.1).
                suggestion["route_feasibility"] = [
                    {
                        "route_index": i,
                        "terminal_precursor": route.terminal_precursor_id,
                        "feasibility": check_feasibility(route.terminal_precursor_id),
                    }
                    for i, route in enumerate(route_result.routes)
                ]
            else:
                # Fall back to legacy behavior: BRENDA on KEGG's enzymes for the target only.
                ec_numbers = kegg.get("enzymes", []) if kegg.get("status") == "success" else []
                if ec_numbers:
                    suggestion["brenda_data"] = query_enantioselectivity_batch(ec_numbers[:5])
                else:
                    suggestion["brenda_data"] = {"status": "no_ec_numbers"}
                suggestion["route_feasibility"] = [{
                    "route_index": 0,
                    "terminal_precursor": compound_id,
                    "feasibility": check_feasibility(compound_id),
                }]

        suggestion["scoring"] = score_suggestion(suggestion)

    # Step 3: Save and display results
    csv_file, json_file = save_suggestions_to_csv(suggestions, out_dir=args.out_dir)
    print(f"Results saved to {csv_file} and {json_file}")

if __name__ == "__main__":
    main()
