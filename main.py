import sys
import json
import os

# Import our organized modules
from ChiraLLM.query_handler import ask_gpt_chirality
from ChiraLLM.feasibility import analyze_feasibility
from utils.file_saver import save_data_to_csv
from utils.display import display_results_in_terminal # <-- Import our new display function

def main():
    """The main execution pipeline for the ChiralAI MVP."""
    print("Welcome to the ChiralAI MVP (Modules 1 & 2)")
    print("-" * 40)

    # --- MODULE 1: DISCOVERY ENGINE ---
    query = input("Enter your natural-language query (e.g., 'a chiral precursor for a biodegradable polymer'):\n> ")
    if not query:
        print("Query cannot be empty. Exiting.")
        sys.exit(1)

    print("\n🔬 Querying Discovery Engine (LLM)...")
    llm_output_str = ask_gpt_chirality(query)
    
    try:
        suggestions = json.loads(llm_output_str)
        if "error" in suggestions:
            print(f"❌ Error from LLM: {suggestions['error']}")
            sys.exit(1)
        if isinstance(suggestions, dict):
            suggestions = [suggestions]

    except json.JSONDecodeError:
        print("❌ Critical Error: Failed to parse the LLM's response.")
        sys.exit(1)

    # --- FIX: Standardize all keys from the LLM to lowercase ---
    suggestions = [{k.lower(): v for k, v in s.items()} for s in suggestions]

    raw_suggestions_path = save_data_to_csv(suggestions, "llm_suggestions")
    print(f"✅ Raw LLM suggestions saved to: {raw_suggestions_path}")


    # --- MODULE 2: FEASIBILITY FILTER ---
    host = input("\nEnter host for feasibility check (e.g., 'E. coli' or 'yeast'):\n> ").strip()
    if not host:
        print("Host cannot be empty. Exiting.")
        sys.exit(1)
        
    print(f"\n⚙️  Running Feasibility Analysis for host '{host}'...")
    
    all_results = []
    for suggestion in suggestions:
        # Now this check will work correctly with lowercase 'smiles'
        if suggestion.get("name") and suggestion.get("smiles"):
            result = analyze_feasibility(suggestion["name"], suggestion["smiles"], host)
            all_results.append(result)
            # We can remove the per-line print to keep the log clean
        else:
            # FIX: Truncate long values in the error message to keep the terminal clean
            truncated_suggestion = {k: (str(v)[:70] + '...' if isinstance(v, str) and len(str(v)) > 70 else v) for k, v in suggestion.items()}
            print(f"  - Skipping invalid entry: {truncated_suggestion}")

    # --- NEW: Display the formatted report in the terminal ---
    display_results_in_terminal(all_results)

    # Save the final, annotated results to a file as before
    if all_results:
        host_str = host.lower().replace(" ", "_")
        final_results_path = save_data_to_csv(all_results, f"feasibility_analysis_{host_str}")
    else:
        print("\nNo valid molecules were processed to save.")

if __name__ == "__main__":
    main()