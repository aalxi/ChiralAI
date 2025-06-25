from utils.file_saver import save_suggestions_to_csv
from ChiraLLM.query_handler import ask_gpt_chirality
from filter import process_csv_and_analyze

if __name__ == "__main__":
    # Step 1: Get user query and generate molecule suggestions
    user_query = input("Enter your query for molecule suggestions: ")
    suggestions = ask_gpt_chirality(user_query)
    
    # Step 2: Save suggestions to a CSV file
    csv_file_path = save_suggestions_to_csv(suggestions)
    print(f"Molecule suggestions saved to {csv_file_path}")
    
    # Step 3: Ask user for host organism
    host = input("Enter the host organism (e.g., 'E. coli', 'yeast'): ")
    
    # Step 4: Process the CSV file and analyze feasibility
    analysis_results = process_csv_and_analyze(csv_file_path, host)
    
    # Step 5: Display the results
    print("Feasibility analysis results:")
    for result in analysis_results:
        print(result)