import rdkit
from rdkit import Chem
import requests
import openai
import os
import re
import json
import pandas as pd
from dotenv import load_dotenv
from datetime import datetime

class ChiralAnalyzer:
    def __init__(self):
        load_dotenv()
        openai.api_key  = os.getenv('OPENAI_API_KEY')

    def validate_chirality(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return {"valid": False, "error": "Invalid SMILES"}
        
        chiral_centers = Chem.FindMolChiralCenters(molecule, includeUnassigned=True)
        return {"valid": len(chiral_centers) > 0, "chiral_centers": chiral_centers}
    
    def query_kegg(self, compound_id):
        url = f"http://rest.kegg.jp/get/{compound_id}"
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        else:
            return {"error": f"Failed to retrieve data for {compound_id}"}
    def ask_gpt_chirality(self, query):
        system_prompt = """You are an expert in chemistry. For each  molecule you suggest, provide:
        1. Molecule name
        2. SMILES notation
        3. KEGG ID (if known)
        4. Brief description of applications
        Format your response  as JSON with key: name, SMILES, KEGG_ID, applications"""

        try:
            response = openai.chat.completions.create(
                model="gpt-4",
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": query}
                ],
                temperature=0.7
            )
            return response.choices[0].message.content
        except Exception as e:
            return f"Error in GPT query: {str(e)}"
    
    def parse_response(self, response):
        try:
            if isinstance(response, str):
                json_match =  re.search(r'\{.*\}', response.replace('\n', ''), re.DOTALL)
                if  json_match:
                    return [json.loads(json_match.group())]
            return json.loads(response)
        except:
            suggestions = []
            lines = response.split('\n')
            current_molecule = {}

            for line in lines:
                if 'name:' in line.lower():
                    if current_molecule:
                        suggestions.append(current_molecule)
                    current_molecule = {}
                    current_molecule['name'] = line.split(':', 1)[1].strip()
                elif 'smiles:' in  line.lower():
                    current_molecule['SMILES'] = line.split(':', 1)[1].strip()
                elif 'kegg' in line.lower():
                    current_molecule['KEGG_ID'] = line.split(':', 1)[1].strip()
                elif 'applications' in line.lower():
                    current_molecule['applications'] = line.split(':', 1)[1].strip()
            if current_molecule:
                suggestions.append(current_molecule)

            return suggestions

    def process_query(self, query):
        response_from_gpt = self.ask_gpt_chirality(query)
        print(f"Response from gpt {response_from_gpt}")
        suggestions = self.parse_response(response_from_gpt)

        for suggestion in suggestions:
            smiles = suggestion.get("SMILES")
            if smiles:
                chirality_data = self.validate_chirality(smiles)
                suggestion["chirality_validation"] = chirality_data
        
        for suggestion  in suggestions:
            compound_id = suggestion.get("KEGG_ID")
            if compound_id:
                kegg_data = self.query_kegg(compound_id)
                suggestion["kegg_data"] = kegg_data

        return  suggestions
    
    def save_suggestions_to_csv(self, suggestions):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"suggestions_{timestamp}.csv"

        flattened_data = []
        for  suggestion in suggestions:
            flat_dict =  {}
            for key, value in suggestion.items():
                if isinstance(value, dict):
                    for sub_key, sub_value in value.items():
                        flat_dict[f"{key}_{sub_key}"] = str(sub_value)
                else:
                    flat_dict[key] = str(value)
            flattened_data.append(flat_dict)
        
        df = pd.DataFrame(flattened_data)
        df.to_csv(filename, index=False)
        return filename



def main():
    analyzer = ChiralAnalyzer()
    # test  validate chirality
    smiles = "CC[C@H](N)C(=O)O"  # L-alanine
    # result = analyzer.validate_chirality(smiles)
    # print(result)

    # test  query  kegg
    compound_id = "C00022"  # L-glutamate
    # kegg_data = analyzer.query_kegg(compound_id)
    # print(kegg_data)

    query = "suggest a biodegradable polymer precursor"
    print(f"Processing query: {query}")

    results  = analyzer.process_query(query)
    print(f"Results: {results}")

    filename = analyzer.save_suggestions_to_csv(results)

    print(f"\n Results saved to: {filename}")




main()