import openai
import os
from dotenv import load_dotenv

load_dotenv()

def ask_gpt_chirality(query):
    """
    Queries GPT to suggest chiral molecules based on a user's input.
    """
    system_prompt = """You are an expert in chemistry. For each molecule you suggest, provide:
    1. Molecule name
    2. SMILES notation
    3. KEGG ID (if known)
    4. Brief description of applications
    Format your response as JSON with keys: name, SMILES, KEGG_ID, applications."""
    
    openai.api_key = os.getenv('OPENAI_API_KEY')
    try:
        response = openai.ChatCompletion.create(
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
