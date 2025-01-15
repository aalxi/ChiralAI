import os

from openai import OpenAI
from openai import ChatCompletion
from openai import Client
from dotenv import load_dotenv

load_dotenv()

client = Client(api_key=os.getenv('OPENAI_API_KEY'))

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
    
    OpenAI.api_key = os.getenv('OPENAI_API_KEY')
    try:
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": query}
            ],
            temperature=0.7
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Error in GPT query: {str(e)}"


