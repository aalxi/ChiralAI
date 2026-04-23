import os
from openai import OpenAI
from dotenv import load_dotenv

load_dotenv()

client = OpenAI(api_key=os.getenv('OPENAI_API_KEY'))

def ask_gpt_chirality(query):
    """
    Queries GPT to suggest chiral molecules based on a user's input.
    """
    system_prompt = """You are an expert in biocatalysis and chiral synthesis. Given a user query, suggest exactly 5 ranked chiral molecule candidates that are biosynthetically accessible.

For each candidate you MUST provide:
1. Molecule name with full stereodescriptor (e.g. (S)-ibuprofen, not just ibuprofen)
2. Stereospecific SMILES using @ and @@ to encode defined chirality. If you cannot commit to a single enantiomer, flag it explicitly — do not return achiral SMILES.
3. R/S configuration at each chiral center
4. KEGG compound ID in format C##### (e.g. C00025) — use null if unknown
5. The enantioselective enzyme class that produces this enantiomer (e.g. ketoreductase (KRED), transaminase (TA), lipase, cytochrome P450, epoxide hydrolase, lyase)
6. Known enantiomeric excess (ee%) from BRENDA or literature — use "unknown" if not established
7. Brief note on pharmaceutical or industrial applications

Return your response as a JSON object with a single key "suggestions" containing an array of exactly 5 objects. Each object must have these keys:
  name, SMILES, R_S_config, KEGG_ID, enzyme_class, known_ee, applications

Rank candidates from most to least well-characterized (best BRENDA/literature data first)."""

    try:
        response = client.chat.completions.create(
            model="gpt-4.1",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": query}
            ],
            response_format={"type": "json_object"},
            temperature=0.7
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Error in GPT query: {str(e)}"
