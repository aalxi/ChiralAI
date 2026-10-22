"""Stage A: ask a language model for candidate chiral molecules.

Two fixes here.

* ``from openai import ChatCompletion`` was a leftover from the pre-1.0 SDK.
  That symbol does not exist in openai>=1.0, so importing this module raised
  ImportError before any code ran.
* The prompt asked for the JSON key ``SMILES`` while the consumer read
  ``smiles``. Every row was silently discarded -- the recorded failure in
  ``suggestions_20250702_155147.csv``. The prompt now specifies lower-case keys
  and a list, and the consumer normalises keys case-insensitively anyway, so
  neither side alone can resurrect the bug.

Nothing this function returns is trusted. A language model produces a correct
SMILES and a degenerate string with equal fluency, so every structure goes
through ``ChiraLLM.structure_gate`` before it is analysed.
"""
from __future__ import annotations
import json, os

from dotenv import load_dotenv
from openai import OpenAI

load_dotenv()

MODEL = os.getenv("CHIRALAI_MODEL", "gpt-4o-mini")

SYSTEM_PROMPT = """You are an expert synthetic chemist. Suggest chiral molecules \
matching the user's request.

Return JSON of the form:
  {"molecules": [{"name": ..., "smiles": ..., "kegg_id": ..., "applications": ...}]}

Rules:
- Use lower-case keys exactly as shown.
- "smiles" must be a valid SMILES string for the named molecule, with
  stereochemistry specified using @ / @@ where the molecule is chiral. A
  structure that does not match the name will be rejected downstream.
- "kegg_id" is a KEGG COMPOUND identifier such as C00186, or "" if unknown.
  Do not invent one; an incorrect id is worse than an absent one.
- Prefer molecules whose stereochemistry matters to their application."""


def _client() -> OpenAI:
    key = os.getenv("OPENAI_API_KEY")
    if not key:
        raise RuntimeError("OPENAI_API_KEY is not set (put it in .env)")
    return OpenAI(api_key=key)


def ask_gpt_chirality(query: str, *, temperature: float = 0.2) -> str:
    """Raw JSON text from the model. Raises rather than returning an error string.

    The old version returned ``f"Error in GPT query: {e}"`` on failure, which the
    caller then tried to parse as JSON -- turning an API error into a confusing
    parse error. Temperature is 0.2, not 0.7: this is structure retrieval from
    the model's memory, not creative writing, and sampling noise here shows up
    as invented SMILES.
    """
    response = _client().chat.completions.create(
        model=MODEL,
        messages=[{"role": "system", "content": SYSTEM_PROMPT},
                  {"role": "user", "content": query}],
        response_format={"type": "json_object"},
        temperature=temperature,
    )
    return response.choices[0].message.content


def parse_suggestions(raw: str) -> list[dict]:
    """Coerce whatever shape came back into a list of dicts."""
    parsed = json.loads(raw)
    if isinstance(parsed, list):
        return parsed
    for key in ("molecules", "results", "suggestions", "data"):
        if isinstance(parsed.get(key), list):
            return parsed[key]
    return [parsed]
