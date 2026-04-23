import os
import hashlib
from dotenv import load_dotenv

load_dotenv()

BRENDA_WSDL = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"

def _auth():
    email = os.getenv("BRENDA_EMAIL")
    password = os.getenv("BRENDA_PASSWORD")
    if not email or not password:
        raise EnvironmentError("BRENDA_EMAIL and BRENDA_PASSWORD must be set in .env")
    pw_hash = hashlib.sha256(password.encode("utf-8")).hexdigest()
    return email, pw_hash

def _get_client():
    from zeep import Client
    return Client(BRENDA_WSDL)

def query_enantioselectivity(ec_number, organism=""):
    """
    Query BRENDA for enantioselectivity (ee%) data by EC number.

    Returns a dict with status and a list of entries, each containing
    substrate name, ee value, organism, and commentary. If BRENDA has
    no data for this EC number, status is 'no_data' with an empty list —
    callers should surface this to the researcher rather than silently skipping.
    """
    try:
        client = _get_client()
        email, pw_hash = _auth()

        param_str = f"ecNumber*{ec_number}#organism*{organism}#"
        results = client.service.getEnantioselectivity(email, pw_hash, param_str)

        if not results:
            return {"status": "no_data", "ec_number": ec_number, "entries": []}

        entries = [
            {
                "ec_number": ec_number,
                "substrate": getattr(r, "substrate", None),
                "enantioselectivity": getattr(r, "enantioselectivity", None),
                "organism": getattr(r, "organism", None),
                "commentary": getattr(r, "commentary", None),
                "literature": getattr(r, "literature", None),
            }
            for r in results
        ]
        return {"status": "success", "ec_number": ec_number, "entries": entries}

    except EnvironmentError as e:
        return {"status": "no_credentials", "message": str(e)}
    except Exception as e:
        return {"status": "error", "ec_number": ec_number, "message": str(e)}


def query_enantioselectivity_batch(ec_numbers, organism=""):
    """
    Query BRENDA for a list of EC numbers. Returns dict keyed by EC number.
    Used to enrich all enzyme candidates from a KEGG lookup in one call.
    """
    return {ec: query_enantioselectivity(ec, organism) for ec in ec_numbers}
