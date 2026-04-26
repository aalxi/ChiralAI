import os
import re
import hashlib
import logging
from dotenv import load_dotenv

load_dotenv()

logger = logging.getLogger(__name__)

BRENDA_WSDL = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"

# Matches: "99% ee", "98% R-enantiomeric excess", "ee 95%", ">99% ee", "99.5% ee", "ee: 97%"
_EE_PATTERN = re.compile(
    r'(\d+(?:\.\d+)?)\s*%\s*(?:[RS]-)?enantiomeric excess'   # "99% R-enantiomeric excess"
    r'|ee[:\s]*>?\s*(\d+(?:\.\d+)?)\s*%'                     # "ee 99%", "ee: >99%"
    r'|(\d+(?:\.\d+)?)\s*%\s*ee'                             # "99% ee"
    r'|enantiomeric excess\D{0,40}?(\d+(?:\.\d+)?)\s*%'      # "enantiomeric excess ... 98%"
    r'|(\d+(?:\.\d+)?)\s*%\D{0,25}?enantiomeric excess',     # "98% ... enantiomeric excess"
    re.IGNORECASE,
)
# Matches "(R)-selective", "(S)-selective", "R-selective", "99% R-enantiomeric excess", etc.
_STEREO_PATTERN = re.compile(
    r'\(([RS])\)[- ]selective'          # "(R)-selective"
    r'|\b([RS])[- ]selective'           # "R-selective"
    r'|\b([RS])[- ]enantiomeric'        # "R-enantiomeric excess"
    r'|\b([RS])[- ]enantio',            # "R-enantio..."
    re.IGNORECASE,
)


def _auth():
    email = os.getenv("BRENDA_EMAIL")
    password = os.getenv("BRENDA_PASSWORD")
    if not email or not password:
        raise EnvironmentError("BRENDA_EMAIL and BRENDA_PASSWORD must be set in .env")
    pw_hash = hashlib.sha256(password.encode("utf-8")).hexdigest()
    return email, pw_hash


def _get_client():
    from zeep import Client, Settings
    # strict=False: BRENDA uses legacy SOAP-ENC arrayType attributes that
    # strict zeep rejects; permissive mode is required for correct parsing.
    return Client(BRENDA_WSDL, settings=Settings(strict=False))


def _extract_ee_from_commentary(commentary: str) -> tuple[float | None, str | None]:
    """
    Parse ee% and enantiopreference (R/S) from a BRENDA commentary string.

    BRENDA stores enantioselectivity as free text in commentarySubstrates /
    commentaryProducts fields — there is no structured ee field in the SOAP schema.
    Returns (ee_value_float, stereo_string) where stereo_string is 'R', 'S', or None.
    """
    if not commentary:
        return None, None
    m = _EE_PATTERN.search(commentary)
    ee_value = None
    if m:
        # One of the five capture groups will be non-None
        raw = next((g for g in m.groups() if g is not None), None)
        if raw is not None:
            ee_value = float(raw)

    stereo = None
    sm = _STEREO_PATTERN.search(commentary)
    if sm:
        raw_stereo = next((g for g in sm.groups() if g is not None), None)
        if raw_stereo:
            stereo = raw_stereo.upper()

    return ee_value, stereo


def query_enantioselectivity(ec_number: str, organism: str = "") -> dict:
    """
    Query BRENDA for enantioselectivity data by EC number.

    Uses getSubstratesProducts — the only BRENDA SOAP method that contains
    ee% data (in commentary fields). There is no standalone getEnantioselectivity
    method in the BRENDA WSDL.

    BRENDA's SOAP API requires positional param-string arguments in the form
    "paramName*filterValue" — keyword arguments silently return 0 results.

    Returns a dict with:
      status: 'success' | 'no_data' | 'no_credentials' | 'error'
      entries: list of dicts with substrate, ee_value, stereo, organism, commentary
    """
    try:
        client = _get_client()
        email, pw_hash = _auth()

        # Positional param-string format (BRENDA SOAP convention):
        # "paramName*filterValue" — empty string after * means no filter.
        results = client.service.getSubstratesProducts(
            email, pw_hash,
            f"ecNumber*{ec_number}",
            f"organism*{organism}",
            "substrates*",
            "commentarySubstrates*",
            "literatureSubstrates*",
            "organismSubstrates*",
            "products*",
            "commentaryProducts*",
            "literatureProducts*",
            "organismProducts*",
            "reversibility*",
        )

        if not results:
            return {"status": "no_data", "ec_number": ec_number, "entries": []}

        entries = []
        for r in results:
            comm_sub = getattr(r, "commentarySubstrates", None) or ""
            comm_prod = getattr(r, "commentaryProducts", None) or ""
            combined = comm_sub + " " + comm_prod

            ee_value, stereo = _extract_ee_from_commentary(combined)
            # Only include entries that carry ee% information — substrate records
            # without commentary are substrate coverage data, not enantioselectivity data.
            if ee_value is None:
                continue

            entries.append({
                "ec_number": ec_number,
                "substrate": getattr(r, "substrates", None),
                "enantioselectivity": ee_value,
                "stereo": stereo,
                "organism": getattr(r, "organism", None),
                "commentary": combined.strip() or None,
            })

        if not entries:
            return {"status": "no_data", "ec_number": ec_number, "entries": []}

        return {"status": "success", "ec_number": ec_number, "entries": entries}

    except EnvironmentError as e:
        return {"status": "no_credentials", "message": str(e)}
    except Exception as e:
        logger.warning("BRENDA query failed for EC %s: %s", ec_number, e)
        return {"status": "error", "ec_number": ec_number, "message": str(e)}


def query_enantioselectivity_batch(ec_numbers: list[str], organism: str = "") -> dict:
    """
    Query BRENDA for a list of EC numbers. Returns dict keyed by EC number.
    Used to enrich all enzyme candidates from a KEGG lookup in one call.
    """
    return {ec: query_enantioselectivity(ec, organism) for ec in ec_numbers}
