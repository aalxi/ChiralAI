"""Tier 1 route predictor: best-first backward search through the KEGG reaction
graph from a target compound to curated central metabolites.

See docs/superpowers/specs/2026-04-27-route-predictor-tier1-design.md for the
full design rationale.

Public API: predict_route(compound_id, mode='top_n', n=3, budget=500) -> RouteResult
"""

import functools
import logging
import os
import re
import shutil
import time
from pathlib import Path

import requests
from rdkit import Chem

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# SCIENTIFIC CONTRACT — these two constants are the module's domain commitments.
# A wet-lab reviewer should be able to read and audit them in the first 100 lines.
# ---------------------------------------------------------------------------

CENTRAL_METABOLITES: dict[str, str] = {
    # TCA cycle
    "C00022": "pyruvate",
    "C00024": "acetyl-CoA",
    "C00036": "oxaloacetate",
    "C00149": "(S)-malate",
    "C00122": "fumarate",
    "C00042": "succinate",
    "C00091": "succinyl-CoA",
    "C00026": "alpha-ketoglutarate",
    "C00311": "isocitrate",
    "C00158": "citrate",
    # Glycolysis / PPP
    "C00031": "D-glucose",
    "C00092": "glucose-6-phosphate",
    "C00085": "fructose-6-phosphate",
    "C00354": "fructose-1,6-bisphosphate",
    "C00111": "DHAP",
    "C00118": "G3P",
    "C00197": "3-phosphoglycerate",
    "C00074": "PEP",
    "C00117": "ribose-5-phosphate",
    "C00199": "ribulose-5-phosphate",
    # 20 proteinogenic amino acids
    "C00041": "L-alanine",
    "C00037": "glycine",
    "C00065": "L-serine",
    "C00188": "L-threonine",
    "C00097": "L-cysteine",
    "C00073": "L-methionine",
    "C00407": "L-isoleucine",
    "C00123": "L-leucine",
    "C00183": "L-valine",
    "C00079": "L-phenylalanine",
    "C00082": "L-tyrosine",
    "C00078": "L-tryptophan",
    "C00135": "L-histidine",
    "C00148": "L-proline",
    "C00064": "L-glutamine",
    "C00025": "L-glutamate",
    "C00049": "L-aspartate",
    "C00152": "L-asparagine",
    "C00047": "L-lysine",
    "C00062": "L-arginine",
    # Branched-chain amino acid intermediates (defensibly central — produced from pyruvate
    # via the BCAA biosynthesis pathway; nodes for valine/leucine/isoleucine biosynthesis)
    "C00141": "alpha-ketoisovalerate",
}


# Entries ending in '.' match any EC under that class prefix (e.g., '1.1.1.' matches all KREDs).
# Entries WITHOUT trailing '.' match the EC exactly (e.g., '1.6.99.1' matches only that one EC,
# not '1.6.99.10' or '1.6.99.12'). The dual-mode matching is in _is_industrially_reversible.
INDUSTRIAL_REVERSIBLE_EC_PREFIXES: list[str] = [
    "1.1.1.",      # KREDs / aldo-keto reductases (class)
    "2.6.1.",      # transaminases (class)
    "1.5.1.",      # IREDs (imine reductases) (class)
    "1.6.99.1",    # Old Yellow Enzyme — exact EC, NOT a class prefix
    "3.1.1.",      # lipases (class)
    "1.14.13.22",  # cyclohexanone monooxygenase (BVMO archetype) — exact EC, NOT a class prefix
]


# ---------------------------------------------------------------------------
# Search defaults — tunable, but with sane starting values.
# ---------------------------------------------------------------------------

DEFAULT_BUDGET = 500
DEFAULT_DEPTH_CAP = 8
DEFAULT_MAX_ROUTES = 3
THERMO_PENALTY_PER_KJ = 0.05
FALLBACK_DELTA_G_KJ = 5.0

# v1 STARTING GUESS — NOT VALIDATED. h(n) Tanimoto distance is scaled by this weight before
# being added to g(n) accumulated edge cost. Weight=2.0 means structural similarity to a
# central metabolite matters 2x as much as accumulated thermo+directional cost (g(n) is on
# the order of ~1 per step). Re-tune against integration fixtures.
# Re-tuning trigger: if top-3 routes for (R)-pantolactone (C00599) do not include the
# KIV-via-ketopantoate-hydroxymethyltransferase route, this constant is too high or too low.
TANIMOTO_HEURISTIC_WEIGHT = 2.0


# ---------------------------------------------------------------------------
# Direction constants for reaction equations
# ---------------------------------------------------------------------------

DIRECTION_REVERSIBLE = "reversible"
DIRECTION_FORWARD_ONLY = "forward_only"

_COMPOUND_TOKEN_RE = re.compile(r"^(?:(\d+)\s+)?(C\d{5})$")


# ---------------------------------------------------------------------------
# Reaction equation parsing
# ---------------------------------------------------------------------------


def _parse_reaction_equation(equation: str) -> tuple[list[tuple[int, str]], list[tuple[int, str]], str]:
    """Splits a KEGG reaction equation into substrates, products, and direction.

    KEGG equation format: 'C00033 + C00010 <=> C00024 + C00011' (reversible),
    'C00033 => C00024' (irreversible). Coefficients written as '2 C00006'.

    Returns (substrates, products, direction). direction is 'reversible' or 'forward_only'.
    Raises ValueError if equation has no arrow or contains an unparseable token.
    """
    if "<=>" in equation:
        direction = DIRECTION_REVERSIBLE
        sides = equation.split("<=>", 1)
    elif "=>" in equation:
        direction = DIRECTION_FORWARD_ONLY
        sides = equation.split("=>", 1)
    else:
        raise ValueError(f"No reaction arrow in equation: {equation!r}")

    if len(sides) != 2:
        raise ValueError(f"Could not split equation into two sides: {equation!r}")

    def _parse_side(side: str) -> list[tuple[int, str]]:
        results = []
        for token in side.split("+"):
            token = token.strip()
            if not token:
                continue
            m = _COMPOUND_TOKEN_RE.match(token)
            if m is None:
                raise ValueError(f"Unparseable token {token!r} in side {side!r}")
            coef_str, cid = m.groups()
            coef = int(coef_str) if coef_str else 1
            results.append((coef, cid))
        return results

    return _parse_side(sides[0]), _parse_side(sides[1]), direction


# ---------------------------------------------------------------------------
# Enzyme and reaction directional logic
# ---------------------------------------------------------------------------


def _is_industrially_reversible(ec_numbers: list[str]) -> bool:
    """Returns True if any EC number matches an entry in INDUSTRIAL_REVERSIBLE_EC_PREFIXES.

    Match semantics depend on the entry shape:
      - entries ending in '.' (e.g., '1.1.1.') match any EC starting with that class prefix
      - entries WITHOUT trailing '.' (e.g., '1.6.99.1') match the EC exactly

    This dual-mode handling prevents '1.6.99.1' from over-matching '1.6.99.12'
    (which would falsely classify EC 1.6.99.12 as Old Yellow Enzyme).

    The override list is the wet-lab domain-knowledge contract: enzyme classes that are
    routinely run in the non-physiological direction in industrial biocatalysis (KREDs,
    transaminases, IREDs, EREDs, lipases, BVMOs). For these, the reverse-direction penalty
    in _compute_edge_cost is dropped to ~0.
    """
    for ec in ec_numbers:
        for entry in INDUSTRIAL_REVERSIBLE_EC_PREFIXES:
            if entry.endswith("."):
                if ec.startswith(entry):
                    return True
            else:
                if ec == entry:
                    return True
    return False


# ---------------------------------------------------------------------------
# Disk cache helpers
# ---------------------------------------------------------------------------


def _cache_root() -> Path:
    """Returns the disk cache root, honoring CHIRALAI_CACHE_ROOT env override."""
    override = os.environ.get("CHIRALAI_CACHE_ROOT")
    if override:
        return Path(override)
    return Path.home() / ".cache" / "chiralai"


def _cache_ttl_seconds() -> int:
    """Returns the cache TTL in seconds; CHIRALAI_CACHE_TTL_DAYS override available.

    Raises ValueError with a clear message if the env var is set to a non-numeric value.
    """
    raw = os.environ.get("CHIRALAI_CACHE_TTL_DAYS", "30")
    try:
        return int(raw) * 86400
    except ValueError:
        raise ValueError(
            f"CHIRALAI_CACHE_TTL_DAYS must be an integer number of days; got {raw!r}"
        )


def _disk_cache_get(category: str, key: str) -> str | None:
    """Returns cached content as a string, or None if missing/expired.

    category: one of 'kegg', 'kegg_mol', 'equilibrator'.
    key: the resource identifier (compound ID, reaction ID).

    Tolerant of races: if the file is deleted between the existence check and the
    stat/read, returns None rather than raising FileNotFoundError.
    """
    path = _cache_root() / category / f"{key}.cache"
    try:
        age_seconds = time.time() - path.stat().st_mtime
    except FileNotFoundError:
        return None
    if age_seconds > _cache_ttl_seconds():
        return None
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return None


def _disk_cache_set(category: str, key: str, content: str) -> None:
    """Writes content to disk cache, creating the subdirectory if needed."""
    path = _cache_root() / category / f"{key}.cache"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _clear_disk_cache() -> int:
    """Removes all .cache files under the cache root and returns the count of
    cache files removed.

    Used by `python -m ChiraLLM.route_predictor --clear-cache`.

    Note: only files matching `*.cache` are removed; any other files a user has
    placed under the cache root are left in place. Empty subdirectories are also
    left in place — they cost nothing and avoid surprising the user.
    """
    root = _cache_root()
    if not root.exists():
        return 0
    count = 0
    for cache_file in root.rglob("*.cache"):
        try:
            cache_file.unlink()
            count += 1
        except FileNotFoundError:
            # Concurrent removal — count it as already-cleared
            count += 1
    return count


# ---------------------------------------------------------------------------
# KEGG REST fetchers
# ---------------------------------------------------------------------------

KEGG_REST_BASE = "http://rest.kegg.jp/get"
KEGG_RETRY_DELAY_SECONDS = 1.0


def _parse_kegg_reaction_flat(text: str, rxn_id: str) -> dict | None:
    """Parses a KEGG reaction flat-file response into the standard reaction dict.

    Returns None if the EQUATION field is absent or the equation cannot be parsed.
    """
    equation = None
    ec_numbers: list[str] = []
    in_enzyme_field = False
    for line in text.splitlines():
        if line.startswith("EQUATION "):
            equation = line[12:].strip()
            in_enzyme_field = False
        elif line.startswith("ENZYME "):
            in_enzyme_field = True
            ec_numbers.extend(line[12:].split())
        elif in_enzyme_field and line.startswith(" "):
            ec_numbers.extend(line.strip().split())
        elif line.startswith("///"):
            break
        else:
            in_enzyme_field = False

    if equation is None:
        return None
    try:
        substrates, products, direction = _parse_reaction_equation(equation)
    except ValueError as e:
        logger.warning("Could not parse equation for %s: %s", rxn_id, e)
        return None
    return {
        "rxn_id": rxn_id,
        "equation": equation,
        "substrates": substrates,
        "products": products,
        "ec_numbers": ec_numbers,
        "direction": direction,
    }


def _parse_kegg_compound_reactions(text: str) -> list[str]:
    """Extracts the REACTION field from a KEGG compound flat-file as a list of R##### IDs.

    Handles multi-line REACTION fields (continuation lines start with whitespace).
    """
    reactions: list[str] = []
    in_reaction_field = False
    for line in text.splitlines():
        if line.startswith("REACTION "):
            in_reaction_field = True
            reactions.extend(line[12:].split())
        elif in_reaction_field and line.startswith(" "):
            reactions.extend(line.strip().split())
        elif line.startswith("///"):
            break
        else:
            in_reaction_field = False
    return reactions


@functools.lru_cache(maxsize=4096)
def _fetch_kegg_reaction(rxn_id: str) -> dict | None:
    """Fetches and parses a KEGG reaction by ID.

    Returns a dict with keys rxn_id, equation, substrates, products, ec_numbers,
    direction — or None on 404 or an unparseable response.

    Caching: checks disk cache first; on a network hit writes back to disk.
    LRU cache prevents repeated disk reads within a process.

    Retry: one retry after KEGG_RETRY_DELAY_SECONDS on RequestException or 5xx.
    404 returns None immediately (not transient).
    """
    cached = _disk_cache_get("kegg", rxn_id)
    if cached is not None:
        return _parse_kegg_reaction_flat(cached, rxn_id)

    for attempt in range(2):
        try:
            resp = requests.get(f"{KEGG_REST_BASE}/{rxn_id}", timeout=10)
        except requests.RequestException as e:
            logger.warning("KEGG network error for %s (attempt %d): %s", rxn_id, attempt + 1, e)
            if attempt == 0:
                time.sleep(KEGG_RETRY_DELAY_SECONDS)
                continue
            return None
        if resp.status_code == 200:
            _disk_cache_set("kegg", rxn_id, resp.text)
            return _parse_kegg_reaction_flat(resp.text, rxn_id)
        if resp.status_code == 404:
            return None
        # 5xx or unexpected: retry once then give up
        if attempt == 0:
            time.sleep(KEGG_RETRY_DELAY_SECONDS)
            continue
        logger.warning("KEGG returned %d for %s after retry", resp.status_code, rxn_id)
        return None
    return None


@functools.lru_cache(maxsize=4096)
def _fetch_compound_reactions(compound_id: str) -> list[str]:
    """Fetches the REACTION field of a KEGG compound entry as a list of R##### IDs.

    Returns an empty list if the compound has no reactions or the fetch fails.
    Uses disk cache; single network attempt (no retry — compound lookups are cheap).
    """
    cached = _disk_cache_get("kegg", f"compound_{compound_id}")
    if cached is not None:
        return _parse_kegg_compound_reactions(cached)

    try:
        resp = requests.get(f"{KEGG_REST_BASE}/{compound_id}", timeout=10)
    except requests.RequestException as e:
        logger.warning("KEGG network error for %s: %s", compound_id, e)
        return []
    if resp.status_code != 200:
        return []
    _disk_cache_set("kegg", f"compound_{compound_id}", resp.text)
    return _parse_kegg_compound_reactions(resp.text)


@functools.lru_cache(maxsize=4096)
def _fetch_kegg_mol(compound_id: str):
    """Fetches a compound's MOL file from KEGG and parses to an RDKit Mol.

    Returns None if compound has no MOL file, KEGG returns 404, or RDKit parse fails.
    Caching: checks disk cache first; on a network hit writes back to disk cache.
    LRU cache prevents repeated disk reads within a process.
    """
    cached = _disk_cache_get("kegg_mol", compound_id)
    if cached is not None:
        if not cached.strip():
            return None
        try:
            mol = Chem.MolFromMolBlock(cached)
            if mol is None:
                logger.warning("RDKit could not parse cached MOL for %s", compound_id)
            return mol
        except Exception as e:
            logger.warning("RDKit cached MOL parse exception for %s: %s", compound_id, e)
            return None

    try:
        resp = requests.get(f"{KEGG_REST_BASE}/{compound_id}/mol", timeout=10)
    except requests.RequestException as e:
        logger.warning("KEGG MOL network error for %s: %s", compound_id, e)
        return None
    if resp.status_code != 200 or not resp.text.strip():
        return None
    _disk_cache_set("kegg_mol", compound_id, resp.text)
    try:
        mol = Chem.MolFromMolBlock(resp.text)
        if mol is None:
            logger.warning("RDKit could not parse MOL for %s", compound_id)
        return mol
    except Exception as e:
        logger.warning("RDKit MOL parse exception for %s: %s", compound_id, e)
        return None


# ---------------------------------------------------------------------------
# eQuilibrator REST fetcher
# ---------------------------------------------------------------------------

EQUILIBRATOR_REST_BASE = "https://equilibrator.weizmann.ac.il/api/v1/reaction"


def _fetch_delta_g_kj_per_mol(rxn_id: str) -> float | None:
    """Fetches standard ΔrG' (kJ/mol, forward direction) from eQuilibrator REST.

    Disk-cached. Returns None on network failure, non-200 response, or unparseable JSON.
    Caller is expected to fall back to FALLBACK_DELTA_G_KJ when None is returned.
    """
    cached = _disk_cache_get("equilibrator", rxn_id)
    if cached is not None:
        try:
            return float(cached)
        except ValueError:
            return None

    try:
        resp = requests.get(f"{EQUILIBRATOR_REST_BASE}/{rxn_id}", timeout=10)
    except requests.RequestException as e:
        logger.warning("eQuilibrator network error for %s: %s", rxn_id, e)
        return None
    if resp.status_code != 200:
        return None
    try:
        payload = resp.json()
        dg = float(payload["standard_dg_prime"])
    except (ValueError, KeyError, TypeError) as e:
        logger.warning("eQuilibrator response parse error for %s: %s", rxn_id, e)
        return None
    _disk_cache_set("equilibrator", rxn_id, str(dg))
    return dg
