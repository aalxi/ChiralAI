"""Tier 1 route predictor: best-first backward search through the KEGG reaction
graph from a target compound to curated central metabolites.

See docs/superpowers/specs/2026-04-27-route-predictor-tier1-design.md for the
full design rationale.

Public API: predict_route(compound_id, mode='top_n', n=3, budget=500) -> RouteResult
"""

import heapq
import logging
import os
import re
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

import requests
from rdkit import Chem
from rdkit.Chem import DataStructs, rdFingerprintGenerator

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


# Cofactors and inorganics excluded from search graph expansion in _astar_search.
# Carbon-bearing intermediates (acetyl-CoA, malonyl-CoA, succinyl-CoA) are NOT in this set —
# they carry real carbon flow and must remain traversable. Free CoA (C00010) IS skipped
# because it's the cofactor moiety, not a biosynthetic carbon source.
#
# Why this set is large: the 2026-05-16 chiral-route benchmark showed that without skipping
# redox/energy/methyl cofactors, the A* search exploits them as graph-connectivity hubs and
# terminates at biologically wrong central metabolites via chemically valid but real-world
# nonsensical 2-hop shortcuts (e.g., mandelate → NADH → glutamate). See benchmarks/REPORT.md.
COFACTOR_SKIP_IDS: frozenset[str] = frozenset({
    # Inorganics — H+, water, O2
    "C00080", "C00001", "C00007",
    # Redox carriers — NAD/NADH, NADP/NADPH, FAD/FADH2
    "C00003", "C00004", "C00005", "C00006", "C00016", "C01352",
    # Energy / phosphate — ATP, ADP, AMP, Pi, PPi
    "C00002", "C00008", "C00020", "C00009", "C00013",
    # Methyl donors — SAM, SAH
    "C00019", "C00021",
    # Free CoA (NOT acyl-CoA species, which carry real carbon)
    "C00010",
})


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

_BASE_EDGE_COST = 1.0
_REVERSIBLE_REVERSE_PENALTY = 0.5  # reverse traversal of a KEGG-marked-reversible reaction
_IRREVERSIBLE_REVERSE_PENALTY = 2.0  # reverse traversal of a KEGG-marked-forward-only reaction

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

KEGG_REST_BASE = "https://rest.kegg.jp/get"
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


def _success_only_cache(maxsize: int = 4096):
    """In-process cache that stores only truthy results.

    Replaces functools.lru_cache for KEGG fetchers. The previous use of
    @functools.lru_cache poisoned the cache when a transient network failure
    returned [] or None: the falsy result was memoized for the entire process,
    so a later call (after the network recovered) silently returned the empty
    failure instead of retrying. Caching only truthy results preserves the
    in-process speedup for successful lookups while letting transient
    failures be retried.

    Provides .cache_clear() for test compatibility with the previous API.
    """
    def decorator(func):
        cache: dict = {}
        def wrapper(arg):
            if arg in cache:
                return cache[arg]
            result = func(arg)
            if result:
                if len(cache) >= maxsize:
                    cache.pop(next(iter(cache)))  # FIFO eviction
                cache[arg] = result
            return result
        wrapper.cache_clear = cache.clear
        return wrapper
    return decorator


@_success_only_cache(maxsize=4096)
def _fetch_kegg_reaction(rxn_id: str) -> dict | None:
    """Fetches and parses a KEGG reaction by ID.

    Returns a dict with keys rxn_id, equation, substrates, products, ec_numbers,
    direction — or None on 404 or an unparseable response.

    Caching: checks disk cache first; on a network hit writes back to disk.
    No in-process LRU: it would memoize transient failures (which return None
    without writing to disk) and poison the result for the rest of the process.

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


@_success_only_cache(maxsize=4096)
def _fetch_compound_reactions(compound_id: str) -> list[str]:
    """Fetches the REACTION field of a KEGG compound entry as a list of R##### IDs.

    Returns an empty list if the compound has no reactions or the fetch fails.
    Uses disk cache; single network attempt (no retry — compound lookups are cheap).
    Uses _success_only_cache (not lru_cache) so a transient empty result does not
    poison the cache for the rest of the process.
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


@_success_only_cache(maxsize=4096)
def _fetch_kegg_mol(compound_id: str):
    """Fetches a compound's MOL file from KEGG and parses to an RDKit Mol.

    Returns None if compound has no MOL file, KEGG returns 404, or RDKit parse fails.
    Caching: checks disk cache first; on a network hit writes back to disk cache.
    Uses _success_only_cache (not lru_cache) so transient None returns are not
    memoized for the rest of the process.
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


# ---------------------------------------------------------------------------
# A* heuristic: Morgan-fingerprint Tanimoto distance to central metabolites
# ---------------------------------------------------------------------------

_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(
    radius=2, fpSize=2048, includeChirality=False
)
_central_fingerprints_cache: dict | None = None
_TANIMOTO_FALLBACK_DISTANCE = 0.5


def _get_central_fingerprints() -> dict:
    """Lazy-initialized Morgan fingerprints for every compound in CENTRAL_METABOLITES.

    Returns a dict mapping compound_id -> ExplicitBitVect. Computed once per module
    lifetime. Lazy init defers ~50ms × N RDKit calls until first search, and supports
    future catalog swaps without module-import cost.
    """
    global _central_fingerprints_cache
    if _central_fingerprints_cache is not None:
        return _central_fingerprints_cache
    fps = {}
    for cid in CENTRAL_METABOLITES:
        mol = _fetch_kegg_mol(cid)
        if mol is not None:
            fps[cid] = _morgan_gen.GetFingerprint(mol)
    _central_fingerprints_cache = fps
    logger.info("Computed Morgan fingerprints for %d central metabolites", len(fps))
    return fps


def _tanimoto_to_central(compound_id: str) -> float:
    """Returns 1 - max(Tanimoto similarity to any central metabolite). Range [0, 1].

    0 means the compound is (essentially) identical to a central metabolite.
    1 means maximally distant from all central metabolites.

    Falls back to uniform 0.5 if the compound's MOL file is unavailable (per
    spec §5.1 case #11): heuristic degrades to a constant rather than excluding
    the node from search.

    includeChirality=False is intentional: KEGG MOL files lack stereochemistry
    annotations for many compounds; including chirality bits would unfairly inflate
    distance for compounds KEGG doesn't represent stereochemistry for.
    """
    mol = _fetch_kegg_mol(compound_id)
    if mol is None:
        return _TANIMOTO_FALLBACK_DISTANCE
    fp = _morgan_gen.GetFingerprint(mol)
    centrals = _get_central_fingerprints()
    if not centrals:
        return _TANIMOTO_FALLBACK_DISTANCE
    max_sim = max(
        DataStructs.TanimotoSimilarity(fp, central_fp) for central_fp in centrals.values()
    )
    return 1.0 - max_sim


# ---------------------------------------------------------------------------
# Edge cost computation for A* search
# ---------------------------------------------------------------------------


def _compute_edge_cost(reaction: dict, traversed_direction: str) -> dict:
    """Computes the cost breakdown for traversing one reaction edge in the search.

    Components:
      base               — flat per-step cost (always 1.0)
      thermodynamic      — proportional to |ΔG| when traversing reverse, 0 when forward
      directionality     — KEGG reversibility penalty (0 forward, varies by direction marker reverse)
      industrial_override — negative discount for industrially-reversible EC families when reverse

    traversed_direction: 'forward' (we're walking in the reaction's natural direction)
                        or 'reverse' (we're walking against it).
    """
    base = _BASE_EDGE_COST
    thermodynamic = 0.0
    directionality = 0.0
    industrial_override = 0.0

    if traversed_direction == "reverse":
        # Thermodynamic penalty: proportional to |ΔG| in the unfavorable direction.
        dg = _fetch_delta_g_kj_per_mol(reaction["rxn_id"])
        dg_magnitude = abs(dg) if dg is not None else FALLBACK_DELTA_G_KJ
        thermodynamic = THERMO_PENALTY_PER_KJ * dg_magnitude

        # Directionality penalty: depends on KEGG's reversibility annotation.
        if reaction["direction"] == DIRECTION_FORWARD_ONLY:
            directionality = _IRREVERSIBLE_REVERSE_PENALTY
        else:
            directionality = _REVERSIBLE_REVERSE_PENALTY

        # Industrial-override discount: applied only on reverse traversal.
        if _is_industrially_reversible(reaction["ec_numbers"]):
            # Discount equal to the directionality penalty (i.e., zero-out the reverse cost
            # for these EC families, but keep the base + thermodynamic components).
            industrial_override = -directionality

    total = base + thermodynamic + directionality + industrial_override
    return {
        "base": base,
        "thermodynamic": thermodynamic,
        "directionality": directionality,
        "industrial_override": industrial_override,
        "total": total,
    }


# ---------------------------------------------------------------------------
# Route dataclasses — typed containers for search results
# ---------------------------------------------------------------------------


@dataclass
class RouteStep:
    reaction_id: str
    ec_numbers: list[str]
    precursor_id: str           # the upstream compound (one step closer to terminal)
    intermediate_id: str        # the downstream compound (we just came from this)
    edge_cost_breakdown: dict[str, float]
    traversed_direction: str    # 'forward' or 'reverse'
    # Naming note: "substrate"/"product" would be ambiguous in retrosynthesis traversal;
    # precursor_id (upstream) and intermediate_id (downstream) name graph position directly.


@dataclass
class Route:
    target_id: str
    steps: list[RouteStep]                 # ordered target → terminal precursor
    terminal_precursor_id: str
    terminal_precursor_name: str
    total_cost: float
    cost_breakdown: dict[str, float]       # summed components across all steps
    warnings: list[str] = field(default_factory=list)


@dataclass
class RouteResult:
    target_id: str
    mode: str
    routes: list[Route]
    nodes_explored: int
    budget_exhausted: bool
    warnings: list[str] = field(default_factory=list)
    status: str = "success"


# ---------------------------------------------------------------------------
# A* backward search
# ---------------------------------------------------------------------------


def _astar_search(target_id: str, budget: int = DEFAULT_BUDGET, depth_cap: int = DEFAULT_DEPTH_CAP) -> dict:
    """Best-first backward search from target. Returns a DAG and search metadata.

    Algorithm: weighted A* / greedy best-first.
      f(n) = g(n) + h(n)
      g(n) = cumulative edge cost from target to current node
      h(n) = TANIMOTO_HEURISTIC_WEIGHT * Tanimoto distance to nearest central metabolite

    Termination: priority queue empty, OR nodes_explored >= budget, OR (per-node) depth > depth_cap.
    Search continues past the first central-metabolite hit so diverse alternates can be found.
    """
    counter = 0  # heap tie-breaker
    initial_h = TANIMOTO_HEURISTIC_WEIGHT * _tanimoto_to_central(target_id)
    queue: list = [(initial_h, counter, 0, target_id, 0.0)]
    heapq.heapify(queue)

    visited_dag: dict[str, list] = {}
    leaf_ids: list[str] = []
    nodes_explored = 0
    budget_exhausted = False

    best_g_seen: dict[str, float] = {target_id: 0.0}

    while queue:
        if nodes_explored >= budget:
            budget_exhausted = True
            break

        f_score, _, depth, compound_id, g_score = heapq.heappop(queue)

        if g_score > best_g_seen.get(compound_id, float("inf")):
            continue

        nodes_explored += 1

        if compound_id in CENTRAL_METABOLITES:
            if compound_id not in leaf_ids:
                leaf_ids.append(compound_id)
            continue

        if depth >= depth_cap:
            continue

        reaction_ids = _fetch_compound_reactions(compound_id)

        for rxn_id in reaction_ids:
            reaction = _fetch_kegg_reaction(rxn_id)
            if reaction is None:
                continue

            substrate_ids = [cid for _, cid in reaction["substrates"]]
            product_ids = [cid for _, cid in reaction["products"]]

            if compound_id in product_ids:
                traversed_direction = "forward"
                precursor_candidates = substrate_ids
            elif compound_id in substrate_ids:
                traversed_direction = "reverse"
                precursor_candidates = product_ids
            else:
                continue

            edge_cost = _compute_edge_cost(reaction, traversed_direction)
            new_g = g_score + edge_cost["total"]

            for precursor_id in precursor_candidates:
                if precursor_id == compound_id:
                    continue
                if precursor_id in COFACTOR_SKIP_IDS:
                    continue

                if new_g >= best_g_seen.get(precursor_id, float("inf")):
                    continue

                best_g_seen[precursor_id] = new_g
                visited_dag.setdefault(precursor_id, []).append({
                    "parent_id": compound_id,
                    "reaction": reaction,
                    "edge_cost": edge_cost,
                    "depth": depth + 1,
                    "g_score": new_g,
                    "traversed_direction": traversed_direction,
                })

                h_new = TANIMOTO_HEURISTIC_WEIGHT * _tanimoto_to_central(precursor_id)
                f_new = new_g + h_new
                counter += 1
                heapq.heappush(queue, (f_new, counter, depth + 1, precursor_id, new_g))

    return {
        "target_id": target_id,
        "visited_dag": visited_dag,
        "leaf_ids": leaf_ids,
        "nodes_explored": nodes_explored,
        "budget_exhausted": budget_exhausted,
    }


# ---------------------------------------------------------------------------
# Output formatters
# ---------------------------------------------------------------------------

def _backtrack_route(dag_result: dict, leaf_id: str) -> "Route | None":
    """Walks the DAG from a leaf back to the target, building an ordered Route.
    Picks the cheapest parent at each level (greedy backtrack).
    Returns None if no path exists."""
    visited_dag = dag_result["visited_dag"]
    target_id = dag_result["target_id"]

    if leaf_id not in visited_dag and leaf_id != target_id:
        return None

    steps_reverse: list[RouteStep] = []
    current = leaf_id
    seen = {leaf_id}

    while current != target_id:
        edges = visited_dag.get(current, [])
        if not edges:
            return None  # broken path
        edge = min(edges, key=lambda e: e["g_score"])
        parent = edge["parent_id"]
        if parent in seen:
            return None  # cycle
        seen.add(parent)

        step = RouteStep(
            reaction_id=edge["reaction"]["rxn_id"],
            ec_numbers=edge["reaction"]["ec_numbers"],
            precursor_id=current,           # upstream
            intermediate_id=parent,         # downstream (closer to target)
            edge_cost_breakdown=dict(edge["edge_cost"]),
            traversed_direction=edge.get("traversed_direction", "forward"),
        )
        steps_reverse.append(step)
        current = parent

    # steps_reverse is leaf → target; reverse so steps go target → terminal precursor
    steps = list(reversed(steps_reverse))

    total_cost = sum(s.edge_cost_breakdown["total"] for s in steps)
    cost_breakdown = {
        "base": sum(s.edge_cost_breakdown["base"] for s in steps),
        "thermodynamic": sum(s.edge_cost_breakdown["thermodynamic"] for s in steps),
        "directionality": sum(s.edge_cost_breakdown["directionality"] for s in steps),
        "industrial_override": sum(s.edge_cost_breakdown["industrial_override"] for s in steps),
    }

    return Route(
        target_id=target_id,
        steps=steps,
        terminal_precursor_id=leaf_id,
        terminal_precursor_name=CENTRAL_METABOLITES.get(leaf_id, leaf_id),
        total_cost=total_cost,
        cost_breakdown=cost_breakdown,
        warnings=[],
    )


def _extract_top_n(dag_result: dict, n: int = DEFAULT_MAX_ROUTES) -> list["Route"]:
    """Returns up to n routes, sorted by total_cost ascending."""
    routes = []
    for leaf_id in dag_result["leaf_ids"]:
        route = _backtrack_route(dag_result, leaf_id)
        if route is not None:
            routes.append(route)
    routes.sort(key=lambda r: r.total_cost)
    return routes[:n]


def _extract_full_tree(dag_result: dict) -> list["Route"]:
    """Returns a single Route whose 'steps' encode the canonical (cheapest) path.
    The full DAG is preserved in JSON sidecar; warnings list other leaves."""
    if not dag_result["leaf_ids"]:
        return []
    cheapest_leaf = min(
        dag_result["leaf_ids"],
        key=lambda lid: min(
            (e["g_score"] for e in dag_result["visited_dag"].get(lid, [])),
            default=float("inf"),
        ),
    )
    canonical = _backtrack_route(dag_result, cheapest_leaf)
    if canonical is None:
        return []
    other_leaves = [lid for lid in dag_result["leaf_ids"] if lid != cheapest_leaf]
    if other_leaves:
        canonical.warnings.append(
            f"full_tree mode: {len(other_leaves)} additional leaves not shown in steps "
            f"(see JSON sidecar for full DAG): {', '.join(other_leaves)}"
        )
    return [canonical]


def _extract_shortest_plus_diverse(dag_result: dict, n_diverse: int = 2) -> list["Route"]:
    """Returns the shortest-hop route plus n_diverse maximally-different alternates.
    Diversity is measured by terminal_precursor_id; alternates with the same terminal as the
    shortest are skipped.
    """
    all_routes = []
    for leaf_id in dag_result["leaf_ids"]:
        route = _backtrack_route(dag_result, leaf_id)
        if route is not None:
            all_routes.append(route)
    if not all_routes:
        return []

    shortest = min(all_routes, key=lambda r: (len(r.steps), r.total_cost))
    selected = [shortest]
    seen_terminals = {shortest.terminal_precursor_id}

    remaining = [r for r in all_routes if r.terminal_precursor_id not in seen_terminals]
    remaining.sort(key=lambda r: r.total_cost)
    for r in remaining[:n_diverse]:
        selected.append(r)
        seen_terminals.add(r.terminal_precursor_id)

    return selected


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

_KEGG_ID_RE = re.compile(r"^C\d{5}$")
_VALID_MODES = {"top_n", "full_tree", "shortest_plus_diverse"}


def predict_route(
    compound_id: "str | None",
    mode: str = "top_n",
    n: int = DEFAULT_MAX_ROUTES,
    budget: int = DEFAULT_BUDGET,
) -> RouteResult:
    """Predicts biosynthetic routes from a target KEGG compound to central metabolites.

    Returns a RouteResult with status indicating success or specific failure mode.
    Never raises — all failures encoded in result.status and result.warnings.

    Modes:
      'top_n'                    — up to n routes ranked by total cost (default)
      'full_tree'                — single canonical route + warning listing other leaves
      'shortest_plus_diverse'    — shortest hop count + n diverse alternates
    """
    # Input validation
    if not compound_id:
        return RouteResult(target_id="", mode=mode, routes=[], nodes_explored=0,
                           budget_exhausted=False, status="no_kegg_id")
    if not _KEGG_ID_RE.match(compound_id):
        # Allow synthetic test IDs (C_*) through with a different status — but production
        # tests use _KEGG_ID_RE which requires C\d{5}. For non-matching IDs, return invalid.
        # Note: tests that need to bypass the regex check pre-populate the disk cache and
        # use synthetic IDs — predict_route's regex catches malformed user input.
        if not (compound_id.startswith("C_") and "_" in compound_id):
            return RouteResult(target_id=compound_id, mode=mode, routes=[],
                               nodes_explored=0, budget_exhausted=False,
                               status="invalid_kegg_id",
                               warnings=[f"compound_id {compound_id!r} does not match KEGG format C#####"])
    if mode not in _VALID_MODES:
        return RouteResult(target_id=compound_id, mode=mode, routes=[],
                           nodes_explored=0, budget_exhausted=False,
                           status="invalid_mode",
                           warnings=[f"mode {mode!r} must be one of: {sorted(_VALID_MODES)}"])

    # Pre-flight: does the target compound have any reactions?
    target_reactions = _fetch_compound_reactions(compound_id)
    if not target_reactions:
        cached_compound = _disk_cache_get("kegg", f"compound_{compound_id}")
        if cached_compound is None:
            return RouteResult(target_id=compound_id, mode=mode, routes=[],
                               nodes_explored=0, budget_exhausted=False,
                               status="target_not_in_kegg",
                               warnings=[f"KEGG returned no data for {compound_id}; compound may have been deprecated or merged"])
        return RouteResult(target_id=compound_id, mode=mode, routes=[],
                           nodes_explored=0, budget_exhausted=False,
                           status="target_has_no_reactions",
                           warnings=[f"KEGG entry for {compound_id} lists no reactions; compound may be a leaf metabolite or unconnected"])

    logger.info("Predicting routes for %s, mode=%s, budget=%d", compound_id, mode, budget)
    dag_result = _astar_search(compound_id, budget=budget)

    if mode == "top_n":
        routes = _extract_top_n(dag_result, n=n)
    elif mode == "full_tree":
        routes = _extract_full_tree(dag_result)
    else:  # shortest_plus_diverse
        routes = _extract_shortest_plus_diverse(dag_result, n_diverse=n)

    warnings: list[str] = []
    if dag_result["budget_exhausted"]:
        warnings.append(f"Budget of {budget} nodes exhausted; some routes may be missing")

    if not routes:
        return RouteResult(target_id=compound_id, mode=mode, routes=[],
                           nodes_explored=dag_result["nodes_explored"],
                           budget_exhausted=dag_result["budget_exhausted"],
                           status="no_route_found",
                           warnings=warnings + [f"Search exhausted {dag_result['nodes_explored']} nodes without reaching a central metabolite"])

    logger.info("Found %d routes for %s, %d nodes explored", len(routes), compound_id, dag_result["nodes_explored"])
    return RouteResult(
        target_id=compound_id,
        mode=mode,
        routes=routes,
        nodes_explored=dag_result["nodes_explored"],
        budget_exhausted=dag_result["budget_exhausted"],
        warnings=warnings,
        status="success",
    )


# Module CLI entry point for cache management
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--clear-cache":
        n_cleared = _clear_disk_cache()
        print(f"Cleared {n_cleared} cached files from {_cache_root()}")
        sys.exit(0)
    print("Usage: python -m ChiraLLM.route_predictor --clear-cache")
    sys.exit(1)
