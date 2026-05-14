# Tier 1 Route Predictor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `ChiraLLM/route_predictor.py` — a best-first KEGG-graph backward search from a target compound to curated central metabolites, wired into the existing pipeline so every suggestion gets ranked routes with per-step ee composition.

**Architecture:** Single-file module (~450 LOC) holding KEGG fetchers, eQuilibrator-grounded edge cost model, Morgan-fingerprint Tanimoto heuristic, weighted A\* search, and 3 output-mode formatters. Wired between `database_validator.py` and `brenda_client.py` in `main.py`; CSV/JSON output extended; smoke test fixture added. **Spec:** `docs/superpowers/specs/2026-04-27-route-predictor-tier1-design.md`.

**Tech Stack:** Python 3.12, RDKit (Morgan fingerprints + MOL parsing), `requests` (KEGG REST + eQuilibrator REST), `functools.lru_cache` (in-process), filesystem disk cache, `dataclasses`. Tests: `pytest` + `pytest-mock`. No new pip deps for the runtime module — only `pytest` / `pytest-mock` added as dev deps.

**Parallelism map for subagent-driven-development:**
- Tasks 1–12 are sequential (all edit `ChiraLLM/route_predictor.py`).
- Tasks 13, 14, 15 are independent and can run in parallel (different files: `main.py`, `utils/file_saver.py`, `ChiraLLM/enantioselectivity_scorer.py`).
- Tasks 16 and 17 each depend on a Wave-2 task (16 depends on 13, 17 depends on 12) — schedule after their dependencies.
- Tasks 18 and 19 are independent of the dependency graph after Task 12 — safe to dispatch in parallel with anything in Wave 2/3.

**Worktree note:** This is a substantial feature touching 7 files. Recommend creating an isolated git worktree via `superpowers:using-git-worktrees` at execution start. The plan is written assuming all paths are relative to the project root regardless of worktree location.

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `ChiraLLM/route_predictor.py` | Create | The module — constants, fetchers, cost model, heuristic, A\* search, formatters, public API |
| `tests/__init__.py` | Create | Make `tests/` a package |
| `tests/conftest.py` | Create | Shared fixtures: `tmp_cache_dir`, `mock_kegg`, `mock_equilibrator` |
| `tests/test_route_predictor_unit.py` | Create | Layer 1 unit tests, mocked external calls |
| `tests/test_route_predictor_integration.py` | Create | Layer 2 integration tests, real KEGG, marked `@pytest.mark.integration` |
| `tests/fixtures/synthetic_kegg.json` | Create | Hand-built synthetic KEGG reaction graph for unit tests |
| `pyproject.toml` | Create | Pytest config: markers (integration, live_equilibrator), default `addopts` |
| `requirements-dev.txt` | Create | `pytest`, `pytest-mock` |
| `main.py` | Modify | Add `predict_route` call after `query_kegg`; dedupe ECs across routes; per-route `check_feasibility` |
| `utils/file_saver.py` | Modify | Add 5 flat columns for `route_top1_*`; preserve full `route_prediction` in JSON sidecar |
| `ChiraLLM/enantioselectivity_scorer.py` | Modify | Add `_compose_route_ee` helper; extend `score_suggestion` to compute per-route composed ee |
| `smoke_test.py` | Modify | Add `route_prediction` fixture per Sprint 1 acceptance criterion #7 |
| `scripts/warm_equilibrator_cache.py` | Create | One-shot script: enumerate iJO1366 reactions + KEGG most-cited, pre-fetch ΔG to disk cache |
| `CLAUDE.md` | Modify | Add `route_predictor.py` to module status table |
| `README.md` | Modify | Move "retrosynthetic route prediction" from Roadmap to Working Features |

---

## Task 0: Test infrastructure setup

**Files:**
- Create: `tests/__init__.py`
- Create: `tests/conftest.py`
- Create: `tests/fixtures/synthetic_kegg.json`
- Create: `pyproject.toml`
- Create: `requirements-dev.txt`

- [ ] **Step 1: Create test package directory**

```bash
mkdir -p tests/fixtures
touch tests/__init__.py
```

- [ ] **Step 2: Add dev dependencies**

Create `requirements-dev.txt`:

```
pytest>=8.0
pytest-mock>=3.12
pytest-cov>=4.1
```

Install:

```bash
.venv/bin/pip install -r requirements-dev.txt
```

Expected output: `Successfully installed pytest-X.X.X pytest-mock-X.X.X pytest-cov-X.X.X`

- [ ] **Step 3: Create pyproject.toml with pytest config**

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
markers = [
    "integration: requires live KEGG access (run with -m integration)",
    "live_equilibrator: requires live eQuilibrator API (run with -m live_equilibrator)",
]
addopts = "-m 'not integration and not live_equilibrator'"
```

- [ ] **Step 4: Create the synthetic KEGG fixture**

Create `tests/fixtures/synthetic_kegg.json`:

```json
{
  "compounds": {
    "C_TARGET": {"name": "synthetic target", "reactions": ["R_S1", "R_S2"], "mol": "synthetic_target.mol"},
    "C_INTERMEDIATE_A": {"name": "intermediate A", "reactions": ["R_S2", "R_S3"], "mol": "intermediate_a.mol"},
    "C_INTERMEDIATE_B": {"name": "intermediate B", "reactions": ["R_S4"], "mol": "intermediate_b.mol"},
    "C_PRECURSOR_A": {"name": "precursor A (central)", "reactions": ["R_S3"], "mol": "precursor_a.mol"},
    "C_PRECURSOR_B": {"name": "precursor B (central)", "reactions": ["R_S4"], "mol": "precursor_b.mol"},
    "C_DEAD_END": {"name": "dead-end compound", "reactions": [], "mol": "dead_end.mol"}
  },
  "reactions": {
    "R_S1": {
      "equation": "C_INTERMEDIATE_A + C00006 <=> C_TARGET + C00005",
      "ec_numbers": ["1.1.1.184"],
      "delta_g_kj_per_mol": -8.0
    },
    "R_S2": {
      "equation": "C_INTERMEDIATE_B + C00006 => C_TARGET + C00005",
      "ec_numbers": ["1.1.1.999"],
      "delta_g_kj_per_mol": -25.0
    },
    "R_S3": {
      "equation": "C_PRECURSOR_A <=> C_INTERMEDIATE_A",
      "ec_numbers": ["5.4.99.1"],
      "delta_g_kj_per_mol": -2.0
    },
    "R_S4": {
      "equation": "C_PRECURSOR_B <=> C_INTERMEDIATE_B",
      "ec_numbers": ["1.1.1.184"],
      "delta_g_kj_per_mol": -1.0
    }
  },
  "central_metabolites_for_test": ["C_PRECURSOR_A", "C_PRECURSOR_B"]
}
```

- [ ] **Step 5: Create conftest.py with shared fixtures**

```python
import json
from pathlib import Path
import pytest


@pytest.fixture
def tmp_cache_dir(tmp_path, monkeypatch):
    """Redirect ChiralAI's disk cache to a tmp directory for the test."""
    cache_dir = tmp_path / "chiralai_cache"
    cache_dir.mkdir()
    monkeypatch.setenv("CHIRALAI_CACHE_ROOT", str(cache_dir))
    return cache_dir


@pytest.fixture
def synthetic_kegg():
    """Returns the parsed synthetic KEGG fixture dict."""
    fixture_path = Path(__file__).parent / "fixtures" / "synthetic_kegg.json"
    return json.loads(fixture_path.read_text(encoding="utf-8"))


@pytest.fixture
def mock_kegg(synthetic_kegg, monkeypatch):
    """Patches _fetch_kegg_reaction to return synthetic data.
    Tests that need different fixtures can extend by mutating synthetic_kegg before use."""
    from ChiraLLM import route_predictor

    def fake_fetch_reaction(rxn_id):
        rxn = synthetic_kegg["reactions"].get(rxn_id)
        if rxn is None:
            return None
        substrates, products, direction = route_predictor._parse_reaction_equation(rxn["equation"])
        return {
            "rxn_id": rxn_id,
            "equation": rxn["equation"],
            "substrates": substrates,
            "products": products,
            "ec_numbers": rxn["ec_numbers"],
            "direction": direction,
        }

    def fake_fetch_compound_reactions(cid):
        compound = synthetic_kegg["compounds"].get(cid)
        return compound["reactions"] if compound else []

    monkeypatch.setattr(route_predictor, "_fetch_kegg_reaction", fake_fetch_reaction)
    monkeypatch.setattr(route_predictor, "_fetch_compound_reactions", fake_fetch_compound_reactions)
    return synthetic_kegg


@pytest.fixture
def mock_equilibrator(synthetic_kegg, monkeypatch):
    """Patches _fetch_delta_g_kj_per_mol to return synthetic ΔG values."""
    from ChiraLLM import route_predictor

    def fake_fetch_dg(rxn_id):
        rxn = synthetic_kegg["reactions"].get(rxn_id)
        return rxn["delta_g_kj_per_mol"] if rxn else None

    monkeypatch.setattr(route_predictor, "_fetch_delta_g_kj_per_mol", fake_fetch_dg)
```

- [ ] **Step 6: Verify pytest discovery works**

```bash
.venv/bin/pytest --collect-only 2>&1 | tail -5
```

Expected: `no tests collected` (no test files yet) without errors.

- [ ] **Step 7: Commit**

```bash
git add tests/ pyproject.toml requirements-dev.txt
git commit -m "add pytest scaffolding and synthetic KEGG fixture for route predictor tests

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 1: Module skeleton + constants

**Files:**
- Create: `ChiraLLM/route_predictor.py`
- Test: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_route_predictor_unit.py`:

```python
"""Layer 1 unit tests for ChiraLLM/route_predictor.py.
All external calls (KEGG, eQuilibrator, RDKit MOL parsing) are mocked.
Coverage target: ≥85%."""

from ChiraLLM import route_predictor


def test_constants_present_and_well_formed():
    """The two scientific-contract constants must exist and have expected shape."""
    assert isinstance(route_predictor.CENTRAL_METABOLITES, dict)
    assert "C00022" in route_predictor.CENTRAL_METABOLITES, "pyruvate must be central"
    assert "C00024" in route_predictor.CENTRAL_METABOLITES, "acetyl-CoA must be central"
    assert "C00141" in route_predictor.CENTRAL_METABOLITES, "alpha-ketoisovalerate must be central"
    assert all(k.startswith("C") and len(k) == 6 for k in route_predictor.CENTRAL_METABOLITES)
    assert all(isinstance(v, str) and v for v in route_predictor.CENTRAL_METABOLITES.values())

    assert isinstance(route_predictor.INDUSTRIAL_REVERSIBLE_EC_PREFIXES, list)
    assert "1.1.1." in route_predictor.INDUSTRIAL_REVERSIBLE_EC_PREFIXES, "KREDs must be industrially reversible"
    assert "2.6.1." in route_predictor.INDUSTRIAL_REVERSIBLE_EC_PREFIXES, "transaminases too"


def test_default_constants_have_sane_values():
    assert route_predictor.DEFAULT_BUDGET == 500
    assert route_predictor.DEFAULT_DEPTH_CAP == 8
    assert route_predictor.DEFAULT_MAX_ROUTES == 3
    assert 0.0 < route_predictor.THERMO_PENALTY_PER_KJ < 1.0
    assert 0.0 < route_predictor.FALLBACK_DELTA_G_KJ < 100.0
    assert 0.0 < route_predictor.TANIMOTO_HEURISTIC_WEIGHT < 10.0
```

- [ ] **Step 2: Run test to verify it fails**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py -v
```

Expected: `ModuleNotFoundError: No module named 'ChiralLLM.route_predictor'` or `ImportError`.

- [ ] **Step 3: Create the module skeleton with constants**

```python
"""Tier 1 route predictor: best-first backward search through the KEGG reaction
graph from a target compound to curated central metabolites.

See docs/superpowers/specs/2026-04-27-route-predictor-tier1-design.md for the
full design rationale.

Public API: predict_route(compound_id, mode='top_n', n=3, budget=500) -> RouteResult
"""

import logging

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


INDUSTRIAL_REVERSIBLE_EC_PREFIXES: list[str] = [
    "1.1.1.",      # KREDs / aldo-keto reductases
    "2.6.1.",      # transaminases
    "1.5.1.",      # IREDs (imine reductases)
    "1.6.99.1",    # Old Yellow Enzyme (ene-reductases)
    "3.1.1.",      # lipases
    "1.14.13.22",  # cyclohexanone monooxygenase (Baeyer-Villiger archetype)
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py -v
```

Expected: 2 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "scaffold route_predictor module with scientific-contract constants

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 2: `_parse_reaction_equation`

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
import pytest


class TestParseReactionEquation:
    def test_simple_reversible(self):
        s, p, d = route_predictor._parse_reaction_equation("C00033 + C00010 <=> C00024 + C00011")
        assert s == [(1, "C00033"), (1, "C00010")]
        assert p == [(1, "C00024"), (1, "C00011")]
        assert d == "reversible"

    def test_simple_forward_only(self):
        s, p, d = route_predictor._parse_reaction_equation("C00033 => C00024")
        assert s == [(1, "C00033")]
        assert p == [(1, "C00024")]
        assert d == "forward_only"

    def test_with_coefficients(self):
        s, p, d = route_predictor._parse_reaction_equation("2 C00006 + C00149 <=> 2 C00005 + C00026")
        assert s == [(2, "C00006"), (1, "C00149")]
        assert p == [(2, "C00005"), (1, "C00026")]
        assert d == "reversible"

    def test_extra_whitespace_tolerated(self):
        s, p, d = route_predictor._parse_reaction_equation("  C00033   +   C00010   <=>   C00024  ")
        assert s == [(1, "C00033"), (1, "C00010")]
        assert p == [(1, "C00024")]

    def test_missing_arrow_raises(self):
        with pytest.raises(ValueError, match="No reaction arrow"):
            route_predictor._parse_reaction_equation("C00033 + C00010 C00024")

    def test_unparseable_token_raises(self):
        with pytest.raises(ValueError, match="Unparseable token"):
            route_predictor._parse_reaction_equation("C00033 + foo <=> C00024")

    def test_empty_equation_raises(self):
        with pytest.raises(ValueError):
            route_predictor._parse_reaction_equation("")
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestParseReactionEquation -v
```

Expected: `AttributeError: module 'ChiraLLM.route_predictor' has no attribute '_parse_reaction_equation'` (or all 7 tests fail with AttributeError).

- [ ] **Step 3: Implement `_parse_reaction_equation`**

Append to `ChiraLLM/route_predictor.py` (after the constants block):

```python
import re

DIRECTION_REVERSIBLE = "reversible"
DIRECTION_FORWARD_ONLY = "forward_only"
_COMPOUND_TOKEN_RE = re.compile(r"^(?:(\d+)\s+)?(C\d{5})$")


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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestParseReactionEquation -v
```

Expected: 7 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add _parse_reaction_equation with full KEGG equation grammar coverage

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 3: `_is_industrially_reversible`

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestIsIndustriallyReversible:
    def test_kred_yes(self):
        assert route_predictor._is_industrially_reversible(["1.1.1.184"]) is True

    def test_transaminase_yes(self):
        assert route_predictor._is_industrially_reversible(["2.6.1.5"]) is True

    def test_phosphatase_no(self):
        assert route_predictor._is_industrially_reversible(["3.1.3.1"]) is False

    def test_empty_list_no(self):
        assert route_predictor._is_industrially_reversible([]) is False

    def test_any_match_returns_true(self):
        assert route_predictor._is_industrially_reversible(["3.1.3.1", "1.1.1.184"]) is True

    def test_specific_bvm_match(self):
        assert route_predictor._is_industrially_reversible(["1.14.13.22"]) is True

    def test_close_but_no_match(self):
        # 1.14.13.21 is NOT in the override list (only .22 is)
        assert route_predictor._is_industrially_reversible(["1.14.13.21"]) is False
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestIsIndustriallyReversible -v
```

Expected: 7 AttributeError failures.

- [ ] **Step 3: Implement `_is_industrially_reversible`**

Append to `ChiraLLM/route_predictor.py`:

```python
def _is_industrially_reversible(ec_numbers: list[str]) -> bool:
    """Returns True if any EC number matches a prefix in INDUSTRIAL_REVERSIBLE_EC_PREFIXES.

    The override list is the wet-lab domain-knowledge contract: enzyme classes that are
    routinely run in the non-physiological direction in industrial biocatalysis (KREDs,
    transaminases, IREDs, EREDs, lipases, BVMOs). For these, the reverse-direction penalty
    in _compute_edge_cost is dropped to ~0.
    """
    return any(
        ec.startswith(prefix)
        for ec in ec_numbers
        for prefix in INDUSTRIAL_REVERSIBLE_EC_PREFIXES
    )
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestIsIndustriallyReversible -v
```

Expected: 7 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add _is_industrially_reversible EC-prefix lookup

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 4: Disk cache helpers

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
import os
import time


class TestDiskCache:
    def test_get_returns_none_for_missing(self, tmp_cache_dir):
        # _disk_cache_get reads CHIRALAI_CACHE_ROOT from env
        result = route_predictor._disk_cache_get("kegg", "R12345")
        assert result is None

    def test_set_then_get_roundtrip(self, tmp_cache_dir):
        route_predictor._disk_cache_set("kegg", "R12345", "raw kegg response")
        assert route_predictor._disk_cache_get("kegg", "R12345") == "raw kegg response"

    def test_creates_subdirectory(self, tmp_cache_dir):
        route_predictor._disk_cache_set("equilibrator", "R99999", "data")
        assert (tmp_cache_dir / "equilibrator" / "R99999.cache").exists()

    def test_expired_returns_none(self, tmp_cache_dir, monkeypatch):
        # Force TTL to 1 second for the test
        monkeypatch.setenv("CHIRALAI_CACHE_TTL_DAYS", "0")  # 0 days = always expired
        route_predictor._disk_cache_set("kegg", "R12345", "stale")
        assert route_predictor._disk_cache_get("kegg", "R12345") is None

    def test_unicode_content_preserved(self, tmp_cache_dir):
        route_predictor._disk_cache_set("kegg", "R12345", "alpha-α-ketoglutarate")
        assert route_predictor._disk_cache_get("kegg", "R12345") == "alpha-α-ketoglutarate"
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestDiskCache -v
```

Expected: 5 AttributeError failures.

- [ ] **Step 3: Implement disk cache helpers**

Append to `ChiraLLM/route_predictor.py`:

```python
import os
import time
from pathlib import Path


def _cache_root() -> Path:
    """Returns the disk cache root, honoring CHIRALAI_CACHE_ROOT env override."""
    override = os.environ.get("CHIRALAI_CACHE_ROOT")
    if override:
        return Path(override)
    return Path.home() / ".cache" / "chiralai"


def _cache_ttl_seconds() -> int:
    """Returns the cache TTL in seconds; CHIRALAI_CACHE_TTL_DAYS override available."""
    return int(os.environ.get("CHIRALAI_CACHE_TTL_DAYS", "30")) * 86400


def _disk_cache_get(category: str, key: str) -> str | None:
    """Returns cached content as a string, or None if missing/expired.

    category: one of 'kegg', 'kegg_mol', 'equilibrator'.
    key: the resource identifier (compound ID, reaction ID).
    """
    path = _cache_root() / category / f"{key}.cache"
    if not path.exists():
        return None
    age_seconds = time.time() - path.stat().st_mtime
    if age_seconds > _cache_ttl_seconds():
        return None
    return path.read_text(encoding="utf-8")


def _disk_cache_set(category: str, key: str, content: str) -> None:
    """Writes content to disk cache, creating the subdirectory if needed."""
    path = _cache_root() / category / f"{key}.cache"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _clear_disk_cache() -> int:
    """Removes all cached content. Returns count of files removed.
    Used by `python -m ChiraLLM.route_predictor --clear-cache`."""
    import shutil
    root = _cache_root()
    if not root.exists():
        return 0
    count = sum(1 for _ in root.rglob("*.cache"))
    shutil.rmtree(root)
    return count
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestDiskCache -v
```

Expected: 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add disk cache helpers with TTL and env-override

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 5: `_fetch_kegg_reaction` and `_fetch_compound_reactions`

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestFetchKeggReaction:
    SAMPLE_KEGG_REACTION_FLAT = """ENTRY       R02472                      Reaction
NAME        2-dehydropantoate:NADP+ 2-oxidoreductase
DEFINITION  (R)-Pantoate + NADP+ <=> 2-Dehydropantoate + NADPH + H+
EQUATION    C00966 + C00006 <=> C00890 + C00005 + C00080
ENZYME      1.1.1.169
///
"""

    def test_parses_well_formed_response(self, tmp_cache_dir, mocker):
        mock_resp = mocker.Mock(status_code=200, text=self.SAMPLE_KEGG_REACTION_FLAT)
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)
        # Clear lru_cache before the test
        route_predictor._fetch_kegg_reaction.cache_clear()

        result = route_predictor._fetch_kegg_reaction("R02472")

        assert result is not None
        assert result["rxn_id"] == "R02472"
        assert result["ec_numbers"] == ["1.1.1.169"]
        assert result["direction"] == "reversible"
        # C00966 is on substrate side per KEGG; (R)-pantoate is upstream
        assert (1, "C00966") in result["substrates"]
        assert (1, "C00890") in result["products"]

    def test_returns_none_on_404(self, tmp_cache_dir, mocker):
        mock_resp = mocker.Mock(status_code=404, text="")
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)
        route_predictor._fetch_kegg_reaction.cache_clear()

        result = route_predictor._fetch_kegg_reaction("R99999")
        assert result is None

    def test_disk_cache_hit_skips_network(self, tmp_cache_dir, mocker):
        # Pre-populate disk cache
        route_predictor._disk_cache_set("kegg", "R02472", self.SAMPLE_KEGG_REACTION_FLAT)
        # Network call should not happen
        get_mock = mocker.patch("ChiraLLM.route_predictor.requests.get")
        route_predictor._fetch_kegg_reaction.cache_clear()

        result = route_predictor._fetch_kegg_reaction("R02472")

        assert result is not None
        assert result["ec_numbers"] == ["1.1.1.169"]
        get_mock.assert_not_called()

    def test_network_failure_retries_then_returns_none(self, tmp_cache_dir, mocker):
        import requests
        mocker.patch(
            "ChiraLLM.route_predictor.requests.get",
            side_effect=requests.RequestException("connection refused"),
        )
        sleep_mock = mocker.patch("ChiraLLM.route_predictor.time.sleep")
        route_predictor._fetch_kegg_reaction.cache_clear()

        result = route_predictor._fetch_kegg_reaction("R02472")

        assert result is None
        sleep_mock.assert_called_once_with(1.0)  # single retry after 1s


class TestFetchCompoundReactions:
    SAMPLE_KEGG_COMPOUND_FLAT = """ENTRY       C00599                      Compound
NAME        (R)-Pantolactone
FORMULA     C6H10O3
REACTION    R02472 R09096
ENZYME      1.1.1.169
///
"""

    def test_extracts_reaction_ids(self, tmp_cache_dir, mocker):
        mock_resp = mocker.Mock(status_code=200, text=self.SAMPLE_KEGG_COMPOUND_FLAT)
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)
        route_predictor._fetch_compound_reactions.cache_clear()

        result = route_predictor._fetch_compound_reactions("C00599")

        assert result == ["R02472", "R09096"]

    def test_empty_reaction_field_returns_empty_list(self, tmp_cache_dir, mocker):
        no_reactions = "ENTRY       C99999                      Compound\nNAME        Foo\n///\n"
        mock_resp = mocker.Mock(status_code=200, text=no_reactions)
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)
        route_predictor._fetch_compound_reactions.cache_clear()

        assert route_predictor._fetch_compound_reactions("C99999") == []
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestFetchKeggReaction tests/test_route_predictor_unit.py::TestFetchCompoundReactions -v
```

Expected: AttributeError on every test.

- [ ] **Step 3: Implement `_fetch_kegg_reaction` and `_fetch_compound_reactions`**

Append to `ChiraLLM/route_predictor.py`:

```python
import functools
import requests

KEGG_REST_BASE = "http://rest.kegg.jp/get"
KEGG_RETRY_DELAY_SECONDS = 1.0


def _parse_kegg_reaction_flat(text: str, rxn_id: str) -> dict | None:
    """Parses a KEGG reaction flat-file response into the standard reaction dict."""
    equation = None
    ec_numbers: list[str] = []
    for line in text.splitlines():
        if line.startswith("EQUATION"):
            equation = line[12:].strip()
        elif line.startswith("ENZYME"):
            ec_numbers.extend(line[12:].split())
        elif line.startswith("///"):
            break
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
    """Extracts the REACTION field as a list of R##### IDs."""
    reactions: list[str] = []
    in_reaction_field = False
    for line in text.splitlines():
        if line.startswith("REACTION"):
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
    """Fetches and parses a KEGG reaction. Returns None on 404 or unparseable response.
    Uses LRU + disk cache. Single retry with 1s delay on network failure, then None.
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
        # 5xx: retry once
        if attempt == 0:
            time.sleep(KEGG_RETRY_DELAY_SECONDS)
            continue
        logger.warning("KEGG returned %d for %s after retry", resp.status_code, rxn_id)
        return None
    return None


@functools.lru_cache(maxsize=4096)
def _fetch_compound_reactions(compound_id: str) -> list[str]:
    """Fetches a KEGG compound's REACTION field as a list of R##### IDs.
    Returns empty list if the compound has no reactions or fetch fails.
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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestFetchKeggReaction tests/test_route_predictor_unit.py::TestFetchCompoundReactions -v
```

Expected: 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add KEGG reaction and compound fetchers with retry + cache

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 6: `_fetch_kegg_mol`

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestFetchKeggMol:
    # Minimal valid MOL block for methane (a real RDKit-parseable example)
    SAMPLE_MOL = """methane
  Mrv0541 01010100002D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""

    def test_parses_valid_mol(self, tmp_cache_dir, mocker):
        mock_resp = mocker.Mock(status_code=200, text=self.SAMPLE_MOL)
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)
        route_predictor._fetch_kegg_mol.cache_clear()

        mol = route_predictor._fetch_kegg_mol("C01438")

        assert mol is not None
        assert mol.GetNumAtoms() == 1

    def test_returns_none_on_404(self, tmp_cache_dir, mocker):
        mock_resp = mocker.Mock(status_code=404, text="")
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)
        route_predictor._fetch_kegg_mol.cache_clear()

        assert route_predictor._fetch_kegg_mol("C99999") is None

    def test_returns_none_on_empty_body(self, tmp_cache_dir, mocker):
        mock_resp = mocker.Mock(status_code=200, text="")
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)
        route_predictor._fetch_kegg_mol.cache_clear()

        assert route_predictor._fetch_kegg_mol("C00001") is None

    def test_returns_none_on_unparseable_mol(self, tmp_cache_dir, mocker):
        mock_resp = mocker.Mock(status_code=200, text="not a mol file")
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)
        route_predictor._fetch_kegg_mol.cache_clear()

        assert route_predictor._fetch_kegg_mol("C00001") is None
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestFetchKeggMol -v
```

Expected: 4 AttributeError failures.

- [ ] **Step 3: Implement `_fetch_kegg_mol`**

Append to `ChiraLLM/route_predictor.py`:

```python
from rdkit import Chem


@functools.lru_cache(maxsize=4096)
def _fetch_kegg_mol(compound_id: str):
    """Fetches a compound's MOL file from KEGG and parses to an RDKit Mol.
    Returns None if compound has no MOL file, KEGG returns 404, or RDKit parse fails.
    """
    cached = _disk_cache_get("kegg_mol", compound_id)
    if cached is not None:
        if not cached.strip():
            return None
        try:
            mol = Chem.MolFromMolBlock(cached)
            return mol
        except Exception:
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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestFetchKeggMol -v
```

Expected: 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add _fetch_kegg_mol with RDKit parse + disk cache

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 7: `_fetch_delta_g_kj_per_mol` (eQuilibrator)

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestFetchDeltaG:
    def test_parses_well_formed_response(self, tmp_cache_dir, mocker):
        mock_resp = mocker.Mock(
            status_code=200,
            json=lambda: {"standard_dg_prime": -29.4, "units": "kJ/mol"},
        )
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)

        result = route_predictor._fetch_delta_g_kj_per_mol("R02472")

        assert result == -29.4

    def test_returns_none_on_network_failure(self, tmp_cache_dir, mocker):
        import requests
        mocker.patch(
            "ChiraLLM.route_predictor.requests.get",
            side_effect=requests.RequestException("eQuilibrator down"),
        )

        assert route_predictor._fetch_delta_g_kj_per_mol("R02472") is None

    def test_returns_none_on_5xx(self, tmp_cache_dir, mocker):
        mock_resp = mocker.Mock(status_code=503, text="")
        mocker.patch("ChiraLLM.route_predictor.requests.get", return_value=mock_resp)

        assert route_predictor._fetch_delta_g_kj_per_mol("R02472") is None

    def test_disk_cache_hit_skips_network(self, tmp_cache_dir, mocker):
        route_predictor._disk_cache_set("equilibrator", "R02472", "-29.4")
        get_mock = mocker.patch("ChiraLLM.route_predictor.requests.get")

        result = route_predictor._fetch_delta_g_kj_per_mol("R02472")

        assert result == -29.4
        get_mock.assert_not_called()
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestFetchDeltaG -v
```

Expected: 4 AttributeError failures.

- [ ] **Step 3: Implement `_fetch_delta_g_kj_per_mol`**

Append to `ChiraLLM/route_predictor.py`:

```python
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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestFetchDeltaG -v
```

Expected: 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add _fetch_delta_g_kj_per_mol with graceful network-failure handling

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 8: Heuristic — `_get_central_fingerprints` + `_tanimoto_to_central`

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestHeuristic:
    def test_central_fingerprints_lazy_init(self, mocker):
        # Reset module state
        route_predictor._central_fingerprints_cache = None

        # Mock _fetch_kegg_mol to return a real RDKit mol for any compound id
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CC(=O)C(=O)O")  # pyruvate
        mocker.patch("ChiraLLM.route_predictor._fetch_kegg_mol", return_value=mol)

        fps = route_predictor._get_central_fingerprints()

        assert isinstance(fps, dict)
        assert len(fps) > 30  # CENTRAL_METABOLITES has ~40 entries

    def test_tanimoto_to_central_zero_for_central_compound(self, mocker):
        from rdkit import Chem
        # Compound IS pyruvate; nearest central should be pyruvate itself, distance ≈ 0
        pyruvate = Chem.MolFromSmiles("CC(=O)C(=O)O")
        mocker.patch("ChiraLLM.route_predictor._fetch_kegg_mol", return_value=pyruvate)
        route_predictor._central_fingerprints_cache = None

        distance = route_predictor._tanimoto_to_central("C00022")

        assert distance < 0.1  # very close

    def test_tanimoto_to_central_positive_for_distant_compound(self, mocker):
        from rdkit import Chem
        # Use a complex molecule (caffeine) very different from central metabolites
        caffeine = Chem.MolFromSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")

        def fake_fetch(cid):
            if cid in route_predictor.CENTRAL_METABOLITES:
                return Chem.MolFromSmiles("CC(=O)C(=O)O")  # pyruvate stand-in
            return caffeine

        mocker.patch("ChiraLLM.route_predictor._fetch_kegg_mol", side_effect=fake_fetch)
        route_predictor._central_fingerprints_cache = None

        distance = route_predictor._tanimoto_to_central("C07481")  # caffeine KEGG ID

        assert distance > 0.5  # very distant

    def test_tanimoto_falls_back_to_uniform_when_mol_missing(self, mocker):
        # Mock _fetch_kegg_mol to return None for the query, real mol for centrals
        from rdkit import Chem

        def fake_fetch(cid):
            if cid == "C99999":
                return None
            return Chem.MolFromSmiles("CC(=O)C(=O)O")

        mocker.patch("ChiraLLM.route_predictor._fetch_kegg_mol", side_effect=fake_fetch)
        route_predictor._central_fingerprints_cache = None

        distance = route_predictor._tanimoto_to_central("C99999")
        assert distance == 0.5  # uniform fallback per spec §5.1 case #11
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestHeuristic -v
```

Expected: 4 AttributeError failures.

- [ ] **Step 3: Implement the heuristic**

Append to `ChiraLLM/route_predictor.py`:

```python
from rdkit.Chem import rdFingerprintGenerator, DataStructs

_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=False)
_central_fingerprints_cache: dict | None = None
_TANIMOTO_FALLBACK_DISTANCE = 0.5


def _get_central_fingerprints() -> dict:
    """Lazy-initialized Morgan fingerprints for every compound in CENTRAL_METABOLITES.
    Returns dict mapping compound_id -> ExplicitBitVect. Computed once per module lifetime.
    Lazy init defers ~50ms × N RDKit calls until first search, supports future catalog swaps.
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
    0 means 'already at a central metabolite', 1 means 'maximally distant'.
    Falls back to uniform 0.5 if the compound's MOL file is unavailable (per spec §5.1 #11).
    """
    mol = _fetch_kegg_mol(compound_id)
    if mol is None:
        return _TANIMOTO_FALLBACK_DISTANCE
    fp = _morgan_gen.GetFingerprint(mol)
    centrals = _get_central_fingerprints()
    if not centrals:
        return _TANIMOTO_FALLBACK_DISTANCE
    max_sim = max(DataStructs.TanimotoSimilarity(fp, central_fp) for central_fp in centrals.values())
    return 1.0 - max_sim
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestHeuristic -v
```

Expected: 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add Morgan-fingerprint Tanimoto heuristic with lazy init + fallback

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 9: `_compute_edge_cost`

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestComputeEdgeCost:
    REVERSIBLE_RXN = {
        "rxn_id": "R_test_rev",
        "ec_numbers": ["1.1.1.184"],  # KRED — industrially reversible
        "direction": "reversible",
    }

    FORWARD_ONLY_RXN = {
        "rxn_id": "R_test_fwd",
        "ec_numbers": ["3.1.3.1"],  # phosphatase, NOT industrially reversible
        "direction": "forward_only",
    }

    def test_forward_traversal_no_directionality_penalty(self, mocker):
        mocker.patch(
            "ChiraLLM.route_predictor._fetch_delta_g_kj_per_mol", return_value=-10.0
        )
        cost = route_predictor._compute_edge_cost(self.REVERSIBLE_RXN, "forward")

        assert cost["directionality"] == 0.0
        assert cost["industrial_override"] == 0.0
        assert cost["base"] == 1.0
        assert cost["total"] == cost["base"] + cost["thermodynamic"] + cost["directionality"] + cost["industrial_override"]

    def test_reverse_traversal_pays_directionality_for_irreversible(self, mocker):
        mocker.patch(
            "ChiraLLM.route_predictor._fetch_delta_g_kj_per_mol", return_value=-30.0
        )
        cost = route_predictor._compute_edge_cost(self.FORWARD_ONLY_RXN, "reverse")

        assert cost["directionality"] > 0
        assert cost["thermodynamic"] > 0
        assert cost["industrial_override"] == 0.0  # phosphatase not in override list
        assert cost["total"] > 1.0

    def test_industrial_override_zeroes_directionality_for_kred_reverse(self, mocker):
        mocker.patch(
            "ChiraLLM.route_predictor._fetch_delta_g_kj_per_mol", return_value=-10.0
        )
        cost = route_predictor._compute_edge_cost(self.REVERSIBLE_RXN, "reverse")

        # KRED in reverse: industrial_override should reduce/zero the directionality penalty
        assert cost["industrial_override"] < 0  # discount is negative
        assert cost["directionality"] >= 0
        # net effect: total should be close to forward cost
        forward_cost = route_predictor._compute_edge_cost(self.REVERSIBLE_RXN, "forward")
        assert abs(cost["total"] - forward_cost["total"]) < 0.5

    def test_fallback_delta_g_when_unreachable(self, mocker):
        mocker.patch("ChiraLLM.route_predictor._fetch_delta_g_kj_per_mol", return_value=None)
        cost = route_predictor._compute_edge_cost(self.FORWARD_ONLY_RXN, "reverse")

        # Should not crash; thermodynamic component computed from FALLBACK_DELTA_G_KJ
        assert cost["thermodynamic"] > 0
        assert "total" in cost

    def test_breakdown_components_sum_to_total(self, mocker):
        mocker.patch(
            "ChiraLLM.route_predictor._fetch_delta_g_kj_per_mol", return_value=-15.0
        )
        cost = route_predictor._compute_edge_cost(self.REVERSIBLE_RXN, "reverse")

        component_sum = (
            cost["base"] + cost["thermodynamic"] + cost["directionality"] + cost["industrial_override"]
        )
        assert abs(cost["total"] - component_sum) < 1e-9
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestComputeEdgeCost -v
```

Expected: 5 AttributeError failures.

- [ ] **Step 3: Implement `_compute_edge_cost`**

Append to `ChiraLLM/route_predictor.py`:

```python
_BASE_EDGE_COST = 1.0
_REVERSIBLE_REVERSE_PENALTY = 0.5  # reverse traversal of a KEGG-marked-reversible reaction
_IRREVERSIBLE_REVERSE_PENALTY = 2.0  # reverse traversal of a KEGG-marked-forward-only reaction


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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestComputeEdgeCost -v
```

Expected: 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add _compute_edge_cost with transparent thermo + directionality + override breakdown

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 10: Dataclasses + `_astar_search`

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestAstarSearch:
    def test_finds_route_in_synthetic_graph(self, mock_kegg, mock_equilibrator, mocker):
        # Override CENTRAL_METABOLITES for this test to use synthetic precursors
        mocker.patch.object(
            route_predictor,
            "CENTRAL_METABOLITES",
            {"C_PRECURSOR_A": "synthetic A", "C_PRECURSOR_B": "synthetic B"},
        )
        # Bypass real Tanimoto computation (no real MOL files in fixture)
        mocker.patch("ChiraLLM.route_predictor._tanimoto_to_central", return_value=0.5)

        result = route_predictor._astar_search("C_TARGET", budget=50, depth_cap=5)

        assert isinstance(result, dict)
        assert result["nodes_explored"] > 0
        assert result["budget_exhausted"] is False
        # At least one route should reach a synthetic central metabolite
        leaf_ids = result["leaf_ids"]
        assert any(leaf in {"C_PRECURSOR_A", "C_PRECURSOR_B"} for leaf in leaf_ids)

    def test_budget_exhaustion(self, mock_kegg, mock_equilibrator, mocker):
        mocker.patch.object(
            route_predictor,
            "CENTRAL_METABOLITES",
            {"C_NEVER_REACHED": "unreachable"},
        )
        mocker.patch("ChiraLLM.route_predictor._tanimoto_to_central", return_value=0.5)

        result = route_predictor._astar_search("C_TARGET", budget=3, depth_cap=5)

        assert result["budget_exhausted"] is True
        assert result["nodes_explored"] == 3

    def test_target_with_no_reactions_returns_empty(self, mock_kegg, mock_equilibrator, mocker):
        mocker.patch("ChiraLLM.route_predictor._tanimoto_to_central", return_value=0.5)
        result = route_predictor._astar_search("C_DEAD_END", budget=50, depth_cap=5)

        assert result["leaf_ids"] == []
        assert result["nodes_explored"] >= 1


class TestRouteDataclasses:
    def test_route_step_construction(self):
        step = route_predictor.RouteStep(
            reaction_id="R_test",
            ec_numbers=["1.1.1.184"],
            precursor_id="C_A",
            intermediate_id="C_B",
            edge_cost_breakdown={"base": 1.0, "thermodynamic": 0.5, "directionality": 0.0, "industrial_override": 0.0, "total": 1.5},
            traversed_direction="forward",
        )
        assert step.precursor_id == "C_A"
        assert step.intermediate_id == "C_B"

    def test_route_construction(self):
        step = route_predictor.RouteStep(
            reaction_id="R_test", ec_numbers=["1.1.1.1"],
            precursor_id="C_A", intermediate_id="C_B",
            edge_cost_breakdown={"base": 1.0, "thermodynamic": 0.0, "directionality": 0.0, "industrial_override": 0.0, "total": 1.0},
            traversed_direction="forward",
        )
        route = route_predictor.Route(
            target_id="C_B", steps=[step],
            terminal_precursor_id="C_A", terminal_precursor_name="alpha",
            total_cost=1.0,
            cost_breakdown={"base": 1.0, "thermodynamic": 0.0, "directionality": 0.0, "industrial_override": 0.0},
            warnings=[],
        )
        assert route.terminal_precursor_id == "C_A"
        assert len(route.steps) == 1
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestAstarSearch tests/test_route_predictor_unit.py::TestRouteDataclasses -v
```

Expected: 5 AttributeError failures.

- [ ] **Step 3: Implement dataclasses + `_astar_search`**

Append to `ChiraLLM/route_predictor.py`:

```python
import heapq
from dataclasses import dataclass, field


@dataclass
class RouteStep:
    reaction_id: str
    ec_numbers: list[str]
    precursor_id: str           # the upstream compound (one step closer to terminal)
    intermediate_id: str        # the downstream compound (we just came from this)
    edge_cost_breakdown: dict[str, float]
    traversed_direction: str    # 'forward' or 'reverse'
    # Naming note: "substrate" and "product" would be ambiguous about graph direction in
    # retrosynthesis. precursor_id (upstream) and intermediate_id (downstream) are
    # unambiguous about graph position regardless of how the chemical reaction is described.


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


def _astar_search(target_id: str, budget: int = DEFAULT_BUDGET, depth_cap: int = DEFAULT_DEPTH_CAP) -> dict:
    """Best-first backward search from target. Returns a DAG and search metadata.

    Algorithm: weighted A* / greedy best-first.
      f(n) = g(n) + h(n)
      g(n) = cumulative edge cost from target to current node
      h(n) = TANIMOTO_HEURISTIC_WEIGHT * Tanimoto distance to nearest central metabolite

    Termination: priority queue empty, OR nodes_explored >= budget, OR (per-node) depth > depth_cap.
    Search continues past the first central-metabolite hit so that diverse alternates can be found.
    """
    # Priority queue items: (f_score, tie_breaker, depth, compound_id, g_score)
    counter = 0  # tie-breaker for heap stability
    initial_h = TANIMOTO_HEURISTIC_WEIGHT * _tanimoto_to_central(target_id)
    queue: list = [(initial_h, counter, 0, target_id, 0.0)]
    heapq.heapify(queue)

    # visited_dag[compound_id] = list of incoming edges
    # Each edge: {"parent_id": str, "reaction": dict, "edge_cost": dict, "depth": int, "g_score": float}
    visited_dag: dict[str, list] = {}
    leaf_ids: list[str] = []
    nodes_explored = 0
    budget_exhausted = False

    # Track best known g_score per compound to prune higher-cost revisits
    best_g_seen: dict[str, float] = {target_id: 0.0}

    while queue:
        if nodes_explored >= budget:
            budget_exhausted = True
            break

        f_score, _, depth, compound_id, g_score = heapq.heappop(queue)

        # Skip if a better path to this node has been found since enqueue
        if g_score > best_g_seen.get(compound_id, float("inf")):
            continue

        nodes_explored += 1

        # Terminal check: central metabolite reached
        if compound_id in CENTRAL_METABOLITES:
            if compound_id not in leaf_ids:
                leaf_ids.append(compound_id)
            continue  # don't expand past central metabolites

        if depth >= depth_cap:
            continue

        # Expand: fetch reactions for this compound
        reaction_ids = _fetch_compound_reactions(compound_id)

        for rxn_id in reaction_ids:
            reaction = _fetch_kegg_reaction(rxn_id)
            if reaction is None:
                continue

            # Determine which side compound_id is on, and what the precursors are.
            substrate_ids = [cid for _, cid in reaction["substrates"]]
            product_ids = [cid for _, cid in reaction["products"]]

            if compound_id in product_ids:
                # We came from the product side; precursors are the substrates.
                traversed_direction = "forward"
                precursor_candidates = substrate_ids
            elif compound_id in substrate_ids:
                # We came from the substrate side; precursors are the products (reverse traversal).
                traversed_direction = "reverse"
                precursor_candidates = product_ids
            else:
                # Compound not actually in this reaction (shouldn't happen with KEGG data).
                continue

            edge_cost = _compute_edge_cost(reaction, traversed_direction)
            new_g = g_score + edge_cost["total"]

            for precursor_id in precursor_candidates:
                if precursor_id == compound_id:
                    continue  # self-loop (cofactor on both sides)
                # Skip cofactors that are common throughout metabolism — keep the search focused
                if precursor_id in {"C00080", "C00001", "C00007"}:  # H+, H2O, O2
                    continue

                if new_g >= best_g_seen.get(precursor_id, float("inf")):
                    continue  # known better path

                best_g_seen[precursor_id] = new_g
                visited_dag.setdefault(precursor_id, []).append({
                    "parent_id": compound_id,
                    "reaction": reaction,
                    "edge_cost": edge_cost,
                    "depth": depth + 1,
                    "g_score": new_g,
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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestAstarSearch tests/test_route_predictor_unit.py::TestRouteDataclasses -v
```

Expected: 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add Route/RouteStep/RouteResult dataclasses + weighted A* search

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 11: Output formatters

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestFormatters:
    @pytest.fixture
    def hand_built_dag(self, mocker):
        """A 2-route DAG: target → A → CENTRAL_1, target → B → CENTRAL_2."""
        mocker.patch.object(
            route_predictor,
            "CENTRAL_METABOLITES",
            {"C_CENTRAL_1": "alpha", "C_CENTRAL_2": "beta"},
        )
        rxn_1 = {"rxn_id": "R1", "ec_numbers": ["1.1.1.1"], "direction": "reversible"}
        rxn_2 = {"rxn_id": "R2", "ec_numbers": ["1.1.1.2"], "direction": "reversible"}
        rxn_3 = {"rxn_id": "R3", "ec_numbers": ["2.6.1.1"], "direction": "reversible"}
        rxn_4 = {"rxn_id": "R4", "ec_numbers": ["3.1.3.1"], "direction": "reversible"}
        cost_low = {"base": 1.0, "thermodynamic": 0.0, "directionality": 0.0, "industrial_override": 0.0, "total": 1.0}
        cost_high = {"base": 1.0, "thermodynamic": 1.0, "directionality": 0.5, "industrial_override": 0.0, "total": 2.5}
        return {
            "target_id": "C_TARGET",
            "visited_dag": {
                "C_A": [{"parent_id": "C_TARGET", "reaction": rxn_1, "edge_cost": cost_low, "depth": 1, "g_score": 1.0}],
                "C_B": [{"parent_id": "C_TARGET", "reaction": rxn_2, "edge_cost": cost_high, "depth": 1, "g_score": 2.5}],
                "C_CENTRAL_1": [{"parent_id": "C_A", "reaction": rxn_3, "edge_cost": cost_low, "depth": 2, "g_score": 2.0}],
                "C_CENTRAL_2": [{"parent_id": "C_B", "reaction": rxn_4, "edge_cost": cost_low, "depth": 2, "g_score": 3.5}],
            },
            "leaf_ids": ["C_CENTRAL_1", "C_CENTRAL_2"],
            "nodes_explored": 5,
            "budget_exhausted": False,
        }

    def test_extract_top_n_returns_routes_sorted_by_cost(self, hand_built_dag):
        routes = route_predictor._extract_top_n(hand_built_dag, n=2)

        assert len(routes) == 2
        assert routes[0].total_cost <= routes[1].total_cost
        assert routes[0].terminal_precursor_id == "C_CENTRAL_1"  # cheaper

    def test_extract_top_n_respects_n(self, hand_built_dag):
        routes = route_predictor._extract_top_n(hand_built_dag, n=1)
        assert len(routes) == 1

    def test_extract_full_tree_returns_one_route_with_full_dag(self, hand_built_dag):
        routes = route_predictor._extract_full_tree(hand_built_dag)
        assert len(routes) == 1
        # full_tree mode packs all leaves into a single route's notes/structure
        assert routes[0].target_id == "C_TARGET"

    def test_extract_shortest_plus_diverse_returns_unique_terminals(self, hand_built_dag):
        routes = route_predictor._extract_shortest_plus_diverse(hand_built_dag, n_diverse=1)
        terminals = {r.terminal_precursor_id for r in routes}
        assert len(routes) >= 1
        assert len(terminals) == len(routes)  # all unique

    def test_route_steps_in_target_to_precursor_order(self, hand_built_dag):
        routes = route_predictor._extract_top_n(hand_built_dag, n=1)
        route = routes[0]
        # First step's intermediate_id should be the target
        assert route.steps[0].intermediate_id == "C_TARGET"
        # Last step's precursor_id should be the terminal
        assert route.steps[-1].precursor_id == route.terminal_precursor_id
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestFormatters -v
```

Expected: 5 AttributeError failures.

- [ ] **Step 3: Implement the three formatters**

Append to `ChiraLLM/route_predictor.py`:

```python
def _backtrack_route(dag_result: dict, leaf_id: str) -> Route | None:
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
        # Pick the cheapest incoming edge (lowest g_score)
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
            traversed_direction="forward",  # direction info is in reaction; orient relative to graph
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


def _extract_top_n(dag_result: dict, n: int = DEFAULT_MAX_ROUTES) -> list[Route]:
    """Returns up to n routes, sorted by total_cost ascending."""
    routes = []
    for leaf_id in dag_result["leaf_ids"]:
        route = _backtrack_route(dag_result, leaf_id)
        if route is not None:
            routes.append(route)
    routes.sort(key=lambda r: r.total_cost)
    return routes[:n]


def _extract_full_tree(dag_result: dict) -> list[Route]:
    """Returns a single Route whose 'steps' encode the complete DAG.
    The first leaf is the canonical terminal; the full DAG is preserved in warnings as JSON.
    For consumers who need the actual tree, use the JSON sidecar route_prediction field.
    """
    if not dag_result["leaf_ids"]:
        return []
    # Pick the cheapest leaf as the canonical route
    cheapest_leaf = min(
        dag_result["leaf_ids"],
        key=lambda lid: min(e["g_score"] for e in dag_result["visited_dag"].get(lid, [{"g_score": float("inf")}])),
    )
    canonical = _backtrack_route(dag_result, cheapest_leaf)
    if canonical is None:
        return []
    # Annotate that the full DAG has more leaves
    other_leaves = [lid for lid in dag_result["leaf_ids"] if lid != cheapest_leaf]
    if other_leaves:
        canonical.warnings.append(
            f"full_tree mode: {len(other_leaves)} additional leaves not shown in steps "
            f"(see JSON sidecar for full DAG): {', '.join(other_leaves)}"
        )
    return [canonical]


def _extract_shortest_plus_diverse(dag_result: dict, n_diverse: int = 2) -> list[Route]:
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

    # Shortest by step count (tie-break by cost)
    shortest = min(all_routes, key=lambda r: (len(r.steps), r.total_cost))
    selected = [shortest]
    seen_terminals = {shortest.terminal_precursor_id}

    # Add maximally-different alternates: prioritize different terminal precursors
    remaining = [r for r in all_routes if r.terminal_precursor_id not in seen_terminals]
    remaining.sort(key=lambda r: r.total_cost)
    for r in remaining[:n_diverse]:
        selected.append(r)
        seen_terminals.add(r.terminal_precursor_id)

    return selected
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestFormatters -v
```

Expected: 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add three output formatters (top_n, full_tree, shortest_plus_diverse)

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 12: `predict_route` public API + failure modes

**Files:**
- Modify: `ChiraLLM/route_predictor.py`
- Modify: `tests/test_route_predictor_unit.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestPredictRoute:
    def test_invalid_kegg_id_format(self):
        result = route_predictor.predict_route("not_a_kegg_id")
        assert result.status == "invalid_kegg_id"
        assert result.routes == []

    def test_none_compound_id(self):
        result = route_predictor.predict_route(None)
        assert result.status == "no_kegg_id"

    def test_empty_compound_id(self):
        result = route_predictor.predict_route("")
        assert result.status == "no_kegg_id"

    def test_invalid_mode(self):
        result = route_predictor.predict_route("C00599", mode="bogus")
        assert result.status == "invalid_mode"

    def test_target_with_no_reactions(self, mock_kegg, mock_equilibrator, mocker):
        mocker.patch("ChiraLLM.route_predictor._tanimoto_to_central", return_value=0.5)
        mocker.patch.object(
            route_predictor, "CENTRAL_METABOLITES",
            {"C_PRECURSOR_A": "alpha"},
        )

        result = route_predictor.predict_route("C_DEAD_END", mode="top_n")
        assert result.status == "target_has_no_reactions"

    def test_successful_top_n_search(self, mock_kegg, mock_equilibrator, mocker):
        mocker.patch.object(
            route_predictor, "CENTRAL_METABOLITES",
            {"C_PRECURSOR_A": "alpha", "C_PRECURSOR_B": "beta"},
        )
        mocker.patch("ChiraLLM.route_predictor._tanimoto_to_central", return_value=0.5)

        result = route_predictor.predict_route("C_TARGET", mode="top_n", n=2, budget=50)

        assert result.status == "success"
        assert len(result.routes) >= 1
        assert all(r.terminal_precursor_id in {"C_PRECURSOR_A", "C_PRECURSOR_B"} for r in result.routes)
        assert result.nodes_explored > 0

    def test_full_tree_mode(self, mock_kegg, mock_equilibrator, mocker):
        mocker.patch.object(
            route_predictor, "CENTRAL_METABOLITES",
            {"C_PRECURSOR_A": "alpha", "C_PRECURSOR_B": "beta"},
        )
        mocker.patch("ChiraLLM.route_predictor._tanimoto_to_central", return_value=0.5)

        result = route_predictor.predict_route("C_TARGET", mode="full_tree", budget=50)

        assert result.status == "success"
        assert result.mode == "full_tree"
        assert len(result.routes) == 1  # full_tree always returns one canonical route

    def test_no_route_found_when_all_branches_dead_end(self, mock_kegg, mock_equilibrator, mocker):
        mocker.patch.object(
            route_predictor, "CENTRAL_METABOLITES",
            {"C_NEVER_REACHED": "unreachable"},
        )
        mocker.patch("ChiraLLM.route_predictor._tanimoto_to_central", return_value=0.5)

        result = route_predictor.predict_route("C_TARGET", mode="top_n", budget=50)

        assert result.status == "no_route_found"
        assert result.routes == []
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestPredictRoute -v
```

Expected: 8 AttributeError failures.

- [ ] **Step 3: Implement `predict_route`**

Append to `ChiraLLM/route_predictor.py`:

```python
_KEGG_ID_RE = re.compile(r"^C\d{5}$")
_VALID_MODES = {"top_n", "full_tree", "shortest_plus_diverse"}


def predict_route(
    compound_id: str | None,
    mode: str = "top_n",
    n: int = DEFAULT_MAX_ROUTES,
    budget: int = DEFAULT_BUDGET,
) -> RouteResult:
    """Predicts biosynthetic routes from a target KEGG compound to central metabolites.

    Returns a RouteResult with status indicating success or specific failure mode.
    Never raises — all failures encoded in result.status and result.warnings.
    """
    # Input validation
    if not compound_id:
        return RouteResult(target_id="", mode=mode, routes=[], nodes_explored=0,
                           budget_exhausted=False, status="no_kegg_id")
    if not _KEGG_ID_RE.match(compound_id):
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
        # Distinguish 'compound exists but has no reactions' from 'compound not in KEGG'
        # by checking whether the compound itself has a KEGG entry.
        # _fetch_compound_reactions already returned [], so KEGG either had no REACTION field
        # or returned non-200. Use a 1-extra HEAD-style check via cache hit.
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

    # Dispatch to formatter
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
        n = _clear_disk_cache()
        print(f"Cleared {n} cached files from {_cache_root()}")
        sys.exit(0)
    print("Usage: python -m ChiraLLM.route_predictor --clear-cache")
    sys.exit(1)
```

- [ ] **Step 4: Run all unit tests + check coverage**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py -v --cov=ChiraLLM.route_predictor --cov-report=term-missing
```

Expected: All tests PASS. Coverage ≥85%.

- [ ] **Step 5: Commit**

```bash
git add ChiraLLM/route_predictor.py tests/test_route_predictor_unit.py
git commit -m "add predict_route public API + cache CLI entry point

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 13: Wire `predict_route` into `main.py` (PARALLELIZABLE with 14 and 15)

**Files:**
- Modify: `main.py`

- [ ] **Step 1: Read the existing main.py suggestion loop**

```bash
.venv/bin/python -c "import ast, sys; src = open('main.py').read(); print(src)" | head -80
```

Confirm the loop structure starts at the `for suggestion in suggestions:` line.

- [ ] **Step 2: Replace the inner block with the route-aware version**

Replace the existing block in `main.py` from `for suggestion in suggestions:` through the end of the loop with:

```python
    from ChiraLLM.route_predictor import predict_route
    from dataclasses import asdict

    for suggestion in suggestions:
        smiles = suggestion.get("SMILES")
        if smiles:
            suggestion["chirality_validation"] = validate_chirality(smiles)

        compound_id = suggestion.get("KEGG_ID")
        if compound_id:
            kegg = query_kegg(compound_id)
            suggestion["kegg_data"] = kegg

            # NEW: route prediction runs after KEGG validation, before BRENDA enrichment.
            # The full set of ECs across all routes' steps becomes the BRENDA query input.
            route_result = predict_route(compound_id, mode="top_n", n=3)
            suggestion["route_prediction"] = asdict(route_result)

            if route_result.status == "success":
                # Deduplicate ECs across all routes' steps before the BRENDA batch call.
                all_ecs = sorted({
                    ec
                    for route in route_result.routes
                    for step in route.steps
                    for ec in step.ec_numbers
                })
                if all_ecs:
                    suggestion["brenda_data"] = query_enantioselectivity_batch(all_ecs[:20])
                else:
                    suggestion["brenda_data"] = {"status": "no_ec_numbers"}

                # Per-route terminal-precursor feasibility (call site change per spec §2.1).
                suggestion["route_feasibility"] = [
                    {
                        "route_index": i,
                        "terminal_precursor": route.terminal_precursor_id,
                        "feasibility": check_feasibility(route.terminal_precursor_id),
                    }
                    for i, route in enumerate(route_result.routes)
                ]
            else:
                # Fall back to legacy behavior: BRENDA on KEGG's enzymes for the target only.
                ec_numbers = kegg.get("enzymes", []) if kegg.get("status") == "success" else []
                if ec_numbers:
                    suggestion["brenda_data"] = query_enantioselectivity_batch(ec_numbers[:5])
                else:
                    suggestion["brenda_data"] = {"status": "no_ec_numbers"}
                suggestion["route_feasibility"] = [{
                    "route_index": 0,
                    "terminal_precursor": compound_id,
                    "feasibility": check_feasibility(compound_id),
                }]

        suggestion["scoring"] = score_suggestion(suggestion)
```

Also delete the duplicate top-of-file `from ChiraLLM.route_predictor import predict_route` if Step 2 caused it to be imported twice; the canonical import must be at the top of the file alongside the others. Final imports section should look like:

```python
import argparse
import json
import sys
from dataclasses import asdict
from ChiraLLM.query_handler import ask_gpt_chirality
from ChiraLLM.database_validator import query_kegg
from ChiraLLM.chirality_checker import validate_chirality
from ChiraLLM.brenda_client import query_enantioselectivity_batch
from ChiraLLM.feasibility_checker import check_feasibility
from ChiraLLM.enantioselectivity_scorer import score_suggestion
from ChiraLLM.route_predictor import predict_route
from utils.file_saver import save_suggestions_to_csv
```

Remove the in-loop `from ... import ...` lines.

- [ ] **Step 3: Smoke-check that main.py imports cleanly**

```bash
.venv/bin/python -c "import main"
```

Expected: silent success (no import errors).

- [ ] **Step 4: Commit**

```bash
git add main.py
git commit -m "wire predict_route into main pipeline; per-route terminal-precursor feasibility

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 14: Extend `utils/file_saver.py` (PARALLELIZABLE with 13 and 15)

**Files:**
- Modify: `utils/file_saver.py`

- [ ] **Step 1: Add `_flatten_route_prediction` helper and wire it in**

Open `utils/file_saver.py` and find the `_flatten_scoring` function. Add this new helper after it:

```python
def _flatten_route_prediction(route_pred: dict) -> dict:
    """Flattens a RouteResult dict into 5 CSV columns covering the top route.
    Full route data (all routes, all step breakdowns) preserved in JSON sidecar.
    """
    routes = route_pred.get("routes") or []
    if not routes:
        return {
            "route_top1_step_count": None,
            "route_top1_terminal_precursor": None,
            "route_top1_total_cost": None,
            "route_top1_composed_ee": None,
            "route_top1_warnings": "; ".join(route_pred.get("warnings") or []),
        }
    top = routes[0]
    return {
        "route_top1_step_count": len(top.get("steps") or []),
        "route_top1_terminal_precursor": top.get("terminal_precursor_name"),
        "route_top1_total_cost": top.get("total_cost"),
        "route_top1_composed_ee": top.get("composed_ee"),  # populated by scorer if available
        "route_top1_warnings": "; ".join(top.get("warnings") or []),
    }
```

- [ ] **Step 2: Update the main flattening loop in `save_suggestions_to_csv`**

Find the inner `for key, value in suggestion.items():` loop. Replace its body with:

```python
        for key, value in suggestion.items():
            if key == "scoring" and isinstance(value, dict):
                flat_dict.update(_flatten_scoring(value))
            elif key == "route_prediction" and isinstance(value, dict):
                flat_dict.update(_flatten_route_prediction(value))
            elif isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    flat_dict[f"{key}_{sub_key}"] = str(sub_value)
            else:
                flat_dict[key] = str(value)
```

- [ ] **Step 3: Verify the file is syntactically valid**

```bash
.venv/bin/python -c "from utils.file_saver import save_suggestions_to_csv, _flatten_route_prediction"
```

Expected: silent success.

- [ ] **Step 4: Smoke-check with a synthetic suggestion**

```bash
.venv/bin/python -c "
from utils.file_saver import save_suggestions_to_csv
import tempfile, os
sample = [{
    'name': 'test',
    'SMILES': 'CCO',
    'route_prediction': {
        'routes': [{'steps': [{}], 'terminal_precursor_name': 'pyruvate', 'total_cost': 1.5, 'composed_ee': 95.0, 'warnings': []}],
        'warnings': [],
    },
}]
with tempfile.TemporaryDirectory() as tmp:
    csv, jpath = save_suggestions_to_csv(sample, out_dir=tmp)
    print('Wrote:', csv)
    print(open(csv).read())
"
```

Expected: CSV output includes `route_top1_step_count,route_top1_terminal_precursor,route_top1_total_cost,...` columns.

- [ ] **Step 5: Commit**

```bash
git add utils/file_saver.py
git commit -m "extend file_saver to flatten route_prediction into 5 CSV columns

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 15: Add `_compose_route_ee` to scorer (PARALLELIZABLE with 13 and 14)

**Files:**
- Modify: `ChiraLLM/enantioselectivity_scorer.py`
- Modify: `tests/test_route_predictor_unit.py` (we add scorer-helper tests here for cohesion)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_route_predictor_unit.py`:

```python
class TestComposeRouteEe:
    def test_multiplicative_composition(self):
        from ChiraLLM.enantioselectivity_scorer import _compose_route_ee

        route = {
            "steps": [
                {"ec_numbers": ["1.1.1.184"]},
                {"ec_numbers": ["2.6.1.5"]},
            ]
        }
        brenda_data = {
            "1.1.1.184": {"status": "success", "entries": [{"enantioselectivity": 99.0}]},
            "2.6.1.5":   {"status": "success", "entries": [{"enantioselectivity": 95.0}]},
        }

        composed = _compose_route_ee(route, brenda_data)

        # 0.99 * 0.95 * 100 = 94.05
        assert composed is not None
        assert abs(composed - 94.05) < 0.01

    def test_returns_none_when_any_step_missing_ee(self):
        from ChiraLLM.enantioselectivity_scorer import _compose_route_ee

        route = {"steps": [{"ec_numbers": ["1.1.1.184"]}, {"ec_numbers": ["9.9.9.9"]}]}
        brenda_data = {"1.1.1.184": {"status": "success", "entries": [{"enantioselectivity": 99.0}]}}

        assert _compose_route_ee(route, brenda_data) is None

    def test_empty_route_returns_none(self):
        from ChiraLLM.enantioselectivity_scorer import _compose_route_ee
        assert _compose_route_ee({"steps": []}, {}) is None

    def test_picks_best_ee_across_multiple_ecs_per_step(self):
        from ChiraLLM.enantioselectivity_scorer import _compose_route_ee
        route = {"steps": [{"ec_numbers": ["1.1.1.184", "1.1.1.999"]}]}
        brenda_data = {
            "1.1.1.184": {"status": "success", "entries": [{"enantioselectivity": 50.0}]},
            "1.1.1.999": {"status": "success", "entries": [{"enantioselectivity": 95.0}]},
        }
        # Should pick 95% (the best)
        composed = _compose_route_ee(route, brenda_data)
        assert abs(composed - 95.0) < 0.01
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestComposeRouteEe -v
```

Expected: 4 ImportError or AttributeError failures.

- [ ] **Step 3: Implement `_compose_route_ee`**

Append to `ChiraLLM/enantioselectivity_scorer.py`:

```python
import math


def _compose_route_ee(route: dict, brenda_data: dict) -> Optional[float]:
    """Computes whole-route ee composition for a multi-step biosynthetic route.

    Formula: ee_overall = ∏(ee_step_i / 100) * 100, taking the best per-step ee
    across all ECs assigned to that step.

    Returns None if any step has no BRENDA-verified ee, or if the route has no steps.

    Known limitation: this assumes step independence. Routes with dynamic kinetic
    resolution (DKR) — where an upstream racemization step combined with a downstream
    enantioselective step produces an artificially-high terminal ee — are underestimated.
    See the route_predictor design spec §7.1.
    """
    steps = route.get("steps") or []
    if not steps:
        return None

    step_ees: list[float] = []
    for step in steps:
        ec_numbers = step.get("ec_numbers", [])
        best_ee = None
        for ec in ec_numbers:
            ec_data = brenda_data.get(ec)
            if not isinstance(ec_data, dict) or ec_data.get("status") != "success":
                continue
            for entry in ec_data.get("entries", []):
                ee_val = entry.get("enantioselectivity")
                if isinstance(ee_val, (int, float)):
                    if best_ee is None or ee_val > best_ee:
                        best_ee = float(ee_val)
        if best_ee is None:
            return None
        step_ees.append(best_ee)

    return math.prod(e / 100.0 for e in step_ees) * 100.0
```

- [ ] **Step 4: Wire `_compose_route_ee` into `score_suggestion`**

In `ChiraLLM/enantioselectivity_scorer.py`, find the end of `score_suggestion` (where the return dict is built). Just before the `return {` statement, insert:

```python
    # Whole-route ee composition (per-route, attached back into the route_prediction structure)
    route_pred = suggestion.get("route_prediction") or {}
    routes = route_pred.get("routes") or []
    for route in routes:
        composed = _compose_route_ee(route, brenda_data)
        route["composed_ee"] = composed
```

- [ ] **Step 5: Run scorer tests to verify integration**

```bash
.venv/bin/pytest tests/test_route_predictor_unit.py::TestComposeRouteEe -v
```

Expected: 4 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add ChiraLLM/enantioselectivity_scorer.py tests/test_route_predictor_unit.py
git commit -m "add _compose_route_ee multiplicative composition + wire into score_suggestion

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 16: Extend `smoke_test.py` with `route_prediction` fixture

**Files:**
- Modify: `smoke_test.py`

**Depends on:** Task 13 (main.py wiring must exist for the smoke test to exercise the new pipeline).

- [ ] **Step 1: Read the existing smoke test to identify the fixture pattern**

```bash
grep -n "MOCK_SUGGESTIONS\|SEGMENTS\|def evaluate" smoke_test.py | head -20
```

- [ ] **Step 2: Add a route_prediction-specific evaluate helper**

Append to `smoke_test.py` (before the `if __name__ == "__main__":` block):

```python
def evaluate_route_prediction(suggestions: list[dict]) -> dict:
    """Acceptance check for route_prediction sprint output.
    Per spec §10 #7: 3 KEGG-covered segments must succeed; 3 others must return
    appropriate error status without raising.
    """
    KEGG_COVERED_SEGMENTS = {"codexis", "pharma", "academic"}
    NON_COVERED_SEGMENTS = {"arnold", "ginkgo", "adversarial"}
    EXPECTED_NON_COVERED_STATUSES = {
        "no_kegg_id", "invalid_kegg_id", "target_not_in_kegg",
        "target_has_no_reactions", "no_route_found",
    }

    results = {"passes": [], "fails": []}
    for s in suggestions:
        seg = s.get("fixture_key", "unknown")
        rp = s.get("route_prediction", {})
        status = rp.get("status")
        n_routes = len(rp.get("routes", []))

        if seg in KEGG_COVERED_SEGMENTS:
            if status == "success" and n_routes >= 1:
                results["passes"].append(f"{seg}/{s.get('name')}: success, {n_routes} routes")
            else:
                results["fails"].append(f"{seg}/{s.get('name')}: expected success, got status={status} routes={n_routes}")
        elif seg in NON_COVERED_SEGMENTS:
            if status in EXPECTED_NON_COVERED_STATUSES:
                results["passes"].append(f"{seg}/{s.get('name')}: clean error status={status}")
            else:
                results["fails"].append(f"{seg}/{s.get('name')}: unexpected status={status}")
        else:
            results["fails"].append(f"{seg}/{s.get('name')}: unknown segment")
    return results
```

- [ ] **Step 3: Add `route_prediction` to the SEGMENTS list and CLI dispatch**

Find the line in `smoke_test.py` that reads:
```python
SEGMENTS = [
    ("Codexis — directed evolution / biocatalysis CRO",                "codexis"),
    ...
]
```

Add a new entry at the end:
```python
    ("Route prediction — Tier 1 acceptance check (all 6 fixture segments)", "route_prediction"),
```

In the CLI dispatch (around `argparse` handling), add:

```python
    if args.fixture == "route_prediction":
        # Run all 6 segments, then evaluate against route-prediction acceptance criteria
        all_results = []
        for label, key in SEGMENTS:
            if key == "route_prediction":
                continue  # don't recurse
            print(f"\n=== Running {label} ===")
            r = run_query(label, args.out_dir, mock=args.mock, fixture_key=key)
            all_results.extend(r.get("suggestions", []))
        verdict = evaluate_route_prediction(all_results)
        print("\n=== Route Prediction Acceptance ===")
        for p in verdict["passes"]:
            print(f"  PASS  {p}")
        for f in verdict["fails"]:
            print(f"  FAIL  {f}")
        if verdict["fails"]:
            sys.exit(1)
        print("Overall: ROUTE PREDICTION ACCEPTANCE PASSED")
        sys.exit(0)
```

- [ ] **Step 4: Run the route_prediction fixture (live KEGG; mocked GPT)**

```bash
.venv/bin/python smoke_test.py --fixture route_prediction --mock
```

Expected: ROUTE PREDICTION ACCEPTANCE PASSED.

- [ ] **Step 5: Commit**

```bash
git add smoke_test.py
git commit -m "add route_prediction fixture to smoke test per Sprint 1 acceptance #7

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 17: Integration tests against live KEGG

**Files:**
- Create: `tests/test_route_predictor_integration.py`

**Depends on:** Task 12 (`predict_route` public API must exist).

- [ ] **Step 1: Create the integration test file**

```python
"""Layer 2 integration tests for route_predictor.
Hits live KEGG. Skipped by default — run with: pytest -m integration

Tests assert structural properties only, not exact route content (KEGG drift)."""

import pytest
from ChiraLLM import route_predictor


pytestmark = pytest.mark.integration


def test_pantolactone_returns_real_routes():
    """(R)-pantolactone (C00599) — Codexis-segment classic; 3-step KEGG route to KIV."""
    result = route_predictor.predict_route("C00599", mode="top_n", n=3, budget=200)

    assert result.status == "success", f"Expected success, got {result.status}: {result.warnings}"
    assert len(result.routes) >= 1
    assert result.nodes_explored > 0
    assert result.nodes_explored < 200, "should not exhaust budget on this well-connected target"

    # Each route should have a terminal in CENTRAL_METABOLITES
    for route in result.routes:
        assert route.terminal_precursor_id in route_predictor.CENTRAL_METABOLITES
        assert route.terminal_precursor_name  # non-empty string
        assert route.total_cost > 0
        assert all(isinstance(s.reaction_id, str) and s.reaction_id.startswith("R") for s in route.steps)


def test_lactic_acid_finds_short_route():
    """(S)-lactic acid (C00186) — should find 1-step route to pyruvate via L-LDH (1.1.1.27)."""
    result = route_predictor.predict_route("C00186", mode="top_n", n=3, budget=100)

    assert result.status == "success"
    assert len(result.routes) >= 1
    # The shortest route should be 1 step (lactate → pyruvate)
    shortest = min(result.routes, key=lambda r: len(r.steps))
    assert len(shortest.steps) == 1
    assert shortest.terminal_precursor_id == "C00022"  # pyruvate


def test_norcoclaurine_handles_poor_kegg_coverage():
    """(S)-norcoclaurine (C09136) has poor KEGG enzyme annotation.
    The predictor should return a clean error or partial-success status without raising."""
    result = route_predictor.predict_route("C09136", mode="top_n", n=3, budget=100)

    # Acceptable outcomes: partial success, no_route_found, target_has_no_reactions
    assert result.status in {"success", "no_route_found", "target_has_no_reactions"}
    # Whatever happens, no exceptions should propagate
    assert isinstance(result.warnings, list)
```

- [ ] **Step 2: Run the integration tests against live KEGG**

```bash
.venv/bin/pytest tests/test_route_predictor_integration.py -v -m integration
```

Expected: 3 tests PASS. (May take 30-90 seconds — first run populates the disk cache.)

- [ ] **Step 3: Commit**

```bash
git add tests/test_route_predictor_integration.py
git commit -m "add integration tests for route predictor against live KEGG

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 18: Pre-warm script for eQuilibrator cache

**Files:**
- Create: `scripts/warm_equilibrator_cache.py`

**Depends on:** Task 7 (`_fetch_delta_g_kj_per_mol` must exist).

- [ ] **Step 1: Create the script directory if needed**

```bash
mkdir -p /Users/alexei/Documents/GitHub/ChiralAI/scripts
```

- [ ] **Step 2: Write the pre-warm script**

Create `scripts/warm_equilibrator_cache.py`:

```python
"""Pre-fetch ΔrG' for KEGG reactions in iJO1366 + a curated short list,
populating the disk cache so cold-start route prediction runs are fast.

Run once after install: python scripts/warm_equilibrator_cache.py
"""

import sys
import time
from ChiraLLM import route_predictor

# Curated reactions worth pre-warming even outside iJO1366 — common biocatalysis benchmarks.
ADDITIONAL_REACTIONS = [
    "R02472",  # 2-dehydropantoate reductase (pantolactone)
    "R00754",  # alcohol dehydrogenase
    "R00277",  # L-lactate dehydrogenase
    "R00351",  # citrate synthase
]


def collect_ijo1366_reactions() -> list[str]:
    """Returns the list of KEGG reaction IDs annotated in iJO1366 metabolites' reactions."""
    try:
        from cobra.io import load_model
    except ImportError:
        print("cobra not installed; skipping iJO1366 reactions")
        return []
    print("Loading iJO1366...", flush=True)
    model = load_model("iJO1366")
    kegg_rxns = set()
    for rxn in model.reactions:
        kegg_ids = rxn.annotation.get("kegg.reaction", [])
        if isinstance(kegg_ids, str):
            kegg_ids = [kegg_ids]
        kegg_rxns.update(kid for kid in kegg_ids if kid.startswith("R"))
    return sorted(kegg_rxns)


def main():
    rxns = collect_ijo1366_reactions() + ADDITIONAL_REACTIONS
    rxns = sorted(set(rxns))
    print(f"Pre-warming eQuilibrator cache for {len(rxns)} reactions...", flush=True)

    n_success = 0
    n_failed = 0
    for i, rxn_id in enumerate(rxns):
        if i % 50 == 0:
            print(f"  [{i}/{len(rxns)}] success={n_success} failed={n_failed}", flush=True)
        dg = route_predictor._fetch_delta_g_kj_per_mol(rxn_id)
        if dg is not None:
            n_success += 1
        else:
            n_failed += 1
        # Be polite to the eQuilibrator REST endpoint (no documented rate limit, but it's free)
        time.sleep(0.1)

    print(f"\nDone: {n_success} cached, {n_failed} failed (likely no eQuilibrator data for those)")
    print(f"Cache location: {route_predictor._cache_root() / 'equilibrator'}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Verify the script imports cleanly (do not run it — it would take minutes)**

```bash
.venv/bin/python -c "import scripts.warm_equilibrator_cache; print('OK')"
```

If `scripts` is not a package: add an empty `scripts/__init__.py` then re-test.

- [ ] **Step 4: Commit**

```bash
git add scripts/
git commit -m "add eQuilibrator cache pre-warm script for iJO1366 reactions

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 19: Documentation updates

**Files:**
- Modify: `CLAUDE.md`
- Modify: `README.md`

- [ ] **Step 1: Update CLAUDE.md module status table**

In `CLAUDE.md`, find the Architecture section's module list. Add the new module after `enantioselectivity_scorer.py`:

```
ChiraLLM/route_predictor.py           Tier 1 biosynthesis route predictor. predict_route(compound_id,
                                      mode='top_n'|'full_tree'|'shortest_plus_diverse', n, budget)
                                      runs weighted A* backward through KEGG reaction graph from
                                      target to curated central metabolites. Edge cost is a
                                      transparent breakdown: thermodynamic (eQuilibrator ΔG),
                                      directionality (KEGG <=> vs =>), industrial-reversibility
                                      override (KREDs/TAs/IREDs/lipases/BVMOs). Heuristic is
                                      Morgan-fingerprint Tanimoto distance to nearest central.
                                      Disk-cached at ~/.cache/chiralai/. Returns RouteResult
                                      with status field; never raises.
```

In the module status table, add:

```
| `route_predictor.py` | Working | Tier 1 (KEGG traversal); Tier 2 RetroRules SMARTS deferred |
```

In Known Gaps, remove "Retrosynthetic route prediction: Not implemented" and replace with:

```
- **Tier 2 retrobiosynthesis**: Tier 1 covers compounds KEGG already knows. Novel-target retrobiosynthesis (target SMILES → enzymatic disconnection via RetroRules SMARTS + RDKit `RunReactants`) is planned for Sprint 2.
```

- [ ] **Step 2: Update README.md**

In `README.md`, find the pipeline diagram. Add a new step after KEGG and before BRENDA:

```
    ↓
Route predictor — backward search through KEGG reactions from target to central metabolites
    ↓
BRENDA — retrieves known ee values for ECs across all route steps (deduplicated)
```

In the Architecture section, add `ChiraLLM/route_predictor.py` to the file list with one-line description.

In the Known Limitations section, change the bullet on retrosynthetic route planning to:

```
- **Tier 1 only.** Route prediction works for compounds KEGG already covers (~12k reactions). Novel targets require Tier 2 (RetroRules SMARTS retrobiosynthesis), planned for the next sprint.
```

In the Roadmap, replace the "Retrosynthetic route prediction" bullet with:

```
- [ ] Tier 2 — novel-target retrobiosynthesis via RetroRules + RDKit RunReactants
```

- [ ] **Step 3: Verify markdown rendering (manual check)**

```bash
head -30 README.md && echo "---" && head -30 CLAUDE.md
```

Expected: Files render as readable markdown without syntax errors.

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md README.md
git commit -m "document Tier 1 route predictor; move retrobiosynthesis from roadmap to working

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Self-Review Pass

After completing all 20 tasks (0–19), verify:

**Spec coverage:**

| Spec section | Implementing task |
|--------------|-------------------|
| §1 Context | All tasks |
| §2 Architecture (insertion point) | Task 13 |
| §2 Single-file module layout | Tasks 1–12 |
| §3.1 Constants | Task 1 |
| §3.2 Function inventory (all 12) | Tasks 2, 3, 5, 6, 7, 8, 9, 10, 11, 12 |
| §3.3 Dataclasses | Task 10 |
| §4 Data flow worked example | Tasks 12, 13 |
| §5.1 Failure mode matrix (14 rows) | Task 12 (most), Tasks 5, 6, 7 (network failures), Task 8 (MOL fallback) |
| §6.1 Layer 1 unit tests | Tasks 0–12 (test-first throughout) |
| §6.1 Layer 2 integration tests | Task 17 |
| §6.1 Layer 3 smoke | Task 16 |
| §10 AC #1 module exists | Task 12 |
| §10 AC #2 unit coverage ≥85% | Task 12 verification step |
| §10 AC #3 integration tests | Task 17 |
| §10 AC #4 main.py wiring | Task 13 |
| §10 AC #5 _compose_route_ee in scorer | Task 15 |
| §10 AC #6 file_saver flat columns | Task 14 |
| §10 AC #7 smoke test fixture | Task 16 |
| §10 AC #8 no new pip deps | Tasks 1–12 use only existing deps |
| §10 AC #9 cache + pre-warm script | Task 4 (cache), Task 18 (script) |
| §10 AC #10 docs updated | Task 19 |

**Type consistency check:**
- `RouteStep.precursor_id` and `RouteStep.intermediate_id` used consistently across Tasks 10, 11, 13, 14, 15.
- `Route.terminal_precursor_id` used consistently across Tasks 10, 13, 14.
- `RouteResult.status` values match the failure-mode matrix and the AC #7 expected statuses (`success`, `no_kegg_id`, `invalid_kegg_id`, `target_not_in_kegg`, `target_has_no_reactions`, `no_route_found`, `invalid_mode`).

**Placeholder scan:**
- No "TBD", "TODO", or "implement later" anywhere in plan steps.
- Every code block is complete and runnable as written.
- Every test has full assertion code, not "// add assertions".
- No "similar to Task N" — each task's code is self-contained.

**Acceptance criteria satisfied per task:**
- AC #8 (no new pip deps) — Task 0 only adds dev deps (`pytest`, `pytest-mock`, `pytest-cov`); runtime `requirements.txt` unchanged.
- AC #2 (≥85% coverage) — verified explicitly in Task 12 Step 4.
- AC #7 (smoke acceptance) — verified explicitly in Task 16 Step 4.
