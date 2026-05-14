"""Layer 2 integration tests for route_predictor.
Hits live KEGG. Skipped by default — run with: pytest -m integration

Tests assert structural properties only, not exact route content (KEGG drift)."""

import pytest
from ChiraLLM import route_predictor


pytestmark = pytest.mark.integration


def test_pantolactone_returns_real_routes():
    """(R)-pantolactone (C00599) — Codexis-segment classic; multi-step KEGG route to a central metabolite."""
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
    """(S)-lactic acid (C00186) — should find a 1-step route to pyruvate via L-LDH (1.1.1.27)."""
    result = route_predictor.predict_route("C00186", mode="top_n", n=3, budget=100)

    assert result.status == "success"
    assert len(result.routes) >= 1
    # The shortest route should be 1 step (lactate → pyruvate)
    shortest = min(result.routes, key=lambda r: len(r.steps))
    assert len(shortest.steps) == 1
    assert shortest.terminal_precursor_id == "C00022"  # pyruvate


def test_norcoclaurine_handles_poor_kegg_coverage():
    """(S)-norcoclaurine (C09136) has poor KEGG enzyme annotation.
    The predictor should return a clean status without raising."""
    result = route_predictor.predict_route("C09136", mode="top_n", n=3, budget=100)

    # Acceptable outcomes: success, no_route_found, target_has_no_reactions
    assert result.status in {"success", "no_route_found", "target_has_no_reactions"}
    # Whatever happens, no exceptions should propagate
    assert isinstance(result.warnings, list)
