"""Layer 1 unit tests for ChiraLLM/route_predictor.py.
All external calls (KEGG, eQuilibrator, RDKit MOL parsing) are mocked.
Coverage target: ≥85%."""

import pytest

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
