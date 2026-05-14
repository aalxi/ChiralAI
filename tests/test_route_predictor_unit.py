"""Layer 1 unit tests for ChiraLLM/route_predictor.py.
All external calls (KEGG, eQuilibrator, RDKit MOL parsing) are mocked.
Coverage target: ≥85%."""

import os
import time

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

    def test_oye_full_ec_does_not_overmatch(self):
        """1.6.99.1 (Old Yellow Enzyme) should match ONLY EC 1.6.99.1, not 1.6.99.10+."""
        assert route_predictor._is_industrially_reversible(["1.6.99.1"]) is True
        assert route_predictor._is_industrially_reversible(["1.6.99.10"]) is False
        assert route_predictor._is_industrially_reversible(["1.6.99.12"]) is False

    def test_bvm_full_ec_does_not_overmatch(self):
        """1.14.13.22 (BVMO) should match ONLY EC 1.14.13.22, not hypothetical 1.14.13.220."""
        assert route_predictor._is_industrially_reversible(["1.14.13.22"]) is True
        assert route_predictor._is_industrially_reversible(["1.14.13.220"]) is False


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
        # Force TTL to 0 days = always expired
        monkeypatch.setenv("CHIRALAI_CACHE_TTL_DAYS", "0")
        route_predictor._disk_cache_set("kegg", "R12345", "stale")
        assert route_predictor._disk_cache_get("kegg", "R12345") is None

    def test_unicode_content_preserved(self, tmp_cache_dir):
        route_predictor._disk_cache_set("kegg", "R12345", "alpha-α-ketoglutarate")
        assert route_predictor._disk_cache_get("kegg", "R12345") == "alpha-α-ketoglutarate"

    def test_invalid_ttl_env_raises_clear_error(self, tmp_cache_dir, monkeypatch):
        monkeypatch.setenv("CHIRALAI_CACHE_TTL_DAYS", "abc")
        route_predictor._disk_cache_set("kegg", "R12345", "data")
        with pytest.raises(ValueError, match="CHIRALAI_CACHE_TTL_DAYS must be an integer"):
            route_predictor._disk_cache_get("kegg", "R12345")

    def test_get_tolerates_concurrent_deletion(self, tmp_cache_dir, monkeypatch):
        # Simulate the file being deleted between the existence check and stat
        path = tmp_cache_dir / "kegg" / "R12345.cache"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("data")
        original_stat = type(path).stat

        def fake_stat(self, *a, **kw):
            if self.name == "R12345.cache":
                raise FileNotFoundError(self)
            return original_stat(self, *a, **kw)

        monkeypatch.setattr(type(path), "stat", fake_stat)
        assert route_predictor._disk_cache_get("kegg", "R12345") is None

    def test_clear_disk_cache_only_removes_cache_files(self, tmp_cache_dir):
        # Set up some cache files and a non-cache file
        route_predictor._disk_cache_set("kegg", "R1", "a")
        route_predictor._disk_cache_set("kegg", "R2", "b")
        route_predictor._disk_cache_set("equilibrator", "R3", "c")
        unrelated = tmp_cache_dir / "unrelated.txt"
        unrelated.write_text("don't delete me")

        n = route_predictor._clear_disk_cache()

        assert n == 3
        assert unrelated.exists()  # non-cache file preserved
        assert route_predictor._disk_cache_get("kegg", "R1") is None  # cache cleared
