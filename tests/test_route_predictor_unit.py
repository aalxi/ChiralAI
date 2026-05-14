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
        assert abs(cost["total"] - forward_cost["total"]) <= 0.5

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
