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
