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
    """Patches _fetch_kegg_reaction and _fetch_compound_reactions to return synthetic data.
    Tests that need different fixtures can extend by mutating synthetic_kegg before use."""
    from ChiraLLM import route_predictor

    def fake_fetch_reaction(rxn_id):
        rxn = synthetic_kegg["reactions"].get(rxn_id)
        if rxn is None:
            return None
        eq = rxn["equation"]
        # Parse the synthetic equation without calling _parse_reaction_equation, because
        # synthetic compound IDs (e.g. "C_INTERMEDIATE_A") don't match KEGG's C\d{5} regex.
        if "<=>" in eq:
            direction, sep = "reversible", "<=>"
        else:
            direction, sep = "forward_only", "=>"
        left, right = eq.split(sep, 1)

        def parse_side(side):
            result = []
            for tok in side.split("+"):
                parts = tok.strip().split()
                if not parts:
                    continue
                if len(parts) == 2:
                    result.append((int(parts[0]), parts[1]))
                else:
                    result.append((1, parts[0]))
            return result

        return {
            "rxn_id": rxn_id,
            "equation": eq,
            "substrates": parse_side(left),
            "products": parse_side(right),
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
    return synthetic_kegg
