"""Thermodynamics tests, written against the real contract.

The test this replaces defined a ``DummyCC`` that accepted a bare float for
``p_h`` and returned ``.to(unit).m``. Both are wrong: eQuilibrator needs
unit-bearing quantities and returns a pint ``Measurement``. The mock encoded
the bug, so the suite stayed green through three live defects and 299 calls
that produced no number. These tests assert on the *inputs we send* and the
*shape we emit*, which is what actually broke.
"""
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utils import thermo


def test_pathway_map_id_is_rejected_as_input():
    """The original defect: a pathway map arriving where a reaction belongs."""
    r = thermo.delta_g_kegg("path:map00260")
    assert r["status"] == "error"
    assert r["delta_g"] is None
    assert "not a reaction equation" in r["reason"]


def test_input_is_validated_before_engine_availability():
    """Otherwise the bug is invisible whenever the engine happens to be absent."""
    assert thermo.delta_g_kegg("")["status"] == "error"


def test_kegg_ids_are_namespaced(monkeypatch):
    """eQuilibrator resolves 'kegg:C00186', not a bare 'C00186'."""
    seen = {}

    class CC:
        def __setattr__(self, k, v):
            seen[k] = v

        def parse_reaction_formula(self, formula):
            seen["formula"] = formula
            raise RuntimeError("stop here; we only care about the formula")

    monkeypatch.setattr(thermo, "_engine", lambda: (CC(), None))
    thermo.delta_g_kegg("C00186 + C00003 <=> C00022 + C00004")
    assert seen["formula"] == "kegg:C00186 + kegg:C00003 <=> kegg:C00022 + kegg:C00004"


def test_conditions_are_sent_as_quantities_not_floats(monkeypatch):
    """A bare float for p_h raises deep inside pint. That was bug number two."""
    seen = {}

    class CC:
        def __setattr__(self, k, v):
            seen[k] = v

        def parse_reaction_formula(self, formula):
            raise RuntimeError("stop")

    monkeypatch.setattr(thermo, "_engine", lambda: (CC(), None))
    thermo.delta_g_kegg("C00186 <=> C00022", ph=7.5, ionic_strength=0.25)
    assert not isinstance(seen["p_h"], float), "p_h must carry units"
    assert not isinstance(seen["ionic_strength"], float)


def test_huge_uncertainty_is_reported_as_not_estimable(monkeypatch):
    """A central value with +/-100000 kJ/mol is eQuilibrator saying 'I cannot'."""
    class Q:
        def __init__(self, v): self.v = v
        def m_as(self, _u): return self.v

    class M:
        value, error = Q(-87.4), Q(100000.0)

    class CC:
        def __setattr__(self, k, v): pass
        def parse_reaction_formula(self, f): return object()
        def standard_dg_prime(self, r): return M()

    monkeypatch.setattr(thermo, "_engine", lambda: (CC(), None))
    r = thermo.delta_g_kegg("C00186 <=> C00022")
    assert r["status"] == "not-estimable"
    assert r["delta_g"] is None, "a meaningless central value must not be reported"


def test_failure_carries_a_reason():
    """The blanket 'except Exception: return None' is why this was unfalsifiable."""
    assert thermo.delta_g_kegg("nonsense").get("reason")
