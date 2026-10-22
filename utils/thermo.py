"""Reaction thermodynamics via eQuilibrator, honestly reported.

Three defects are fixed here.

1. **It was being fed the wrong kind of thing.** The caller passed
   ``path:map00260`` -- a pathway *map* identifier -- where a stoichiometric
   reaction equation was expected. Benchmarked: 299 calls, 0 numbers returned.
   The entry point now takes a KEGG reaction equation and namespaces the
   compound ids the way eQuilibrator requires.

2. **A blanket ``except Exception: return None``** turned every such failure
   into a silent, confident "no value available", so the defect was invisible
   for the lifetime of the repo. Failures now come back with a reason attached.

3. **``ComponentContribution()`` ran at import time**, so merely importing this
   module tried to build a multi-hundred-megabyte local cache. It is now
   constructed on first use.
"""
from __future__ import annotations
import functools, re

_CID = re.compile(r"\bC\d{5}\b")

# Above this, eQuilibrator is signalling "no estimate", not a wide estimate.
NOT_ESTIMABLE_KJ = 200.0


@functools.lru_cache(maxsize=1)
def _engine():
    """(ComponentContribution | None, reason_unavailable)."""
    try:
        import warnings
        with warnings.catch_warnings():
            # equilibrator_cache's ORM models trip a SQLAlchemy naming warning
            # on import. It is theirs, not ours, and not actionable here.
            warnings.simplefilter("ignore")
            from equilibrator_api import ComponentContribution
    except Exception as exc:
        return None, (f"equilibrator-api is not installed ({type(exc).__name__}). "
                      f"pip install equilibrator-api")
    try:
        return ComponentContribution(), None
    except Exception as exc:
        return None, (f"eQuilibrator could not initialise its thermodynamic cache: "
                      f"{type(exc).__name__}: {str(exc)[:160]}")


def available() -> tuple[bool, str | None]:
    cc, reason = _engine()
    return cc is not None, reason


def delta_g_kegg(equation: str, ph: float = 7.0, ionic_strength: float = 0.1) -> dict:
    """ΔG'° in kJ/mol for a KEGG reaction equation.

    ``equation`` is the EQUATION line of a KEGG reaction record, for example
    ``C00186 + 2 C00125 <=> C00022 + 2 C00126 + 2 C00080``. Bare KEGG ids are
    namespaced to ``kegg:C00186`` because that is what eQuilibrator resolves.
    """
    # Validate the INPUT before checking the engine. The original defect was a
    # pathway map id arriving here; if availability were tested first, that bug
    # would still be invisible whenever the engine happens to be missing.
    if not equation or "=" not in equation:
        return {"status": "error", "delta_g": None,
                "reason": f"not a reaction equation: {equation!r}. Expected a KEGG "
                          f"EQUATION line such as 'C00186 + C00003 <=> C00022 + C00004'."}
    cc, reason = _engine()
    if cc is None:
        return {"status": "unavailable", "delta_g": None, "reason": reason}
    formula = _CID.sub(lambda m: f"kegg:{m.group(0)}", equation)
    try:
        # These are unit-bearing quantities, not floats. Assigning a bare float
        # raises AttributeError deep inside pint, far from the assignment -- and
        # the blanket except that used to live here swallowed it, which is how a
        # one-line bug stayed invisible for the lifetime of the repo.
        from equilibrator_api import Q_
        cc.p_h = Q_(ph)
        cc.ionic_strength = Q_(f"{ionic_strength}M")
        rxn = cc.parse_reaction_formula(formula)
        dg = cc.standard_dg_prime(rxn)
        value = float(dg.value.m_as("kJ/mol"))
        error = float(dg.error.m_as("kJ/mol"))
        # eQuilibrator returns a central value with a huge uncertainty when it
        # cannot actually estimate the reaction -- a participant outside its
        # training set, e.g. cytochrome c. The central value is then meaningless,
        # and reporting it alone manufactures precision that does not exist.
        if error > NOT_ESTIMABLE_KJ:
            return {"status": "not-estimable", "delta_g": None,
                    "uncertainty": round(error, 1), "equation": formula,
                    "reason": f"uncertainty +/-{error:.0f} kJ/mol: at least one "
                              f"participant is outside eQuilibrator's coverage, so "
                              f"the central value ({value:.1f}) carries no information"}
        return {"status": "ok", "delta_g": round(value, 2),
                "uncertainty": round(error, 2),
                "equation": formula, "ph": ph, "ionic_strength": ionic_strength}
    except Exception as exc:
        return {"status": "error", "delta_g": None, "equation": formula,
                "reason": f"{type(exc).__name__}: {str(exc)[:160]}"}


def delta_g(equation: str, ph: float = 7.0, ionic_strength: float = 0.1) -> float | None:
    """Backwards-compatible shim: the number, or None. Prefer ``delta_g_kegg``."""
    return delta_g_kegg(equation, ph, ionic_strength)["delta_g"]
