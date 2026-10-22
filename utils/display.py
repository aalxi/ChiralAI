"""Terminal report.

The old report printed ``Confidence: 88%``. That number took two distinct
values across a 36-compound benchmark, and 0.88 was emitted both for "too large
to make" and for "no route found" -- two opposite diagnoses sharing one figure.
It is gone. What replaces it is the evidence itself: which checks passed, which
layers were never consulted, and what the stereochemical answer was.

A reader should be able to see at a glance which claims are grounded in a
curated database and which are absent.
"""
from __future__ import annotations

VERDICT_LINE = {
    "route-documented-in-host": "ROUTE DOCUMENTED IN HOST",
    "route-documented-but-not-in-host": "ROUTE DOCUMENTED, NOT IN THIS HOST",
    "route-documented-host-unknown": "ROUTE DOCUMENTED, HOST NOT CHECKED",
    "route-documented-for-the-wrong-stereoisomer": "ROUTE DOCUMENTED, BUT FOR THE WRONG ENANTIOMER",
    "no-documented-biosynthetic-route": "NO DOCUMENTED BIOSYNTHETIC ROUTE",
    "no-kegg-record": "NO KEGG RECORD FOR THIS STRUCTURE",
    "rejected-invalid-structure": "REJECTED - STRUCTURE DOES NOT MATCH THE NAME",
}

STEREO_LINE = {
    "EXACT": "the requested stereoisomer is enzymatically accessible",
    "WRONG-FORM": "only the OTHER stereoisomer is curated - wrong enantiomer",
    "SKELETON": "stereochemistry unspecified in the query",
    "NONE": "no curated reaction for this connectivity",
    "NOT-ASSESSED": "not assessed",
}


def display_results_in_terminal(results: list[dict]) -> None:
    if not results:
        print("\nNo results to display.")
        return

    print("\n" + "=" * 72)
    print("   CHIRALAI ASSESSMENT")
    print("=" * 72)

    for i, r in enumerate(results, 1):
        if r.get("error"):
            print(f"\n--- {i}. error: {r['error']} ---")
            continue

        print(f"\n--- {i}. {r.get('name', '?')} ---")
        print(f"  verdict      : {VERDICT_LINE.get(r.get('verdict'), r.get('verdict'))}")

        for reason in r.get("reasons", []):
            for j, chunk in enumerate(_wrap(reason, 64)):
                print(f"      {'* ' if j == 0 else '  '}{chunk}")

        st = r.get("stereochemistry", {})
        out = st.get("outcome", "NOT-ASSESSED")
        print(f"\n  stereochemistry: {out} - {STEREO_LINE.get(out, '')}")
        if st.get("enzymes_for_requested_form"):
            print(f"      makes your form : {', '.join(st['enzymes_for_requested_form'][:8])}")
        risk = st.get("selectivity_risk")
        if risk:
            print(f"      SELECTIVITY RISK: {risk['enantiomer']}")
            print(f"        made by       : {', '.join(risk['made_by_ec'])}")
        for e in st.get("enzymes_for_other_forms", [])[:4]:
            print(f"      other isomer    : {e['compound'][:42]:42s} EC {', '.join(e['ec'])[:24]}")

        ret = r.get("retrieval", {})
        if ret:
            print(f"\n  retrieval    : KEGG {ret.get('kegg_id') or '--'} "
                  f"via {ret.get('method')} ({ret.get('how', '')[:46]})")

        thermo = [t for t in r.get("thermodynamics", []) if t.get("status") == "ok"]
        if thermo:
            print("  thermodynamics (ΔG'° kJ/mol, pH 7):")
            for t in thermo[:5]:
                print(f"      {t['reaction']:10s} {t['delta_g']:>8.1f} +/- {t['uncertainty']:.1f}")

        ev = r.get("evidence", {})
        if ev:
            print("\n  evidence     : " + ", ".join(
                f"{k}={v}" for k, v in ev.items() if k != "stereochemistry_outcome"))

        for layer in r.get("layers_unavailable", []):
            print(f"  NOT CONSULTED: {layer['layer']} - {layer['reason'][:58]}")

        if r.get("next_steps"):
            print("  next steps   : " + "; ".join(r["next_steps"])[:150])

    print("\n" + "=" * 72)
    print("No field above is a prediction. Each is either a curated database")
    print("record or an explicit statement that a layer was not consulted.")
    print("=" * 72)


def _wrap(text: str, width: int) -> list[str]:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return lines
