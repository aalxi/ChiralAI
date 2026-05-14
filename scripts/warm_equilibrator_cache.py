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
        # Be polite to the eQuilibrator REST endpoint
        time.sleep(0.1)

    print(f"\nDone: {n_success} cached, {n_failed} failed (likely no eQuilibrator data for those)")
    print(f"Cache location: {route_predictor._cache_root() / 'equilibrator'}")


if __name__ == "__main__":
    main()
