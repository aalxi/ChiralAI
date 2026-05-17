"""Run the chiral-route recall benchmark.

Scope: validates ChiraLLM.route_predictor.predict_route in isolation against
literature-documented biosynthesis routes. Does NOT exercise the full main.py
pipeline (LLM is non-deterministic; BRENDA returns no_credentials in this
environment, which would give degenerate ee scoring).

Per-compound metrics:
  route_found              — at least one route returned
  terminal_match           — terminal_precursor_id in expected_terminals
  ec_overlap_count         — number of expected ECs found anywhere in the top route
  step_count_predicted     — steps in top-1 route (or None)
  step_within_2            — abs(predicted - expected) ≤ 2

Aggregate metrics:
  routes_found_pct
  terminal_match_pct
  ec_overlap_mean / median
  step_within_2_pct

Usage:
  python3 -m benchmarks.run_benchmark
"""

import json
import logging
import statistics
import sys
import time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

# Allow running as `python3 benchmarks/run_benchmark.py` from repo root
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from ChiraLLM.route_predictor import predict_route  # noqa: E402
from benchmarks.targets import BENCHMARK_TARGETS  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("benchmark")

# Quiet route_predictor's per-step INFO chatter; keep our own progress visible
logging.getLogger("ChiraLLM.route_predictor").setLevel(logging.WARNING)


def score_target(target: dict, budget: int = 200) -> dict:
    """Run predict_route and score against the target's expected ground truth."""
    t0 = time.perf_counter()
    result = predict_route(target["kegg_id"], mode="top_n", n=3, budget=budget)
    elapsed = time.perf_counter() - t0

    base = {
        "kegg_id": target["kegg_id"],
        "name": target["name"],
        "category": target["category"],
        "status": result.status,
        "elapsed_sec": round(elapsed, 2),
        "nodes_explored": result.nodes_explored,
        "budget_exhausted": result.budget_exhausted,
        "expected_terminals": target["expected_terminals"],
        "expected_steps": target["expected_steps"],
        "expected_ec_numbers": target["expected_ec_numbers"],
        "warnings": list(result.warnings),
    }

    if result.status != "success" or not result.routes:
        return {
            **base,
            "route_found": False,
            "terminal_match": False,
            "ec_overlap_count": 0,
            "step_count_predicted": None,
            "step_within_2": False,
            "predicted_terminal": None,
            "predicted_ecs": [],
        }

    top_route = result.routes[0]
    predicted_ecs = sorted({ec for step in top_route.steps for ec in (step.ec_numbers or [])})
    expected_ecs = set(target["expected_ec_numbers"])
    ec_overlap = sorted(expected_ecs & set(predicted_ecs))

    return {
        **base,
        "route_found": True,
        "terminal_match": top_route.terminal_precursor_id in target["expected_terminals"],
        "ec_overlap_count": len(ec_overlap),
        "ec_overlap": ec_overlap,
        "step_count_predicted": len(top_route.steps),
        "step_within_2": abs(len(top_route.steps) - target["expected_steps"]) <= 2,
        "predicted_terminal": top_route.terminal_precursor_id,
        "predicted_terminal_name": top_route.terminal_precursor_name,
        "predicted_ecs": predicted_ecs,
        "top_route_total_cost": round(top_route.total_cost, 3),
        "top_route_step_summary": [
            {
                "rxn": s.reaction_id,
                "ec": s.ec_numbers,
                "from": s.precursor_id,
                "to": s.intermediate_id,
            }
            for s in top_route.steps
        ],
    }


def aggregate(per_target: list[dict]) -> dict:
    n = len(per_target)
    found = [r for r in per_target if r["route_found"]]
    matched = [r for r in found if r["terminal_match"]]
    within_2 = [r for r in found if r["step_within_2"]]
    ec_overlaps = [r["ec_overlap_count"] for r in found]
    ec_overlaps_norm = [
        r["ec_overlap_count"] / max(1, len(r["expected_ec_numbers"])) for r in found
    ]

    return {
        "n_total": n,
        "n_route_found": len(found),
        "n_terminal_match": len(matched),
        "n_step_within_2": len(within_2),
        "routes_found_pct": round(100 * len(found) / n, 1),
        "terminal_match_pct": round(100 * len(matched) / n, 1),
        "step_within_2_pct": round(100 * len(within_2) / n, 1) if n else 0,
        "ec_overlap_mean": round(statistics.mean(ec_overlaps), 2) if ec_overlaps else 0,
        "ec_overlap_median": statistics.median(ec_overlaps) if ec_overlaps else 0,
        "ec_overlap_normalized_mean": (
            round(statistics.mean(ec_overlaps_norm), 3) if ec_overlaps_norm else 0
        ),
    }


def main() -> int:
    logger.info("Running benchmark on %d chiral targets", len(BENCHMARK_TARGETS))

    results = []
    for i, target in enumerate(BENCHMARK_TARGETS, 1):
        logger.info("[%d/%d] %s (%s)", i, len(BENCHMARK_TARGETS), target["name"], target["kegg_id"])
        try:
            scored = score_target(target)
        except Exception as e:
            logger.exception("scoring failed for %s", target["kegg_id"])
            scored = {
                "kegg_id": target["kegg_id"],
                "name": target["name"],
                "category": target["category"],
                "status": "harness_error",
                "error": str(e),
                "route_found": False,
                "terminal_match": False,
                "ec_overlap_count": 0,
                "step_count_predicted": None,
                "step_within_2": False,
            }
        status = "OK" if scored.get("terminal_match") else (
            "PARTIAL" if scored.get("route_found") else "MISS"
        )
        logger.info("  %s — status=%s elapsed=%.1fs nodes=%d",
                    status, scored.get("status"), scored.get("elapsed_sec", 0),
                    scored.get("nodes_explored", 0))
        results.append(scored)

    summary = aggregate(results)
    logger.info("Aggregate: routes_found=%s%% terminal_match=%s%% step±2=%s%% ec_overlap_mean=%s",
                summary["routes_found_pct"], summary["terminal_match_pct"],
                summary["step_within_2_pct"], summary["ec_overlap_mean"])

    timestamp = datetime.now().strftime("%Y-%m-%d")
    out_path = REPO_ROOT / "benchmarks" / "results" / f"{timestamp}.json"
    out_path.write_text(json.dumps(
        {"summary": summary, "results": results, "run_at": datetime.now().isoformat()},
        indent=2, default=str,
    ))
    logger.info("Wrote results to %s", out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
