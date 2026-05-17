# Tier 1 Route Predictor — Benchmark Report

**Run:** 2026-05-16 · **N:** 10 curated chiral targets with literature-documented biosynthesis · **Scope:** `ChiraLLM.route_predictor.predict_route` in isolation (LLM and BRENDA excluded for reproducibility)

## Headline numbers

| Metric | Value |
|---|---|
| Routes found | **90%** (9/10) |
| Terminal precursor matches literature | **10%** (1/10) |
| Step count within ±2 of literature | **70%** (7/10) |
| EC overlap with literature route (mean) | 0.11 enzymes |
| EC overlap normalized by expected (mean) | 5.6% |

**Read this honestly:** the predictor reliably finds *some* route through KEGG, but most predicted routes terminate at the wrong central metabolite. The one clean hit (L-DOPA → L-tyrosine via EC 1.14.16.2) is the literature pathway.

## Per-target results

| Target | Status | Found | Terminal match | Predicted terminal | Expected terminal(s) | EC overlap | Steps (pred/lit) |
|---|---|---|---|---|---|---|---|
| L-Lactate | ✓ | yes | no | (S)-malate | pyruvate | 0/2 | 1 / 1 |
| (R)-Pantolactone | ✓ | yes | no | L-lysine | KIV / pantoate | 0/2 | 2 / 3 |
| L-DOPA | ✓ | yes | **yes** | L-tyrosine | L-tyrosine | **1/2** | 1 / 1 |
| (R)-3-Hydroxybutyrate | ✓ | yes | no | α-ketoglutarate | acetyl-CoA | 0/3 | 2 / 3 |
| (R)-Mandelate | ✓ | yes | no | L-glutamate | phenylglyoxylate | 0/2 | 2 / 2 |
| (R)-Mevalonate | ✓ | yes | no | succinyl-CoA | acetyl-CoA | 0/3 | 2 / 3 |
| (R)-Noradrenaline | ✓ | yes | no | glycine | L-tyrosine | 0/3 | 2 / 3 |
| (S)-Naringenin | ✓ | yes | no | succinate | L-Phe / 4-coumarate | 0/5 | 1 / 5 |
| (S)-Norcoclaurine | KEGG gap | no | — | — | L-tyrosine | — | — / 4 |
| Coniferyl alcohol | ✓ | yes | no | L-methionine | L-Phe | 0/4 | 2 / 6 |

## Diagnosis: the "cofactor shortcut" failure mode

Inspecting the actual predicted step graphs reveals a consistent pattern. Most failures terminate by hopping through KEGG edges where the connectivity is mediated by a **shared cofactor**, not by real carbon flow:

| Target | Predicted shortcut | What's actually happening |
|---|---|---|
| (R)-Mandelate ← L-glutamate | C00025 → C00004 (NADH) → C01984 | Glutamate dehydrogenase produces NADH; mandelate dehydrogenase consumes NADH. Graph says "connected"; biology says these are two independent reactions sharing a redox cofactor. |
| (R)-Noradrenaline ← glycine | C00037 → C00021 (SAH) → C00547 | Glycine methyltransferase produces S-adenosyl-homocysteine; phenylethanolamine N-methyltransferase consumes the corresponding SAM/SAH pair. Methyl-donor connectivity, not carbon flow. |
| Coniferyl alcohol ← L-methionine | C00073 → C00019 (SAM) → C00590 | Same pattern — SAM acts as a methyl-donor hub. |
| (S)-Naringenin ← succinate | C00042 → C00509 via naringenin dioxygenase | Succinate is a *byproduct* of the 2-OG-dependent dioxygenase reaction, not a precursor. The graph edge is directionally misleading. |
| L-Lactate ← (S)-malate | C00149 → C00186 via malolactic enzyme | This one is defensible — malolactic fermentation in lactic acid bacteria is a real pathway, just not the textbook one. |

**Root cause confirmed in code.** `route_predictor.COFACTOR_SKIP_IDS` currently excludes only H+, H2O, O2. NADH/NADPH, ATP/ADP, CoA, SAM, FAD, and S-adenosyl-homocysteine are all valid graph edges in the current implementation, which lets the A* search use them as connectivity hubs.

## What this means

The Tier 1 algorithm (weighted A* over KEGG with a Tanimoto chemical-distance heuristic) is doing exactly what was specified, but the search space has an unstated assumption: **carbon-flow connectivity ≠ shared-cofactor connectivity**. With only the gas-and-water cofactor exclusion, KEGG's reaction graph is densely connected through cofactor recycling, and the shortest-cost path is often a 2-hop shortcut through one of those hubs.

This is a real, quantifiable calibration gap surfaced by empirical testing — exactly the kind of finding the benchmark exists to catch. The fix is well-understood (extend `COFACTOR_SKIP_IDS` to include the major redox/energy/methyl cofactors, and/or weight cofactor edges down rather than skipping them).

## What's defensible to claim today

- **Route discovery is reliable** — 90% of targets yield a route through KEGG; step counts are roughly biological (70% within ±2 of literature).
- **The system works correctly on tightly-coupled redox steps** — L-DOPA recovers cleanly because the literature path (L-Tyr → L-DOPA via tyrosine hydroxylase) is a single, cofactor-clean enzymatic step.
- **The cost decomposition is transparent** — a chemist can read the predicted route and immediately see *why* a wrong terminal was chosen (the cofactor edge is visible in the EC list).

## What's NOT defensible to claim today

- That Tier 1 predicts biologically natural routes by default — for multi-step targets that depend on natural-product carbon-flow rules, it doesn't yet.
- That the EC list of the predicted route is reliable for downstream BRENDA enantioselectivity scoring — most predicted EC numbers were NOT in the literature route.

## Next calibration target

Extend `COFACTOR_SKIP_IDS` (or introduce a cofactor-edge penalty) to break the connectivity hubs:
- Redox: NADH/NAD+ (C00004/C00003), NADPH/NADP+ (C00005/C00006), FADH2/FAD (C01352/C00016)
- Energy: ATP/ADP (C00002/C00008), PPi/Pi (C00013/C00009), AMP (C00020)
- Carriers: CoA (C00010), acetyl-CoA carrier carbons are real (don't skip C00024), but free CoA released as a byproduct should be ignored
- Methyl: SAM (C00019), SAH (C00021)

This is the obvious next sprint and is well-scoped — re-running this same benchmark after the change will give a clean before/after comparison.

## Reproducing this report

```bash
python3 -m benchmarks.run_benchmark
# Writes benchmarks/results/YYYY-MM-DD.json
```

Each per-target dict includes the full step-by-step predicted route, so any number in this report can be drilled into. Re-running on a fully warm cache takes ~30 sec; cold cache takes ~30 min (KEGG fetches dominate).
