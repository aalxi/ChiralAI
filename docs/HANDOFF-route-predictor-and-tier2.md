# Hand-off: Route Predictor (Tier 1) → Tier 2 Retrobiosynthesis

**Audience:** the next model or developer continuing route-prediction work — either (a) further Tier 1 calibration, or (b) building Tier 2 novel-target retrobiosynthesis.

**Last updated:** 2026-05-23. Read [`CLAUDE.md`](../CLAUDE.md) first for the project brief; this doc is the route-predictor-specific deep context that doesn't belong inline there.

---

## 1. Where things stand

Tier 1 route prediction is **working and empirically benchmarked**. The full pipeline (GPT → chirality → KEGG → route → BRENDA → FBA → scorer → CSV) is wired in `main.py`. The route predictor is the only module with a standalone validation harness (`benchmarks/`).

Two calibration passes have landed, both driven by the benchmark surfacing its own next bug:

1. **Cofactor-hub fix (commit `d025bfd`, 2026-05-17):** extended `COFACTOR_SKIP_IDS` from 3 → 17 entries (redox carriers, energy/phosphate, methyl donors, free CoA). Eliminated the failure mode where A* exploited NADH/SAM/ATP as graph-connectivity hubs to reach biologically wrong central metabolites.
2. **Pseudo-compound fix (this hand-off, 2026-05-23):** extended `COFACTOR_SKIP_IDS` to 30 entries — added generic redox placeholders ("Acceptor" C00028 / "Reduced acceptor" C00030), electron carriers (ferredoxin, thioredoxin, cytochrome c, ubiquinone/ubiquinol, e-), and the generic "Protein" placeholder C00017. Targets the KEGG pseudo-compound hub failure mode named in the 05-17 report (mandelate routing through C00030).

**Benchmark numbers: see `benchmarks/REPORT.md` for the authoritative, current table.** Do not trust any percentage hard-coded elsewhere; REPORT.md is regenerated from the latest `benchmarks/results/*.json`.

---

## 2. Tier 1 architecture you are inheriting

`ChiraLLM/route_predictor.py`, ~1000 lines. Weighted A* / greedy best-first **backward** search through the KEGG reaction graph, from a target compound to curated central metabolites.

### The scientific contract (first ~140 lines — auditable by a wet-lab reviewer)

Two domain-commitment constants, deliberately at the top of the file:

- **`CENTRAL_METABOLITES`** (~41 nodes): the stop set. Full TCA, glycolysis, PPP, the 20 amino acids, KIV. Backward search terminates here — these are the biological anchors a route should bottom out at.
- **`COFACTOR_SKIP_IDS`** (30 nodes): compounds the traversal refuses to step *through* as if they were carbon-carrying intermediates. Cofactors, redox/electron carriers, and pseudo-compounds. **Acyl-CoA species (acetyl-CoA, succinyl-CoA) are deliberately NOT in this set** — they carry real biosynthetic carbon. Free CoA (C00010) is.

Two more tunables:
- **`INDUSTRIAL_REVERSIBLE_EC_PREFIXES`**: EC classes (KREDs, transaminases, IREDs, OYE, lipases, BVMOs) whose reverse-direction penalty is dropped to ~0, because they are routinely run "uphill" in industrial biocatalysis. Dual-mode matching: trailing `.` = class prefix; no trailing `.` = exact EC.
- **`TANIMOTO_HEURISTIC_WEIGHT = 2.0`**: **NOT validated.** A v1 starting guess. The A* heuristic `h(n)` is `WEIGHT × Tanimoto distance to nearest central metabolite`. This is the single most impactful un-tuned knob — see residual failure mode #3.

### The cost model (transparent, per-step, auditable)

Each `RouteStep.edge_cost_breakdown` decomposes into:
`base (1.0) + thermodynamic (eQuilibrator ΔG × 0.05/kJ) + directionality (KEGG ⇌ vs → reverse penalty) + industrial-reversibility override`.
Total route cost is the sum. A user can read the breakdown and disagree with any term.

### The data contract — Tier 2 MUST reproduce this

```
RouteStep(reaction_id, ec_numbers, precursor_id, intermediate_id,
          edge_cost_breakdown, traversed_direction)
Route(target_id, steps, terminal_precursor_id, terminal_precursor_name,
      total_cost, cost_breakdown, warnings)
RouteResult(target_id, mode, routes, nodes_explored, budget_exhausted,
            warnings, status)
```

`predict_route()` **never raises** — every failure path returns a `RouteResult` with a `status` string (`success`, `no_kegg_id`, `invalid_kegg_id`, `target_not_in_kegg`, `target_has_no_reactions`, `no_route_found`, `invalid_mode`). The scorer and `file_saver` consume these dataclasses. **If Tier 2 emits the same `Route`/`RouteStep` objects, the entire downstream pipeline — BRENDA ee annotation, FBA, scoring, CSV — works unchanged.** That is the integration design goal.

### Caching

`_success_only_cache` decorator (custom, not `functools.lru_cache`) on the three KEGG fetchers — only memoizes truthy results, so transient network failures don't poison the cache. Disk cache for KEGG reactions, MOL files, and eQuilibrator ΔG at `~/.cache/chiralai/` (30-day TTL). The eQuilibrator pre-warm script (`scripts/`, commit `a4b908a`) seeds iJO1366 reaction ΔGs into this cache.

---

## 3. Tier 2 — what "ready for B" means

**Goal:** given a target **SMILES** that KEGG does not know, enumerate enzymatic disconnections backward to central metabolites via reaction-template (SMARTS) matching. This is the genuinely novel differentiator — no OSS tool does chiral-aware biocatalytic retrosynthesis with provenance.

### The hook point — already exists

`predict_route()` returns `status="target_not_in_kegg"` when KEGG has no reactions for the compound. **That is the exact fallback trigger for Tier 2.** Orchestration logic (in `main.py`, kept thin): if Tier 1 returns `target_not_in_kegg` or `no_route_found`, hand the target SMILES to Tier 2. Both tiers return `RouteResult`; the caller does not care which produced it.

### Recommended build (from the retrobiosynthesis-landscape research)

- **RetroRules SQLite** (~50MB, CC-licensed, ~6M reaction SMARTS indexed by EC + radius) is the raw asset. Download once, query locally — do NOT take the RetroPath2/KNIME 400MB dependency unless Tier 2's own search proves insufficient (that's Tier 3).
- **RDKit `RunReactants`** applies a reaction SMARTS backward to a target mol → candidate precursor sets. Recurse to the `CENTRAL_METABOLITES` stop set.
- Per step: RDKit CIP check (reuse `chirality_checker`), BRENDA ee (reuse `brenda_client`), terminal feasibility (reuse `feasibility_checker`).
- Rank routes by ∑ee × terminal feasibility flux.

### What Tier 2 reuses from Tier 1 (and the refactor this implies)

Tier 2 should import, not duplicate: `CENTRAL_METABOLITES`, `COFACTOR_SKIP_IDS`, the Morgan-fingerprint Tanimoto heuristic, the `_success_only_cache` pattern, and the `Route`/`RouteStep`/`RouteResult` dataclasses.

**Refactor recommendation:** before Tier 2 grows, lift these shared pieces out of `route_predictor.py` into a small `ChiraLLM/route_common.py` (constants + dataclasses + heuristic). This is the one place the two tiers would otherwise collide on the same file. Doing it first makes the two work-streams cleanly file-disjoint (see §5).

### COBRApy's role does NOT change

FBA validates the terminal precursor's producibility *after* a route is proposed. It cannot enumerate routes. Do not try to use it for search.

---

## 4. Environment caveats — READ BEFORE RE-RUNNING THE BENCHMARK

- **eQuilibrator SSL failure (active as of 2026-05-23).** `equilibrator.weizmann.ac.il` is returning `CERTIFICATE_VERIFY_FAILED` (CA-bundle / proxy issue in this environment). Reactions already in the disk cache (pre-warmed iJO1366 set) still resolve; **uncached reactions fall back to `FALLBACK_DELTA_G_KJ = 5.0`**, i.e. the thermodynamic cost term goes inert for them. Consequence: benchmark runs are slow (per-reaction read timeouts) AND the thermo term is partially degraded. When comparing benchmark runs across dates, confirm thermo state is consistent (compare `top_route_total_cost` on a target whose route is identical between runs) before attributing any delta to a code change. This bit us once — the 05-17 → 05-23 comparison required this exact check.
- **BRENDA: no credentials in this environment.** Every ee falls back to `llm_claim`. The benchmark deliberately excludes BRENDA (route_predictor in isolation) so this doesn't affect route numbers — but it does mean the full-pipeline ee column is unverified here.

---

## 5. Git workflow for a parallel hand-off

**The rule that matters:** two models must not commit to the same branch concurrently — that collides like two people editing one doc (it's what previously left a worktree 32 commits behind).

- `main` is current at the pseudo-compound-fix commit.
- The existing `worktree-tier1-route-predictor` branch is **stale** (32 commits behind). Don't reuse it as-is; reset it to current `main` first, or cut a fresh branch.

**Recommended split:**
- **You / Tier 2** work on `main` (or a `tier2-retrobiosynthesis` branch).
- **Hand off continued Tier 1 calibration** on its own fresh branch cut from current `main` (e.g. a new worktree `tier1-calibration`).
- Because Tier 1 (`route_predictor.py`) and Tier 2 (new module) touch different files, the branches **merge back cleanly via PR** — the only shared surface is the constants/dataclasses, which the `route_common.py` refactor (§3) removes as a conflict point.

If you do the refactor first and commit it to `main`, both downstream branches start from a clean shared base.

---

## 6. How to run things

```bash
# Unit tests (no network) — 75 tests, ~2s
python3 -m pytest tests/test_route_predictor_unit.py -q

# Integration tests (live KEGG)
python3 -m pytest tests/test_route_predictor_integration.py -q

# Benchmark — writes benchmarks/results/YYYY-MM-DD.json
# Cold/SSL-degraded: many minutes. Warm cache, eQuilibrator reachable: ~30s.
python3 -m benchmarks.run_benchmark

# Then regenerate the narrative report from the new JSON (manual edit of REPORT.md).
```

The benchmark scores each target on: route_found, terminal_match (terminal in `expected_terminals`), ec_overlap_count, step_within_2. Ground truth is in `benchmarks/targets.py` (10 curated chiral targets with literature citations).

---

## 7. The three residual failure modes (Tier 1's next calibration targets)

Characterized in `benchmarks/REPORT.md`, ranked by tractability:

1. **KEGG pseudo-compound hubs** — *partially addressed by the 05-23 fix.* Re-check whether mandelate now routes correctly; if other placeholder IDs surface, extend `COFACTOR_SKIP_IDS` further.
2. **Reaction-direction semantics** — byproducts treated as precursors (e.g. succinate from a 2-OG-dependent dioxygenase cosubstrate cycle read as a carbon source). Harder: needs curated per-reaction direction annotations or a cosubstrate-byproduct heuristic.
3. **Tanimoto heuristic preferring chemically-similar over biologically-natural terminals** — e.g. mevalonate via the shorter succinyl-CoA path instead of the literature acetyl-CoA path. Retune `TANIMOTO_HEURISTIC_WEIGHT` (currently 2.0, unvalidated) against the benchmark, or downweight it when alternate routes reach the same intermediate.
