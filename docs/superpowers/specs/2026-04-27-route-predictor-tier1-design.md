# Route Predictor Tier 1 — Design Spec

**Status:** Draft for review
**Date:** 2026-04-27
**Author:** Brainstormed in collaboration with the user via `superpowers:brainstorming`
**Sprint:** 1 of the 4-sprint roadmap (see `MEMORY.md` → `feedback_assistant_role.md` for full roadmap)
**Next skill:** `superpowers:writing-plans` upon approval

---

## 1. Context — Why this module exists

ChiralAI today is a **target-validation pipeline**: GPT proposes named molecules, KEGG/BRENDA/RDKit/COBRApy annotate them, the scorer ranks. It is not a biosynthesis predictor — it does not reason backward from a target to enumerate enzymatic steps.

This module, `ChiraLLM/route_predictor.py`, adds the missing capability for compounds KEGG already covers: given a target KEGG compound ID, walk the reaction graph backward via a best-first search, terminating at curated central-carbon metabolites, returning ranked routes. Every step gets enzyme EC numbers that flow into the existing BRENDA enrichment, and the scorer is extended to compute **whole-route ee composition** as a multiplicative composite — a per-route metric no other tool produces.

This is **Tier 1 of the route-predictor roadmap.** Tier 2 (RetroRules SMARTS search via RDKit `RunReactants`) and Tier 3 (RetroPath2 wrapper integration) are deferred to subsequent sprints. Tier 1's value is end-to-end demonstration and zero new pip dependencies; Tier 2's is novel-target retrobiosynthesis when KEGG doesn't already cover the target.

---

## 2. Architecture

### 2.1 Pipeline insertion point

```
main.py
   │
   ├── ask_gpt_chirality()                       (existing — proposes 5 candidates)
   ├── validate_chirality()                      (existing — RDKit R/S check)
   ├── query_kegg(compound_id)                   (existing — flat-file parser)
   │
   ├── predict_route(compound_id, mode='top_n')  ← NEW: ChiraLLM/route_predictor.py
   │       returns RouteResult dict
   │       ↓
   ├── query_enantioselectivity_batch(           (existing — but now called with all ECs
   │       all ECs across all routes — deduped)   collected across every step in every route)
   │
   ├── for each route:                            (existing primitive, NEW call site)
   │       check_feasibility(route.terminal_precursor_id)
   │
   ├── score_suggestion(...)                     (existing — extended for whole-route ee)
   └── save_suggestions_to_csv(...)              (existing — extended for route columns)
```

**Why predict_route runs before BRENDA, not after KEGG-only:** The predictor's output is the *full set of ECs across every step of every route*. Without the predictor, BRENDA was only queried for the target compound's ECs — missing every intermediate's enzymes. Inserting the predictor here converts BRENDA from "ee for the final enzyme" into "ee at every step of the biosynthesis," which is what the whole-route ee composition needs.

**Why `check_feasibility` keeps its current signature (`compound_id` only):** It's a reusable primitive — FBA on a single compound. Renaming it to know about routes would leak the predictor's domain model into a module that should remain route-agnostic. The change is at the **call site** in `main.py`: instead of one call on the target, iterate `for route in routes: check_feasibility(route.terminal_precursor_id)`.

### 2.2 Single-file module layout (Approach A)

Approach A (single rich module) was selected over Approach B (decomposition by responsibility) because the genuinely-swappable seam — the cost function — is one of many internal pieces, none of which have demonstrated reuse pressure today. Function-level seams inside one file give us refactor flexibility without paying the cross-file abstraction tax now.

```
ChiraLLM/route_predictor.py   (~450 LOC including docstrings)
   ├── Constants (top of file — the scientific contract)
   ├── KEGG layer:        _fetch_kegg_reaction, _fetch_kegg_mol, _parse_reaction_equation
   ├── Cost model:        _fetch_delta_g_kj_per_mol, _is_industrially_reversible,
   │                      _compute_edge_cost
   ├── Heuristic:         _get_central_fingerprints (lazy), _tanimoto_to_central
   ├── Search:            _astar_search
   ├── Output formatters: _extract_top_n, _extract_full_tree,
   │                      _extract_shortest_plus_diverse
   ├── Dataclasses:       Route, RouteStep, RouteResult
   └── Public API:        predict_route
```

### 2.3 Dependencies

**Net new pip dependencies: zero.** Existing `rdkit`, `requests`, `functools` (stdlib), `dataclasses` (stdlib), `pathlib` (stdlib) cover everything.

**External services:**
- KEGG REST (`rest.kegg.jp`) — already used by `database_validator.py`, free, 3 req/sec, 30-day disk cache.
- eQuilibrator REST (`equilibrator.weizmann.ac.il/api/`) — academic free service, no SLA. Disk-cached. Graceful degradation to constant fallback when unreachable.

**Caching:**
- In-process LRU (`functools.lru_cache(maxsize=4096)`) on `_fetch_kegg_reaction` and `_fetch_kegg_mol`.
- Disk cache at `~/.cache/chiralai/{kegg,equilibrator,kegg_mol}/{id}.{json,mol}` with 30-day TTL.
- TTL overridable via `CHIRALAI_CACHE_TTL_DAYS` env var.
- Manual invalidation via `python -m ChiraLLM.route_predictor --clear-cache`.

---

## 3. Components

### 3.1 Constants (the scientific contract — top of file)

```python
CENTRAL_METABOLITES: dict[str, str] = {
    # TCA cycle
    "C00022": "pyruvate", "C00024": "acetyl-CoA",
    "C00036": "oxaloacetate", "C00149": "(S)-malate",
    "C00122": "fumarate", "C00042": "succinate",
    "C00091": "succinyl-CoA", "C00026": "alpha-ketoglutarate",
    "C00311": "isocitrate", "C00158": "citrate",
    # Glycolysis / PPP
    "C00031": "D-glucose", "C00092": "glucose-6-phosphate",
    "C00085": "fructose-6-phosphate", "C00354": "fructose-1,6-bisphosphate",
    "C00111": "DHAP", "C00118": "G3P",
    "C00197": "3-phosphoglycerate", "C00074": "PEP",
    "C00117": "ribose-5-phosphate", "C00199": "ribulose-5-phosphate",
    # Amino acids (20 proteinogenic)
    "C00041": "L-alanine", "C00037": "glycine", "C00065": "L-serine",
    "C00188": "L-threonine", "C00097": "L-cysteine", "C00073": "L-methionine",
    "C00407": "L-isoleucine", "C00123": "L-leucine", "C00183": "L-valine",
    "C00079": "L-phenylalanine", "C00082": "L-tyrosine", "C00078": "L-tryptophan",
    "C00135": "L-histidine", "C00148": "L-proline", "C00064": "L-glutamine",
    "C00025": "L-glutamate", "C00049": "L-aspartate", "C00152": "L-asparagine",
    "C00047": "L-lysine", "C00062": "L-arginine",
    # Branched-chain amino acid intermediates (defensibly central — produced from pyruvate
    # via the BCAA biosynthesis pathway, used as nodes for valine/leucine/isoleucine and
    # for downstream secondary metabolism)
    "C00141": "alpha-ketoisovalerate",
}

INDUSTRIAL_REVERSIBLE_EC_PREFIXES: list[str] = [
    "1.1.1.",      # KREDs / aldo-keto reductases
    "2.6.1.",      # transaminases
    "1.5.1.",      # IREDs (imine reductases)
    "1.6.99.1",    # Old Yellow Enzyme (ene-reductases)
    "3.1.1.",      # lipases
    "1.14.13.22",  # cyclohexanone monooxygenase (Baeyer-Villiger archetype)
]

DEFAULT_BUDGET = 500           # max nodes A* will expand
DEFAULT_DEPTH_CAP = 8          # safety net only — budget terminates first in practice
DEFAULT_MAX_ROUTES = 3
THERMO_PENALTY_PER_KJ = 0.05   # cost units per kJ/mol when traversing reverse
FALLBACK_DELTA_G_KJ = 5.0      # used when eQuilibrator unreachable (median |ΔG| magnitude)

# v1 STARTING GUESS — NOT VALIDATED. h(n) Tanimoto distance is scaled by this weight before
# being added to g(n) accumulated edge cost. Weight=2.0 means structural similarity to a
# central metabolite matters 2x as much as accumulated thermo+directional cost (g(n) is on
# the order of ~1 per step). Re-tune against integration fixtures.
# Re-tuning trigger: if top-3 routes for (R)-pantolactone (C00599) do not include the
# KIV-via-ketopantoate-hydroxymethyltransferase route, this constant is too high or too low.
TANIMOTO_HEURISTIC_WEIGHT = 2.0
```

These two dicts (`CENTRAL_METABOLITES`, `INDUSTRIAL_REVERSIBLE_EC_PREFIXES`) are the scientific contract a wet-lab reviewer reads first. Every entry is a domain decision a chemist would defend; both are extensible by editing the constants.

### 3.2 Function inventory

| Function | Signature | Role |
|----------|-----------|------|
| `_fetch_kegg_reaction(rxn_id)` → `dict\|None` | LRU + disk cached. Parsed reaction with substrates/products/EC/direction. |
| `_fetch_kegg_mol(compound_id)` → `Chem.Mol\|None` | LRU + disk cached. RDKit Mol from KEGG MOL file. |
| `_parse_reaction_equation(equation)` → `(substrates, products, direction)` | Pure function. Splits KEGG equation strings. |
| `_fetch_delta_g_kj_per_mol(rxn_id)` → `float\|None` | eQuilibrator REST + disk cache. None on network failure. |
| `_is_industrially_reversible(ec_numbers)` → `bool` | Lookup against `INDUSTRIAL_REVERSIBLE_EC_PREFIXES`. |
| `_compute_edge_cost(reaction, traversed_direction)` → `dict` | Returns breakdown: `{base, thermodynamic, directionality, industrial_override, total}`. |
| `_get_central_fingerprints()` → `dict[str, ExplicitBitVect]` | Lazy-init Morgan fingerprints for `CENTRAL_METABOLITES`. |
| `_tanimoto_to_central(compound_id)` → `float` | Returns `1 - max(Tanimoto)` against central set. Range [0,1]. |
| `_astar_search(target_id, budget, depth_cap)` → `dict` | Best-first search. Returns DAG. |
| `_extract_top_n(dag, target_id, n)` → `list[Route]` | Backtrack + sort by cost ascending. |
| `_extract_full_tree(dag, target_id)` → `Route` | Complete DAG as serializable nested structure. |
| `_extract_shortest_plus_diverse(dag, target_id, n_diverse)` → `list[Route]` | Shortest hop + N maximally-different alternates. |
| `predict_route(compound_id, mode, n, budget)` → `RouteResult` | Public entrypoint. Never raises. |

### 3.3 Dataclasses

```python
@dataclass
class RouteStep:
    reaction_id: str
    ec_numbers: list[str]                      # all ECs from KEGG
    precursor_id: str                          # the upstream compound (one step closer to terminal)
    intermediate_id: str                       # the downstream compound (we just came from this)
    edge_cost_breakdown: dict[str, float]      # base, thermo, directionality, industrial, total
    traversed_direction: str                   # 'forward' or 'reverse'
    # Naming note: in retrosynthesis we walk product→substrate, so "substrate" and "product"
    # would be ambiguous about which direction we mean. precursor_id (upstream) and
    # intermediate_id (downstream) are unambiguous about graph position regardless of how
    # the chemical reaction itself is described in KEGG.

@dataclass
class Route:
    target_id: str
    steps: list[RouteStep]                    # ordered target → terminal precursor
    terminal_precursor_id: str
    terminal_precursor_name: str               # from CENTRAL_METABOLITES lookup
    total_cost: float
    cost_breakdown: dict[str, float]           # summed components across all steps
    warnings: list[str]

@dataclass
class RouteResult:
    target_id: str
    mode: str                                  # 'top_n' | 'full_tree' | 'shortest_plus_diverse'
    routes: list[Route]                        # always a list, even for full_tree
    nodes_explored: int
    budget_exhausted: bool
    warnings: list[str]                        # global (not per-route) warnings
    status: str                                # see Section 5 status taxonomy
```

---

## 4. Data flow (worked example)

Target: **(R)-pantolactone** (KEGG `C00599`), `mode="top_n"`, `n=3`, `budget=500`.

1. **Entry validation**: `C00599` matches `C\d{5}` ✓. KEGG ID format check passes.
2. **A\* init**: priority queue seeded with `(f=0.0, depth=0, "C00599", path=[])`. Empty `visited_dag`.
3. **Main loop iteration**: pop cheapest. Not in `CENTRAL_METABOLITES` → expand.
   - Fetch reactions for `C00599` from KEGG (cache hit, since `query_kegg` ran for the target earlier in `main.py`).
   - For each reaction: parse equation, determine which side `C00599` is on. If on product side → `traversed_direction = "forward"`; if substrate side → `traversed_direction = "reverse"`.
   - For each precursor on the relevant side: compute `edge_cost`, compute `f = g + h`, push to queue, record edge in `visited_dag`.
4. **Termination**: continues until `nodes_explored >= 500`, or queue empty, or all leaves popped. **Search continues past the first central-metabolite hit** — this is intentional, to surface non-obvious industrial routes that thermodynamically-cheap routes would mask.
5. **Formatter dispatch**: `mode="top_n"` → `_extract_top_n(dag, "C00599", n=3)` backtracks from each terminating central-metabolite leaf to target, sorts by `total_cost`, returns top 3.
6. **RouteResult assembly**: dataclass with routes, nodes explored, budget exhaustion flag, warnings, status.
7. **Downstream consumption in `main.py`**:
   - `suggestion["route_prediction"] = asdict(result)`
   - `all_ecs = {ec for route in result.routes for step in route.steps for ec in step.ec_numbers}` (deduplicated)
   - `suggestion["brenda_data"] = query_enantioselectivity_batch(list(all_ecs))`
   - `suggestion["route_feasibility"] = [check_feasibility(r.terminal_precursor_id) for r in result.routes]`
8. **Scorer extension** — multiplicative whole-route ee:
   ```python
   for route in routes:
       step_ees = [best_ee_for_ecs(step.ec_numbers, brenda_data) for step in route.steps]
       if all(e is not None for e in step_ees):
           composed_ee = math.prod(e / 100.0 for e in step_ees) * 100.0
   ```
   This computed per-route. **Assumes step independence** (see §7 known limitations).
9. **Persistence**: `file_saver.py` extension adds flat columns `route_top1_step_count`, `route_top1_terminal_precursor`, `route_top1_total_cost`, `route_top1_composed_ee`, `route_top1_warnings`. Full route trees go to JSON sidecar.

---

## 5. Error handling

**Contract: `predict_route` never raises.** All failures encode as `RouteResult.status` + structured warnings.

### 5.1 Failure mode matrix

| # | Failure | Detection | Returned status | Warning |
|---|---------|-----------|-----------------|---------|
| 1 | `compound_id` is None/empty | input validation | `no_kegg_id` | none |
| 2 | Malformed KEGG ID (not `C\d{5}`) | regex check | `invalid_kegg_id` | format violation message |
| 3 | KEGG REST 404 on target | HTTP status | `target_not_in_kegg` | "compound may have been deprecated or merged" |
| 4 | KEGG REST 5xx / network timeout (per reaction) | requests exception | `success` (partial) | "KEGG unreachable for reactions: [...]" — single retry after 1s delay first; if both fail, mark reaction unfetched and skip |
| 5 | Target compound has empty REACTION field | `kegg_entry["reactions"] == []` | `target_has_no_reactions` | "compound may be a leaf metabolite or unconnected" |
| 6 | Reaction equation parse fails | regex no match | `success` | "Could not parse equation for {rxn_id}" |
| 7 | All branches exhausted, no central reached | queue empty, no leaves | `no_route_found` | "Search exhausted N nodes without reaching a central metabolite" |
| 8 | Budget exhausted, ≥1 leaf found | `nodes_explored >= budget` | `success` (partial) | "Budget of {budget} nodes exhausted; some routes may be missing" |
| 9 | Depth cap hit | `depth > DEFAULT_DEPTH_CAP` | `success` | aggregated if ≥10 occurrences |
| 10 | eQuilibrator unreachable | HTTP timeout / 5xx | `success` | "eQuilibrator unreachable for N reactions; fallback ΔG used" |
| 11 | KEGG MOL file missing | empty response body | `success` | "MOL unavailable for N compounds; Tanimoto heuristic defaulted" — heuristic degrades to uniform 0.5 |
| 12 | RDKit fails to parse MOL | `MolFromMolBlock` returns None | `success` | "RDKit could not parse MOL for {cid}" |
| 13 | Cycle in DAG (revisited compound) | `visited_dag` lookup | (silent — expected) | none |
| 14 | Unrecognized `mode` argument | string match | `invalid_mode` | "mode must be one of: top_n, full_tree, shortest_plus_diverse" |

### 5.2 Differential treatment rationale

- **Two-attempt retry for KEGG, single-shot for eQuilibrator.** KEGG is the only source of truth for the reaction graph; losing a reaction loses an entire branch. eQuilibrator is one input among three to a single edge cost; the fallback constant is acceptable.
- **MOL unavailability degrades the heuristic, doesn't exclude nodes.** A\* with an uninformative heuristic degenerates to Dijkstra (still correct, just less efficient). Excluding nodes would cripple search through any compound class with poor PubChem cross-referencing.

---

## 6. Testing strategy

### 6.1 Three test layers

**Layer 1 — Unit (no network, deterministic, ~10s).**
Lives in `tests/test_route_predictor_unit.py`. All external calls mocked. Covers every public function, every status branch, every edge cost component. Coverage target: ≥85%. Mocks use synthetic KEGG IDs (`C_A`, `R_synth_001`) to decouple test correctness from KEGG release cycle.

**Layer 2 — Integration (real KEGG, mocked eQuilibrator).**
Lives in `tests/test_route_predictor_integration.py`. Marked `@pytest.mark.integration` — skipped by default. Three real-target fixtures:
- `C00599` (R-pantolactone): happy path with industrial relevance.
- `C00186` (S-lactic acid): trivial route to pyruvate via L-LDH.
- `C09136` (S-norcoclaurine): KEGG has poor enzyme annotation — tests `target_has_no_reactions` / `no_route_found` against real data.

Asserts only structural properties (status, `nodes_explored < budget`, `terminal_precursor_id ∈ CENTRAL_METABOLITES`). Does **not** assert exact route content (KEGG drift).

**Layer 3 — End-to-end smoke (live everything, manual).**
Extends existing `smoke_test.py` with a `route_prediction` fixture. Run via `.venv/bin/python smoke_test.py --fixture route_prediction`.

### 6.2 CI configuration

Default `pytest` runs Layer 1 only (fast, deterministic, no network). Pytest markers control Layer 2 / Layer 3.

### 6.3 Test data location

```
tests/
├── test_route_predictor_unit.py
├── test_route_predictor_integration.py
└── fixtures/
    └── synthetic_kegg.json     # Layer 1 mock data, version-controlled
```

---

## 7. Known limitations (v1)

These are explicit scope decisions, not bugs.

1. **Multiplicative ee composition assumes step independence.** Real biocatalysis sometimes has correlated stereochemistry across steps — most notably **dynamic kinetic resolution (DKR)**, where an upstream racemization step combined with a downstream enantioselective step produces an artificially-high terminal ee. The multiplicative formula `∏(ee_i/100) * 100` underestimates DKR routes' true ee. Documented; v2 should add a correlated-step ee model.

2. **`_astar_search` is greedy best-first, not provable A\*.** Tanimoto-to-central is a structural-similarity heuristic, not a true cost lower bound — it can either overestimate or underestimate true cost-to-goal depending on the metabolic landscape. An overestimating heuristic loses the optimality guarantee A\* gives with admissible heuristics, but is often used deliberately because it's faster. We trade optimality for speed within the 500-node budget. Within that budget the algorithm finds good routes; outside it, the routes returned are best-effort, not provably best.

3. **eQuilibrator REST has no SLA.** Mitigations baked in: (a) graceful fallback to constant ΔG, (b) disk caching with 30-day TTL, (c) optional pre-warm script populating ~2700 KEGG reactions in iJO1366 + most-cited reactions. But on a cold cache with eQuilibrator down, route quality degrades.

4. **Tier 1 cannot find routes for compounds outside KEGG.** This is by design — Tier 2 (RetroRules + RDKit `RunReactants`) is the answer for novel-target retrobiosynthesis.

5. **`_compute_edge_cost` returns flat scalar weighting of components.** A wet-lab researcher who wants to bias toward thermodynamic feasibility (or away from it, for engineered routes) can edit constants but cannot supply a custom cost function. Not painful in v1; revisit in v2 if pre-tuning becomes needed.

---

## 8. Out of scope (deferred to subsequent sprints)

- **Tier 2 — RetroRules SMARTS search via RDKit `RunReactants`.** Sprint 2.
- **Tier 3 — RetroPath2 wrapper integration.** Sprint 3 if Tier 2 has correctness gaps.
- **Multi-host FBA.** `feasibility_checker.py` keeps iJO1366; iMM904 (yeast) and iJN1463 (P. putida) deferred.
- **AlphaFold-grounded engineering candidates.** Deferred to Sprint 3.
- **OpenTrons OT-2 wet-lab protocol export.** Deferred to Sprint 4.
- **Custom skills (CIP validator, BRENDA SOAP wrapper, literature grounding).** Authored separately via `superpowers:writing-skills` in Sprint 2.

---

## 9. Design decisions log

The user's substantive design contributions, captured for the record:

- **2026-04-27 — Curated central-carbon stop set** (over iJO1366-broad, depth-only, or iJO1366+depth-cap). Rationale: scientifically meaningful termination — chemist wants "starts from glucose," not "starts from random metabolite X."
- **2026-04-27 — Three output modes** (top-N, full tree, shortest+diverse) chosen over single-best. Rationale: surface non-obvious routes; downstream consumer chooses interpretation.
- **2026-04-27 — Reversibility is graded, not binary.** KEGG `<=>` vs `=>` markers parsed; reverse-traversal of `=>` gets a baseline directional cost. Cross-referenced with eQuilibrator ΔrG′ data when available — penalty proportional to `|ΔG|`. Industrial override set (KREDs, transaminases, IREDs, EREDs, lipases, BVMOs) drops reverse-direction penalty for these EC families. Output is a transparent feasibility breakdown, not a black-box scalar.
- **2026-04-27 — A\* with global node budget + Tanimoto heuristic** (over per-node branching limits). Rationale: best-first search with informed heuristic naturally concentrates compute on promising branches; single global budget gives easy targets fast finishes and hard targets bounded effort.
- **2026-04-27 — All ECs per step query BRENDA**, ee chosen per step. Rationale: most thorough; caching makes the call volume acceptable.
- **2026-04-27 — Search continues past first central-metabolite hit.** Rationale: the whole point of ChiralAI is non-obvious routes; first-hit termination defaults to the most thermodynamically favorable route, which is rarely the most industrially interesting.
- **2026-04-27 — EC deduplication across routes before BRENDA batch.** Rationale: BRENDA returns the same data regardless of which route asked; per-route queries waste API budget without adding information.
- **2026-04-27 — DKR limitation documented upfront.** Rationale: a wet-lab reviewer will catch this; acknowledging the simplification is more credible than having it discovered later.
- **2026-04-27 — Single-file module (Approach A)** over decomposed-by-responsibility (Approach B). Rationale: the genuinely-swappable seam (cost function) is one of many internal pieces; function-level seams inside one file give refactor flexibility without paying cross-file abstraction tax now.
- **2026-04-27 — eQuilibrator REST direct, not `equilibrator-api` package.** Rationale: package ships ~1GB of training data; REST is lightweight with disk-cache absorbing repeat hits.
- **2026-04-27 — `check_feasibility` keeps single-compound primitive signature.** Call site in `main.py` iterates over routes' terminal precursors. Rationale: avoid leaking route concept into the FBA module.
- **2026-04-27 — Lazy fingerprint init via `_get_central_fingerprints()` accessor.** Rationale: future-proofs against catalog swap (e.g., to Sigma compounds with hundreds of thousands of entries) where eager module-import cost would become seconds.
- **2026-04-27 (review pass)** — Field naming `precursor_id` / `intermediate_id` chosen over `substrate_id` / `product_id` in `RouteStep`. Rationale: in retrosynthesis we walk product→substrate, so chemistry-conventional names create directional ambiguity for code readers. Graph-position names are unambiguous regardless of how KEGG describes the underlying reaction.
- **2026-04-27 (review pass)** — KIV (C00141) re-justified as a branched-chain amino acid intermediate (defensibly central), not a pantolactone-specific shortcut. Comment in `CENTRAL_METABOLITES` rewritten to reflect this. Rationale: avoids the slippery slope of target-specific stop entries.
- **2026-04-27 (review pass)** — `TANIMOTO_HEURISTIC_WEIGHT = 2.0` flagged in source as v1 starting guess with a concrete re-tuning trigger (KIV-via-KPHMT route presence in pantolactone top-3). Rationale: prevents the constant from rotting as an "unvalidated forever" magic number.
- **2026-04-27 (review pass)** — DKR detection removed from acceptance criterion #5. Stays in §7.1 as documented limitation. Rationale: racemase detection has real edge cases (substrate-specific racemases, mutase confusion) that v1 cannot honestly solve. Better to ship limited and documented than half-shipped.

---

## 10. Acceptance criteria for Sprint 1

The sprint is complete when:

1. `ChiraLLM/route_predictor.py` exists with all 12 functions, 3 dataclasses, 2 constant dicts, public `predict_route` API.
2. `tests/test_route_predictor_unit.py` exists, ≥85% coverage on Layer 1, all tests passing under `pytest -m 'not integration'`.
3. `tests/test_route_predictor_integration.py` exists with the 3 real-target fixtures, all passing under `pytest -m integration` against live KEGG.
4. `main.py` calls `predict_route` after `query_kegg`, deduplicates ECs across routes, calls `query_enantioselectivity_batch` once with the full set, calls `check_feasibility` per route's terminal precursor.
5. `enantioselectivity_scorer.py` has a `_compose_route_ee(route, brenda_data)` helper computing multiplicative ee per route. DKR-suspect detection (e.g., upstream racemase) is **not** in scope for v1 — see §7.1 known limitation.
6. `utils/file_saver.py` adds the 5 flat columns for `route_top1_*`. JSON sidecar contains full `route_prediction` nested structure.
7. `smoke_test.py` extended with a `route_prediction` fixture. For the 3 KEGG-covered segments (codexis, pharma, academic), `predict_route` returns `status="success"` with at least one route. For the other 3 segments (arnold, ginkgo, adversarial — targets with no KEGG ID or no KEGG reactions), `predict_route` returns the appropriate error status (`no_kegg_id`, `target_not_in_kegg`, or `target_has_no_reactions`) without raising. All 6 segments must complete the fixture run cleanly.
8. `requirements.txt` is **unchanged** — no new pip dependencies.
9. Disk cache directory `~/.cache/chiralai/` is created at module import; at least one warm-cache pre-fetch script exists at `scripts/warm_equilibrator_cache.py`.
10. `CLAUDE.md` and `README.md` updated to reflect that route prediction is now a working module (not a roadmap item).

---

## 11. Spec self-review (inline)

**Placeholder scan:** No "TBD" or "TODO" tags in this spec. The v2-deferred items (DKR, custom cost functions, multi-host FBA, AlphaFold, OT-2) are explicitly scoped in §7 and §8, not handwaved.

**Internal consistency:** Section 2's pipeline ordering (`predict_route` between KEGG and BRENDA) matches Section 4's worked-example data flow and Section 10's acceptance criterion #4. Section 3's function inventory matches Section 4's worked example. Section 5's failure-matrix statuses match Section 3's `RouteResult.status` field type.

**Scope check:** This spec covers a single subsystem (route prediction) wired into existing modules with bounded edits (1 new module + 4 modules touched: `main.py`, `enantioselectivity_scorer.py`, `file_saver.py`, `smoke_test.py`). It does NOT bundle Tier 2 (separate sprint), custom skills (separate sprint via `writing-skills`), or pipeline-wide refactors. Implementable as a single plan.

**Ambiguity check:**
- "Multiplicative ee composition" is precisely defined in §4 step 8 with formula.
- "Industrial override" is precisely defined: drops directionality penalty for ECs matching prefixes in `INDUSTRIAL_REVERSIBLE_EC_PREFIXES`.
- "DKR detection in scorer" (acceptance #5) clarified inline as heuristic only — upstream racemization-class detection, not full DKR analysis.
- "Pre-warm script populates ~2700 reactions" (§7.3) — clarified as iJO1366 reactions (a concrete, knowable set) plus optionally others.

No contradictions found. No ambiguity remains.
