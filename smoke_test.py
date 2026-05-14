"""
ChiralAI Smoke Test
===================
Runs queries through the full pipeline and evaluates output quality against
scientific pass criteria. Five fixture sets model real target user segments:
  codexis   — directed evolution / biocatalysis CRO sanity-check molecules
  arnold    — new-to-nature enzyme reconnaissance + graceful failure
  pharma    — Merck/BMS/Pfizer route-scouting molecules
  ginkgo    — full-stack strain engineering (where COBRApy matters most)
  academic  — teaching-grade "hello world" biomanufacturing molecules
  adversarial — exotic/non-natural: tests graceful fallback to LLM-claim

Usage:
    .venv/bin/python smoke_test.py --mock                  # Codexis fixture (default)
    .venv/bin/python smoke_test.py --mock --all            # all five user-segment fixtures
    .venv/bin/python smoke_test.py --mock --fixture arnold # single named fixture
"""
import argparse
import json
import os
import sys
from datetime import datetime
from dotenv import load_dotenv

load_dotenv()

# ── Pass/fail criteria ────────────────────────────────────────────────────────

def check_stereospecific_smiles(s):
    smiles = s.get("SMILES", "")
    ok = "@" in smiles
    return ok, f"SMILES={'...' + smiles[-20:] if len(smiles)>20 else smiles!r}"

def check_chirality_assigned(s):
    centers = (s.get("chirality_validation") or {}).get("chiral_centers", [])
    if not centers:
        return False, "no chiral_centers reported"
    unassigned = [c for c in centers if c[1] == "?"]
    ok = len(unassigned) == 0
    return ok, f"centers={centers}, unassigned={unassigned}"

def check_score_not_zero(s):
    score = (s.get("scoring") or {}).get("composite_score", 0)
    ok = score > 0.0
    return ok, f"composite_score={score}"

def check_ec_number_present(s):
    ec = (s.get("scoring") or {}).get("top_enzyme", {}).get("ec_number")
    ok = ec is not None and ec != "unknown"
    return ok, f"top_enzyme.ec_number={ec!r}"


PER_SUGGESTION_CHECKS = [
    ("stereospecific_SMILES",   check_stereospecific_smiles),
    ("chirality_assigned",      check_chirality_assigned),
    ("score_nonzero",           check_score_not_zero),
    ("ec_number_present",       check_ec_number_present),
]

# ── Fixtures ──────────────────────────────────────────────────────────────────
# Five segments from the target user analysis + one adversarial regression set.
# All SMILES CIP assignments verified with RDKit. KEGG IDs and ECs live-confirmed.
MOCK_SUGGESTIONS = {

    # ── Segment 1: Codexis — directed evolution / biocatalysis CRO ───────────
    # Sanity-check molecules every industrial biocatalysis shop would test first.
    # (R)-3-hydroxybutyrate: KEGG C01089, EC 1.1.1.30 — canonical KRED benchmark.
    # (S)-DMPEA: sitagliptin/montelukast structural family — tests TA retrieval.
    "codexis": [
        {
            "name": "(R)-3-hydroxybutyrate",
            "SMILES": "C[C@@H](O)CC(=O)[O-]",
            "R_S_config": "R",
            "KEGG_ID": "C01089",
            "enzyme_class": "3-hydroxybutyrate dehydrogenase",
            "known_ee": ">98% ee",
            "applications": "PHA polymer precursor; canonical KRED sanity-check molecule"
        },
        {
            "name": "(S)-1-(3,4-dimethoxyphenyl)ethanamine",
            "SMILES": "COc1ccc([C@H](C)N)cc1OC",
            "R_S_config": "S",
            "KEGG_ID": "C05587",
            "enzyme_class": "monoamine oxidase / omega-transaminase",
            "known_ee": ">99% ee",
            "applications": "Sitagliptin/montelukast structural family; TA engineering target"
        },
        {
            "name": "(R)-pantolactone",
            "SMILES": "O=C1OC[C@@](O)(C)C1",
            "R_S_config": "R",
            "KEGG_ID": "C01012",
            "enzyme_class": "carbonyl reductase",
            "known_ee": "99.5% ee",
            "applications": "Pantothenic acid (vitamin B5) synthesis intermediate"
        },
        {
            "name": "(S)-1-phenylethanol",
            "SMILES": "[C@@H](O)(c1ccccc1)C",
            "R_S_config": "S",
            "KEGG_ID": "C02934",
            "enzyme_class": "ketoreductase (KRED)",
            "known_ee": ">99% ee",
            "applications": "Chiral building block; Lactobacillus kefir ADH benchmark"
        },
        {
            "name": "(R)-1-(4-fluorophenyl)ethanol",
            "SMILES": "C[C@@H](O)c1ccc(F)cc1",
            "R_S_config": "R",
            "KEGG_ID": None,
            "enzyme_class": "ketoreductase (KRED)",
            "known_ee": "97% ee",
            "applications": "Paroxetine intermediate; KRED directed evolution benchmark"
        },
    ],

    # ── Segment 2: Frances Arnold Lab — new-to-nature enzyme chemistry ────────
    # (S)-styrene oxide: tests P450/monooxygenase retrieval (EC 5.3.3.8 isomerase via KEGG,
    #   but biosynthetically it's styrene monooxygenase EC 1.14.14.11).
    # Non-natural cyclopropane: no KEGG — should flag as outside natural pathway space.
    "arnold": [
        {
            "name": "(S)-styrene oxide",
            "SMILES": "[C@H]1(c2ccccc2)CO1",
            "R_S_config": "S",
            "KEGG_ID": "C02944",
            "enzyme_class": "styrene monooxygenase",
            "known_ee": ">99% ee",
            "applications": "Asymmetric epoxidation; CYP450 / styrene monooxygenase"
        },
        {
            "name": "(1R,2S)-cyclopropane-1-carboxylic acid ethyl ester from styrene",
            "SMILES": "CCOC(=O)[C@@H]1C[C@@H]1c1ccccc1",
            "R_S_config": "1R,2S",
            "KEGG_ID": None,
            "enzyme_class": "engineered cytochrome P450 (carbene transferase)",
            "known_ee": "97% ee",
            "applications": "Non-natural cyclopropanation; Arnold lab P450 BM3 directed evolution"
        },
        {
            "name": "(R)-mandelic acid",
            "SMILES": "O[C@@H](C(=O)O)c1ccccc1",
            "R_S_config": "R",
            "KEGG_ID": "C01983",
            "enzyme_class": "mandelate dehydrogenase (EC 1.1.1.379)",
            "known_ee": ">99% ee",
            "applications": "Pharmaceutical building block; hydroxymandelate synthase route"
        },
        {
            "name": "(S)-4-hydroxyphenylglycine",
            "SMILES": "N[C@@H](C(=O)O)c1ccc(O)cc1",
            "R_S_config": "S",
            "KEGG_ID": "C03618",
            "enzyme_class": "4-hydroxyphenylglycine aminotransferase",
            "known_ee": ">99% ee",
            "applications": "Vancomycin antibiotic biosynthetic intermediate"
        },
        {
            "name": "(S)-2-aminobutyric acid",
            "SMILES": "CC[C@@H](N)C(=O)O",
            "R_S_config": "S",
            "KEGG_ID": "C02261",
            "enzyme_class": "threonine deaminase / transaminase",
            "known_ee": ">99% ee",
            "applications": "Non-proteinogenic amino acid; levetiracetam precursor"
        },
    ],

    # ── Segment 3: Pharma process chemistry — Merck / BMS / Pfizer ───────────
    # Route-scouting molecules: established biocatalytic precedent, BRENDA coverage.
    # (S)-lactic acid: iJO1366 native metabolite — feasibility positive control.
    "pharma": [
        {
            "name": "(S)-lactic acid",
            "SMILES": "C[C@H](O)C(=O)O",
            "R_S_config": "S",
            "KEGG_ID": "C00186",
            "enzyme_class": "L-lactate dehydrogenase",
            "known_ee": ">99.9% ee",
            "applications": "PLA bioplastic monomer; Lactobacillus fermentation; iJO1366 native"
        },
        {
            "name": "(S)-2-chloropropanoic acid",
            "SMILES": "Cl[C@@H](C)C(=O)O",
            "R_S_config": "S",
            "KEGG_ID": "C02231",
            "enzyme_class": "haloalkanoic acid dehalogenase",
            "known_ee": ">98% ee",
            "applications": "Herbicide (S)-dichlorprop intermediate; dehalogenase kinetic resolution"
        },
        {
            "name": "(S)-1-phenylalanine",
            "SMILES": "N[C@@H](Cc1ccccc1)C(=O)O",
            "R_S_config": "S",
            "KEGG_ID": "C00079",
            "enzyme_class": "phenylalanine transaminase",
            "known_ee": ">99% ee",
            "applications": "Aspartame precursor; well-mapped pharma API pathway"
        },
        {
            "name": "(R)-1,3-butanediol",
            "SMILES": "C[C@@H](O)CCO",
            "R_S_config": "R",
            "KEGG_ID": "C04483",
            "enzyme_class": "3-hydroxybutyrate dehydrogenase / KRED",
            "known_ee": ">98% ee",
            "applications": "Bioplastic precursor; Genomatica target; cofactor-dependent reduction"
        },
        {
            "name": "(R)-3-amino-1-butanol",
            "SMILES": "C[C@@H](N)CCO",
            "R_S_config": "R",
            "KEGG_ID": None,
            "enzyme_class": "omega-transaminase (ATA)",
            "known_ee": "95% ee",
            "applications": "Prochiral ketoamine; sitagliptin-analog route-scouting"
        },
    ],

    # ── Segment 4: Ginkgo/Zymergen — full-stack strain engineering ────────────
    # Where the COBRApy layer becomes the differentiating feature.
    # (S)-reticuline: KEGG C01516, BIA landmark, well-expressed in S. cerevisiae.
    # (S)-norcoclaurine: KEGG C09136, no KEGG enzymes — tests pathway-gap handling.
    "ginkgo": [
        {
            "name": "(S)-reticuline",
            "SMILES": "COc1ccc([C@@H]2c3cc(OC)c(O)cc3CCN2C)cc1O",
            "R_S_config": "S",
            "KEGG_ID": "C01516",
            "enzyme_class": "berberine bridge enzyme / CYP80B1",
            "known_ee": ">99% ee",
            "applications": "Central BIA biosynthesis intermediate; opioid pathway landmark. "
                           "Note: IUPAC CIP assigns R at C1 of the THIQ ring due to phantom-atom "
                           "priority inflation; 'S' is the historical alkaloid nomenclature convention."
        },
        {
            "name": "(S)-norcoclaurine",
            "SMILES": "Oc1ccc(C[C@@H]2NCCc3cc(O)c(O)cc32)cc1",
            "R_S_config": "S",
            "KEGG_ID": "C09136",
            "enzyme_class": "norcoclaurine synthase (NCS, EC 4.2.1.78)",
            "known_ee": ">99% ee",
            "applications": "First committed step BIA biosynthesis; Pictet-Speraserase product"
        },
        {
            "name": "(R)-mandelic acid",
            "SMILES": "O[C@@H](C(=O)O)c1ccccc1",
            "R_S_config": "R",
            "KEGG_ID": "C01983",
            "enzyme_class": "mandelate dehydrogenase (EC 1.1.1.379)",
            "known_ee": ">99% ee",
            "applications": "Strain engineering target; E. coli heterologous production via HmaS"
        },
        {
            "name": "(S)-lactic acid",
            "SMILES": "C[C@H](O)C(=O)O",
            "R_S_config": "S",
            "KEGG_ID": "C00186",
            "enzyme_class": "L-lactate dehydrogenase",
            "known_ee": ">99.9% ee",
            "applications": "E. coli native metabolite; iJO1366 feasibility positive control"
        },
        {
            "name": "(R)-1,3-butanediol",
            "SMILES": "C[C@@H](O)CCO",
            "R_S_config": "R",
            "KEGG_ID": "C04483",
            "enzyme_class": "3-hydroxybutyrate dehydrogenase",
            "known_ee": ">98% ee",
            "applications": "Bioplastic precursor; Genomatica commercial strain target"
        },
    ],

    # ── Segment 5: Academic synbio — Drew Endy / Pam Silver labs ─────────────
    # (S)-lactic acid is the "hello world" — must produce clean stereo + feasibility.
    # (S)-norcoclaurine: no KEGG enzymes, but good for teaching pathway gap discussion.
    "academic": [
        {
            "name": "(S)-lactic acid",
            "SMILES": "C[C@H](O)C(=O)O",
            "R_S_config": "S",
            "KEGG_ID": "C00186",
            "enzyme_class": "L-lactate dehydrogenase",
            "known_ee": ">99.9% ee",
            "applications": "Hello-world biomanufacturing; Lactobacillus fermentation"
        },
        {
            "name": "(R)-1,3-butanediol",
            "SMILES": "C[C@@H](O)CCO",
            "R_S_config": "R",
            "KEGG_ID": "C04483",
            "enzyme_class": "3-hydroxybutyrate dehydrogenase / KRED",
            "known_ee": ">98% ee",
            "applications": "Bioplastic precursor; cofactor balance teaching case (NAD+/NADH)"
        },
        {
            "name": "(S)-1-phenylalanine",
            "SMILES": "N[C@@H](Cc1ccccc1)C(=O)O",
            "R_S_config": "S",
            "KEGG_ID": "C00079",
            "enzyme_class": "phenylalanine transaminase",
            "known_ee": ">99% ee",
            "applications": "Aromatic amino acid biosynthesis; shikimate pathway teaching"
        },
        {
            "name": "(S)-norcoclaurine",
            "SMILES": "Oc1ccc(C[C@@H]2NCCc3cc(O)c(O)cc32)cc1",
            "R_S_config": "S",
            "KEGG_ID": "C09136",
            "enzyme_class": "norcoclaurine synthase (NCS)",
            "known_ee": ">99% ee",
            "applications": "BIA pathway entry; Pictet-Speraserase mechanistic teaching case"
        },
        {
            "name": "(R)-3-hydroxybutyrate",
            "SMILES": "C[C@@H](O)CC(=O)[O-]",
            "R_S_config": "R",
            "KEGG_ID": "C01089",
            "enzyme_class": "3-hydroxybutyrate dehydrogenase",
            "known_ee": ">98% ee",
            "applications": "PHA polymer precursor; green chemistry teaching molecule"
        },
    ],

    # ── Adversarial: exotic/non-natural — tests graceful fallback ─────────────
    "adversarial": [
        {
            "name": "(1R,2S)-2-fluorocyclopropane-1-carboxylic acid",
            "SMILES": "F[C@@H]1C[C@H]1C(=O)O",
            "R_S_config": "1R,2S",
            "KEGG_ID": None,
            "enzyme_class": "cyclopropane synthase",
            "known_ee": "unknown",
            "applications": "Novel chiral scaffold for agrochemical design"
        },
        {
            "name": "(1R,2S)-2-fluorocyclopropylamine",
            "SMILES": "N[C@H]1C[C@@H]1F",
            "R_S_config": "1R,2S",
            "KEGG_ID": None,
            "enzyme_class": "unknown",
            "known_ee": "unknown",
            "applications": "Mechanism-based MAO inhibitor scaffold"
        },
        {
            "name": "(S)-1-trifluoromethylethanol",
            "SMILES": "FC(F)(F)[C@@H](O)C",
            "R_S_config": "S",
            "KEGG_ID": None,
            "enzyme_class": "ketoreductase (KRED)",
            "known_ee": "85% ee",
            "applications": "Fluorinated building block; expensive asymmetric synthesis"
        },
        {
            "name": "(R)-hexafluoroleucine",
            "SMILES": "N[C@H](CC(C(F)(F)F)C(F)(F)F)C(=O)O",
            "R_S_config": "R",
            "KEGG_ID": None,
            "enzyme_class": "unknown",
            "known_ee": "unknown",
            "applications": "Fluorinated amino acid for protein engineering"
        },
        {
            "name": "(1R,2R)-2-methylcyclopropane-1-carbaldehyde",
            "SMILES": "O=C[C@H]1C[C@@H]1C",
            "R_S_config": "1R,2R",
            "KEGG_ID": None,
            "enzyme_class": "cyclopropane carboxaldehyde synthase",
            "known_ee": "unknown",
            "applications": "Exotic chiral aldehyde; no biosynthetic route established"
        },
    ],
}

# ── Batch quality checks ───────────────────────────────────────────────────────

def check_score_discrimination(suggestions):
    scores = [(s.get("scoring") or {}).get("composite_score", 0) for s in suggestions]
    if not scores:
        return False, "no suggestions"
    spread = max(scores) - min(scores)
    ok = spread >= 0.05
    return ok, f"score range={min(scores):.4f}–{max(scores):.4f}, spread={spread:.4f}"

def check_any_feasible(suggestions):
    statuses = [(s.get("feasibility") or {}).get("status") for s in suggestions]
    # "not_in_model" is a valid scientific answer. Only fail if every status is "error".
    ok = not all(st == "error" for st in statuses)
    return ok, f"feasibility statuses={statuses}"

def check_confidence_not_all_low(suggestions):
    confs = [(s.get("scoring") or {}).get("confidence") for s in suggestions]
    ok = not all(c == "low" for c in confs)
    return ok, f"confidence tiers={confs}"


BATCH_CHECKS = [
    ("score_discrimination",    check_score_discrimination),
    ("feasibility_resolves",    check_any_feasible),
    ("confidence_not_all_low",  check_confidence_not_all_low),
]

# ── Runner ────────────────────────────────────────────────────────────────────

# Ordered list of (query_label, fixture_key) pairs for --all
SEGMENTS = [
    ("Codexis — directed evolution / biocatalysis CRO",                "codexis"),
    ("Frances Arnold Lab — new-to-nature enzyme reconnaissance",        "arnold"),
    ("Merck/BMS/Pfizer — pharma route-scouting",                       "pharma"),
    ("Ginkgo/Zymergen — full-stack strain engineering",                 "ginkgo"),
    ("Academic synbio — teaching-grade hello-world molecules",          "academic"),
    ("Adversarial — exotic/non-natural, graceful fallback stress test", "adversarial"),
    ("Route prediction — Tier 1 acceptance check (all 6 fixture segments)", "route_prediction"),
]

# Map query label → fixture key (used by evaluate() to detect adversarial)
QUERY_TO_FIXTURE = {label: key for label, key in SEGMENTS}


def run_query(query: str, out_dir: str, mock: bool = False,
              fixture_key: str = None) -> dict:
    """Run full pipeline for one query. Returns results dict.

    If mock=True, bypasses the GPT call and injects fixture suggestions.
    fixture_key overrides automatic fixture selection when provided.
    """
    import json as _json

    from ChiraLLM.database_validator import query_kegg
    from ChiraLLM.chirality_checker import validate_chirality
    from ChiraLLM.brenda_client import query_enantioselectivity_batch
    from ChiraLLM.feasibility_checker import check_feasibility
    from ChiraLLM.enantioselectivity_scorer import score_suggestion
    from ChiraLLM.route_predictor import predict_route
    from dataclasses import asdict
    from utils.file_saver import save_suggestions_to_csv

    print(f"\n{'='*70}")
    print(f"QUERY: {query}{' [MOCK GPT]' if mock else ''}")
    print('='*70)

    if mock:
        import copy
        key = fixture_key or QUERY_TO_FIXTURE.get(query, "codexis")
        suggestions = copy.deepcopy(MOCK_SUGGESTIONS[key])
        print(f"  Using {len(suggestions)} fixture suggestions (fixture={key!r})")
    else:
        from ChiraLLM.query_handler import ask_gpt_chirality
        raw = ask_gpt_chirality(query)
        try:
            parsed = _json.loads(raw)
            suggestions = parsed.get("suggestions", [parsed]) if isinstance(parsed, dict) else parsed
        except Exception as e:
            print(f"  [ERROR] JSON parse failed: {e}")
            print(f"  GPT raw response: {raw[:300]}")
            return {"query": query, "suggestions": [], "error": str(e)}

    print(f"  Pipeline processing {len(suggestions)} suggestions")

    for suggestion in suggestions:
        smiles = suggestion.get("SMILES")
        if smiles:
            suggestion["chirality_validation"] = validate_chirality(smiles)

        compound_id = suggestion.get("KEGG_ID")
        if compound_id:
            print(f"  KEGG lookup: {compound_id} ({suggestion.get('name', '?')})")
            kegg = query_kegg(compound_id)
            suggestion["kegg_data"] = kegg
            ec_numbers = kegg.get("enzymes", []) if kegg.get("status") == "success" else []
            if ec_numbers:
                print(f"    BRENDA query for ECs: {ec_numbers[:5]}")
                suggestion["brenda_data"] = query_enantioselectivity_batch(ec_numbers[:5])
            else:
                suggestion["brenda_data"] = {"status": "no_ec_numbers"}
            suggestion["feasibility"] = check_feasibility(compound_id)
            route_result = predict_route(compound_id, mode="top_n", n=3)
            suggestion["route_prediction"] = asdict(route_result)
        else:
            route_result = predict_route(None, mode="top_n", n=3)
            suggestion["route_prediction"] = asdict(route_result)

        suggestion["scoring"] = score_suggestion(suggestion)

    csv_path, json_path = save_suggestions_to_csv(suggestions, out_dir=out_dir)
    print(f"  Saved: {csv_path}")
    print(f"         {json_path}")

    return {"query": query, "suggestions": suggestions, "csv": csv_path, "json": json_path,
            "fixture_key": fixture_key or QUERY_TO_FIXTURE.get(query)}


def evaluate(result: dict) -> bool:
    suggestions = result.get("suggestions", [])
    query = result.get("query", "")
    fkey = result.get("fixture_key") or QUERY_TO_FIXTURE.get(query)
    is_adversarial = fkey == "adversarial"
    all_passed = True

    print(f"\n── Quality evaluation: {query[:65]} ──")
    if is_adversarial:
        print("   [adversarial fixture: relaxed ec_number + confidence checks]")

    # Per-suggestion checks
    for i, s in enumerate(suggestions):
        name = s.get("name", f"suggestion_{i+1}")
        print(f"\n  [{i+1}] {name}")
        for check_name, check_fn in PER_SUGGESTION_CHECKS:
            if is_adversarial and check_name == "ec_number_present":
                ec = (s.get("scoring") or {}).get("top_enzyme", {}).get("ec_number")
                if ec == "unknown":
                    print(f"      ~ ec_number_present: 'unknown' [expected for adversarial]")
                    continue
            passed, detail = check_fn(s)
            icon = "✓" if passed else "✗"
            print(f"      {icon} {check_name}: {detail}")
            if not passed:
                all_passed = False

    # Batch checks
    print(f"\n  [batch]")
    for check_name, check_fn in BATCH_CHECKS:
        if is_adversarial and check_name == "confidence_not_all_low":
            _, detail = check_fn(suggestions)
            print(f"      ~ confidence_not_all_low: {detail} [all-low expected for adversarial]")
            continue
        passed, detail = check_fn(suggestions)
        icon = "✓" if passed else "✗"
        print(f"    {icon} {check_name}: {detail}")
        if not passed:
            all_passed = False

    # Best candidate spotlight
    if suggestions:
        best = max(suggestions, key=lambda s: (s.get("scoring") or {}).get("composite_score", 0))
        sc = best.get("scoring", {})
        te = sc.get("top_enzyme", {})
        print(f"\n  Best candidate: {best.get('name')}")
        print(f"    composite_score:  {sc.get('composite_score')}")
        print(f"    confidence:       {sc.get('confidence')}")
        print(f"    ee_source:        {te.get('ee_source')}")
        print(f"    ee_value:         {te.get('ee_value')}")
        print(f"    tanimoto:         {te.get('tanimoto_similarity')}")
        print(f"    stereo_confirmed: {sc.get('stereo_confirmed')}")
        print(f"    feasibility_flux: {sc.get('feasibility_flux')}")
        if sc.get("scoring_notes"):
            print(f"    notes:")
            for note in sc["scoring_notes"]:
                print(f"      - {note}")

    return all_passed


def evaluate_route_prediction(suggestions: list) -> dict:
    """Acceptance check for the route_prediction Sprint 1 deliverable.

    Per spec §10 #7: KEGG-covered segments must succeed; non-covered segments must
    return appropriate error status without raising.
    """
    KEGG_COVERED_SEGMENTS = {"codexis", "pharma", "academic"}
    NON_COVERED_SEGMENTS = {"arnold", "ginkgo", "adversarial"}
    EXPECTED_NON_COVERED_STATUSES = {
        "no_kegg_id", "invalid_kegg_id", "target_not_in_kegg",
        "target_has_no_reactions", "no_route_found",
    }

    results = {"passes": [], "fails": []}
    for s in suggestions:
        seg = s.get("_fixture_key", "unknown")
        rp = s.get("route_prediction") or {}
        status = rp.get("status")
        n_routes = len(rp.get("routes", []))

        if seg in KEGG_COVERED_SEGMENTS:
            if status == "success" and n_routes >= 1:
                results["passes"].append(f"{seg}/{s.get('name')}: success, {n_routes} routes")
            else:
                results["fails"].append(
                    f"{seg}/{s.get('name')}: expected success, got status={status!r} routes={n_routes}"
                )
        elif seg in NON_COVERED_SEGMENTS:
            if status in EXPECTED_NON_COVERED_STATUSES:
                results["passes"].append(f"{seg}/{s.get('name')}: clean error status={status!r}")
            else:
                results["fails"].append(
                    f"{seg}/{s.get('name')}: unexpected status={status!r}"
                )
        else:
            results["fails"].append(f"{seg}/{s.get('name')}: unknown segment {seg!r}")
    return results


def main():
    parser = argparse.ArgumentParser(description="ChiralAI smoke test — user-segment fixtures")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--query", "-q", type=str,
                      help="Single natural-language query (live GPT call unless --mock)")
    mode.add_argument("--all", action="store_true",
                      help="Run all six fixture segments")
    mode.add_argument("--fixture", type=str,
                      choices=list(MOCK_SUGGESTIONS.keys()) + ["route_prediction"],
                      help="Run a single named fixture segment in mock mode")
    parser.add_argument("--out-dir", type=str, default="smoke_test_output")
    parser.add_argument("--mock", action="store_true",
                        help="Skip GPT; inject fixture suggestions (no API key needed)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if args.fixture == "route_prediction":
        # Meta-fixture: run all 6 base segments, tag suggestions with their segment,
        # then evaluate route_prediction outputs against bipartite acceptance criteria.
        args.mock = True
        all_suggestions = []
        for label, key in SEGMENTS:
            if key == "route_prediction":
                continue  # don't recurse
            print(f"\n=== Running {label} ===")
            result = run_query(label, args.out_dir, mock=args.mock, fixture_key=key)
            for s in result.get("suggestions", []):
                s["_fixture_key"] = key
            all_suggestions.extend(result.get("suggestions", []))
        verdict = evaluate_route_prediction(all_suggestions)
        print("\n=== Route Prediction Acceptance ===")
        for p in verdict["passes"]:
            print(f"  PASS  {p}")
        for f in verdict["fails"]:
            print(f"  FAIL  {f}")
        if verdict["fails"]:
            print(f"\nOverall: ROUTE PREDICTION ACCEPTANCE FAILED "
                  f"({len(verdict['fails'])} failures, {len(verdict['passes'])} passes)")
            sys.exit(1)
        print(f"\nOverall: ROUTE PREDICTION ACCEPTANCE PASSED "
              f"({len(verdict['passes'])} passes)")
        sys.exit(0)
    elif args.fixture:
        # Single named fixture
        label = next((l for l, k in SEGMENTS if k == args.fixture), args.fixture)
        run_pairs = [(label, args.fixture)]
        args.mock = True
    elif args.all:
        run_pairs = SEGMENTS
        args.mock = True
    elif args.query:
        run_pairs = [(args.query, None)]
    else:
        # Default: Codexis fixture
        label, key = SEGMENTS[0]
        run_pairs = [(label, key)]
        args.mock = True

    overall_pass = True
    results_summary = []

    for query_label, fkey in run_pairs:
        result = run_query(query_label, args.out_dir, mock=args.mock, fixture_key=fkey)
        passed = evaluate(result)
        results_summary.append({"query": query_label, "fixture": fkey, "passed": passed})
        if not passed:
            overall_pass = False

    print(f"\n{'='*70}")
    print("SMOKE TEST SUMMARY")
    print('='*70)
    for r in results_summary:
        icon = "PASS" if r["passed"] else "FAIL"
        fixture_tag = f" [{r['fixture']}]" if r["fixture"] else ""
        print(f"  [{icon}]{fixture_tag} {r['query'][:60]}")
    print(f"\nOverall: {'ALL PASS' if overall_pass else 'SOME FAILURES'}")
    print('='*70)

    sys.exit(0 if overall_pass else 1)


if __name__ == "__main__":
    main()
