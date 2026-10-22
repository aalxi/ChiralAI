"""Regression tests for every defect the benchmark found.

Each test names the defect it pins down. These are the checks that would have
failed on the shipped code, so if one of them ever fails again the regression is
named rather than merely red.

Tests marked ``network`` hit KEGG. Responses are cached on disk under
~/.cache/chiralai/kegg, so a second run is offline. Skip them with
``-m "not network"``.
"""
import os, sys
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from ChiraLLM import feasibility, retrieval, structure_gate
from stereo import api as stereo_api
from utils import display, kegg

S_LACTIC = "C[C@H](O)C(=O)O"
R_LACTIC = "C[C@@H](O)C(=O)O"
R_LIMONENE = "CC1=CC[C@@H](CC1)C(=C)C"
FLAT_LIMONENE = "CC1=CCC(CC1)C(=C)C"


# ── the CLI used to be unable to start at all ────────────────────────────────

def test_entry_point_imports():
    """main.py:4 imported save_suggestions_to_csv, which never existed."""
    import main
    assert callable(main.main)


def test_llm_stage_imports():
    """query_handler imported ChatCompletion, absent from openai>=1.0."""
    import ChiraLLM.query_handler as qh
    assert callable(qh.ask_gpt_chirality)


def test_no_stray_diff_markers():
    src = open(os.path.join(os.path.dirname(__file__), "..", "ChiraLLM", "feasibility.py")).read()
    assert not src.lstrip().startswith(("+", "-"))


# ── the fabricated confidence score ──────────────────────────────────────────

def test_no_confidence_field_anywhere():
    """It took 2 values over 36 compounds; 0.88 meant both 'too big' and 'no route'."""
    rec = feasibility.assess("L-lactic acid", S_LACTIC, "E. coli", cross_check=False)
    assert "confidence" not in rec
    assert "evidence" in rec and rec["evidence"]


def test_report_does_not_print_a_confidence_percentage(capsys):
    rec = feasibility.assess("L-lactic acid", S_LACTIC, "E. coli", cross_check=False)
    display.display_results_in_terminal([rec])
    assert "confidence" not in capsys.readouterr().out.lower()


# ── the key-case mismatch that discarded every LLM row ───────────────────────

@pytest.mark.parametrize("key", ["smiles", "SMILES", "Smiles", "structure"])
def test_csv_keys_are_read_case_insensitively(key):
    """The prompt asked for 'SMILES'; the consumer read 'smiles'. Silent data loss."""
    assert feasibility._row_get({"Name": "x", key: S_LACTIC}, "smiles", "structure") == S_LACTIC


# ── the structure gate: the LLM's output is not trusted ──────────────────────

def test_gate_rejects_unparseable_structure():
    assert not structure_gate.check("x", "C1CC(", cross_check=False)["ok"]


def test_gate_warns_when_stereochemistry_is_unspecified():
    r = structure_gate.check("limonene", FLAT_LIMONENE, cross_check=False)
    assert r["ok"] and any("none specified" in w for w in r["warnings"])


@pytest.mark.network
@pytest.mark.parametrize("name,smiles,why", [
    ("Limonene", open(os.path.join(os.path.dirname(__file__), "..", "output",
                                   "llm_suggestions_20250630_174747.csv")).read()
                 .split("Limonene,")[1].split(",")[0].strip('"'),
     "896-character degenerate polyterpene"),
    ("L-Lactic Acid", "C(C(=O)O)C(C(=O)O)C(=O)O", "a tricarboxylic acid"),
    ("L-DOPA", "C1=CC(=C(C=C1)O)C(C(=O)O)N", "m-tyrosine: one hydroxyl short of DOPA"),
])
def test_gate_rejects_this_repos_own_recorded_llm_failures(name, smiles, why):
    """Replays the actual bad structures saved in output/. Every one of these was
    analysed as if correct by the shipped pipeline, because every downstream
    check keyed on the name."""
    r = structure_gate.check(name, smiles)
    assert not r["ok"], f"{why} was accepted"
    assert "formula mismatch" in r["rejected_because"]


@pytest.mark.network
def test_gate_flags_an_unverifiable_name_rather_than_blessing_it():
    """'Poly(lactic acid)' -> malonic acid cannot be refuted: KEGG has no compound
    by that name. The gate must say the structure is unverified, not pass silently."""
    r = structure_gate.check("Poly(lactic acid)", "C(C(=O)O)C(=O)O")
    assert any("unverified" in w.lower() for w in r["warnings"])


@pytest.mark.network
def test_gate_rejects_name_structure_mismatch():
    """A tricarboxylic acid labelled 'lactic acid' is in this repo's own output."""
    r = structure_gate.check("lactic acid", "OC(=O)C(O)(C(=O)O)C(=O)O")
    assert not r["ok"] and "formula mismatch" in r["rejected_because"]


@pytest.mark.network
def test_gate_does_not_judge_against_the_first_name_hit():
    """KEGG returns UDP-glucose first for 'glucose'. First-hit comparison rejects
    correct glucose -- a false rejection I introduced and had to fix."""
    r = structure_gate.check("glucose", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
    assert r["ok"], r.get("rejected_because")


# ── retrieval by structure rather than name ──────────────────────────────────

@pytest.mark.network
def test_retrieval_distinguishes_enantiomers():
    """Name search cannot do this; it is why the pipeline was stereo-blind."""
    l = retrieval.resolve("L-lactic acid", S_LACTIC)
    d = retrieval.resolve("D-lactic acid", R_LACTIC)
    assert l["kegg_id"] == "C00186"
    assert d["kegg_id"] == "C00256"
    assert l["kegg_id"] != d["kegg_id"]


# ── pathway classification ───────────────────────────────────────────────────

def test_degradation_is_not_counted_as_biosynthesis():
    """'any pathway' is what passed citalopram via a P450 drug-degradation map."""
    assert kegg.classify_pathway("path:map00930", "Caprolactam degradation") == "degradation"
    assert kegg.classify_pathway("path:map00982", "Drug metabolism - cytochrome P450") == "xenobiotic_degradation"
    assert kegg.classify_pathway("path:map01100", "Metabolic pathways") == "global_overview"
    assert kegg.classify_pathway("path:map00030", "Pentose phosphate pathway") == "other"
    assert kegg.classify_pathway("path:map00400", "Phenylalanine biosynthesis") == "biosynthesis"


def test_degradation_test_precedes_metabolism_test():
    """Ordering bug risk: 'Drug metabolism' contains the word metabolism."""
    assert kegg.classify_pathway("path:map00983", "Drug metabolism - other enzymes") != "metabolism"


# ── stereochemistry: the layer that did not exist ────────────────────────────

def test_wrong_form_is_not_reported_as_feasible():
    """L-glucose inherits D-glucose's KEGG maps by connectivity. Without a
    stereochemistry-aware verdict the headline reads ROUTE DOCUMENTED IN HOST
    for a compound with no documented route to the requested enantiomer."""
    rec = feasibility.assess(
        "L-glucose", "OC[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O",
        "E. coli", cross_check=False)
    assert rec["stereochemistry"]["outcome"] == "WRONG-FORM"
    assert rec["verdict"] == "route-documented-for-the-wrong-stereoisomer"
    assert rec["biosynthetic_feasibility"] is False


def test_all_four_stereo_outcomes_are_reachable():
    assert stereo_api.lookup(S_LACTIC)["outcome"] == "EXACT"
    assert stereo_api.lookup(FLAT_LIMONENE)["outcome"] == "SKELETON"
    assert stereo_api.lookup("C1CC2(CCC1)CCCCN2CCCF")["outcome"] == "NONE"


def test_enantiomers_get_different_answers():
    """9/9 enantiomer pairs previously received byte-identical analysis."""
    a, b = stereo_api.lookup(S_LACTIC), stereo_api.lookup(R_LACTIC)
    assert a["enzymes_for_requested_form"] != b["enzymes_for_requested_form"]


def test_selectivity_risk_is_the_mirror_of_the_requested_form():
    """(S)-lactate's risk ECs must be (R)-lactate's producing ECs, and vice versa."""
    s, r = stereo_api.lookup(S_LACTIC), stereo_api.lookup(R_LACTIC)
    assert set(s["selectivity_risk"]["made_by_ec"]) <= set(r["enzymes_for_requested_form"])
    assert set(r["selectivity_risk"]["made_by_ec"]) <= set(s["enzymes_for_requested_form"])


def test_limonene_case_this_repo_originally_got_wrong():
    res = stereo_api.lookup(R_LIMONENE)
    assert res["outcome"] == "EXACT"
    assert "EC:4.2.3.20" in res["enzymes_for_requested_form"]
    assert "EC:4.2.3.16" in res["selectivity_risk"]["made_by_ec"]


def test_protonation_state_does_not_cause_a_false_wrong_form():
    """Neutral (S)-lactic acid vs ChEBI's (S)-lactate anion differ in InChIKey
    block 3. Matching the full key returned a confident, wrong WRONG-FORM."""
    for smi in ("C[C@H](O)C(=O)O", "C[C@H](O)C(=O)[O-]"):
        assert stereo_api.lookup(smi)["outcome"] == "EXACT"


# ── the host argument now does something ─────────────────────────────────────

@pytest.mark.network
def test_host_changes_the_verdict():
    """Host changed 0/36 outputs before; eco/sce were defined and never used.

    Ergosterol is the clean case and the biology is unambiguous: yeast encodes
    the steroid biosynthesis map, E. coli makes no sterols at all.
    """
    # Canonicalised from KEGG's own MOL record for C01694. A hand-drawn version
    # that omitted the trans double-bond geometry was read as a different
    # stereoisomer and came back WRONG-FORM -- the stereo index is strict about
    # double bonds too, not only tetrahedral centres.
    ergosterol = ("CC(C)[C@@H](C)/C=C/[C@@H](C)[C@H]1CC[C@H]2C3=CC=C4C[C@@H](O)"
                  "CC[C@]4(C)[C@H]3CC[C@]12C")
    eco = feasibility.assess("ergosterol", ergosterol, "E. coli", cross_check=False)
    sce = feasibility.assess("ergosterol", ergosterol, "yeast", cross_check=False)
    assert eco["host_kegg_code"] == "eco" and sce["host_kegg_code"] == "sce"
    assert not eco["host_biosynthetic_pathways"], "E. coli does not make sterols"
    assert sce["host_biosynthetic_pathways"], "yeast does make ergosterol"
    assert eco["verdict"] != sce["verdict"]


@pytest.mark.network
def test_host_pathway_sets_actually_differ():
    """If these sets were identical the host layer would be decorative."""
    eco, sce = kegg.host_pathway_ids("eco"), kegg.host_pathway_ids("sce")
    assert len(eco) > 100 and len(sce) > 100
    assert eco != sce and (eco - sce) and (sce - eco)


@pytest.mark.network
def test_unavailable_layers_are_declared_not_silently_empty():
    rec = feasibility.assess("L-lactic acid", S_LACTIC, "E. coli", cross_check=False)
    layers = {l["layer"] for l in rec["layers_unavailable"]}
    assert "BRENDA" in layers
    assert all(l["reason"] for l in rec["layers_unavailable"])
