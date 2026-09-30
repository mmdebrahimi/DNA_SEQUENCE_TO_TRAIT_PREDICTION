"""Guards for the all-RT GLM headroom generalization (2026-09-30)."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "glm_headroom_allrt", REPO / "scripts" / "glm_headroom_allrt.py")
A = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(A)

ART = REPO / "wiki" / "glm_headroom_allrt_2026-09-30.json"


def _art():
    if not ART.exists():
        pytest.skip("artifact absent (gitignored data or ESM cache missing at build time)")
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_bar_is_IMPORTED_not_restated():
    """One definition across both GLM scripts. A second copy is a second thing to drift."""
    G = A._bar()
    assert G.PREREGISTERED["min_auroc"] == 0.65
    assert "hiv_esm_vs_catalog_2026-07-09" in G.PREREGISTERED["source"]
    src = (REPO / "scripts" / "glm_headroom_allrt.py").read_text(encoding="utf-8")
    assert "min_auroc\": 0.65" not in src and "PREREGISTERED = {" not in src, \
        "the bar must be imported from glm_alphabet_headroom, never re-declared here"


def test_the_underpowered_guard_matches_the_committed_run():
    """Reusing the committed run's own <10-per-class guard is what makes the two drug SETS agree; a
    different guard would silently compare different drugs."""
    assert A.MIN_PER_CLASS == 10
    from scripts.hiv_esm_vs_catalog_allrt import CUTOFF
    assert CUTOFF == 3.0


def test_nominal_and_genuine_are_BOTH_reported_and_small_subsets_are_named():
    """THE DOCUMENTED TRAP. The committed zero-shot run reports nominal 1 / genuine 0 -- the nominal pass
    is DOR on a 37-isolate subset with burden 0.783. Reporting a bare pass-count launders that."""
    s = A.summarize([
        {"drug": "BIG", "subset_n": 500, "zero_shot_esm2_auroc": 0.45, "supervised_auroc": 0.82,
         "verdict_supervised": "PASS", "verdict_zero_shot": "FAIL_BELOW_MIN_AUROC"},
        {"drug": "TINY", "subset_n": 37, "zero_shot_esm2_auroc": 0.80, "supervised_auroc": 0.86,
         "verdict_supervised": "PASS", "verdict_zero_shot": "PASS"},
    ])
    assert s["zero_shot_pass_nominal"] == 1 and s["zero_shot_pass_genuine"] == 0
    assert s["supervised_pass_nominal"] == 2 and s["supervised_pass_genuine"] == 1
    assert "TINY(n=37)" in s["small_subsets"], "an excluded small subset must be NAMED"


def test_self_check_passes():
    A.self_check()


def test_every_per_drug_anchor_reproduced_the_committed_run():
    """Six anchors, not one. A drug whose subset disagrees with the committed run cannot be compared
    against it, so its supervised number is WITHHELD rather than reported."""
    d = _art()
    n_ok, n_tot = d["anchors_reproducing"].split("/")
    assert n_ok == n_tot and int(n_tot) >= 6, f"anchors: {d['anchors_reproducing']}"
    for c in d["cells"]:
        if c.get("status") == "ANCHOR_MISMATCH":
            assert "supervised_auroc" not in c, "a mismatched anchor must WITHHOLD the supervised number"


def test_the_committed_baseline_is_cited_and_its_median_reproduces():
    """Our own median zero-shot on the anchored (non-Full) variant must equal the committed
    median_esm_subset_auroc -- an independent check that the whole drug set was rebuilt faithfully."""
    d = _art()
    base = d["committed_zero_shot_baseline"]
    assert (REPO / base["artifact"]).exists()
    assert base["n_drugs_passing_bar_genuine"] == 0
    anchor = d["summary_by_variant"]["anchor"]
    assert abs(anchor["median_zero_shot"] - base["median_esm_subset_auroc"]) <= 0.02


def test_the_honest_headline_is_the_GROUPED_split_and_zero_shot_never_passes_genuinely():
    d = _art()
    full = d["summary_by_variant"]["full"]
    assert full["zero_shot_pass_genuine"] == 0, "zero-shot must not acquire a genuine pass"
    assert full["supervised_pass_genuine"] >= 6
    assert full["median_supervised"] > full["median_zero_shot"] + 0.2
    grouped = [c for c in d["cells"]
               if c["variant"] == "full" and c.get("status") == "SCORED"]
    assert grouped and all(c["split"] == "leave_one_study_out" for c in grouped), \
        "the .Full variant must be scored leave-one-STUDY-out, not ungrouped"


def test_the_ungrouped_variant_is_labelled_optimistic():
    """The non-Full files carry no RefID, so that arm cannot be grouped -- it must say so in its own
    split name rather than being quoted as if it were held out by study."""
    d = _art()
    for c in d["cells"]:
        if c["variant"] == "anchor" and c.get("status") == "SCORED":
            assert "OPTIMISTIC" in c["split"]


def test_a_failing_drug_is_reported_not_dropped():
    """The bar must be shown to BITE. D4T fails on the shuffled null (unstable at 11 resistant of ~900);
    a run where every drug passes would be a bar that cannot discriminate."""
    d = _art()
    full = [c for c in d["cells"] if c["variant"] == "full" and c.get("status") == "SCORED"]
    fails = [c for c in full if c["verdict_supervised"] != "PASS"]
    assert fails, "no drug failed -- check the bar still discriminates"
    for c in fails:
        assert "supervised_auroc" in c and c["verdict_supervised"].startswith("FAIL")


def test_honest_scope_names_the_single_gene_limit():
    """Every drug here shares ONE gene (RT), so this generalizes across DRUGS and not across GENES.
    Losing that caveat would overstate the result substantially."""
    d = _art()
    s = d["honest_scope"].lower()
    assert "one gene" in s or "one gene (rt)" in s
    assert "weakest" in s and "floor" in s
    assert "underpowered" in s
