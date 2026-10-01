"""Guards for the per-gene HIV deployability measurement + the shipped-number correction (2026-09-30)."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "hiv_deployability_per_gene", REPO / "scripts" / "hiv_deployability_per_gene.py")
G = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(G)

ART = REPO / "wiki" / "hiv_deployability_per_gene_2026-09-30.json"
MODELS = {c: REPO / "data" / "hiv_ref" / f"hiv_{c.lower()}_supervised_complement.json"
          for c in ("NNRTI", "PI", "INSTI")}


def _art():
    if not ART.exists():
        pytest.skip("artifact absent (gitignored Stanford data missing at build time)")
    return json.loads(ART.read_text(encoding="utf-8"))


def _model(cls):
    return json.loads(MODELS[cls].read_text(encoding="utf-8"))


def test_the_split_and_metric_are_IMPORTED_from_the_nnrti_harness():
    """Re-writing the split or the pass rule would let this run disagree with the committed RT artifact
    for reasons unrelated to the gene being measured."""
    G.self_check()
    src = (REPO / "scripts" / "hiv_deployability_per_gene.py").read_text(encoding="utf-8")
    assert "D._fit_predict_grouped" in src and "D._blindspot_metrics" in src
    assert "def _blindspot_metrics" not in src, "the metric must be imported, never re-declared here"
    # Prose may DESCRIBE the split; what must not appear is a second implementation of it.
    assert "GroupKFold(" not in src and "import GroupKFold" not in src, \
        "the split must come from the harness, not be re-instantiated here"
    assert "cross_val_predict" not in src, "the CV call belongs to the harness"


def test_the_rt_anchor_reproduces_the_committed_artifact():
    """Six shipped numbers rest on this harness; if it cannot reproduce a known result it is not a
    measurement, and the script exits 3 rather than reporting the other genes."""
    d = _art()
    a = d["anchor"]
    assert a["reproduces"] is True, a
    assert a["committed_n"] == a["measured_n"] and a["committed_R"] == a["measured_R"]
    assert abs(a["measured_auroc"] - a["committed_auroc"]) <= G.ANCHOR_TOL
    assert (REPO / a["artifact"]).exists()


def test_the_deployed_comparator_agrees_with_the_original_majors_set():
    """The original harness scored the NNRTI majors set directly while this uses the deployed
    `call_hiv_observed`; the anchor is only interpretable because they agree on every isolate."""
    d = _art()
    agree = d["cells"]["NNRTI"]["deployed_vs_majors_agreement"]
    n, tot = (int(x) for x in agree.split("/"))
    assert n == tot and tot > 1000, agree


def test_pi_is_REFUSED_not_scored_and_carries_no_number():
    """THE DEFECT THIS RUN FOUND. PI's own position-based catalog over-calls, which empties its own blind
    spot of resistant isolates -- so a 0.89 was not merely unsourced, it is not measurable here."""
    d = _art()
    pi = d["cells"]["PI"]
    assert pi["status"] == "BLIND_SPOT_UNSCOREABLE"
    assert "leave_study_out_blindspot_auroc_MEASURED" not in pi, \
        "a refused cell must WITHHOLD the number, not report an unpowered one"
    assert pi["blind_spot_R"] < G.MIN_PER_CLASS
    assert pi["catalog_called_positive_fraction"] > 0.5, "the refusal is caused by an over-calling catalog"
    assert "refusal" in pi and "over-calling" in pi["refusal"]


def test_insti_is_a_SECOND_GENE_and_passes_the_same_bar():
    """The 2026-09-30 GLM result generalized across DRUGS but every drug shared ONE gene (RT). Integrase
    is a different gene with its own catalog and its own dataset."""
    d = _art()
    insti = d["cells"]["INSTI"]
    assert insti["status"] == "SCORED"
    assert insti["gene"] == "IN" and d["cells"]["NNRTI"]["gene"] == "RT"
    bs = insti["blind_spot"]
    assert bs["pass"] is True
    assert bs["auroc"] >= 0.65 and bs["auroc"] > bs["burden"] and bs["null"] < 0.55
    assert insti["n_studies"] >= 5, "leave-one-STUDY-out needs several studies to mean anything"


def test_every_shipped_model_carries_measured_provenance_or_an_explicit_refusal():
    """The failure mode being fixed: a number in a shipped, git-tracked file with no artifact behind it."""
    for cls, path in MODELS.items():
        dep = _model(cls)["deployability"]
        assert "measured" in dep, f"{cls}: deployability must state whether it was measured"
        assert dep.get("artifact"), f"{cls}: deployability must cite its artifact"
        assert (REPO / dep["artifact"]).exists(), f"{cls}: cited artifact must resolve"
        if dep["measured"]:
            assert isinstance(dep["leave_study_out_blindspot_auroc"], float)
            assert dep.get("passes") is True and dep.get("pass_rule")
        else:
            assert dep["leave_study_out_blindspot_auroc"] is None, \
                f"{cls}: an unmeasured cell must carry None, never a number"
            assert dep.get("blind_spot_status") == "UNSCOREABLE" and dep.get("reason")


def test_the_shipped_numbers_equal_the_measured_ones():
    d = _art()
    for cls in MODELS:
        want = d["cells"][cls].get("leave_study_out_blindspot_auroc_MEASURED")
        assert _model(cls)["deployability"]["leave_study_out_blindspot_auroc"] == want, cls


def test_the_builder_reproduces_the_committed_blocks_so_a_rebuild_cannot_reintroduce_a_literal():
    """THE ROOT-CAUSE GUARD. The values were literals in the builder's CLASSES dict; if the builder still
    carried them, any rebuild would silently restore the unsourced figures."""
    spec = importlib.util.spec_from_file_location(
        "build_hiv_complement_model", REPO / "scripts" / "build_hiv_complement_model.py")
    B = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(B)
    for cls, cfg in B.CLASSES.items():
        assert "leave_study_out" not in cfg, \
            f"{cls}: the builder must DERIVE deployability, not carry a literal"
        assert B.deployability_block(cls) == _model(cls)["deployability"], \
            f"{cls}: a rebuild would change the shipped deployability block"


def test_the_builder_REFUSES_when_the_anchor_did_not_bind(tmp_path, monkeypatch):
    """Non-vacuity: stamping deployability from a harness that cannot reproduce a known result is exactly
    how an unsourced number ships, so the builder must raise rather than fall back."""
    spec = importlib.util.spec_from_file_location(
        "build_hiv_complement_model", REPO / "scripts" / "build_hiv_complement_model.py")
    B = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(B)
    bad = json.loads(ART.read_text(encoding="utf-8"))
    bad["anchor"]["reproduces"] = False
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(bad), encoding="utf-8")
    monkeypatch.setattr(B, "DEPLOYABILITY_ARTIFACT", p)
    with pytest.raises(ValueError, match="anchor did not bind"):
        B.deployability_block("INSTI")
    missing = tmp_path / "nope.json"
    monkeypatch.setattr(B, "DEPLOYABILITY_ARTIFACT", missing)
    with pytest.raises(FileNotFoundError):
        B.deployability_block("INSTI")


# The commit immediately BEFORE this correction. Pinned deliberately: the first version of the guard below
# compared against `HEAD`, which passed while the work was uncommitted and then FAILED the moment it was
# committed -- HEAD had become the corrected state, so "only deployability changed" became "nothing
# changed". A guard that asserts a diff needs a FIXED baseline, not a moving one.
PRE_CORRECTION_COMMIT = "b67ee2f"


def test_the_scorer_weights_are_untouched_by_the_correction():
    """This changed a DISCLOSURE, not a model: every complement must score exactly as before."""
    import subprocess
    checked = 0
    for cls, path in MODELS.items():
        rel = str(path.relative_to(REPO)).replace("\\", "/")
        old = subprocess.run(["git", "show", f"{PRE_CORRECTION_COMMIT}:{rel}"],
                             capture_output=True, text=True, cwd=REPO)
        if old.returncode != 0:
            pytest.skip(f"{PRE_CORRECTION_COMMIT} unreachable (shallow clone?)")
        o, n = json.loads(old.stdout), _model(cls)
        assert o["weights"] == n["weights"], f"{cls}: weights changed"
        assert o["intercept"] == n["intercept"], f"{cls}: intercept changed"
        changed = {k for k in set(o) | set(n) if o.get(k) != n.get(k)}
        assert changed == {"deployability"}, f"{cls}: unexpected fields changed: {changed}"
        # Non-vacuity: the baseline must actually be the PRE-correction state, or this proves nothing.
        assert o["deployability"] == {"leave_study_out_blindspot_auroc": pytest.approx(
            G.SUPERSEDED_LITERALS[cls])}, f"{cls}: {PRE_CORRECTION_COMMIT} is not the pre-correction state"
        checked += 1
    assert checked == 3, "all three complements must be compared"


def test_no_memo_or_module_restates_an_unsourced_089_for_pi():
    """The number was restated in a builder docstring, a doubt.py comment and two memos. A restated value
    drifts from its source -- this repo's recurring defect -- so PI must not carry a bare 0.89 anywhere."""
    targets = [REPO / "scripts" / "build_hiv_complement_model.py",
               REPO / "dna_decode" / "eval" / "doubt.py",
               REPO / "wiki" / "supervised_complement_wired_2026-09-30.md",
               REPO / "wiki" / "supervised_complement_per_class_nonvacuity_2026-09-30.md"]
    markers = ("unsourced", "superseded", "until 2026-09-30", "restated", "0.8923", "read \"",
               "no number", "only the first was measured", "not measurable")
    for t in targets:
        if not t.exists():
            continue
        lines = t.read_text(encoding="utf-8").splitlines()
        for i, line in enumerate(lines):
            if "0.89" not in line:
                continue
            # A correction legitimately spans a sentence, so the marker may sit on a neighbouring line;
            # checking a small window avoids failing an honest multi-line correction while still catching
            # a bare restatement.
            window = " ".join(lines[max(0, i - 3):i + 4]).lower()
            assert any(w in window for w in markers), \
                f"{t.name}: 0.89 restated without marking it superseded:\n  {line.strip()}"


def test_the_derived_nonvacuity_artifact_propagated_the_correction():
    """It reads `model_info()`, so it is the check that the fix flows through rather than being local."""
    p = REPO / "wiki" / "supervised_complement_nonvacuity_2026-09-30.json"
    if not p.exists():
        pytest.skip("non-vacuity artifact absent")
    d = json.loads(p.read_text(encoding="utf-8"))
    for c in d["cells"]:
        got = c.get("leave_study_out_blindspot_auroc")
        if c["class"] == "PI":
            assert got is None, "PI must no longer publish a deployability number"
            assert c["verdict"] == "THRESHOLD_NEVER_FIRES"
        else:
            assert got is not None and got > 0.65
