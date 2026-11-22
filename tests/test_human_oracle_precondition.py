"""Pin the fold-9 re-measurement work of 2026-11-22.

WHY THESE EXIST. Three separate defects, each of which produced or would have produced a number
that looks like a measurement:

  * THE CEILING WAS MEASURED ON A BIASED SUBSAMPLE. `quintiles()` returns each stratum sorted by
    activity ASCENDING and the real_ceiling arm took `[:per_stratum]` -- the n LOWEST-activity
    members of every stratum, truncating the top quintile's tail where most of the activity signal
    lives. Measured: that understated the ceiling by 0.1095 in K562 (-4.14 sd of the random-draw
    spread) and 0.0300 in HepG2 (-2.74 sd), while random-n and all-records agree to 0.04 sd -- so
    sample SIZE was never the issue, the ORDER was. The ceiling is the DENOMINATOR of
    stratified_frac_of_ceiling, so understating it OVERSTATES every generator's fraction of
    achievable signal: the significance-manufacturing direction.
  * THE DECISION-GRADE SET INCLUDED THE RESERVED SHARD. Fold 10 is FINAL_ONLY; any threshold or
    verdict rule informed by a folds-9+10 number has spent it. The held-out fold set is now a
    parameter and the output path is fold-tagged, so a fold-9 result cannot land on the folds-9+10
    record.
  * NOTHING CHECKED THE ORACLE SCORES A KNOWN INPUT TO A KNOWN VALUE. The protocol named this third
    precondition and it did not exist. It RAISES rather than returning a verdict a loop could skip:
    a floor that never computes is not a failing floor, it is no floor.

No torch, no substrate, no network. The end-to-end values are pinned as ARTIFACT records rather
than recomputed -- a recompute is ~5 minutes per script and belongs in a run, not a test suite.
"""
from __future__ import annotations

import ast
import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ORACLE_SRC = REPO / "scripts" / "glm_human_activity_oracle.py"
NULL_SRC = REPO / "scripts" / "glm_human_stratified_null.py"
ORACLE_F9 = REPO / "wiki" / "glm_human_activity_oracle_fold9_2026-11-22.json"
ORACLE_F910 = REPO / "wiki" / "glm_human_activity_oracle_2026-11-22.json"
NULL_F9 = REPO / "wiki" / "glm_human_stratified_null_fold9_2026-11-22.json"
NULL_F910 = REPO / "wiki" / "glm_human_stratified_null_2026-11-22.json"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def oracle():
    return _load("_oracle_precond_under_test", ORACLE_SRC)


def _code_only(src: str) -> str:
    """Source with comment lines and string literals removed.

    A text search over raw source cannot tell an executable reference from a comment that EXPLAINS
    the thing being forbidden -- and a correction necessarily quotes what it corrects, so a guard
    written beside one is unusually likely to trip on its own prose. Two tests in this repo already
    failed that way. Comments are stripped by line; string literals are blanked via the real parser.
    """
    lines = [ln for ln in src.splitlines() if not ln.lstrip().startswith("#")]
    stripped = "\n".join(lines)
    tree = ast.parse(stripped)
    spans = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if node.end_lineno is not None:
                spans.append((node.lineno, node.end_lineno))
    out = stripped.splitlines()
    for lo, hi in spans:
        for i in range(lo - 1, min(hi, len(out))):
            out[i] = ""
    return "\n".join(out)


# ------------------------------------------------------- the third precondition, behaviourally

def test_the_precondition_passes_on_the_shipped_pipeline(oracle):
    rec = oracle.assert_oracle_fixture()
    assert rec["status"] == "PASSED"
    assert rec["fixture_seq"] == "ACGTACGTACGT"
    # It must report WHAT it checked. "The precondition passed" and "the precondition never ran"
    # are different states and a bare boolean cannot tell them apart.
    assert len(rec["checked"]) >= 6
    assert rec["limit"], "a precondition with no stated limit is overclaiming"


def test_the_expected_values_are_hand_derivable_not_merely_current_output(oracle):
    """The fixture is a CORRECTNESS check only if its values come from the definitions.

    kmer4 index arithmetic, itertools.product("ACGT", repeat=4) order, A0 C1 G2 T3:
        ACGT = 0*64 + 1*16 + 2*4 + 3 = 27   CGTA = 1*64 + 2*16 + 3*4 + 0 = 108
        GTAC = 2*64 + 3*16 + 0*4 + 1 = 177  TACG = 3*64 + 0*16 + 1*4 + 2 = 198
    and "ACGTACGTACGT" has 9 4-mer windows: ACGT x3, CGTA x2, GTAC x2, TACG x2.
    """
    assert oracle.FIXTURE_KMER4 == {27: 3 / 9, 108: 2 / 9, 177: 2 / 9, 198: 2 / 9}
    assert sum(oracle.FIXTURE_KMER4.values()) == pytest.approx(1.0)
    # bins=3 over "AAACCCGGGTTT" -> AAAC | CCGG | GTTT, each base count / 4
    assert oracle.FIXTURE_POSITIONAL_ASYM == [0.75, 0.25, 0.0, 0.0,
                                              0.0, 0.50, 0.50, 0.0,
                                              0.0, 0.0, 0.25, 0.75]


def test_the_asymmetric_fixture_is_non_square_which_is_what_makes_layout_checkable(oracle):
    """The symmetric fixture is BLIND to a base-major/bin-major reshape flip, found by perturbing
    the builder and watching the check fail to fire. 4 bins x 4 bases over an ACGT-periodic
    sequence is symmetric, so its transpose is itself. 3 x 4 cannot coincide with its transpose."""
    assert len(oracle.FIXTURE_POSITIONAL) == 16           # 4 bins x 4 bases -- square, symmetric
    assert len(oracle.FIXTURE_POSITIONAL_ASYM) == 12      # 3 bins x 4 bases -- non-square
    assert len(oracle.FIXTURE_SEQ_ASYM) % 3 == 0


@pytest.mark.parametrize("what", ["reindexed", "unnormalised", "n_as_a", "gc_g_only",
                                  "ragged_bins", "transposed"])
def test_the_precondition_refuses_each_plausible_pipeline_change(oracle, what, monkeypatch):
    """Non-vacuity, one perturbation at a time. Every one is a defect a wrong-version dependency or
    a well-meant refactor could plausibly introduce; a guard that cannot fire is decoration."""
    import numpy as np

    def kmer_reindexed(seq, k):
        import itertools
        idx = {"".join(p): i for i, p in enumerate(itertools.product("TGCA", repeat=k))}
        v = np.zeros(len(idx), dtype=np.float32)
        n = 0
        for i in range(len(seq) - k + 1):
            j = idx.get(seq[i:i + k])
            if j is not None:
                v[j] += 1.0
                n += 1
        return (v / n).tolist() if n else v.tolist()

    def kmer_unnormalised(seq, k):
        idx = oracle._kmer_index(k)
        v = np.zeros(len(idx), dtype=np.float32)
        for i in range(len(seq) - k + 1):
            j = idx.get(seq[i:i + k])
            if j is not None:
                v[j] += 1.0
        return v.tolist()

    def kmer_n_as_a(seq, k):
        idx = oracle._kmer_index(k)
        v = np.zeros(len(idx), dtype=np.float32)
        n = 0
        for i in range(len(seq) - k + 1):
            j = idx.get(seq[i:i + k].replace("N", "A"))
            if j is not None:
                v[j] += 1.0
                n += 1
        return (v / n).tolist() if n else v.tolist()

    def positional_ragged(seq, bins):
        w = max(1, len(seq) // bins)
        out = np.zeros((bins, 4), dtype=np.float32)
        for b in range(bins):
            chunk = seq[b * w:(b + 1) * w] or "A"
            for ci, base in enumerate("ACGT"):
                out[b, ci] = chunk.count(base) / len(chunk)
        return out.reshape(-1).tolist()

    def positional_transposed(seq, bins):
        if len(seq) % bins:
            raise ValueError(f"{len(seq)} positions do not divide into {bins} bins")
        w = len(seq) // bins
        out = np.zeros((bins, 4), dtype=np.float32)
        for b in range(bins):
            chunk = seq[b * w:(b + 1) * w]
            for ci, base in enumerate("ACGT"):
                out[b, ci] = chunk.count(base) / w
        return out.T.reshape(-1).tolist()

    target, repl = {
        "reindexed": ("kmer_feats", kmer_reindexed),
        "unnormalised": ("kmer_feats", kmer_unnormalised),
        "n_as_a": ("kmer_feats", kmer_n_as_a),
        "gc_g_only": ("gc", lambda s: [s.count("G") / max(len(s), 1)]),
        "ragged_bins": ("positional_onehot", positional_ragged),
        "transposed": ("positional_onehot", positional_transposed),
    }[what]

    monkeypatch.setattr(oracle, target, repl)
    with pytest.raises(oracle.OraclePreconditionFailed):
        oracle.assert_oracle_fixture()


def test_the_precondition_raises_rather_than_returning_an_ignorable_verdict(oracle):
    """P2's scorer preflight returned None on ImageError through the same channel as a legitimate
    'too small to split', so a loop with a defensive `if vals:` produced a silent empty result.
    This must RAISE: there is no falsy return a caller can skip past."""
    assert issubclass(oracle.OraclePreconditionFailed, Exception)
    assert not issubclass(oracle.OraclePreconditionFailed, (ValueError, KeyError)), (
        "a bespoke precondition failure must not masquerade as a routine ValueError/KeyError that "
        "surrounding code is likely to be catching already")
    src = _code_only(ORACLE_SRC.read_text(encoding="utf-8"))
    body = src[src.index("def assert_oracle_fixture"):]
    body = body[:body.index("\ndef ") if "\ndef " in body[1:] else len(body)]
    assert "return None" not in body and "return False" not in body


# ------------------------------------------------------- the fold parameter + path separation

def test_the_held_out_fold_set_defaults_to_the_committed_measurement(oracle):
    assert oracle.parse_fold_set(None, oracle.HELD_OUT_FOLDS) == {9, 10}
    assert oracle.fold_tag({9, 10}) == "", (
        "the default must keep the committed artifact path, or re-running the shipped command "
        "silently relocates the record")
    assert oracle.fold_tag({9}) == "_fold9"


def test_the_working_and_final_only_shards_are_distinct(oracle):
    assert oracle.WORKING_HELD_OUT_FOLDS == {9}
    assert oracle.FINAL_ONLY_FOLDS == {10}
    assert not (oracle.WORKING_HELD_OUT_FOLDS & oracle.FINAL_ONLY_FOLDS)
    assert oracle.WORKING_HELD_OUT_FOLDS | oracle.FINAL_ONLY_FOLDS == oracle.HELD_OUT_FOLDS


@pytest.mark.parametrize("spec", ["7", "8", "1", "9,7", "0", "11", "x"])
def test_scoring_on_a_training_fold_or_a_nonsense_fold_is_refused(oracle, spec):
    with pytest.raises(ValueError):
        oracle.parse_fold_set(spec, oracle.HELD_OUT_FOLDS)


def test_a_fold9_result_cannot_overwrite_the_fold_9_and_10_record():
    """Both artifacts must exist as SEPARATE files carrying DIFFERENT numbers. They are two
    measurements on two different shards; collapsing them onto one path is the documented
    shared-key silent overwrite."""
    for p in (ORACLE_F9, ORACLE_F910, NULL_F9, NULL_F910):
        assert p.exists(), f"{p.name} missing"
    a = json.loads(ORACLE_F910.read_text(encoding="utf-8"))
    b = json.loads(ORACLE_F9.read_text(encoding="utf-8"))
    assert a["fold_contract"]["held_out"] == [9, 10]
    assert b["fold_contract"]["held_out"] == [9]
    assert b["cells"]["K562"]["n_held_out"] < a["cells"]["K562"]["n_held_out"]


# ------------------------------------------------------- the ceiling fix, as a record

def test_the_ceiling_arm_draws_randomly_rather_than_taking_the_first_n():
    code = _code_only(NULL_SRC.read_text(encoding="utf-8"))
    assert "random.Random(" in code, (
        "the real_ceiling arm must draw a RANDOM matched-n sample. quintiles() returns each "
        "stratum sorted by activity ascending, so a first-n slice samples every stratum's bottom "
        "edge and truncates the top quintile's tail -- measured at -4.14 sd in K562")
    assert "[:a.per_stratum]" not in code, (
        "a first-n slice of an activity-SORTED stratum was the defect; it must not come back")


def test_the_null_artifact_carries_the_ceilings_own_noise_scale():
    d = json.loads(NULL_F9.read_text(encoding="utf-8"))
    for cell, c in d["cells"].items():
        g = c["ceiling_diagnostics"]
        assert len(g["matched_n_rho_by_seed"]) >= 10, f"{cell}: too few draws to estimate a spread"
        assert g["matched_n_sd"] is not None
        # The bar refuses frozen numeric tolerances: a tolerance is a MEASURED spread.
        assert 0.0 < g["matched_n_sd"] < 0.2, f"{cell}: implausible ceiling sd {g['matched_n_sd']}"
        assert g["all_held_out_records_rho"] > 0


def test_the_fraction_is_divided_by_the_all_records_ceiling_not_a_single_draw():
    """A denominator that moves with the seed makes the verdict rule seed-dependent."""
    d = json.loads(NULL_F9.read_text(encoding="utf-8"))
    for cell, c in d["cells"].items():
        assert d["ceiling_relative"]["ceiling_spearman"][cell] == pytest.approx(
            c["ceiling_diagnostics"]["all_held_out_records_rho"]), (
            f"{cell}: the published ceiling must be the all-records estimate")
        frac = d["ceiling_relative"]["stratified_frac_of_ceiling"][cell]
        num = c["arms"]["stratified"]["target_monotonicity_spearman"]
        den = d["ceiling_relative"]["ceiling_spearman"][cell]
        assert frac == pytest.approx(num / den, abs=1e-4)
    assert d["ceiling_relative"]["frac_sd_note"], "the fraction's uncertainty must be stated"


def test_the_null_no_longer_claims_the_held_out_folds_are_untouched():
    """The real_ceiling arm READS the held-out folds. The old key said folds_untouched, which its
    own code refutes, and the reservation had already counted that read as a spent query."""
    d = json.loads(NULL_F9.read_text(encoding="utf-8"))
    assert "folds_untouched" not in d["leak"]
    assert d["leak"]["folds_read_by_real_ceiling_arm_never_fitted"] == [9]
    assert d["leak"]["folds_not_read_at_all"] == [10], (
        "fold 10 is FINAL_ONLY; a decision-grade run must not read it")


def test_both_fold9_artifacts_record_that_the_precondition_ran():
    for p in (ORACLE_F9, NULL_F9):
        d = json.loads(p.read_text(encoding="utf-8"))
        assert d["precondition"]["status"] == "PASSED", f"{p.name} has no passing precondition"
        assert d["precondition"]["precondition"] == "oracle-feature-fixture-v1"


def test_the_oracle_gate_still_passes_on_fold_9_alone():
    """The decision the folds-9+10 measurement informed, re-measured without the reserved shard."""
    d = json.loads(ORACLE_F9.read_text(encoding="utf-8"))
    assert d["verdict"] == "ORACLE_VALID_PROCEED"
    for cell in ("HepG2", "K562"):
        g = d["cells"][cell]["geometries"]["element_200nt_adaptors_stripped"]
        assert g["passes_gate"] is True
        assert g["margin_over_gc"] >= 0.05
        assert g["best_learned"] >= 0.20
