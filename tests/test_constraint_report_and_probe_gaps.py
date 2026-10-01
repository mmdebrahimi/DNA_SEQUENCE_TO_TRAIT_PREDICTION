"""Gap coverage for `constraints/report.py` and the kill-count probe (test-epilogue, 2026-10-01).

The plan's two files cover the headline contracts (three derived states, inert vs satisfied, the
non-vacuous `assert_no_silent_inert`, and the real-data artifact assertions). What they leave:

REPORT
* `note_not_evaluated` and `status()` on an unseen key -- the `not_evaluated` state is only ever read
  off the committed artifact, never constructed directly;
* `not_evaluated` must NOT require an acknowledgement. Only `inert` does, and conflating the two would
  force a false claim about a constraint that never ran;
* an acknowledgement for a key with no counts is INVISIBLE in `as_dict()`, so a stale entry can neither
  help nor be noticed.

PROBE -- and this one is a real finding, not a missing assertion
* the probe's acknowledgement loop is `acknowledgements.get(key, "inert on curated reference input")`.
  The catch-all default means a NEW inert law is auto-acknowledged with a generic string, so
  `assert_no_silent_inert()` can never raise inside this probe and its **exit-1 branch is structurally
  unreachable**. Measured, not argued: injecting an inert law makes the probe print "every inert law
  carries an explicit reason" and return 0.

  That is the same shape as the two guards this repo already retracted (the HIV operator-prefix filter
  over zero censored values; the ResFinder POINT-row exclusion over zero POINT rows) -- a check that
  cannot fire. It is pinned from the test side rather than patched: the tests below require every
  CURRENT inert law to have a SPECIFIC reason, so adding a law that falls through to the default fails
  here even though the probe stays green.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib
import re

import pytest

from dna_decode.constraints.registry import (
    AXIOMATIC,
    INAPPLICABLE,
    REFUTED,
    SATISFIED,
    UNIVERSAL,
    Constraint,
    ConstraintVerdict,
)
from dna_decode.constraints.report import (
    ACTIVE,
    INERT,
    NOT_EVALUATED,
    ConstraintReport,
    SilentInertConstraintError,
)

REPO = pathlib.Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "constraint_kill_count_probe_gaps", REPO / "scripts" / "constraint_kill_count_probe.py")
P = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P)

PROBE_SRC = (REPO / "scripts" / "constraint_kill_count_probe.py").read_text(encoding="utf-8")


# ===================================================================== report.py

def test_status_of_an_unseen_key_is_not_evaluated_and_does_NOT_create_an_entry():
    """Querying must not materialise a row -- otherwise `as_dict` would grow a phantom constraint for
    every key anyone ever asked about."""
    r = ConstraintReport()
    assert r.status("never_heard_of_it") == NOT_EVALUATED
    assert r.counts == {}
    assert r.as_dict()["constraints"] == {}


def test_note_not_evaluated_makes_a_constraint_VISIBLE_with_zero_counts():
    """The point of the method: a law that never ran must appear in the artifact rather than vanish,
    because a missing row reads as 'no problem here'."""
    r = ConstraintReport()
    r.note_not_evaluated("law_x")
    d = r.as_dict()
    assert d["constraints"]["law_x"] == {
        "n_evaluated": 0, "n_refuted": 0, "n_satisfied": 0, "n_inapplicable": 0,
        "status": NOT_EVALUATED, "acknowledged_reason": None}
    assert d["summary"] == {"n_active": 0, "n_inert": 0, "n_not_evaluated": 1}


def test_note_not_evaluated_is_idempotent_and_never_clobbers_real_counts():
    r = ConstraintReport()
    r.record("law_x", ConstraintVerdict(REFUTED, "d"))
    r.note_not_evaluated("law_x")
    r.note_not_evaluated("law_x")
    assert r.status("law_x") == ACTIVE
    assert r.counts["law_x"].n_refuted == 1


def test_a_NOT_EVALUATED_constraint_needs_NO_acknowledgement():
    """Only `inert` is a claim about protection. Requiring a reason for a law that never ran would make
    the author assert something they did not observe -- the blanket-guard failure mode."""
    r = ConstraintReport()
    r.note_not_evaluated("never_ran")
    r.record("only_inapplicable", ConstraintVerdict(INAPPLICABLE, "no precondition"))
    assert r.inert_keys() == []
    r.assert_no_silent_inert()                      # must not raise

    # non-vacuous: an inert one in the same report DOES raise
    r.record("ran_and_found_nothing", ConstraintVerdict(SATISFIED, "clean"))
    assert r.inert_keys() == ["ran_and_found_nothing"]
    with pytest.raises(SilentInertConstraintError, match="ran_and_found_nothing"):
        r.assert_no_silent_inert()


def test_an_acknowledgement_for_a_key_with_no_counts_is_INVISIBLE_in_the_artifact():
    """`as_dict` iterates `counts`, so an acknowledgement keyed on a law that never ran is carried in
    memory and rendered nowhere. It can neither mask a real inert law (that law has counts, so it gets
    its own row) nor be pruned by a reader. Documented as a known limit of the ledger."""
    r = ConstraintReport()
    r.acknowledge_inert("ghost_law", "a reason for a law that was never evaluated")
    assert "ghost_law" in r.acknowledged
    assert "ghost_law" not in r.as_dict()["constraints"]
    assert r.as_dict()["summary"]["n_not_evaluated"] == 0


def test_a_stale_acknowledgement_cannot_silence_a_DIFFERENT_inert_law():
    """The safety consequence of the above, stated positively."""
    r = ConstraintReport()
    r.acknowledge_inert("ghost_law", "stale entry left behind by an earlier run")
    r.record("real_law", ConstraintVerdict(SATISFIED, "clean"))
    with pytest.raises(SilentInertConstraintError, match="real_law"):
        r.assert_no_silent_inert()


def test_acknowledge_inert_overwrites_rather_than_accumulating():
    r = ConstraintReport()
    r.record("k", ConstraintVerdict(SATISFIED))
    r.acknowledge_inert("k", "first reason")
    r.acknowledge_inert("k", "second reason")
    assert r.acknowledged["k"] == "second reason"
    assert r.as_dict()["constraints"]["k"]["acknowledged_reason"] == "second reason"


def test_an_acknowledgement_reason_must_be_non_empty_and_the_check_is_non_vacuous():
    r = ConstraintReport()
    with pytest.raises(ValueError, match="needs a reason"):
        r.acknowledge_inert("k", "")
    r.acknowledge_inert("k", "a real reason")
    assert r.acknowledged == {"k": "a real reason"}


def test_inert_and_active_key_lists_are_SORTED_so_the_artifact_is_stable():
    """An unsorted list would make the committed artifact churn on dict-ordering alone."""
    r = ConstraintReport()
    for key in ("zeta", "alpha", "mu"):
        r.record(key, ConstraintVerdict(SATISFIED))
    for key in ("yankee", "bravo"):
        r.record(key, ConstraintVerdict(REFUTED, "x"))
    assert r.inert_keys() == ["alpha", "mu", "zeta"]
    assert r.active_keys() == ["bravo", "yankee"]
    assert list(r.as_dict()["constraints"]) == sorted(r.as_dict()["constraints"])


def test_the_three_summary_counters_PARTITION_the_constraints():
    """Every row must land in exactly one bucket; a status that counted in two would overstate the
    layer's reach."""
    r = ConstraintReport()
    r.record("a", ConstraintVerdict(REFUTED, "x"))
    r.record("b", ConstraintVerdict(SATISFIED))
    r.record("c", ConstraintVerdict(INAPPLICABLE, "no precondition"))
    r.note_not_evaluated("d")
    s = r.as_dict()["summary"]
    assert (s["n_active"], s["n_inert"], s["n_not_evaluated"]) == (1, 1, 2)
    assert s["n_active"] + s["n_inert"] + s["n_not_evaluated"] == len(r.counts)


def test_a_mixed_constraint_with_one_refutation_among_many_passes_reads_ACTIVE():
    """`status` is `n_refuted > 0`, not a rate -- one kill is enough to stop it posing as inert."""
    r = ConstraintReport()
    for _ in range(50):
        r.record("k", ConstraintVerdict(SATISFIED))
    r.record("k", ConstraintVerdict(REFUTED, "the one"))
    c = r.counts["k"]
    assert (c.n_evaluated, c.n_refuted, c.n_satisfied) == (51, 1, 50)
    assert r.status("k") == ACTIVE
    r.assert_no_silent_inert()


# ===================================================================== the probe's inert gate

def test_the_probes_inert_acknowledgement_loop_has_NO_catch_all_default():
    """FIXED 2026-10-01. The loop used to read
    `acknowledgements.get(key, "inert on curated reference input")`, which auto-acknowledged ANY new
    inert law and made the probe's exit-1 contract decoration. The generic default is gone: only keys
    present in the hand-written dict are acknowledged, so an unlisted inert law falls through to the
    assertion. Source-level pin so the next reader sees the invariant without running anything."""
    assert 'acknowledgements.get(key, "inert on curated reference input")' not in PROBE_SRC, (
        "the catch-all default is back; the probe's exit-1 branch is unreachable again")
    assert "if key in acknowledgements:" in PROBE_SRC, (
        "the acknowledgement loop no longer gates on explicit membership")
    assert "NO GENERIC DEFAULT" in PROBE_SRC, "the reason for the invariant must stay documented"


def test_the_probes_exit_1_branch_is_REACHABLE(tmp_path, monkeypatch):
    """MEASURED, not argued, in both directions. Injecting a brand-new inert law that the probe's
    acknowledgement dict knows nothing about now makes the probe FAIL, which is what makes its exit-1
    contract real rather than decorative. Before the fix this same injection returned 0 and printed
    that every inert law carried an explicit reason -- the exact shape of the two guards this repo has
    already retracted, reproduced by the script written to prevent them."""
    from dna_decode.constraints import universal as U

    extra = Constraint("injected_inert_law", UNIVERSAL, AXIOMATIC,
                       lambda ctx: ConstraintVerdict(SATISFIED, "always clean"))
    monkeypatch.setattr(U, "UNIVERSAL_CONSTRAINTS", U.UNIVERSAL_CONSTRAINTS + (extra,))

    assert P.main(["--out", str(tmp_path / "o.json")]) == 1, (
        "an unacknowledged inert law did not fail the probe; the catch-all default is back")


def test_the_probe_still_exits_0_on_the_real_reference_set(tmp_path):
    """Non-vacuity for the test above: the failure must be caused by the injection, not by the probe
    being broken for everything."""
    assert P.main(["--out", str(tmp_path / "ok.json")]) == 0


def test_every_CURRENT_inert_law_has_a_SPECIFIC_reason_not_the_generic_fallback():
    """The test-side replacement for the unreachable exit-1 branch: run the probe for real and require
    each inert law's reason to be one the author actually wrote for THAT law. Adding a law without an
    acknowledgement entry fails here."""
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        out = pathlib.Path(d) / "o.json"
        assert P.main(["--out", str(out)]) == 0
        rows = json.loads(out.read_text(encoding="utf-8"))["kill_counts"]["constraints"]

    inert = {k: v["acknowledged_reason"] for k, v in rows.items() if v["status"] == INERT}
    assert inert, "no inert law at all -- the probe's shape changed; re-derive this guard"
    for key, reason in inert.items():
        assert reason, f"{key} is inert with no reason"
        assert reason != "inert on curated reference input", (
            f"{key} fell through to the generic fallback. Add a specific entry to the probe's "
            "`acknowledgements` dict saying why THAT law removed nothing on this input.")
        assert len(reason) > 40, f"{key}'s reason is a placeholder: {reason!r}"


def test_the_acknowledgements_dict_covers_EXACTLY_the_shipped_laws():
    """A hand-written dict beside a live tuple of laws is this repo's recorded drift shape. Keys are
    read out of the probe's source and diffed against `UNIVERSAL_CONSTRAINTS`, so both a new law and a
    stale entry are caught."""
    from dna_decode.constraints import universal as U

    block = PROBE_SRC.split("acknowledgements = {", 1)[1].split("\n    }", 1)[0]
    declared = set(re.findall(r'^\s{8}"([a-z_]+)":', block, flags=re.M))
    shipped = {c.key for c in U.UNIVERSAL_CONSTRAINTS}
    assert declared, "could not parse the acknowledgements dict; re-check this guard"
    assert declared == shipped, (
        f"missing acknowledgements for {sorted(shipped - declared)}; "
        f"stale entries for {sorted(declared - shipped)}")


# ===================================================================== the probe's degraded paths

def test_a_MISSING_reference_is_REPORTED_not_silently_dropped(tmp_path, monkeypatch):
    """A partial read still produces a number, so the artifact must say which references were absent --
    otherwise a shrinking cohort looks like a clean run."""
    kept = {k: v for k, v in list(P.REFERENCE_CDS.items())[:2]}
    monkeypatch.setattr(P, "REFERENCE_CDS", {**kept, "ghost": "data/not_a_real_reference.fna"})
    out = tmp_path / "o.json"
    assert P.main(["--out", str(out)]) == 0
    art = json.loads(out.read_text(encoding="utf-8"))
    assert art["missing_references"] == ["data/not_a_real_reference.fna"]
    assert sorted(art["reference_cds_scored"]) == sorted(kept)
    assert "ghost" not in art["per_cds_verdicts"]


def test_the_clade_demo_reads_UNAVAILABLE_rather_than_real_when_the_fixture_is_absent(
        tmp_path, monkeypatch):
    """The honesty hinge of the probe. With the Mycoplasma fixture gone the demo must label itself
    `unavailable` and carry no verdicts -- never default to `real`, which is what the committed
    artifact's `source == 'real'` assertion is trusting."""
    monkeypatch.setattr(P, "MYCOPLASMA_FIXTURE", "data/mycoplasma_ref/absent.fna")
    out = tmp_path / "o.json"
    assert P.main(["--out", str(out)]) == 0
    demo = json.loads(out.read_text(encoding="utf-8"))["clade_demo"]
    assert demo == {"source": "unavailable"}
    assert "verdict_under_clade_table" not in demo


def test_without_the_clade_demo_NO_law_is_active_which_is_what_makes_the_demo_load_bearing(
        tmp_path, monkeypatch):
    """The probe's pre-stated prediction is that the ONLY active law is `no_internal_stop_codon`, and
    only via the Mycoplasma CDS scored under the wrong code. Remove the fixture and the layer should go
    entirely inert -- which is the honest reading of 'curated references already satisfy the laws' and
    proves the active count is not coming from the reference set."""
    monkeypatch.setattr(P, "MYCOPLASMA_FIXTURE", "data/mycoplasma_ref/absent.fna")
    out = tmp_path / "o.json"
    assert P.main(["--out", str(out)]) == 0
    kc = json.loads(out.read_text(encoding="utf-8"))["kill_counts"]
    assert kc["summary"]["n_active"] == 0, "a reference CDS refuted a law; the prediction was wrong"
    assert kc["summary"]["n_inert"] >= 1


def test_read_fasta_tolerates_blank_lines_and_multiple_records():
    """`read_fasta` concatenates every non-header line, so a multi-record FASTA silently becomes one
    sequence. Pinned because the probe's references are single-record and the behaviour would bite a
    caller who pointed it at a CDS set."""
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        f = pathlib.Path(d) / "multi.fna"
        f.write_text(">a\nacg\n\n>b\ntga\n", encoding="utf-8")
        assert P.read_fasta(f) == "ACGTGA"


def test_the_mycoplasma_fixture_is_committed_in_frame_and_clean_under_table_4():
    """The fixture is the whole clade demo. Check the sequence itself, not just the artifact that
    summarises it: in frame, refuted under table 1, satisfied under table 4."""
    from dna_decode.constraints import universal as U

    p = REPO / P.MYCOPLASMA_FIXTURE
    if not p.exists():
        pytest.skip("mycoplasma fixture absent")
    cds = P.read_fasta(p)
    assert len(cds) % 3 == 0, f"fixture length {len(cds)} is not in frame"
    assert set(cds) <= set("ACGT"), "fixture carries ambiguity codes"
    assert U.no_internal_stop_codon({"cds": cds}).verdict == REFUTED
    assert U.no_internal_stop_codon(
        {"cds": cds, "clade": "bacteria.mycoplasma"}).verdict == SATISFIED
    assert U.cds_length_multiple_of_three({"cds": cds}).verdict == SATISFIED
