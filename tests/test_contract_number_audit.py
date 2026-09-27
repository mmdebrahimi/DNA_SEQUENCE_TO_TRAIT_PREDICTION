"""Tests for the contract-number audit.

The load-bearing tests are the NON-VACUITY ones. This audit reached a clean 0/130 once by accepting any
number traceable to the package source -- which accepted 12 of 12 randomly-generated numbers. A check that
passes everything reads exactly like a check that found nothing, so the guards below plant drift and require
it to be caught.
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from scripts.contract_number_audit import (
    ADJUDICATED_BENIGN,
    DECIMAL_RE,
    _cited_artifacts,
    _expand_braces,
    _matches_by_rounding,
    audit,
)

REPO = Path(__file__).resolve().parent.parent
ARTIFACT = REPO / "wiki" / "contract_number_audit_2026-09-24.json"

_REPORT_CACHE: dict = {}


def _cached_report() -> dict:
    """One shared report for the read-only assertions.

    Step 1 replaced a 10-artifact sample with a full-pool scan, which costs ~13 s per `audit()` call.
    Every test calling it independently put this one file at 7m20s. The report is never mutated here, so
    the tests that merely assert ON it share one. Tests that are ABOUT re-running -- determinism, and the
    circularity guard -- deliberately call `audit()` directly instead.
    """
    if "rep" not in _REPORT_CACHE:
        _REPORT_CACHE["rep"] = audit()
    return _REPORT_CACHE["rep"]


# --- citation parsing ----------------------------------------------------------------------------

def test_brace_expansion():
    """The repo cites artifacts in brace form; not expanding it made the audit invent drift on the
    inverse cell's +53.0%, whose source sits in a roundtrip artifact the regex never loaded."""
    assert set(_expand_braces("wiki/a_2026-09-24.{md,json}")) == {
        "wiki/a_2026-09-24.md", "wiki/a_2026-09-24.json"}
    multi = set(_expand_braces("wiki/f_{x,y}_2026-07-1{6,7}.{md,json}"))
    assert len(multi) == 8
    assert "wiki/f_x_2026-07-16.md" in multi
    assert _expand_braces("wiki/plain.json") == ["wiki/plain.json"]


def test_cited_artifacts_tries_both_sidecars_for_a_bare_stem():
    got = _cited_artifacts("see wiki/foo_2026-01-01 for detail")
    assert "wiki/foo_2026-01-01.json" in got and "wiki/foo_2026-01-01.md" in got


def test_cited_artifacts_strips_trailing_punctuation():
    assert "wiki/foo.json" in _cited_artifacts("(wiki/foo.json).")


# --- what counts as a cited NUMBER ---------------------------------------------------------------

def test_identifier_fragments_are_not_numbers():
    """canFam3.1 / GCA_000005845.2 / v1.2 are identifiers. Capturing them produced a false
    CANDIDATE_DRIFT on the dog cell's '3.1', which is an assembly name."""
    for ident in ("canFam3.1", "GCA_000005845.2", "v1.2", "hg19.2"):
        assert DECIMAL_RE.findall(ident) == [], f"{ident} should not yield a measurement"


def test_real_decimals_are_still_captured():
    assert DECIMAL_RE.findall("accuracy 0.8125 and spec 0.939") == ["0.8125", "0.939"]


# --- rounding-aware matching ---------------------------------------------------------------------

def test_rounding_match_accepts_prose_rounding():
    """Prose rounds: a cited 0.398 must match an artifact's 0.3977. Without this the audit flagged
    SLCO1B1's correct figure."""
    assert _matches_by_rounding("0.398", [0.3977]) == pytest.approx(0.3977)
    assert _matches_by_rounding("0.40", [0.3977]) == pytest.approx(0.3977)


def test_rounding_match_accepts_percent_fraction_pairing():
    assert _matches_by_rounding("53.0", [0.5295]) == pytest.approx(0.5295)
    assert _matches_by_rounding("0.29", [29.0]) == pytest.approx(29.0)


def test_rounding_match_rejects_a_genuinely_different_number():
    """It must NOT be a universal acceptor -- that is the vacuity failure."""
    assert _matches_by_rounding("0.79", [0.2128]) is None
    assert _matches_by_rounding("0.926", [0.454]) is None


# --- NON-VACUITY: planted drift must be caught ---------------------------------------------------

def test_planted_drift_is_caught():
    """THE guard. Feed the real matcher a number that is in no artifact and confirm it does not pass."""
    assert _matches_by_rounding("0.7431", [0.1, 0.2, 0.9999]) is None


def test_random_numbers_do_not_pass_against_a_small_artifact_set():
    """The vacuity regression, stated as a property: 12 random 4-decimal numbers must not all match a
    handful of artifact values. The earlier package-source variant accepted 12 of 12."""
    random.seed(7)
    art = [0.3977, 0.8125, 0.2128, 64.0, 88.0]
    toks = [f"0.{random.randint(1000, 9999)}" for _ in range(12)]
    matched = [t for t in toks if _matches_by_rounding(t, art) is not None]
    assert len(matched) <= 1, f"too permissive: {matched}"


def test_adjudicated_benign_is_small_and_every_entry_has_a_reason():
    """Exemptions are named one at a time, with a reason, precisely because the blanket version was
    vacuous. A growing table is a smell, so this pins it small."""
    assert len(ADJUDICATED_BENIGN) <= 5
    for (cell_id, tok), why in ADJUDICATED_BENIGN.items():
        assert cell_id and tok
        assert len(why) > 40, f"{cell_id} {tok} needs a real reason, got {why!r}"


# --- live audit ----------------------------------------------------------------------------------

def test_live_audit_is_clean_and_non_trivial():
    """Clean AND actually checking something -- a zero over zero checks is not a result."""
    rep = _cached_report()
    assert rep["n_numbers_checked"] >= 100
    assert rep["n_cells_audited"] >= 10
    assert rep["n_candidate_drift"] == 0, rep["verdict"]


def test_live_audit_would_flag_an_injected_bad_citation():
    """End-to-end non-vacuity: a cell whose prose cites a number absent from its artifact must be
    reported, using the REAL audit path on a synthetic contract."""
    import scripts.contract_number_audit as m

    class FakeCell:
        cell_id = "synthetic:test:cell"
        claim = "accuracy 0.7431 on the panel (wiki/contract_number_audit_2026-09-24.json)"
        validation_slice = label_provenance = demotion_rule = claim_status = ""

    real = m.audit.__globals__["__builtins__"]  # touch nothing; we monkeypatch the import instead
    assert real is not None
    # drive the internals directly rather than patching the registry import
    blob = FakeCell.claim
    arts = m._cited_artifacts(blob)
    text = "\n".join(t for t in (m._load_artifact_text(a) for a in arts) if t)
    assert text, "fixture artifact must resolve, else this test proves nothing"
    nums = m.DECIMAL_RE.findall(blob)
    assert "0.7431" in nums
    assert not any(v in text for v in m._number_variants("0.7431"))
    assert m._matches_by_rounding("0.7431", m._artifact_numbers(text)) is None


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_committed_artifact_records_its_own_posture():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert "TRIAGE" in d["posture"].upper()
    assert d["n_numbers_checked"] >= 100


# --- coverage disclosure -------------------------------------------------------------------------
# `scope` always said the audit covers "cells whose prose names a wiki/ artifact", so the exclusion
# was disclosed -- but not its SIZE, which is what says whether the audited set is near-complete or a
# fraction. It is a fraction, and the skipped set is LARGER than the checked one.

def test_coverage_block_is_present_and_internally_consistent():
    rep = _cached_report()
    assert rep["n_numbers_unverifiable"] == sum(u["n_numbers"] for u in rep["cells_unverifiable"])
    assert rep["n_cells_with_numbers_but_no_citation"] == len(rep["cells_unverifiable"])
    for u in rep["cells_unverifiable"]:
        assert u["n_numbers"] == len(u["numbers"]) > 0, "a cell with no numbers is not a coverage gap"
        assert u["why"]


def test_coverage_is_non_trivial_so_the_disclosure_is_not_decoration():
    """If nothing were excluded the block would be noise; the point is that plenty is."""
    rep = _cached_report()
    assert rep["n_cells_with_numbers_but_no_citation"] > 0
    assert rep["n_numbers_unverifiable"] > 0


def test_unverifiable_cells_are_disjoint_from_audited_cells():
    """A cell is either checked against an artifact or reported as uncheckable — never both."""
    rep = _cached_report()
    audited = {r["cell_id"] for r in rep["cells"]}
    skipped = {u["cell_id"] for u in rep["cells_unverifiable"]}
    assert not (audited & skipped)


def test_coverage_does_not_contaminate_the_drift_verdict():
    """Unchecked is not wrong. The verdict must key on drift alone, or an uncited number would read
    as a defect and the audit would stop being a triage funnel."""
    rep = _cached_report()
    expected = "NOTHING_TO_ADJUDICATE" if rep["n_candidate_drift"] == 0 else "ADJUDICATION_REQUIRED"
    assert rep["verdict"] == expected
    assert rep["n_numbers_checked"] == sum(r["checked"] for r in rep["cells"])


def test_an_uncited_cell_whose_only_numbers_are_adjudicated_bars_is_not_flagged():
    """Non-vacuity in the OTHER direction: a threshold has no artifact to live in, so a cell carrying
    only adjudicated bars is not a coverage gap. Exercises the real filter on a synthetic cell."""
    import scripts.contract_number_audit as m

    (cell_id, tok) = next(iter(ADJUDICATED_BENIGN))
    blob = f"the bar is {tok}"
    bare = [t for t in dict.fromkeys(m.DECIMAL_RE.findall(blob))
            if not m.NON_MEASUREMENT.match(t) and (cell_id, t) not in ADJUDICATED_BENIGN]
    assert bare == [], f"{tok} should be exempt for {cell_id}"
    # ... while the same number on a DIFFERENT cell is still a coverage gap (the table is keyed by both)
    bare_other = [t for t in dict.fromkeys(m.DECIMAL_RE.findall(blob))
                  if not m.NON_MEASUREMENT.match(t) and ("other:cell:x", t) not in ADJUDICATED_BENIGN]
    assert bare_other == [tok]


# --- decoy discrimination ------------------------------------------------------------------------
# A citation only VERIFIES a number if the number could have failed to match. These tests exist
# because the naive version of this check is vacuous: a cell whose figures are few and round matches
# a numeric-dense JSON about anything, so "found" alone cannot distinguish provenance from luck.

def test_discrimination_is_reported_with_a_legal_value():
    rep = _cached_report()
    for r in rep["cells"]:
        assert r["discrimination"] in {"HIGH", "MODERATE", "LOW", "n/a"}
        assert 0 <= r["decoy_full_match"] <= r["n_decoys"]
        if r["checked"] and r["n_decoys"]:
            assert r["discrimination"] != "n/a"


def test_discrimination_actually_discriminates():
    """Uniformly HIGH would mean the decoys are too weak to ever match; uniformly LOW would mean the
    control is broken. The signal is only useful if it separates cells."""
    rep = _cached_report()
    seen = {r["discrimination"] for r in rep["cells"] if r["checked"]}
    assert "HIGH" in seen, "no cell resists decoys — the control is not exercised"
    assert len(seen) > 1, f"discrimination is constant at {seen} — it separates nothing"


def test_a_cell_with_many_distinctive_numbers_resists_decoys():
    """Concrete anchor: salmserovar carries 30 measured figures; unrelated artifacts should almost never
    contain all of them. If this ever fails, the decoy match has become too permissive.

    EDITED when the 10-artifact sample was replaced by the full-pool scan, and the reason is a
    measurement, not a convenience. This asserted `decoy_full_match == 0`, which was reachable only
    because 10 draws is too few to find the coincidence: against ~600 trials, salmserovar -- the MOST
    distinctive cell in the registry -- still has 1 five-file haystack containing all 30 of its numbers.
    Exactly-zero is therefore not a property any cell has, so the assertion is now on the RATE, which is
    what the grade is defined on. Kept deliberately tight (1%) so it still fails if matching loosens.
    """
    rep = _cached_report()
    rows = {r["cell_id"]: r for r in rep["cells"]}
    r = rows["typing:Salmonella:salmserovar"]
    assert r["checked"] >= 20
    assert r["decoy_match_rate"] < 0.01, r["decoy_match_rate"]
    assert r["discrimination"] == "HIGH"


def test_low_discrimination_does_not_count_as_drift():
    """Weak evidence is not a defect. A cell can be clean AND poorly verified, and the verdict must
    key on drift alone or the audit stops being a triage funnel."""
    rep = _cached_report()
    assert rep["verdict"] == ("NOTHING_TO_ADJUDICATE" if rep["n_candidate_drift"] == 0
                             else "ADJUDICATION_REQUIRED")
    low = [r for r in rep["cells"] if r["discrimination"] == "LOW"]
    for r in low:
        assert r["candidate_drift"] == [], "a LOW-discrimination cell is unverified, not drifted"
    assert rep["n_cells_low_discrimination"] == len(low)
    assert rep["cells_low_discrimination"] == sorted(r["cell_id"] for r in low)


def test_discrimination_is_deterministic():
    """A report whose numbers move run to run is not auditable. Determinism now comes from scanning the
    WHOLE pool rather than from a fixed seed -- docstring corrected; the assertion is unchanged."""
    a, b = audit(), audit()
    assert [(r["cell_id"], r["decoy_full_match"]) for r in a["cells"]] == \
           [(r["cell_id"], r["decoy_full_match"]) for r in b["cells"]]


# --- Step 1: the control is a deterministic, size-matched, FULL-POOL scan -------------------------
# The 10-artifact sample it replaced was underpowered to the point of being wrong in both directions.
# These tests exist to make that un-reintroducible.

def test_the_control_has_no_seed_and_no_sample_size():
    """A seed is only needed by a sampler. Their absence is the structural proof this is a full scan."""
    import scripts.contract_number_audit as mod
    assert not hasattr(mod, "_DECOY_SEED")
    assert not hasattr(mod, "_N_DECOYS")
    assert "random" not in mod._decoy_full_match_count.__code__.co_names


def test_every_graded_cell_is_compared_against_the_whole_pool():
    """`n_decoys` must be the pool size (hundreds), not 10. A cell graded off 10 trials is the defect."""
    rep = _cached_report()
    graded = [r for r in rep["cells"] if r["checked"] and r["n_decoys"]]
    assert graded
    for r in graded:
        # grp N means the pool is consumed in N-sized groups, so trials ~= pool/N
        assert r["n_decoys"] * r["decoy_group_size"] >= 100, r


def test_full_pool_anchors_measured_before_implementation():
    """Anchors measured BEFORE this control was written, so a mismatch means the code is wrong.

    The planning measurement was 6/602, 1/602, 308/602, 355/602. The denominators are 599 here and two
    counts are one lower, because the pool now EXCLUDES this script's own accrued outputs -- artifacts
    that list the very numbers under test. That exclusion is the fix in
    `test_the_pool_excludes_the_audits_own_output`; the rates are unchanged to three decimals.
    """
    rep = _cached_report()
    rows = {r["cell_id"]: r for r in rep["cells"]}
    for cell_id, expect in (("typing:bacteria:mlst", 6),
                            ("finder:bacteria:resfinder", 1),
                            ("typing:Klebsiella:ktype", 307),
                            ("finder:Escherichia_coli:pointfinder", 354)):
        assert rows[cell_id]["decoy_full_match"] == expect, (cell_id, rows[cell_id])
        assert rows[cell_id]["n_decoys"] == 599, (cell_id, rows[cell_id])


def test_the_pool_excludes_the_audits_own_output():
    """Circularity guard. `main()` writes wiki/contract_number_audit_<date>.json, and that artifact
    RECORDS the number tokens it checked -- so an accrued run is a decoy guaranteed to match. Found by
    measurement: `finder:any:forward` moved 74/300 -> 75/300 between two runs, the second having
    inherited the first's artifact."""
    import scripts.contract_number_audit as mod
    own = sorted(p.name for p in (REPO / "wiki").glob("contract_number_audit_*.json"))
    assert own, "no accrued audit artifacts — this guard would be vacuous"
    rep = mod.audit()
    # every graded cell's trial count must be consistent with a pool that dropped those files
    pool_all = len(list((REPO / "wiki").glob("*.json")))
    for r in rep["cells"]:
        if r["checked"] and r["decoy_group_size"] == 1:
            assert r["n_decoys"] <= pool_all - len(own), (r["cell_id"], r["n_decoys"], pool_all, len(own))


def test_essentiality_does_not_grade_high_the_case_that_exposed_the_underpowering():
    """THE load-bearing non-vacuity test of Step 1.

    `essentiality:any:essentiality`'s seven numbers matched 0 of 10 sampled artifacts, which graded HIGH
    -- the strongest verdict available -- while against the full pool 45 of 600 completely unrelated
    artifacts contain ALL seven (a yeast benchmark, two TB baselines, staleness snapshots, a pigment run).
    If this ever grades HIGH again, the control has gone back to being underpowered.

    NOTE the plan for this step predicted LOW; measured, it is MODERATE at 7.5%. The substantive claim is
    "not HIGH", which is what is asserted -- bending the band to make the predicted word come true is the
    failure mode this repo has recorded twice as a mis-specified bar.
    """
    import scripts.contract_number_audit as mod
    from dna_decode.data.cell_registry import cells

    pool = sorted(p.relative_to(REPO).as_posix() for p in (REPO / "wiki").glob("*.json")
                  if not p.name.startswith("contract_number_audit_"))
    cell = next(c for c in cells() if c.cell_id == "essentiality:any:essentiality")
    blob = " ".join(str(getattr(cell, f, "") or "") for f in mod.PROSE_FIELDS)
    nums = [t for t in dict.fromkeys(DECIMAL_RE.findall(blob))
            if not mod.NON_MEASUREMENT.match(t) and (cell.cell_id, t) not in ADJUDICATED_BENIGN]
    # 8, not the 7 of the original measurement: Step 4's citation placement moved a sentence-final
    # period off `AUROC 0.580`, which DECIMAL_RE's (?![\d.]) lookahead had been hiding. The figure is
    # 44/600 either way -- adding a number can only make a full match harder, and it barely moved.
    assert len(nums) == 8, nums
    full, n = mod._decoy_full_match_count(nums, set(), pool, {}, group_size=1)
    assert n >= 500
    assert full >= 20, f"only {full} decoys matched — has the pool or the matcher changed?"
    assert mod._grade(full, n, nums) != "HIGH"


def test_size_matching_removes_the_multi_citation_haystack_advantage():
    """The second defect: a cell citing N artifacts was scored against N concatenated files while each
    decoy was ONE file, so multi-citation cells earned discrimination for free (a 33-file bundle scored
    6/6 against single-file decoys). A single round number must not look distinctive just because its
    own haystack is big."""
    import scripts.contract_number_audit as mod

    pool = sorted(p.relative_to(REPO).as_posix() for p in (REPO / "wiki").glob("*.json")
                  if not p.name.startswith("contract_number_audit_"))
    nums = ["1.0"]
    f1, n1 = mod._decoy_full_match_count(nums, set(), pool, {}, group_size=1)
    f8, n8 = mod._decoy_full_match_count(nums, set(), pool, {}, group_size=8)
    assert n8 < n1, "grouping must consume the pool in N-sized blocks"
    # a bigger haystack can only ever make a full match EASIER, which is exactly the advantage being
    # neutralised: the grade for a many-citation cell is computed against equally-big decoys.
    assert f8 / n8 >= f1 / n1
    assert mod._grade(f8, n8, nums) != "HIGH", "one round number must not grade HIGH at any group size"


def test_grouping_is_non_overlapping_and_drops_a_partial_tail():
    """Every trial must be exactly group_size files, or the control compares unequal haystacks."""
    import scripts.contract_number_audit as mod
    pool = [f"wiki/p{i}.json" for i in range(10)]
    cache = {p: (f"filler {i}", [float(i)]) for i, p in enumerate(pool)}
    _, n = mod._decoy_full_match_count(["9.99"], set(), pool, cache, group_size=3)
    assert n == 3, n          # 10 // 3 == 3 trials; the trailing single file is dropped
    _, n1 = mod._decoy_full_match_count(["9.99"], set(), pool, cache, group_size=1)
    assert n1 == 10


def test_step1_leaves_verdict_and_drift_untouched():
    """A grading change must never read as drift. Baseline recorded at bcfcbc2 before any edit."""
    rep = _cached_report()
    assert rep["verdict"] == "NOTHING_TO_ADJUDICATE"
    assert rep["n_candidate_drift"] == 0
    assert rep["n_numbers_checked"] == 176
    assert rep["n_cells_audited"] == 19


# --- Step 2: provenance kinds ---------------------------------------------------------------------
# "Unverifiable" was doing the work of four statements at once, which made 34 numbers read as 34
# defects. These tests exist so the typing can never become a silent exemption mechanism.

def test_doi_prefix_is_typed_structural_only_when_it_is_actually_a_doi():
    import scripts.contract_number_audit as mod
    real = "wet-lab DpoTropiSearch (Concha-Eloko 2025; Zenodo 10.5281/zenodo.14065540) prophage-LCA"
    assert mod._structural_kind(real, "10.5281") == mod.KIND_STRUCTURAL
    # the same token WITHOUT the slash that makes it a DOI is an ordinary number
    assert mod._structural_kind("a ratio of 10.5281 over the pool", "10.5281") is None
    # and an ordinary decimal is never a DOI
    assert mod._structural_kind("top-1 ~0.45 / top-5 ~0.60", "0.45") is None


def test_kb_length_is_typed_structural_and_a_ratio_is_not():
    import scripts.contract_number_audit as mod
    assert mod._structural_kind("ARHGAP36 5.1-kb intron-1 deletion", "5.1") == mod.KIND_STRUCTURAL
    assert mod._structural_kind("STX17 4.6-kb dup grey", "4.6") == mod.KIND_STRUCTURAL
    assert mod._structural_kind("ASIP 11-bp-del black", "11.0") is None
    assert mod._structural_kind("top-1 ~0.45 over KL-types", "0.45") is None


def test_a_measurement_merely_NEAR_the_word_kb_is_not_typed_a_length():
    """NON-VACUITY IN THE DANGEROUS DIRECTION. An over-broad context match would type a genuine
    measurement as structural and drop it out of the audit entirely -- a check quietly stopping
    checking, which is the exact failure this whole script exists to catch."""
    import scripts.contract_number_audit as mod
    near = "the 5.1-kb deletion cohort reached sens 0.45 across every kb window we tried"
    assert mod._structural_kind(near, "0.45") is None, "a real measurement was exempted as a length"
    assert mod._structural_kind(near, "5.1") == mod.KIND_STRUCTURAL


def test_a_token_used_both_structurally_and_as_a_measurement_is_not_exempted():
    """All-or-nothing on purpose: if one occurrence is a real number, the token stays auditable."""
    import scripts.contract_number_audit as mod
    mixed = "a 4.6-kb duplication, and separately a concordance of 4.6 against the panel"
    assert mod._structural_kind(mixed, "4.6") is None


def test_per_kind_counts_sum_to_every_extracted_number():
    """The invariant that makes a SHRINKING unverifiable count trustworthy -- the residual got smaller
    because numbers were NAMED, not because a typing bug lost them. audit() raises if this breaks, so
    this also pins that the guard is reachable."""
    rep = _cached_report()
    assert sum(rep["n_by_kind"].values()) == rep["n_numbers_extracted"]
    assert rep["n_numbers_extracted"] == 192
    assert set(rep["n_by_kind"]) <= set(rep["kind_definitions"])


def test_the_kind_tally_reconciles_with_the_cited_and_uncited_halves():
    rep = _cached_report()
    cited = sum(r["checked"] for r in rep["cells"])
    uncited = sum(len(u["kinds"]) for u in rep["cells_unverifiable"]) \
        + sum(len(u["kinds"]) for u in rep["cells_all_numbers_typed_no_citation"])
    assert cited + uncited == rep["n_numbers_extracted"], (cited, uncited)
    assert cited == rep["n_numbers_checked"] == 176


def test_unverifiable_now_means_only_the_residual():
    """34 -> 24 in Step 2, because six numbers were NAMED: three structural (a Zenodo DOI prefix and
    two kb lengths), two enforced by an exact-equality test assert, one external reference value.
    Then 24 -> 9 in Step 4, because three cells gained the artifact that actually holds their numbers
    (kleb 3 + essentiality 7 + hla b5701 5 = 15). The nine that remain are the two vacuous-citation
    cells (cyp2d6, pigment) and pneumoserotype, whose ~2.1% no-call rate is derived."""
    rep = _cached_report()
    assert rep["n_by_kind"]["unverifiable"] == 9
    assert rep["n_numbers_unverifiable"] == 9
    assert sum(u["n_numbers"] for u in rep["cells_unverifiable"]) == 9
    for u in rep["cells_unverifiable"]:
        assert all(k == "unverifiable" for k in
                   (u["kinds"][t] for t in u["numbers"])), u["cell_id"]


def test_cells_whose_every_number_is_named_leave_the_residual_but_not_the_report():
    rep = _cached_report()
    named = {u["cell_id"] for u in rep["cells_all_numbers_typed_no_citation"]}
    assert {"pgx:human:nudt15", "typing:Escherichia_coli:pathotype",
            "typing:cat:catcolor", "typing:horse:horsecolor"} <= named
    assert named.isdisjoint({u["cell_id"] for u in rep["cells_unverifiable"]})


def test_the_doi_pattern_generalised_beyond_the_numbers_it_was_written_for():
    """Why a PATTERN and not a list. It was written for the Zenodo prefix in typing:klebsiella:kleb and
    it also caught `doi:10.5061/dryad...` in BOTH dog cells, which the hand triage had missed -- those
    were previously counted as artifact matches by coincidence. A per-number list would have shipped
    that miss, which is the hand-enumerated-exclusion trap this repo has hit five times."""
    rep = _cached_report()
    rows = {r["cell_id"]: r for r in rep["cells"]}
    for cell_id in ("typing:dog:coatcolor", "typing:dog:morphology"):
        assert rows[cell_id]["number_kinds"]["10.5061"] == "structural-non-measurement"
    assert rep["n_by_kind"]["structural-non-measurement"] == 5


def test_every_declared_kind_entry_carries_a_reason_and_the_tables_stay_small():
    """Same discipline as ADJUDICATED_BENIGN: a declaration without a reason is an unaudited exemption."""
    import scripts.contract_number_audit as mod
    for name in ("EXTERNAL_REFERENCE", "ENFORCED_BY_TEST", "SUPERSEDED_VALUE", "DECLARED_THRESHOLD"):
        table = getattr(mod, name)
        assert table, name
        assert len(table) <= 5, f"{name} is growing into a per-number exclusion list ({len(table)})"
        for key, why in table.items():
            assert isinstance(key, tuple) and len(key) == 2, key
            assert isinstance(why, str) and len(why) > 40, (name, key, why)


def test_adjudicated_benign_pin_is_untouched_by_the_new_tables():
    """The new tables are separate precisely so this pre-existing pin does not have to move."""
    assert len(ADJUDICATED_BENIGN) <= 5
    import scripts.contract_number_audit as mod
    for name in ("EXTERNAL_REFERENCE", "ENFORCED_BY_TEST", "SUPERSEDED_VALUE", "DECLARED_THRESHOLD"):
        assert not (set(getattr(mod, name)) & set(ADJUDICATED_BENIGN)), name


def test_a_declared_kind_outranks_an_artifact_match():
    """Precedence matters: a DOI prefix that happens to appear in some JSON is still not a measurement,
    and a threshold that appears in an artifact is still a bar."""
    import scripts.contract_number_audit as mod
    k = mod._classify_kind("typing:klebsiella:kleb", "Zenodo 10.5281/zenodo.1", "10.5281",
                           matched_in_artifact=True, cited=True)
    assert k == mod.KIND_STRUCTURAL
    k2 = mod._classify_kind("typing:Escherichia_coli:pathotype", "each >=0.80", "0.80",
                            matched_in_artifact=True, cited=True)
    assert k2 == mod.KIND_THRESHOLD


def test_cited_but_absent_is_candidate_drift_not_unverifiable():
    """Collapsing the two would hide real drift inside a coverage statistic."""
    import scripts.contract_number_audit as mod
    assert mod._classify_kind("some:cell", "sens 0.777 measured", "0.777",
                              matched_in_artifact=False, cited=True) == mod.KIND_DRIFT
    assert mod._classify_kind("some:cell", "sens 0.777 measured", "0.777",
                              matched_in_artifact=False, cited=False) == mod.KIND_UNVERIFIABLE


def test_step2_leaves_verdict_and_drift_untouched():
    """A typing change must never read as drift."""
    rep = _cached_report()
    assert rep["verdict"] == "NOTHING_TO_ADJUDICATE"
    assert rep["n_candidate_drift"] == 0
    assert rep["n_by_kind"].get("candidate-drift", 0) == 0
    assert rep["n_numbers_checked"] == 176


# --- Step 4: only the citations the fixed rule admitted --------------------------------------------

def test_the_three_admitted_citations_resolve_and_carry_no_drift():
    """Step 3's rule: cite only when the artifact resolves AND every number is present AND full-pool
    discrimination beats LOW. Three cells qualified; each must stay clean."""
    rep = _cached_report()
    rows = {r["cell_id"]: r for r in rep["cells"]}
    for cell_id, artifact in (
        ("typing:klebsiella:kleb", "wiki/klebsiella_topk_ksweep_2026-07-25.json"),
        ("essentiality:any:essentiality", "wiki/essentiality_e3_human_2026-07-28.json"),
        ("hla:human:b5701", "wiki/hla_validation_2026-07-06.md"),
    ):
        r = rows[cell_id]
        assert artifact in r["cited"], (cell_id, r["cited"])
        assert r["candidate_drift"] == [], (cell_id, r["candidate_drift"])
        assert r["unresolved_citations"] == [], (cell_id, r["unresolved_citations"])
        assert r["discrimination"] != "LOW", (cell_id, r["discrimination"], r["decoy_match_rate"])


def test_the_one_strong_citation_is_the_klebsiella_one():
    """kleb is the only admitted cell whose citation genuinely could have failed: 8 of 599 unrelated
    haystacks contain its number set. The other two are MODERATE and are recorded as weak."""
    rep = _cached_report()
    rows = {r["cell_id"]: r for r in rep["cells"]}
    assert rows["typing:klebsiella:kleb"]["discrimination"] == "HIGH"
    assert rows["typing:klebsiella:kleb"]["decoy_match_rate"] < 0.05


def test_the_cells_the_rule_refused_still_have_no_citation():
    """Refusals are the point. cyp2d6 and pigment are LOW-discrimination by measurement (two numbers,
    one of them 1.0, which 730+ wiki files contain); pneumoserotype's ~2.1% no-call rate is DERIVED, so
    citing its cohort artifact would manufacture an adjudication item out of a correct claim."""
    rep = _cached_report()
    refused = {u["cell_id"] for u in rep["cells_unverifiable"]}
    assert refused == {"pgx:human:cyp2d6", "typing:human:pigment",
                       "typing:Streptococcus_pneumoniae:pneumoserotype"}, refused
    cited = {r["cell_id"] for r in rep["cells"]}
    assert not (refused & cited)


def test_step4_raised_coverage_without_creating_drift():
    """The whole point: more numbers verified, nothing newly flagged."""
    rep = _cached_report()
    assert rep["n_cells_audited"] == 19          # was 16
    assert rep["n_numbers_unverifiable"] == 9    # was 24 after Step 2, 34 before Step 2
    assert rep["n_candidate_drift"] == 0
    assert rep["verdict"] == "NOTHING_TO_ADJUDICATE"


def test_the_full_pool_control_reclassified_cells_the_sample_had_flattered():
    """The point of Step 1, pinned so a silent regression to a weaker control is visible: the sample
    reported 2 LOW cells, the full pool reports 6. Three of them (ktype 51%, pointfinder 59%, slco1b1
    65%) were graded MODERATE or better by 10 draws."""
    rep = _cached_report()
    rows = {r["cell_id"]: r for r in rep["cells"]}
    assert rep["n_cells_low_discrimination"] >= 6
    for cell_id in ("typing:Klebsiella:ktype", "finder:Escherichia_coli:pointfinder",
                    "pgx:human:slco1b1"):
        assert rows[cell_id]["discrimination"] == "LOW", (cell_id, rows[cell_id])
        assert rows[cell_id]["decoy_match_rate"] >= 0.50
