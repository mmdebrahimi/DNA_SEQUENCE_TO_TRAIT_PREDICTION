"""Solver-status audit for FBA deletion sweeps — pure, synthetic frames, no solver and no feba.db."""
from __future__ import annotations

import ast
import json
from pathlib import Path

import pandas as pd
import pytest

from dna_decode.fba.solver_audit import (
    audit_deletion_frame,
    merge_audits,
    suspect_cell_set,
)

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"

# Every script that runs a cobrapy gene-deletion sweep. A new one added here without an audit is the
# regression this file exists to stop.
DELETION_SCRIPTS = (
    "fba_conditional_carbon_validate.py",
    "fba_conditional_essentiality_validate.py",
    "fba_regulatory_conditional_test.py",
    "fba_gapfill_carbon_recheck.py",
    "fba_gapfill_conditional_test.py",
)


def frame(rows):
    """Build a cobrapy-shaped deletion frame: ids is a FROZENSET, per the verified 0.31.1 contract."""
    return pd.DataFrame(
        [{"ids": frozenset([g]), "growth": growth, "status": status} for g, growth, status in rows],
        columns=["ids", "growth", "status"],
    )


def test_an_all_optimal_sweep_reports_nothing_suspect():
    a = audit_deletion_frame(frame([("b0001", 0.5, "optimal"), ("b0002", 0.0, "optimal")]), "glucose")
    assert a.n_rows == 2
    assert a.n_nonoptimal == 0 and a.n_nan_growth == 0 and a.n_suspect == 0
    assert a.statuses == {"optimal": 2}


def test_the_audit_records_WHICH_cells_failed_not_just_how_many():
    """The whole reason this module exists: a count cannot be crossed against the commit set."""
    a = audit_deletion_frame(
        frame([("b0001", 0.5, "optimal"), ("b0002", float("nan"), "infeasible")]), "acetate")
    assert a.n_nonoptimal == 1
    assert a.nonoptimal_cells == {"b0002"}          # the id, not the tally
    assert a.suspect_cells() == {"b0002"}


def test_nan_growth_on_an_OPTIMAL_solve_is_tracked_separately():
    """The case a status check alone would miss, which is why the two sets are not collapsed."""
    a = audit_deletion_frame(frame([("b0003", float("nan"), "optimal")]), "glycerol")
    assert a.n_nonoptimal == 0                      # status says fine
    assert a.nan_cells == {"b0003"}                 # growth says otherwise
    assert a.n_suspect == 1                         # an abstention arm still excludes it


def test_a_frame_without_a_status_column_degrades_instead_of_raising():
    """An older cobrapy must not crash a run; the partial audit rides into the artifact as a flag."""
    df = pd.DataFrame(
        [{"ids": frozenset(["b0001"]), "growth": float("nan")}], columns=["ids", "growth"])
    a = audit_deletion_frame(df, "glucose")
    assert a.status_available is False
    assert a.n_nonoptimal == 0                      # cannot know
    assert a.nan_cells == {"b0001"}                 # can still see this


def test_merge_is_json_serialisable_and_totals_across_conditions():
    audits = {
        "glucose": audit_deletion_frame(frame([("b0001", 0.5, "optimal")]), "glucose"),
        "acetate": audit_deletion_frame(
            frame([("b0001", float("nan"), "infeasible"), ("b0002", 0.2, "optimal")]), "acetate"),
    }
    merged = merge_audits(audits)
    assert merged["n_rows_total"] == 3
    assert merged["n_suspect_total"] == 1
    assert merged["n_conditions_with_nonoptimal"] == 1
    assert merged["suspect_cells"] == [["b0001", "acetate"]]
    assert merged["suspect_fraction"] == pytest.approx(1 / 3, abs=1e-6)
    json.dumps(merged)                              # must survive the artifact round-trip


def test_suspect_cell_set_is_keyed_by_gene_AND_condition():
    """A gene can fail in one condition and solve fine in another -- abstention is per CELL."""
    audits = {
        "glucose": audit_deletion_frame(frame([("b0001", 0.5, "optimal")]), "glucose"),
        "acetate": audit_deletion_frame(frame([("b0001", float("nan"), "infeasible")]), "acetate"),
    }
    assert suspect_cell_set(audits) == {("b0001", "acetate")}


# --- static coverage guards (Step 6) -------------------------------------------------------------
# Mirrors the repo's existing "no unmigrated consumer" regression-guard pattern: the wirings themselves
# cannot enforce a coverage rule, so the rule is a test.

@pytest.mark.parametrize("script", DELETION_SCRIPTS)
def test_every_deletion_script_audits_solver_status(script):
    """A sixth deletion script cannot ship unaudited without turning this red."""
    src = (SCRIPTS / script).read_text(encoding="utf-8")
    tree = ast.parse(src)
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    names |= {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
    assert "audit_deletion_frame" in names, (
        f"{script} runs gene deletions but never audits solver status. A non-optimal solve returns NaN "
        f"growth, which every script in this repo codes as ESSENTIAL -- so an unaudited script cannot "
        f"distinguish a failed solve from a real biological result.")


@pytest.mark.parametrize("script", DELETION_SCRIPTS)
def test_no_script_reads_growth_without_also_reading_status(script):
    """The NaN-to-essential coding must never again appear without an audit beside it."""
    src = (SCRIPTS / script).read_text(encoding="utf-8")
    if 'row["growth"]' not in src:
        pytest.skip(f"{script} does not read a deletion growth column directly")
    assert "audit_deletion_frame" in src or 'row["status"]' in src or 'row.get("status")' in src


# --- commit stratification (Step 7) ---------------------------------------------------------------

def _strata_fixture():
    from dna_decode.fba.conditional_essentiality import GeneRecord

    keys = tuple(f"c{i}" for i in range(5))
    recs = [
        GeneRecord("gConst", "gConst", {k: (k == "c0") for k in keys}, {}, True),
        GeneRecord("gOne", "gOne", {k: (k == "c1") for k in keys}, {}, True),
        GeneRecord("gTwo", "gTwo", {k: (k in ("c2", "c3")) for k in keys}, {}, True),
    ]
    calls = {
        "c0": {"gConst": False, "gOne": False, "gTwo": False},
        "c1": {"gConst": False, "gOne": True, "gTwo": False},
        "c2": {"gConst": False, "gOne": False, "gTwo": True},
        "c3": {"gConst": False, "gOne": False, "gTwo": True},
        "c4": {"gConst": False, "gOne": False, "gTwo": False},
    }
    return keys, recs, calls


def test_the_three_strata_partition_the_gene_set_exactly():
    from scripts.fba_conditional_carbon_validate import commit_strata

    keys, recs, calls = _strata_fixture()
    got = commit_strata(recs, calls, keys)
    assert sum(s["n_genes"] for s in got.values()) == len(recs)
    assert got["predicted_constant"]["n_genes"] == 1        # gConst: all-dispensable
    assert got["predicted_1_of_n"]["n_genes"] == 1          # gOne
    assert got["predicted_2plus"]["n_genes"] == 1           # gTwo


def test_a_constant_prediction_can_never_exact_match_a_two_sided_gene():
    """The arithmetic invariant the retracted bug violated: all exact matches must come from the
    non-constant strata, so a run reporting matches alongside 100%-constant is self-contradictory."""
    from scripts.fba_conditional_carbon_validate import commit_strata

    keys, recs, calls = _strata_fixture()
    got = commit_strata(recs, calls, keys)
    assert got["predicted_constant"]["n_exact_set_match"] == 0
    assert got["predicted_1_of_n"]["n_exact_set_match"] == 1
    assert got["predicted_2plus"]["n_exact_set_match"] == 1


def test_strata_report_which_genes_touch_a_suspect_solve():
    """A global suspect rate does not clear a CONCENTRATED subset -- one bad cell can make or break a
    1-of-N exact match."""
    from scripts.fba_conditional_carbon_validate import commit_strata

    keys, recs, calls = _strata_fixture()
    got = commit_strata(recs, calls, keys, suspect={("gOne", "c1")})
    assert got["predicted_1_of_n"]["n_genes_touching_a_nonoptimal_cell"] == 1
    assert got["predicted_2plus"]["n_genes_touching_a_nonoptimal_cell"] == 0


# --- what a non-optimal solve MEANS (the probe that reversed the plan's premise) ------------------

def test_a_deterministic_infeasible_resolve_is_genuine_essentiality_not_a_failure():
    """The finding that reversed this whole remediation. With ATPM lower_bound=6.86, deleting the only
    catabolic route for the sole carbon source leaves NO feasible flux distribution -- so the LP is
    genuinely infeasible rather than feasible-with-zero-growth."""
    from scripts.fba_infeasibility_probe import classify_resolve, verdict_for

    assert classify_resolve(None, "infeasible") == "still_infeasible"
    assert classify_resolve(float("nan"), "infeasible") == "still_infeasible"
    assert verdict_for({"still_infeasible": 39}) == "INFEASIBLE_IS_DETERMINISTIC_GENUINE_ESSENTIALITY"


def test_an_ordinary_lethal_knockout_is_feasible_with_zero_growth():
    """This is how essentiality normally presents, which is why 'infeasible' needed explaining."""
    from scripts.fba_infeasibility_probe import classify_resolve

    assert classify_resolve(0.0, "optimal") == "optimal_zero_growth"
    assert classify_resolve(1e-12, "optimal") == "optimal_zero_growth"


def test_a_cell_that_resolves_to_REAL_growth_would_have_been_a_spurious_call():
    """The branch that would have vindicated the abstention arm. It did not fire: 0 of 39."""
    from scripts.fba_infeasibility_probe import classify_resolve, verdict_for

    assert classify_resolve(0.42, "optimal") == "optimal_real_growth"
    assert verdict_for({"optimal_real_growth": 20, "still_infeasible": 19}) == \
        "INFEASIBLE_SOLVES_ARE_SPURIOUS"


def test_the_probe_refuses_to_guess_when_there_is_nothing_to_resolve():
    from scripts.fba_infeasibility_probe import verdict_for

    assert verdict_for({}) == "NO_SUSPECT_CELLS"


def test_every_fba_wiki_artifact_is_parseable_json():
    """A hand-edited artifact silently stopped being JSON for a day: the scope note was written with
    Python implicit string concatenation across lines, which JSON has no equivalent for. Nothing read
    the file, so nothing complained. An artifact no consumer can parse is not an artifact."""
    import json as _json

    root = Path(__file__).resolve().parent.parent
    broken = []
    for f in sorted((root / "wiki").glob("fba_*.json")):
        try:
            _json.loads(f.read_text(encoding="utf-8"))
        except _json.JSONDecodeError as e:
            broken.append(f"{f.name}: {e}")
    assert not broken, "unparseable FBA artifacts: " + "; ".join(broken)
