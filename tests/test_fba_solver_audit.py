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
