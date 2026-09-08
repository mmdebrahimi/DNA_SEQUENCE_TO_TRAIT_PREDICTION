"""Every disclosure layer `trust_block` attaches must be able to REACH a human, not just the JSON.

THE BUG CLASS. `trust_block` attaches five namespace-separate layers (`DISCLOSURE_LAYERS`). Each exists
because some caveat is decision-relevant at call time. A layer that is attached but never rendered is a
silent half-disclosure: the machine-readable record is honest and the human output is not. The CLI's own
comment names the standard -- "a block carried only in JSON is not a disclosure" -- and two layers had
already violated it: `organism_scope` (fixed 2026-09-05) and `prospective` (fixed here).

`prospective` is the costly one to have missed. It is the project's STRONGEST tier -- an isolate
post-dating the lock is leakage-free BY CONSTRUCTION -- and both E. coli cells currently sit at
`superseded_by_surface_change` after the v2 gentamicin lock revised `amr_rules.py`. Their post-lock
numbers describe the RETIRED rule and are withheld. A reader who remembers "this cell has prospective
validation" was, at the CLI, relying on evidence for a rule that no longer ships.
"""
from __future__ import annotations

import io
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.amr.cli import main as amr_main  # noqa: E402
from dna_decode.data.trust_surface import (DISCLOSURE_LAYERS, prospective_one_line,  # noqa: E402
                                           trust_block)

_HEADER = ("Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\t"
           "Subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t"
           "% Coverage of reference\t% Identity to reference\tAlignment length\t"
           "Closest reference accession\tClosest reference name\tHMM accession\tHMM description")


def _run_dir(tmp: Path, rows) -> Path:
    cells_rows = []
    for sym, cls, sub, meth in rows:
        cells = [""] * 22
        cells[5] = sym; cells[10] = cls; cells[11] = sub; cells[12] = meth
        cells_rows.append("\t".join(cells))
    d = tmp / "run"; d.mkdir()
    (d / "main.tsv").write_text("\n".join([_HEADER, *cells_rows]) + "\n", encoding="utf-8")
    return d


def _invoke(argv):
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = amr_main(argv)
    return rc, buf.getvalue()


_RMT_ROW = [("rmtB", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE", "EXACTX")]


# --- the real surface -----------------------------------------------------------------------------

def test_prospective_status_reaches_the_printed_output():
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                           "--organism", "Escherichia_coli_Shigella"])
    assert rc == 0
    assert "prospective-lock" in out, f"prospective layer never printed:\n{out}"
    assert "SUPERSEDED" in out


def test_the_superseded_line_says_the_numbers_are_withheld_and_the_clock_restarts():
    """Naming the state without its consequence would leave a reader thinking the evidence still counts."""
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        _rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                            "--organism", "Escherichia_coli_Shigella"])
    assert "WITHHELD" in out
    assert "clock restarts" in out
    assert "amr_rules.py" in out          # names WHICH file drifted


def test_it_is_silent_for_a_cell_with_no_prospective_evidence():
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), [("blaKPC-2", "BETA-LACTAM", "CARBAPENEM", "EXACTX")])
        _rc, out = _invoke(["--drug", "meropenem", "--amrfinder-run", str(rd),
                            "--organism", "Klebsiella_pneumoniae"])
    assert "prospective-lock" not in out


def test_the_call_is_unchanged_by_the_disclosure():
    """L2 qualifies a call; it must never alter one."""
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        _rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                            "--organism", "Escherichia_coli_Shigella", "--json-only"])
    import json
    assert json.loads(out)["prediction"] == "R"


# --- the renderer's states ------------------------------------------------------------------------

def test_scored_reports_numbers_and_names_the_construction():
    line = prospective_one_line({"prospective": {
        "status": "scored", "n_scored": 61, "acc": 0.967, "sens": 0.917, "spec": 1.0,
        "powering": "POWERED"}})
    assert "61" in line and "0.967" in line
    assert "BY CONSTRUCTION" in line


def test_an_underpowered_score_is_flagged_not_quoted_flat():
    line = prospective_one_line({"prospective": {
        "status": "scored", "n_scored": 8, "acc": 0.5, "powering": "UNDERPOWERED"}})
    assert "UNDERPOWERED" in line and "indicative" in line


def test_lock_unverified_withholds_the_numbers():
    line = prospective_one_line({"prospective": {"status": "lock_unverified",
                                                 "generated": "2026-08-24"}})
    assert "WITHHELD" in line


def test_not_accrued_is_silent_and_silence_never_means_validated():
    """An armed lock with nothing accrued has no call-time consequence; printing it every call is noise."""
    assert prospective_one_line({"prospective": {"status": "not_accrued"}}) is None
    assert prospective_one_line({}) is None
    assert prospective_one_line({"prospective": None}) is None


# --- the structural guard -------------------------------------------------------------------------

_RENDERED = {"lineage", "source_concentration", "prospective", "organism_scope"}


def test_the_known_layer_set_has_not_grown_unnoticed():
    """A NEW layer added to DISCLOSURE_LAYERS without a renderer is the exact bug this file exists for.

    `doubt_layer` is deliberately NOT in `_RENDERED`: see the tripwire below.
    """
    assert set(DISCLOSURE_LAYERS) == _RENDERED | {"doubt_layer"}, (
        "DISCLOSURE_LAYERS changed -- wire a renderer for the new layer into dna_decode/amr/cli.py "
        "and add it to _RENDERED, or justify it here the way doubt_layer is justified.")


def test_doubt_layer_is_inert_today_and_this_trips_the_moment_it_is_not():
    """TRIPWIRE, not a renderer.

    `doubt_layer` (the per-cell determinant-COMPLETENESS screen -- distinct from the record's top-level
    per-call `doubt`) is attached to the badge and has no renderer. Writing one now would be dead code:
    measured across every deployed cell, NO cell carries a strong completeness signal, so the line could
    never fire. The honest move is to leave it unwritten and make its absence DETECTABLE the moment it
    starts mattering, rather than ship speculative output.

    If this fails, a cell has acquired a strong completeness signal that a CLI reader cannot see. Write
    the renderer then.
    """
    offenders = []
    for drug in ("ciprofloxacin", "gentamicin", "ceftriaxone", "tetracycline", "meropenem"):
        for org in ("Escherichia_coli_Shigella", "Klebsiella_pneumoniae", "Salmonella_enterica"):
            dl = (trust_block(drug, org) or {}).get("doubt_layer") or {}
            if dl.get("n_strong_completeness_signals"):
                offenders.append((drug, org, dl.get("strong_families")))
    assert not offenders, (
        f"doubt_layer now carries a strong completeness signal {offenders} but has NO renderer, so it "
        "reaches the JSON and never a human. Write a doubt_layer one-liner and wire it into the CLI.")
