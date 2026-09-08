"""The organism-scope warning must reach the PRINTED output, not only the JSON.

WHY THIS EXISTS. `dna_decode/data/organism_scope.py` records that the gentamicin `symbol_rescue`
(rmt/npmA) was validated on E. coli but measured to over-call in Klebsiella. `trust_block` attached it
to the record, the CLI put the record in `rec["validation"]` -- and then printed the validation,
lineage and concentration one-liners while silently dropping this one. So a human running the CLI
against a Klebsiella genome never saw the warning that layer exists to give. The CLI's own comment
names that failure mode three lines above the bug: "A block carried only in JSON is not a disclosure."

These run through the REAL CLI print path with a synthetic AMRFinder run directory -- no Docker, no
network, no D: drive.
"""
from __future__ import annotations

import io
import json
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.amr.cli import main as amr_main  # noqa: E402
from dna_decode.data.organism_scope import one_line, overcall_for  # noqa: E402
from dna_decode.data.trust_surface import trust_block  # noqa: E402

_HEADER = ("Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\t"
           "Subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t"
           "% Coverage of reference\t% Identity to reference\tAlignment length\t"
           "Closest reference accession\tClosest reference name\tHMM accession\tHMM description")

WARN = "ORGANISM-SCOPE WARNING"


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


# An rmt carrier: AMRFinder files rmt* under the GENERIC aminoglycoside subclass, which is exactly why
# the frozen GENTAMICIN-subclass rule was blind to it and why the symbol_rescue exists.
_RMT_ROW = [("rmtB", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE", "EXACTX")]


def test_the_warning_is_printed_for_the_organism_it_was_measured_in():
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                           "--organism", "Klebsiella_pneumoniae"])
    assert rc == 0
    assert WARN in out, f"organism-scope warning missing from printed output:\n{out}"


def test_the_printed_warning_carries_both_caveats_not_a_bare_ppv():
    """A bare PPV would read as an archive-level claim. Both caveats are load-bearing prose."""
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        _rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                            "--organism", "Klebsiella_pneumoniae"])
    assert "SINGLE-SOURCE-DOMINATED" in out
    assert "CONTRADICTED by" in out


def test_it_is_silent_in_the_validated_scope():
    """E. coli is the scope the rescue WAS validated on (12/12, 146/146). Warning there would cry wolf."""
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                           "--organism", "Escherichia_coli_Shigella"])
    assert rc == 0
    assert WARN not in out


def test_it_is_silent_for_a_drug_with_no_measured_overcall():
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), [("gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINTX")])
        _rc, out = _invoke(["--drug", "ciprofloxacin", "--amrfinder-run", str(rd),
                            "--organism", "Klebsiella_pneumoniae"])
    assert WARN not in out


def test_the_json_still_carries_the_block_so_this_is_additive():
    """Printing must not have replaced the machine-readable block."""
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        _rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                            "--organism", "Klebsiella_pneumoniae", "--json-only"])
    rec = json.loads(out)
    assert rec["validation"].get("organism_scope")


def test_the_call_itself_is_unchanged_by_the_disclosure():
    """L2 qualifies a call; it must never alter one. The rescue still fires in both organisms."""
    preds = {}
    for org in ("Klebsiella_pneumoniae", "Escherichia_coli_Shigella"):
        with tempfile.TemporaryDirectory() as td:
            rd = _run_dir(Path(td), _RMT_ROW)
            _rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                                "--organism", org, "--json-only"])
        preds[org] = json.loads(out)["prediction"]
    assert preds["Klebsiella_pneumoniae"] == preds["Escherichia_coli_Shigella"] == "R"


# --- non-vacuity ----------------------------------------------------------------------------------

def test_the_underlying_record_is_real_not_a_fixture():
    """If overcall_for returned None for Klebsiella gentamicin, every assertion above would be
    vacuous and the 'silent' tests would pass for the wrong reason."""
    blk = overcall_for("gentamicin", "Klebsiella_pneumoniae")
    assert blk is not None
    assert blk["ppv_in_this_organism"] < blk["ppv_in_validated_scope"]
    assert one_line(blk)


def test_one_line_is_none_when_nothing_is_measured():
    assert one_line(None) is None
    assert one_line((trust_block("ciprofloxacin", "Klebsiella_pneumoniae") or {}).get(
        "organism_scope")) is None
