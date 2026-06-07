"""Tests for the dna-amr CLI (dna_decode/amr/cli) — cached-run mode (no Docker).

Pins the --organism passthrough into provenance + the v1 record schema. Genome mode (Docker) is not
exercised here. Runnable via pytest OR standalone.
"""
import io
import json
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.amr.cli import main as amr_main

_HEADER = ("Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\t"
           "Subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t"
           "% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\t"
           "Closest reference name\tHMM accession\tHMM description")


def _run_dir(tmp, rows):
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


def test_organism_flag_recorded_in_provenance():
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), [("gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINTX"),
                                 ("parC_S80I", "QUINOLONE", "QUINOLONE", "POINTX")])
        rc, out = _invoke(["--drug", "ciprofloxacin", "--amrfinder-run", str(rd),
                           "--organism", "Klebsiella_pneumoniae", "--json-only"])
    assert rc == 0
    rec = json.loads(out)
    assert rec["provenance"]["amrfinder_organism"] == "Klebsiella_pneumoniae"
    assert rec["prediction"] == "R"   # 2 QRDR point mutations → R (cross-organism rule)


def test_organism_defaults_to_escherichia():
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), [("blaKPC-2", "BETA-LACTAM", "CARBAPENEM", "EXACTX")])
        rc, out = _invoke(["--drug", "meropenem", "--amrfinder-run", str(rd), "--json-only"])
    assert rc == 0
    rec = json.loads(out)
    assert rec["provenance"]["amrfinder_organism"] == "Escherichia"
    assert rec["prediction"] == "R"   # carbapenemase → R


if __name__ == "__main__":
    for k, v in sorted(globals().items()):
        if k.startswith("test_") and callable(v):
            v(); print(f"PASS {k}")
    print("done")
