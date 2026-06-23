"""CLI routing tests for the SARS-CoV-2 Mpro branch of dna-amr (mirrors test_amr_cli_hiv)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.amr.cli import main  # noqa: E402


def test_cli_observed_R_record(tmp_path, capsys):
    out = tmp_path / "rec.json"
    rc = main(["--drug", "nirmatrelvir", "--observed", "Mpro:E166V", "--sample-id", "iso1",
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["schema"] == "amr-mechanism-call-v1"
    assert rec["drug"] == "nirmatrelvir" and rec["prediction"] == "R"
    assert rec["provenance"]["organism"] == "SARS-CoV-2"        # bacterial default relabelled
    assert rec["determinants"][0]["symbol"] == "Mpro"
    assert rec["determinants"][0]["class"] == "TARGET_SITE_MUTATION"
    assert rec["caller"]["caller_is_independent_baseline"] is False


def test_cli_observed_S_exit_code(capsys):
    rc = main(["--drug", "ensitrelvir", "--observed", "Mpro:P132H"])   # benign -> S
    assert rc == 0
    assert "CALL: S" in capsys.readouterr().out


def test_cli_amrfinder_run_rejected(capsys):
    rc = main(["--drug", "nirmatrelvir", "--amrfinder-run", "somedir"])
    assert rc == 2
    assert "bacterial-only" in capsys.readouterr().err


def test_cli_requires_observed_or_genome(capsys):
    # mutually-exclusive group requires exactly one source; none -> argparse error (SystemExit 2)
    try:
        main(["--drug", "nirmatrelvir"])
    except SystemExit as e:
        assert e.code == 2
    else:
        assert False, "expected SystemExit for missing source"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
