"""dna-amr CLI tests for the HIV NNRTI/NRTI target-site branch (wheel-only --observed mode)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.amr.cli import main  # noqa: E402


def _json(capsys, argv):
    rc = main(argv)
    out = capsys.readouterr().out
    return rc, json.loads(out)


def test_hiv_nnrti_resistant_uniform_record(capsys):
    rc, rec = _json(capsys, ["--drug", "efavirenz", "--observed", "RT:K103N", "--json-only"])
    assert rc == 0
    assert rec["schema"] == "amr-mechanism-call-v1"        # same record as every other kingdom
    assert rec["drug"] == "efavirenz" and rec["prediction"] == "R"
    assert rec["provenance"]["organism"] == "HIV-1"
    assert any(d["subclass"] == "K103N" for d in rec["determinants"])
    assert rec["caller"]["caller_is_independent_baseline"] is False


def test_hiv_nrti_resistant(capsys):
    rc, rec = _json(capsys, ["--drug", "zidovudine", "--observed", "RT:T215Y", "--json-only"])
    assert rc == 0 and rec["prediction"] == "R"
    assert rec["caller"]["name"] == "dna_decode-hiv-nrti-major-position-v0"


def test_hiv_susceptible_no_drm(capsys):
    rc, rec = _json(capsys, ["--drug", "efavirenz", "--observed", "RT:A98G", "--json-only"])
    assert rec["prediction"] == "S"
    assert rec["undetectable_mechanisms"]  # an S call surfaces the blind spots


def test_hiv_genome_mode_deferred(capsys, tmp_path):
    fa = tmp_path / "g.fna"
    fa.write_text(">c\nACGT\n", encoding="utf-8")
    rc = main(["--drug", "nevirapine", "--genome-fasta", str(fa)])
    assert rc == 2  # genome-FASTA mode is v0.1 (needs the HXB2 RT reference) -> errors, not a fake call
    assert "v0.1" in capsys.readouterr().err


def test_hiv_amrfinder_run_rejected(capsys, tmp_path):
    rc = main(["--drug", "lamivudine", "--amrfinder-run", str(tmp_path)])
    assert rc == 2 and "bacterial-only" in capsys.readouterr().err


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
