"""dna-amr CLI tests for the HIV NNRTI/NRTI target-site branch (--observed wheel mode + genome-FASTA mode)."""
import json
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.amr.cli import main  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))
_REF = Path(__file__).resolve().parent.parent / "data" / "hiv_ref" / "HIV1_RT_HXB2_cds.fna"


def _planted_genome(tmp_path: Path, pos_1based: int, new_codon: str) -> Path:
    """Build a genome from the committed HXB2 RT reference with one codon mutated, wrapped in flanks."""
    seq = "".join(l for l in _REF.read_text().splitlines() if not l.startswith(">"))
    i = 3 * (pos_1based - 1)
    seq = seq[:i] + new_codon + seq[i + 3:]
    fa = tmp_path / "planted.fna"
    fa.write_text(">contig1\n" + "ACGT" * 60 + seq + "ACGT" * 60 + "\n", encoding="utf-8")
    return fa


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


def test_hiv_genome_mode_missing_fasta(capsys, tmp_path):
    rc = main(["--drug", "nevirapine", "--genome-fasta", str(tmp_path / "nope.fna")])
    assert rc == 2 and "not found" in capsys.readouterr().err


def test_hiv_genome_mode_missing_reference(capsys, tmp_path):
    fa = tmp_path / "g.fna"
    fa.write_text(">c\nACGT\n", encoding="utf-8")
    rc = main(["--drug", "nevirapine", "--genome-fasta", str(fa),
               "--hiv-rt-ref", str(tmp_path / "absent_ref.fna")])
    assert rc == 3 and "RT CDS reference" in capsys.readouterr().err


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
@pytest.mark.skipif(not _REF.exists(), reason="HXB2 RT reference fixture absent")
def test_hiv_genome_mode_calls_nnrti_K103N(capsys, tmp_path):
    fa = _planted_genome(tmp_path, 103, "AAC")   # K103N
    rc, rec = _json(capsys, ["--drug", "efavirenz", "--genome-fasta", str(fa), "--json-only"])
    assert rc == 0 and rec["prediction"] == "R"
    assert rec["provenance"]["mode"] == "blast-rt"
    assert any(d["subclass"] == "K103N" for d in rec["determinants"])


def test_hiv_amrfinder_run_rejected(capsys, tmp_path):
    rc = main(["--drug", "lamivudine", "--amrfinder-run", str(tmp_path)])
    assert rc == 2 and "bacterial-only" in capsys.readouterr().err


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
