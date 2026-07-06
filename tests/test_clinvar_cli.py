"""Pins the first-class `dna-clinvar` Mendelian CLI (`dna_decode.clinvar.cli`)."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.clinvar.cli import main  # noqa: E402
from dna_decode.clinvar.decode import carried_alts, decode_vcf  # noqa: E402
from dna_decode.data.clinvar import ClinVarDecoder  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")


def _panel_tsv(tmp_path):
    p = tmp_path / "panel.tsv"
    p.write_text("chrom\tpos\tref\talt\tgene\tsignificance\treview_status\tdisease\tclinvar_id\n"
                 "7\t100\tA\tT\tCFTR\tPathogenic\treviewed_by_expert_panel\tCystic fibrosis\t1\n"
                 "11\t200\tC\tG\tHBB\tBenign\tcriteria_provided,_single_submitter\tbeta-thal\t2\n",
                 encoding="utf-8")
    return p


def test_carried_alts_reexport():
    assert carried_alts("A", "T", "0/1") == ["T"]
    assert carried_alts("A", "T,G", "1/2") == ["T", "G"]


def test_decode_vcf_pathogenic(tmp_path):
    vcf = tmp_path / "s.vcf"
    vcf.write_text(_HEADER + "7\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
    rep = decode_vcf(vcf, ClinVarDecoder.from_tsv(_panel_tsv(tmp_path)), sample="S1")
    assert rep["n_pathogenic"] == 1 and rep["pathogenic_hits"][0]["gene"] == "CFTR"


def test_cli_runs_on_vcf(tmp_path, capsys):
    vcf = tmp_path / "s.vcf"
    vcf.write_text(_HEADER + "11\t200\t.\tC\tG\t.\tPASS\t.\tGT\t1/1\n", encoding="utf-8")
    rc = main([str(vcf), "--sample-id", "S1", "--panel", str(_panel_tsv(tmp_path))])
    assert rc == 0
    out = capsys.readouterr().out
    assert "benign carrier" in out and "NOT a clinical tool" in out


def test_cli_missing_vcf_exits_2(tmp_path, capsys):
    rc = main([str(tmp_path / "nope.vcf"), "--panel", str(_panel_tsv(tmp_path))])
    assert rc == 2


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
