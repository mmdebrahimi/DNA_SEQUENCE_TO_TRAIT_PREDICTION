"""dna-amr antimalarial branch — productizes the P. falciparum K13 decoder into the shipped CLI.

Pure --observed mode (no BLAST, wheel-only) + real-BLAST --genome-fasta mode (skipif no BLAST/fixture) +
routing guards. Mirrors test_amr_cli_fungal.py for the 3rd kingdom (protozoan).
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.amr.cli import main  # noqa: E402

_K13_REF = Path(__file__).resolve().parent.parent / "data" / "antimalarial_ref" / "Pf3D7_K13_cds.fna"
_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))


def test_antimalarial_observed_resistant(tmp_path):
    out = tmp_path / "r.json"
    rc = main(["--drug", "artemisinin", "--observed", "K13:C580Y", "--sample-id", "iso1",
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "R" and rec["drug"] == "artemisinin"
    assert rec["schema"] == "amr-mechanism-call-v1"            # uniform with bacterial + fungal
    assert rec["determinants"][0]["symbol"] == "K13"
    assert rec["caller"]["name"] == "dna_decode-antimalarial-k13-target-mutation-v0"
    assert rec["provenance"]["organism"] == "Plasmodium_falciparum"   # bacterial default relabeled


def test_antimalarial_observed_susceptible_surfaces_blindspots(tmp_path):
    out = tmp_path / "s.json"
    rc = main(["--drug", "artemisinin", "--observed", "", "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "S"
    assert rec["undetectable_mechanisms"]                       # non-K13 / partner-drug blind spots surfaced


def test_antimalarial_rejects_amrfinder_run(capsys):
    rc = main(["--drug", "artemisinin", "--amrfinder-run", "/nonexistent"])
    assert rc == 2
    assert "bacterial-only" in capsys.readouterr().err


def test_chloroquine_observed_resistant(tmp_path):
    out = tmp_path / "cq.json"
    rc = main(["--drug", "chloroquine", "--observed", "pfcrt:K76T", "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "R" and rec["determinants"][0]["symbol"] == "pfcrt"
    assert rec["caller"]["name"] == "dna_decode-antimalarial-k13-target-mutation-v0"


def test_chloroquine_genome_mode_deferred_intron_guard(tmp_path):
    g = tmp_path / "g.fna"; g.write_text(">c\nACGTACGTACGT\n", encoding="utf-8")
    rc = main(["--drug", "chloroquine", "--genome-fasta", str(g)])
    assert rc == 3          # pfcrt is intron-containing -> genome mode deferred (no footgun)


@pytest.mark.skipif(not (_HAS_BLAST and _K13_REF.exists()), reason="BLAST+ or K13 reference absent")
def test_antimalarial_genome_mode_real_blast(tmp_path):
    """Real makeblastdb+blastn through the CLI: a genome = real 3D7 K13 ref with a planted C580Y -> R."""
    ref_seq = "".join(l.strip() for l in _K13_REF.read_text().splitlines() if not l.startswith(">")).upper()
    pos = 3 * 580 - 3
    r_seq = ref_seq[:pos] + "TAC" + ref_seq[pos + 3:]          # C580Y on the real reference
    flank = "ACGT" * 60
    g = tmp_path / "R.fna"
    g.write_text(f">contig1\n{flank}{r_seq}{flank}\n", encoding="utf-8")
    out = tmp_path / "g.json"
    rc = main(["--drug", "artemisinin", "--genome-fasta", str(g), "--sample-id", "isoR",
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "R"
    assert any(d["symbol"] == "K13" for d in rec["determinants"])
    assert rec["provenance"]["mode"] == "blast-k13"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
