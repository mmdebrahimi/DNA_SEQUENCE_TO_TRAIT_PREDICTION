"""dna-amr fungal branch — productizes the G1-validated C. auris decoder into the shipped CLI.

Pure --observed mode (no BLAST, wheel-only) + real-BLAST --genome-fasta mode (skipif no BLAST/fixtures) +
guards that the fungal/bacterial engines don't cross wires. The bacterial path is covered elsewhere; here we
only assert the fungal routing + the uniform amr-mechanism-call-v1 record shape.
"""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.amr.cli import _parse_observed, main  # noqa: E402

_FREF = Path(__file__).resolve().parent.parent / "data" / "fungal_ref"
_REF = _FREF / "Cauris_ERG11_cds.fna"
_Y132F = _FREF / "Cauris_ERG11_PV630305_Y132F.fna"
_WT = _FREF / "Cauris_ERG11_PV630306_WT.fna"
_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))
_HAS_FIX = _REF.exists() and _Y132F.exists() and _WT.exists()


def test_parse_observed():
    assert _parse_observed("ERG11:Y132F,ERG11:K143R,FKS1:S639F") == {
        "ERG11": {"Y132F", "K143R"}, "FKS1": {"S639F"}}
    with pytest.raises(ValueError):
        _parse_observed("ERG11Y132F")          # no GENE:SUB separator


def test_fungal_observed_resistant(tmp_path, capsys):
    out = tmp_path / "r.json"
    rc = main(["--drug", "fluconazole", "--observed", "ERG11:Y132F", "--sample-id", "iso1",
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "R" and rec["drug"] == "fluconazole"
    assert rec["schema"] == "amr-mechanism-call-v1"           # uniform with bacterial
    assert rec["determinants"][0]["symbol"] == "ERG11"
    assert rec["caller"]["name"] == "dna_decode-fungal-target-mutation-v0"


def test_fungal_observed_susceptible_surfaces_blindspots(tmp_path):
    out = tmp_path / "s.json"
    rc = main(["--drug", "fluconazole", "--observed", "", "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "S"
    assert rec["undetectable_mechanisms"]                      # efflux/aneuploidy blind spots surfaced


def test_bacterial_rejects_observed(capsys):
    rc = main(["--drug", "ciprofloxacin", "--observed", "gyrA:S83L"])
    assert rc == 2
    assert "fungal-only" in capsys.readouterr().err


def test_fungal_rejects_amrfinder_run(capsys):
    rc = main(["--drug", "fluconazole", "--amrfinder-run", "/nonexistent"])
    assert rc == 2
    assert "bacterial-only" in capsys.readouterr().err


@pytest.mark.skipif(not (_HAS_BLAST and _HAS_FIX), reason="BLAST+ or fungal_ref fixtures absent")
def test_fungal_genome_mode_real_blast(tmp_path):
    out = tmp_path / "g.json"
    rc = main(["--drug", "fluconazole", "--genome-fasta", str(_Y132F), "--erg11-ref", str(_REF),
               "--sample-id", "PV630305", "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "R"
    assert any(d["symbol"] == "ERG11" for d in rec["determinants"])
    assert rec["provenance"]["mode"] == "blast-erg11"
    # wild-type genome -> S
    rc2 = main(["--drug", "fluconazole", "--genome-fasta", str(_WT), "--erg11-ref", str(_REF),
                "--json-only"])
    assert rc2 == 0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
