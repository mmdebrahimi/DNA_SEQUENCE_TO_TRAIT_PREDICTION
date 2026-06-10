"""dna-amr antiviral branch — productizes the influenza NA inhibitor decoder into the shipped CLI (4th kingdom).

Catalog unit tests + pure --observed mode (no BLAST, wheel-only) + real-BLAST --genome-fasta mode
(skipif no BLAST/fixture) + routing guards. Mirrors test_amr_cli_antimalarial.py for the viral kingdom.
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.amr.cli import main  # noqa: E402
from dna_decode.data.antiviral_amr import (  # noqa: E402
    call_from_observed_substitutions,
    gene_for_drug,
    is_resistance_mutation,
    resistance_mutations_for,
    supported_antiviral_drugs,
)

_NA_REF = Path(__file__).resolve().parent.parent / "data" / "antiviral_ref" / "N1_NA_NC026434_cds.fna"
_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))


# --- catalog unit tests --------------------------------------------------------------------------

def test_catalog_has_three_nai_drugs():
    assert set(supported_antiviral_drugs()) == {"oseltamivir", "peramivir", "zanamivir"}
    assert gene_for_drug("oseltamivir") == "NA"


def test_h275y_is_oseltamivir_and_peramivir_not_zanamivir():
    assert is_resistance_mutation("oseltamivir", "NA", "H275Y")
    assert is_resistance_mutation("peramivir", "NA", "H275Y")
    assert not is_resistance_mutation("zanamivir", "NA", "H275Y")   # H275Y barely affects zanamivir


def test_call_from_observed_resistant_surfaces_no_blindspots():
    call = call_from_observed_substitutions("oseltamivir", {"NA": {"H275Y"}})
    assert call.prediction == "R" and call.determinants == ["NA:H275Y"]
    assert call.undetectable_mechanisms == []


def test_call_from_observed_susceptible_surfaces_blindspots():
    call = call_from_observed_substitutions("oseltamivir", {"NA": set()})
    assert call.prediction == "S"
    assert "PA_PB2_baloxavir_class_resistance" in call.undetectable_mechanisms


def test_unknown_drug_raises():
    with pytest.raises(KeyError):
        resistance_mutations_for("baloxavir")          # cap-endonuclease inhibitor, not an NAI — not catalogued


# --- CLI --observed (wheel-only) -----------------------------------------------------------------

def test_antiviral_observed_resistant(tmp_path):
    out = tmp_path / "r.json"
    rc = main(["--drug", "oseltamivir", "--observed", "NA:H275Y", "--sample-id", "iso1",
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "R" and rec["drug"] == "oseltamivir"
    assert rec["schema"] == "amr-mechanism-call-v1"            # uniform with bacterial + fungal + antimalarial
    assert rec["determinants"][0]["symbol"] == "NA"
    assert rec["caller"]["name"] == "dna_decode-antiviral-na-target-mutation-v0"
    assert rec["provenance"]["organism"] == "Influenza_A_virus"   # bacterial default relabeled


def test_antiviral_observed_susceptible_surfaces_blindspots(tmp_path):
    out = tmp_path / "s.json"
    rc = main(["--drug", "oseltamivir", "--observed", "", "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "S"
    assert rec["undetectable_mechanisms"]


def test_antiviral_rejects_amrfinder_run(capsys):
    rc = main(["--drug", "oseltamivir", "--amrfinder-run", "/nonexistent"])
    assert rc == 2
    assert "bacterial-only" in capsys.readouterr().err


def test_zanamivir_observed_h275y_is_susceptible(tmp_path):
    """H275Y reduces oseltamivir but NOT zanamivir — the per-drug catalog must split them."""
    out = tmp_path / "z.json"
    rc = main(["--drug", "zanamivir", "--observed", "NA:H275Y", "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "S"


# --- real-BLAST genome mode ----------------------------------------------------------------------

@pytest.mark.skipif(not (_HAS_BLAST and _NA_REF.exists()), reason="BLAST+ or NA reference absent")
def test_antiviral_genome_mode_real_blast(tmp_path):
    """Real makeblastdb+blastn through the CLI: a 'genome' = real N1 NA ref with a planted H275Y -> R;
    the unmodified WT reference -> S. Codon 275 nt span = 3*275-2 .. 3*275 (1-based, WT His=CAC)."""
    ref_seq = "".join(l.strip() for l in _NA_REF.read_text().splitlines() if not l.startswith(">")).upper()
    pos = 3 * 275 - 3
    assert ref_seq[pos:pos + 3] == "CAC"                      # WT His at 275 (sanity)
    r_seq = ref_seq[:pos] + "TAC" + ref_seq[pos + 3:]         # H275Y (His->Tyr) on the real reference
    flank = "ACGT" * 60
    g = tmp_path / "R.fna"
    g.write_text(f">segment6\n{flank}{r_seq}{flank}\n", encoding="utf-8")
    out = tmp_path / "g.json"
    rc = main(["--drug", "oseltamivir", "--genome-fasta", str(g), "--sample-id", "isoR",
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["prediction"] == "R"
    assert any(d["symbol"] == "NA" and d["subclass"] == "H275Y" for d in rec["determinants"])
    assert rec["provenance"]["mode"] == "blast-na"
    # WT reference as a 'genome' -> S
    wt = tmp_path / "S.fna"
    wt.write_text(f">segment6\n{flank}{ref_seq}{flank}\n", encoding="utf-8")
    rc2 = main(["--drug", "oseltamivir", "--genome-fasta", str(wt), "--json-only"])
    assert rc2 == 0


_REAL_R = next(iter((Path(__file__).resolve().parent.parent / "data" / "antiviral_ref").glob("Flu_N1_NA_*H275Y.fna")), None)
_REAL_S = next(iter((Path(__file__).resolve().parent.parent / "data" / "antiviral_ref").glob("Flu_N1_NA_*WT.fna")), None)


@pytest.mark.skipif(not (_HAS_BLAST and _REAL_R and _REAL_S), reason="BLAST+ or real NA field fixtures absent")
def test_antiviral_real_field_isolates_R_S_and_drug_specificity(tmp_path):
    """REAL 2008-2009 N1 NA field isolates (committed fixtures): a real H275Y carrier -> oseltamivir R,
    a real WT -> S, and the SAME H275Y isolate -> zanamivir S (H275Y barely affects zanamivir). The last
    is a genuine biological discriminator, not plumbing — the per-drug catalog must split the two NAIs."""
    rc_r = main(["--drug", "oseltamivir", "--genome-fasta", str(_REAL_R), "--out", str(tmp_path / "r.json"), "--json-only"])
    assert rc_r == 0 and json.loads((tmp_path / "r.json").read_text())["prediction"] == "R"
    rc_s = main(["--drug", "oseltamivir", "--genome-fasta", str(_REAL_S), "--out", str(tmp_path / "s.json"), "--json-only"])
    assert rc_s == 0 and json.loads((tmp_path / "s.json").read_text())["prediction"] == "S"
    rc_z = main(["--drug", "zanamivir", "--genome-fasta", str(_REAL_R), "--out", str(tmp_path / "z.json"), "--json-only"])
    assert rc_z == 0 and json.loads((tmp_path / "z.json").read_text())["prediction"] == "S"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
