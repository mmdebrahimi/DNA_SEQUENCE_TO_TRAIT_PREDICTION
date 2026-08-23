"""Offline tests for the genome-level forward edit path (dna_decode/forward/genome_edit) + the ESM-method
wiring in predict_effect (mock table — no torch/transformers needed)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward import (  # noqa: E402
    cds_point_edit,
    parse_hgvs_c,
    predict_effect,
    predict_genome_edit,
    translate_cds,
    translate_codon,
)


def test_translate_codon():
    assert translate_codon("ATG") == "M" and translate_codon("atg") == "M"
    assert translate_codon("TAA") == "*" and translate_codon("TGA") == "*" and translate_codon("TAG") == "*"
    assert translate_codon("GAA") == "E"
    for bad in ("AT", "ATGG", "ATX", ""):
        with pytest.raises(ValueError):
            translate_codon(bad)


def test_cds_point_edit_consequences():
    cds = "ATGAAAGTT"  # M K V
    # nt_pos 4 (codon2 start) A->C : AAA -> CAA = K->Q (missense)
    info = cds_point_edit(cds, 4, "A", "C")
    assert info["aa_pos"] == 2 and info["wt_aa"] == "K" and info["alt_aa"] == "Q" and info["within_codon"] == 0
    # nt_pos 6 (codon2 3rd base) A->G : AAA -> AAG = K->K (synonymous)
    syn = cds_point_edit(cds, 6, "A", "G")
    assert syn["wt_aa"] == "K" and syn["alt_aa"] == "K"
    # REF mismatch -> loud failure
    with pytest.raises(ValueError, match="REF mismatch"):
        cds_point_edit(cds, 4, "G", "C")
    with pytest.raises(ValueError):
        cds_point_edit(cds, 99, "A", "C")


def test_predict_genome_edit_missense_silent_nonsense():
    cds = "ATGGAAGTT"          # M E V
    prot = "MEV"
    # missense: nt_pos5 (codon2 2nd base) A->G : GAA -> GGA = E->G
    mis = predict_genome_edit(cds, 5, "A", "G", protein_seq=prot, protein="toy")
    assert mis.consequence == "missense" and mis.aa_mutation == "E2G"
    assert mis.protein_prediction is not None and mis.protein_prediction.wt == "E"
    # silent: nt_pos6 (codon2 3rd base) A->G : GAA -> GAG = E->E
    sil = predict_genome_edit(cds, 6, "A", "G", protein_seq=prot, protein="toy")
    assert sil.consequence == "silent" and sil.aa_mutation is None and sil.protein_prediction is None
    # nonsense: nt_pos4 (codon2 1st base) G->T : GAA -> TAA = E->* (stop)
    non = predict_genome_edit(cds, 4, "G", "T", protein_seq=prot, protein="toy")
    assert non.consequence == "nonsense" and non.alt_aa == "*" and non.aa_mutation == "E2*"
    assert non.protein_prediction.predicted_effect == "damaging"


def test_predict_genome_edit_double_coordinate_check():
    # translated WT AA must also match protein_seq (predict_effect re-verifies) -> a wrong protein_seq raises
    cds = "ATGGAAGTT"          # M E V
    with pytest.raises(ValueError, match="WT mismatch"):
        predict_genome_edit(cds, 5, "A", "G", protein_seq="MQV", protein="toy")  # says Q at pos2, really E


def test_predict_effect_esm2_mock_table():
    seq = "MKV"
    table = {2: {"K": -1.2, "R": -0.4, "A": -7.0}}    # ESM zero-shot log-probs at position 2
    p = predict_effect(seq, "K2R", method="esm2", esm_table=table)
    assert p.method == "esm2" and abs(p.raw_score - 0.8) < 1e-9 and p.predicted_effect == "preserved"
    pd = predict_effect(seq, "K2A", method="esm2", esm_table=table)
    assert abs(pd.raw_score + 5.8) < 1e-9 and pd.predicted_effect == "damaging"
    with pytest.raises(ValueError, match="requires esm_table"):
        predict_effect(seq, "K2R", method="esm2")


def test_committed_real_blatem_cds_translates_to_286aa():
    """The committed real blaTEM CDS (PZ538321.1) translates to the 286-aa TEM-1 protein (offline, no D:).
    Pins the real coordinate frame the genome demo depends on."""
    from dna_decode.forward import translate_codon
    fasta = Path(__file__).resolve().parent.parent / "data" / "forward_ref" / "blatem_3349172526.fna"
    if not fasta.exists():
        pytest.skip("real blaTEM CDS fixture not present")
    cds = "".join(ln.strip() for ln in fasta.read_text(encoding="utf-8").splitlines()
                  if not ln.startswith(">")).upper()
    assert len(cds) == 861                                   # 286 codons + stop
    prot = "".join(translate_codon(cds[i:i + 3]) for i in range(0, len(cds) - 2, 3)).rstrip("*")
    assert len(prot) == 286 and prot.startswith("MSIQHFRVALIPFFAAFCLPVFA") and prot.endswith("GASLIKHW")


# --------------------------------------------------------------------------------------------------
# The CLI entry point for the genome-edit path (added 2026-08-23).
#
# `predict_genome_edit` shipped complete, tested and VALIDATED (the 0.7611 genome-edit number over 1,715
# variants) but was reachable only as a library import + one demo script -- while `decode_router` already
# advertised `dna-forward ... --genome-fasta cds.fna` to users, a flag that did not exist. These tests pin
# the entry point that makes that advertisement true, and the refusals that are the whole point of it.
# --------------------------------------------------------------------------------------------------

def test_parse_hgvs_c_accepts_canonical_and_bare():
    assert parse_hgvs_c("c.205G>A") == (205, "G", "A")
    assert parse_hgvs_c("205G>A") == (205, "G", "A")          # bare form
    assert parse_hgvs_c(" c.1a>t ") == (1, "A", "T")           # whitespace + lowercase


def test_parse_hgvs_c_refuses_everything_that_is_not_a_cds_substitution():
    """Each of these has a DIFFERENT coordinate meaning; coercing any of them would be a silent
    coordinate error, which is exactly the failure this module exists to make loud."""
    for bad in ("c.205delG", "c.205_207del", "g.205G>A", "p.M69L", "M69L", "c.205G", "c.G205A",
                "c.0G>A", "c.205N>A", "", "205"):
        with pytest.raises(ValueError, match="HGVS|1-based"):
            parse_hgvs_c(bad)


def test_translate_cds_frame_stop_and_ambiguity():
    assert translate_cds("ATGAAAGTGCTGTAA") == "MKVL"          # ONE trailing stop dropped
    assert translate_cds("ATGAAA") == "MK"                     # no stop present -> unchanged
    assert translate_cds("ATGTAAAAA") == "M*K"                 # INTERNAL stop kept (surfaces loudly)
    assert translate_cds("ATGNNNAAA") == "MXK"                 # ambiguity codon -> X, does not raise
    with pytest.raises(ValueError, match="multiple of 3"):
        translate_cds("ATGAAAGTGCTGTA")                        # out of frame -> refuse, never guess


def _forward_cli(args):
    from dna_decode.forward.cli import main
    return main(args)


def test_cli_cds_path_missense_silent_nonsense(capsys):
    cds = "ATGAAAGTGCTGTAA"                                    # M K V L *
    assert _forward_cli(["--cds-seq", cds, "--mutation", "c.5A>G"]) == 0
    out = capsys.readouterr().out
    assert "MISSENSE" in out and "K2R" in out                  # AAA->AGA

    assert _forward_cli(["--cds-seq", cds, "--mutation", "c.6A>G"]) == 0
    out = capsys.readouterr().out
    assert "SILENT" in out and "synonymous" in out
    assert "predicted_effect" not in out                       # a silent edit makes NO effect prediction

    assert _forward_cli(["--cds-seq", cds, "--mutation", "c.4A>T"]) == 0
    out = capsys.readouterr().out
    assert "NONSENSE" in out and "K2*" in out and "damaging" in out


def test_cli_cds_path_refuses_ref_mismatch_and_protein_disagreement(capsys):
    cds = "ATGAAAGTGCTGTAA"                                    # translates to MKVL
    assert _forward_cli(["--cds-seq", cds, "--mutation", "c.5G>A"]) == 2      # CDS has A at 5, not G
    assert "REF mismatch" in capsys.readouterr().err

    # an INDEPENDENTLY supplied protein that disagrees with the translation is refused, not silently
    # preferred -- the double coordinate check. One that AGREES passes through.
    assert _forward_cli(["--cds-seq", cds, "--mutation", "c.5A>G", "--protein-seq", "MQVL"]) == 2
    assert "disagrees" in capsys.readouterr().err
    assert _forward_cli(["--cds-seq", cds, "--mutation", "c.5A>G", "--protein-seq", "MKVL"]) == 0
    assert "K2R" in capsys.readouterr().out


def test_cli_cds_path_on_the_real_committed_blatem_cds(capsys):
    """End-to-end on REAL data: the famous TEM-1 M69L.

    Literature `M69L` is AMBLER numbering; this precursor's catalytic Ambler-S70 is linear residue 68, so
    Ambler 69 is linear 67 -> codon 67 -> nt 199. Pinning the real coordinate keeps the demo honest: a
    naive reading of "position 69" lands on codon 69 (ACT/Thr) and would silently decode a different
    residue entirely.
    """
    fasta = Path(__file__).resolve().parent.parent / "data" / "forward_ref" / "blatem_3349172526.fna"
    if not fasta.exists():
        pytest.skip("real blaTEM CDS fixture not present")
    assert _forward_cli(["--cds-fasta", str(fasta), "--mutation", "c.199A>T", "--protein", "blaTEM"]) == 0
    out = capsys.readouterr().out
    assert "861 nt (286 aa)" in out
    assert "codon 67: ATG -> TTG" in out and "M -> L" in out and "MISSENSE" in out
    assert "M67L" in out
    # the naive-coordinate trap, pinned: codon 69 is a DIFFERENT residue (Thr), not Met
    assert _forward_cli(["--cds-fasta", str(fasta), "--mutation", "c.205A>T", "--protein", "blaTEM"]) == 0
    assert "codon 69: ACT -> TCT" in capsys.readouterr().out


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
