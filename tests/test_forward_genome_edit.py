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


# --------------------------------------------------------------------------------------------------
# Genome + GFF3 input: `--genome-fasta X.fna --annotations Y.gff3 --gene gyrA` (added 2026-08-23).
# Lifts the input from "an isolated CDS" (which a user rarely has) to a real genome + its annotation.
# --------------------------------------------------------------------------------------------------

_REFSEQ = Path(__file__).resolve().parent.parent / "data" / "cache" / "refseq" / "GCF_000005845.2"


def test_select_gene_cds_prefers_gene_symbol_over_gene_id():
    """The priority is load-bearing, not cosmetic: `gene_id` is strain-unique by construction
    (`gene-b0001`), so resolving a user's `gyrA` against it is the documented 0%-overlap trap."""
    from dna_decode.forward import select_gene_cds
    recs = [
        {"type": "CDS", "gene_symbol": "gyrA", "locus_tag": "b2231", "gene_id": "cds-A", "seqid": "c", "start": 1, "end": 9},
        {"type": "CDS", "gene_symbol": "parC", "locus_tag": "gyrA", "gene_id": "cds-B", "seqid": "c", "start": 10, "end": 18},
        {"type": "gene", "gene_symbol": "gyrA", "locus_tag": "", "gene_id": "gene-x", "seqid": "c", "start": 1, "end": 9},
    ]
    assert select_gene_cds(recs, "gyrA")["gene_id"] == "cds-A"       # symbol beats the locus_tag collision
    assert select_gene_cds(recs, "GYRA")["gene_id"] == "cds-A"       # case-insensitive
    assert select_gene_cds(recs, "cds-B")["gene_symbol"] == "parC"   # falls through to gene_id


def test_select_gene_cds_refuses_ambiguity_and_names_the_separating_field():
    """A multi-copy / alternative-product gene must REFUSE, not silently decode an arbitrary copy."""
    from dna_decode.forward import select_gene_cds
    alt = [  # one locus, two protein accessions -- locus_tag is IDENTICAL, so gene_id is the separator
        {"type": "CDS", "gene_symbol": "mrcB", "locus_tag": "b0149", "gene_id": "cds-NP_1", "seqid": "c", "start": 1, "end": 9},
        {"type": "CDS", "gene_symbol": "mrcB", "locus_tag": "b0149", "gene_id": "cds-YP_1", "seqid": "c", "start": 4, "end": 9},
    ]
    with pytest.raises(ValueError) as e:
        select_gene_cds(alt, "mrcB")
    assert "--gene <gene_id>" in str(e.value) and "cds-NP_1" in str(e.value)
    assert "locus_tag so the decoded copy" not in str(e.value)       # the old, useless advice

    joined = [  # ONE accession over two rows: a -1 programmed frameshift, not a contiguous CDS
        {"type": "CDS", "gene_symbol": "dnaX", "locus_tag": "b0470", "gene_id": "cds-YP_2", "seqid": "c", "start": 1, "end": 9},
        {"type": "CDS", "gene_symbol": "dnaX", "locus_tag": "b0470", "gene_id": "cds-YP_2", "seqid": "c", "start": 9, "end": 21},
    ]
    with pytest.raises(ValueError, match="JOINED multi-segment"):
        select_gene_cds(joined, "dnaX")


def test_select_gene_cds_unknown_gene_suggests_neighbours():
    from dna_decode.forward import select_gene_cds
    recs = [{"type": "CDS", "gene_symbol": s, "locus_tag": "", "gene_id": f"cds-{s}",
             "seqid": "c", "start": 1, "end": 9} for s in ("gyrA", "gyrB", "parC")]
    with pytest.raises(ValueError, match="Did you mean: gyrA, gyrB"):
        select_gene_cds(recs, "gyrZ")


def test_cds_record_key_mirrors_the_extractor_rule():
    """If this drifts from annotations.extract_cds_sequences the lookup silently misses."""
    from dna_decode.forward import cds_record_key
    assert cds_record_key({"gene_id": "g", "locus_tag": "l", "seqid": "c", "start": 1, "end": 2}) == "g"
    assert cds_record_key({"gene_id": "", "locus_tag": "l", "seqid": "c", "start": 1, "end": 2}) == "l"
    assert cds_record_key({"gene_id": "", "locus_tag": "", "seqid": "c", "start": 1, "end": 2}) == "c:1-2"


@pytest.mark.parametrize("gene,hgvs,codon,call", [
    ("gyrA", "c.248C>T", "codon 83: TCG -> TTG", "S83L"),   # THE cipro QRDR mutation
    ("parC", "c.239G>T", "codon 80: AGC -> ATC", "S80I"),   # the second QRDR gene
])
def test_genome_plus_gff_decodes_the_real_qrdr_mutations(gene, hgvs, codon, call, capsys):
    """END-TO-END on the REAL E. coli K-12 MG1655 reference genome + its RefSeq GFF3.

    Both genes are on the MINUS strand, so this also proves the reverse-complement is applied: without
    it, codon 83 is not TCG and the REF check would fail outright. These are the two mutations this
    project's whole cipro arc is built on, so a silent coordinate error here would be maximally costly.
    """
    if not (_REFSEQ / "genome.fna").exists() or not (_REFSEQ / "annotations.gff3").exists():
        pytest.skip("MG1655 reference fixture not present")
    from dna_decode.forward.cli import main
    rc = main(["--genome-fasta", str(_REFSEQ / "genome.fna"),
               "--annotations", str(_REFSEQ / "annotations.gff3"),
               "--gene", gene, "--mutation", hgvs])
    out = capsys.readouterr().out
    assert rc == 0
    assert "strand -" in out                       # minus strand => revcomp actually exercised
    assert codon in out and call in out and "MISSENSE" in out


# --------------------------------------------------------------------------------------------------
# GENOMIC coordinates (VCF-style): --genomic-pos / --ref / --alt (added 2026-08-23).
# A user holding a VCF has genome coordinates, not `c.` ones. The strand math is the whole risk.
# --------------------------------------------------------------------------------------------------

def test_complement_base_refuses_ambiguity_codes():
    from dna_decode.forward import complement_base
    assert (complement_base("a"), complement_base("G")) == ("T", "C")
    for bad in ("N", "R", "", "AT"):
        with pytest.raises(ValueError):
            complement_base(bad)


def test_genomic_to_cds_edit_plus_and_minus_strand():
    """On a MINUS-strand gene BOTH the coordinate and the bases flip; on a plus strand neither does."""
    from dna_decode.forward import genomic_to_cds_edit
    plus = {"start": 337, "end": 2799, "strand": "+"}
    assert genomic_to_cds_edit(plus, 340, "C", "T") == (4, "C", "T")
    assert genomic_to_cds_edit(plus, 337, "A", "G") == (1, "A", "G")

    minus = {"start": 2336793, "end": 2339420, "strand": "-"}          # real MG1655 gyrA
    assert genomic_to_cds_edit(minus, 2339173, "G", "A") == (248, "C", "T")
    assert genomic_to_cds_edit(minus, 2339420, "T", "C") == (1, "A", "G")   # first CDS base = last genomic

    with pytest.raises(ValueError, match="outside"):
        genomic_to_cds_edit(plus, 1, "A", "T")


def test_cds_at_genomic_position_refuses_none_and_overlaps():
    from dna_decode.forward import cds_at_genomic_position
    recs = [
        {"type": "CDS", "seqid": "c1", "start": 10, "end": 20, "strand": "+", "gene_symbol": "a"},
        {"type": "CDS", "seqid": "c1", "start": 18, "end": 30, "strand": "-", "gene_symbol": "b"},
        {"type": "gene", "seqid": "c1", "start": 1, "end": 100, "strand": "+", "gene_symbol": "g"},
    ]
    assert cds_at_genomic_position(recs, 12)["gene_symbol"] == "a"
    with pytest.raises(ValueError, match="no CDS feature covers"):
        cds_at_genomic_position(recs, 5)                 # intergenic -> no codon to decode
    with pytest.raises(ValueError, match="2 CDS features cover"):
        cds_at_genomic_position(recs, 19)                # overlapping reading frames -> refuse, don't pick


def test_genomic_coordinate_decode_on_the_real_reference(capsys):
    """END-TO-END from a bare VCF-style coordinate: no gene named, no `c.` coordinate computed by hand.

    gyrA is MINUS strand, so the genome carries G where the CDS carries C -- if the complement were
    skipped the REF check would reject it, which is exactly why this is the decisive test.
    """
    if not (_REFSEQ / "genome.fna").exists():
        pytest.skip("MG1655 reference fixture not present")
    from dna_decode.forward.cli import main
    base = ["--genome-fasta", str(_REFSEQ / "genome.fna"),
            "--annotations", str(_REFSEQ / "annotations.gff3")]

    assert main([*base, "--genomic-pos", "2339173", "--ref", "G", "--alt", "A"]) == 0
    out = capsys.readouterr().out
    assert "gyrA" in out                                  # gene auto-identified from the coordinate alone
    assert "c.248C>T (reverse-complemented onto the minus strand)" in out
    assert "S83L" in out and "MISSENSE" in out

    # plus strand: no complement, and the conversion is a plain offset
    assert main([*base, "--genomic-pos", "340", "--ref", "C", "--alt", "T"]) == 0
    out = capsys.readouterr().out
    assert "thrA" in out and "-> c.4C>T" in out and "reverse-complemented" not in out
    assert "R2*" in out and "NONSENSE" in out

    # a REF the genome does not carry is refused, not silently decoded
    assert main([*base, "--genomic-pos", "2339173", "--ref", "A", "--alt", "T"]) == 2
    assert "REF mismatch" in capsys.readouterr().err


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
