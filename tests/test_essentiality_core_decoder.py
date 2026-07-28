"""Tests for the conserved-core essentiality decoder (v0, label-independent)."""
from dna_decode.essentiality.core_decoder import score_gene, decode_genome


def test_known_essential_core_genes_predicted():
    # canonical essential E. coli genes by function -> essential
    assert score_gene("ftsZ", "cell division protein FtsZ").prediction == "essential"
    assert score_gene("rpsA", "30S ribosomal protein S1").prediction == "essential"
    assert score_gene("ileS", "isoleucine--tRNA ligase").prediction == "essential"
    assert score_gene("rpoB", "DNA-directed RNA polymerase subunit beta").prediction == "essential"


def test_clearly_nonessential_predicted_non():
    # transposase / catabolism / hypothetical -> not essential
    assert score_gene("insH", "IS5 transposase").prediction == "non_essential"
    assert score_gene("lacZ", "beta-galactosidase (lactose catabolism)").prediction == "non_essential"
    assert score_gene("yaaX", "hypothetical protein").prediction == "non_essential"


def test_score_monotonic_and_matched():
    c = score_gene("rplA", "50S ribosomal protein L1")
    assert c.core_score >= 2.0 and c.matched


def test_decode_genome_shape():
    calls = decode_genome([("ftsZ", "cell division protein FtsZ"), ("insH", "IS5 transposase")])
    assert len(calls) == 2 and calls[0].prediction == "essential" and calls[1].prediction == "non_essential"
