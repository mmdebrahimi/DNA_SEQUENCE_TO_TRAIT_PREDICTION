"""Offline tests for the single-SNP -> GeT-RM actionable-allele concordance (pure parsing logic)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.pgx_single_snp_concordance import actionable_dosage, GENE_SNPS  # noqa: E402


def test_dosage_het():
    assert actionable_dosage("*1/*5", {"*5", "*15", "*17"}) == 1


def test_dosage_hom():
    assert actionable_dosage("*15/*5", {"*5", "*15", "*17"}) == 2


def test_dosage_zero():
    assert actionable_dosage("*1/*1B", {"*5", "*15", "*17"}) == 0


def test_compound_haplotype_counts_once():
    # '*28 + *60' is ONE haplotype carrying *28 -> dosage 1 (not 2)
    assert actionable_dosage("*28 + *60/*1", {"*28", "*37"}) == 1


def test_ambiguous_parenthetical_skipped():
    assert actionable_dosage("*1/(*37)", {"*28", "*37"}) is None
    assert actionable_dosage("*5 or *15/*1", {"*5"}) is None


def test_ugt1a1_28_hom():
    assert actionable_dosage("*28/*28", {"*28", "*37"}) == 2


def test_config_grounded():
    assert GENE_SNPS["slco1b1"][0]["actionable"] == {"*5", "*15", "*17"}
    assert GENE_SNPS["cyp4f2"][0]["actionable"] == {"*3"}
    assert {s["rsid"] for s in GENE_SNPS["ugt1a1"]} == {"rs887829", "rs4148323"}


if __name__ == "__main__":
    fns = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for k, fn in fns:
        fn(); print(f"PASS {k}")
    print(f"\n{len(fns)} passed")
