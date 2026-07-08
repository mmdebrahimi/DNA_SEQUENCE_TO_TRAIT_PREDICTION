"""Full ABO O/A/B/AB decoder — sourced 3-variant rule (rs8176719 + rs8176746/47)."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.data.abo_full import call_abo_full  # noqa: E402


def test_type_O_homozygous_deletion():
    # DD -> O regardless of tag background (O carries A-type tags)
    assert call_abo_full("DD", "GG", "CC") == "O"


def test_type_A_homozygous_A_tags():
    assert call_abo_full("II", "GG", "CC") == "A"     # AA
    assert call_abo_full("DI", "GG", "CC") == "A"     # AO


def test_type_B_homozygous_B_tags():
    assert call_abo_full("II", "TT", "GG") == "B"     # BB
    assert call_abo_full("DI", "TT", "GG") == "B"     # BO


def test_type_AB_het_tags_two_functional():
    assert call_abo_full("II", "GT", "CG") == "AB"    # one A + one B allele


def test_di_het_tags_is_B_not_AB():
    # one deletion (O, A-background tags) + one functional B -> het tags but phenotype B
    assert call_abo_full("DI", "GT", "CG") == "B"


def test_indeterminate_on_bad_or_inconsistent_input():
    assert call_abo_full(None, "GG", "CC") == "Indeterminate"       # no O-call
    assert call_abo_full("II", "GG", "GG") == "Indeterminate"       # tags disagree (A vs B)
    assert call_abo_full("II", "GG", None) == "Indeterminate"       # missing tag
    assert call_abo_full("II", "AA", "TT") == "Indeterminate"       # non-A/B alleles
