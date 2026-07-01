"""Pin the ABO O-status decoder (M4 serological cell): rs8176719 deletion -> O; label binner."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data.abo_blood import INDETERMINATE, bin_blood_type, call_abo_o_status  # noqa: E402


def test_o_status_from_rs8176719():
    assert call_abo_o_status("DD") == "O"        # homozygous deletion
    assert call_abo_o_status("DI") == "non-O"    # one functional allele
    assert call_abo_o_status("II") == "non-O"    # two functional
    assert call_abo_o_status("ID") == "non-O"
    assert call_abo_o_status("--") == INDETERMINATE   # no-call, NOT guessed as O
    assert call_abo_o_status("") == INDETERMINATE
    assert call_abo_o_status("D") == INDETERMINATE    # single allele


def test_o_status_tolerates_acgt_indel_coding():
    assert call_abo_o_status("GG") == "non-O"    # functional (insertion) present
    assert call_abo_o_status("-G") == "non-O"


def test_bin_blood_type():
    assert bin_blood_type("O +") == "O"
    assert bin_blood_type("O -") == "O"
    assert bin_blood_type("A +") == "non-O"
    assert bin_blood_type("B -") == "non-O"
    assert bin_blood_type("AB +") == "non-O"
    assert bin_blood_type("Don't know") is None
    assert bin_blood_type("") is None
