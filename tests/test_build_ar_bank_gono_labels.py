"""Pin the INT-based gono label logic (CDC S/I/R -> R/S; INTERMEDIATE excluded; gono drug map)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import scripts.build_ar_bank_gono_labels as g  # noqa: E402
from dna_decode.data.ar_isolate_bank import IsolateDetail  # noqa: E402


def _det(bs, calls):
    return IsolateDetail(isolate_id=1, panel_id=10, ar_number="0001", organism="Neisseria gonorrhoeae",
                         biosample=bs, mics={}, calls=calls)


def test_int_labels_R_S_and_exclude_I():
    dets = [
        _det("SAMN1", {"Ceftriaxone": "R", "Ciprofloxacin": "S"}),
        _det("SAMN2", {"Ceftriaxone": "S", "Ciprofloxacin": "I"}),   # cipro I -> excluded
        _det("SAMN3", {"Ceftriaxone": "I"}),                          # cef I -> excluded
        _det("", {"Ceftriaxone": "R"}),                               # no biosample -> skipped
    ]
    cro = g.int_labels_for_drug(dets, "ceftriaxone")
    assert cro == {"SAMN1": "R", "SAMN2": "S"}                        # SAMN3 (I) + blank excluded
    cip = g.int_labels_for_drug(dets, "ciprofloxacin")
    assert cip == {"SAMN1": "S"}                                      # SAMN2 cipro I excluded


def test_gono_drug_map_covers_six_scorable():
    assert set(g.GONO_DRUGS) == {"azithromycin", "cefixime", "ceftriaxone", "ciprofloxacin",
                                 "penicillin", "tetracycline"}
    assert "gentamicin" not in g.GONO_DRUGS                           # abstains, not label-built
    # the shown names match the AR Bank page column labels
    assert g.GONO_DRUGS["ceftriaxone"] == "Ceftriaxone"


def test_case_insensitive_int():
    dets = [_det("SAMN9", {"Penicillin": "r"}), _det("SAMN10", {"Penicillin": "s"})]
    assert g.int_labels_for_drug(dets, "penicillin") == {"SAMN9": "R", "SAMN10": "S"}
