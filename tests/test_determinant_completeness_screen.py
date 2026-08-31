"""Pin the L2 completeness screen, and the two defects found by reading its output.

Both defects produced plausible-looking output. Neither would have been caught by a passing test suite —
only by asking whether the ranking meant what it claimed.

Offline; the pure-helper tests need no cache, and the real-data regression skips when it is absent.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from determinant_completeness_screen import (determinant_key, rank_candidates,  # noqa: E402
                                             rule_counts_determinant)


def _fam(r, s, n=None):
    return {"n_genomes": n if n is not None else r + s, "r_carriers": r, "s_carriers": s}


# --- determinant identity ---

def test_key_is_symbol_class_subclass_not_the_allele():
    """rmtB and rmtE1 are different alleles of one blind family, but merging on a bare prefix would
    over-merge unrelated genes. Keep them distinct and let the ranking aggregate the evidence."""
    a = determinant_key({"Element symbol": "rmtB", "Class": "AMINOGLYCOSIDE", "Subclass": "AMINOGLYCOSIDE"})
    b = determinant_key({"Element symbol": "rmtE1", "Class": "AMINOGLYCOSIDE", "Subclass": "AMINOGLYCOSIDE"})
    assert a != b and a[1:] == b[1:]


def test_key_separates_same_symbol_under_different_subclass():
    """Subclass is what the rule matches on, so two rows sharing a symbol but differing in Subclass are
    genuinely different cases for this question."""
    a = determinant_key({"Element symbol": "blaOXA", "Class": "BETA-LACTAM", "Subclass": "BETA-LACTAM"})
    b = determinant_key({"Element symbol": "blaOXA", "Class": "BETA-LACTAM", "Subclass": "CARBAPENEM"})
    assert a != b


# --- DEFECT 1: ranking buried the actionable family under prevalent mixed ones ---

def test_pure_family_outranks_a_higher_volume_mixed_one():
    """The real case: rmtE1 (36R/0S — the KNOWN gap) ranked 5th beneath aph/aadA at 62R/28S, which are
    CORRECT exclusions. Purity is what separates a gap from a deliberate exclusion, so it must lead."""
    ranked = rank_candidates({
        ("aph(6)-Id", "AMINOGLYCOSIDE", "STREPTOMYCIN"): _fam(62, 28),
        ("rmtE1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE"): _fam(36, 0),
    })
    assert ranked[0]["symbol"] == "rmtE1"
    assert ranked[0]["signature"] == "rmt_like"
    assert ranked[1]["signature"] == "mixed"


def test_a_family_carried_by_both_classes_is_never_rmt_like():
    """One S carrier is enough to make it a probable correct exclusion, not a gap."""
    ranked = rank_candidates({("x", "C", "S"): _fam(40, 1)})
    assert ranked[0]["signature"] == "mixed"


def test_unlabelled_families_sort_last_and_get_zero_purity():
    """Absence of labels is unassessable, not innocent — it must not outrank real evidence."""
    ranked = rank_candidates({
        ("unlab", "C", "S"): _fam(0, 0, n=500),
        ("real", "C", "S"): _fam(3, 0),
    })
    assert [c["symbol"] for c in ranked] == ["real", "unlab"]
    assert ranked[1]["r_purity"] == 0.0


def test_rmt_like_needs_at_least_three_r_carriers():
    """A single pure carrier is noise, not a signature."""
    assert rank_candidates({("x", "C", "S"): _fam(1, 0)})[0]["signature"] == "mixed"
    assert rank_candidates({("x", "C", "S"): _fam(3, 0)})[0]["signature"] == "rmt_like"


# --- DEFECT 2: a one-row probe can never satisfy a threshold-2 rule ---

def test_probe_replicates_the_row_to_the_rules_threshold():
    """The defect: ciprofloxacin's rule needs TWO QRDR hits, so a one-row probe reported 0 of 51
    determinants counted and flagged every QRDR point mutation (parC_S80I at 60R/0S) as a 'gap' — when
    the rule represents them perfectly. Probing at threshold asks the representable question."""
    from dna_decode.eval.amr_rules import rule_for
    assert rule_for("ciprofloxacin")["threshold"] >= 2, "precondition: cipro is a multi-hit rule"

    # `Method` is the column the point-mutation parser keys on -- NOT `Element subtype`. The screen
    # passes real rows verbatim so it was never affected; only this synthetic fixture was wrong.
    header = ["Element symbol", "Element name", "Class", "Subclass", "% Identity to reference",
              "Element type", "Method"]
    row = {"Element symbol": "gyrA_S83L", "Element name": "GyrA", "Class": "QUINOLONE",
           "Subclass": "QUINOLONE", "% Identity to reference": "100", "Element type": "AMR",
           "Method": "POINTX"}
    # Representable: the rule counts QRDR point mutations, it just needs two of them.
    assert rule_counts_determinant(header, row, "ciprofloxacin", None) is True


def test_probe_still_rejects_a_determinant_the_rule_cannot_represent():
    """The screen must not become a rubber stamp: repeating an out-of-scope determinant to threshold
    must NOT make it count."""
    header = ["Element symbol", "Element name", "Class", "Subclass", "% Identity to reference",
              "Element type", "Method"]
    row = {"Element symbol": "rmtE1", "Element name": "16S methyltransferase",
           "Class": "AMINOGLYCOSIDE", "Subclass": "AMINOGLYCOSIDE", "% Identity to reference": "100",
           "Element type": "AMR", "Method": "EXACTX"}
    assert rule_counts_determinant(header, row, "gentamicin", None) is False


# --- real-data regression: the screen must rediscover the known gap ---

def test_screen_artifact_ranks_rmt_first_for_gentamicin():
    """End-to-end: a general screen that knows nothing about gentamicin must still surface the known
    blind family at the top."""
    hits = sorted(ROOT.glob("wiki/determinant_completeness_screen_*.json"))
    if not hits:
        pytest.skip("screen artifact not generated on this checkout")
    d = json.loads(hits[-1].read_text(encoding="utf-8"))
    gent = [x for x in d["drugs"] if x["drug"] == "gentamicin"]
    if not gent or not gent[0]["candidates"]:
        pytest.skip("gentamicin not in the committed run")
    top = gent[0]["candidates"][0]
    assert top["symbol"].lower().startswith("rmt"), f"expected an rmt family on top, got {top['symbol']}"
    assert top["s_carriers"] == 0
