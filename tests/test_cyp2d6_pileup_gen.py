"""Offline pins for the CYP2D6/CYP2D7 pileup generator (`scripts/cyp2d6_pileup_gen.py`).

The live path is Docker+network (proven on NA12156 2026-07-07: 117/117 PSV coords covered). These pin the
PURE reformatting + span logic that the `_read_pileup` contract depends on.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import cyp2d6_pileup_gen as g  # noqa: E402


def test_mpileup_reformat_to_pos_ref_bases():
    # samtools mpileup: chrom pos ref depth bases quals -> pos<TAB>ref<TAB>bases
    raw = "chr22\t42128945\tC\t44\tCcCtTtT\tIIIIIII\nchr22\t42128946\tA\t40\t....,,,\tIIIIII\n"
    out = g.mpileup_to_pos_lines(raw)
    lines = out.strip().split("\n")
    assert lines[0] == "42128945\tC\tCcCtTtT"        # field0=pos digit, field2=bases (the _read_pileup contract)
    assert lines[1] == "42128946\tA\t....,,,"


def test_mpileup_reformat_skips_malformed():
    raw = "garbage line\nchr22\tNOTNUM\tC\t1\tA\tI\nchr22\t100\tG\t2\tAA\tII\n"
    out = g.mpileup_to_pos_lines(raw)
    assert out.strip() == "100\tG\tAA"  # only the well-formed digit-pos line survives


def test_mpileup_reformat_empty():
    assert g.mpileup_to_pos_lines("") == ""


def test_psv_spans_pad_and_split_paralogs():
    psvs = [
        {"chrom": "chr22", "pos_d6": 42123196, "pos_d7": 42135348},
        {"chrom": "chr22", "pos_d6": 42132017, "pos_d7": 42145724},
    ]
    d6, d7 = g.psv_spans(psvs)
    assert d6 == ("chr22", 42123196 - g._PAD, 42132017 + g._PAD)
    assert d7 == ("chr22", 42135348 - g._PAD, 42145724 + g._PAD)
    # the two paralog windows must not overlap (d6 body vs d7 body are distinct genomic regions)
    assert d6[2] < d7[1]


def test_psv_spans_uses_committed_catalog():
    """Round-trip on the real 117-PSV catalog: spans are chr22, D6 below D7, both non-empty."""
    from cyp2d6_psv_evidence import load_psvs
    d6, d7 = g.psv_spans(load_psvs())
    assert d6[0] == "chr22" and d7[0] == "chr22"
    assert d6[1] < d6[2] < d7[1] < d7[2]
