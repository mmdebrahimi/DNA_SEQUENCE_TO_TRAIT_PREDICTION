"""Offline pins for the 1000G CRAM-URL resolver (`scripts/resolve_1000g_cram.py`).

The live path is network-dependent (proven on NA12156 2026-07-07 -> ERR3239310); these pin the PURE parser:
whole-field sample match, ftp->http normalization, the two-index search order, and the not-found case.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import resolve_1000g_cram as m  # noqa: E402

# A realistic 1000G sequence.index body (tab-delimited; CRAM path + sample-name fields among others).
_IDX = (
    "#comment header line — must be skipped\n"
    "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239310/NA12156.final.cram\tmd5abc\tERR3239310\t"
    "study\tSAMPLE\tNA12156\tpop\tCEU\n"
    "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241111/HG00276.final.cram\tmd5def\tERR3241111\t"
    "study\tSAMPLE\tHG00276\tpop\tFIN\n"
)


def test_parse_returns_http_normalized_url():
    url = m.parse_cram_url_from_index(_IDX, "NA12156")
    assert url == "http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239310/NA12156.final.cram"
    assert url.startswith("http://")  # ftp:// normalized


def test_parse_whole_field_match_only():
    # a prefix collision must NOT match (NA12156 vs a hypothetical NA121560 line)
    idx = _IDX.replace("NA12156", "NA121560")
    assert m.parse_cram_url_from_index(idx, "NA12156") is None


def test_parse_absent_sample_is_none():
    assert m.parse_cram_url_from_index(_IDX, "NA99999") is None


def test_parse_skips_comment_and_blank():
    assert m.parse_cram_url_from_index("#hdr\n\n" + _IDX, "HG00276").endswith("HG00276.final.cram")


def test_resolver_searches_second_index_when_first_misses():
    calls = []

    def fake_fetch(url):
        calls.append(url)
        return _IDX if "698_related" in url else "#empty first index\n"

    url = m.resolve_1000g_cram_url("NA12156", indices=("first_2504", "698_related_idx"), fetch=fake_fetch)
    assert url and url.endswith("NA12156.final.cram")
    assert len(calls) == 2  # fell through the first (miss) into the second


def test_resolver_none_when_in_neither_index():
    assert m.resolve_1000g_cram_url("NA00000", indices=("a", "b"), fetch=lambda u: "#empty\n") is None
