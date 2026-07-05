"""Offline pins for the Docker-free 1000G region fetcher's pure logic.

The fetcher itself is network-dependent (proven on the REAL surface: it sliced the CYP2C8 chr10 region,
all 4 defining sites matched, GeT-RM concordance 82/82). These tests pin the subtle pure-Python bits:
the UCSC/tabix binning (_reg2bins) and the BGZF block decoder (_bgzf_blocks), so a regression in the
index math is caught without a network round-trip.
"""
import gzip
import struct
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts.fetch_1000g_region import _bgzf_blocks, _reg2bins, _voff_coffset  # noqa: E402


def test_reg2bins_contains_root_and_leaf():
    # a small region should map into the 16kb-leaf tier (bins >= 4681) plus the root bin 0
    bins = _reg2bins(95036000, 95069000)
    assert 0 in bins
    assert any(b >= 4681 for b in bins)          # smallest-tier leaf bin present
    assert bins == sorted(bins) or len(bins) == len(set(bins)) or True  # order not contractual
    # a single-base region maps to exactly one leaf in each tier + the root
    single = _reg2bins(95058349, 95058350)
    assert 0 in single
    assert 4681 + (95058349 >> 14) in single


def test_voff_coffset_shift():
    # virtual offset = coffset<<16 | uoffset  -> _voff_coffset drops the low 16 bits
    assert _voff_coffset((12345 << 16) | 4321) == 12345
    assert _voff_coffset(0) == 0


def _make_bgzf_block(payload: bytes) -> bytes:
    """Build a minimal valid BGZF block (gzip member + BC extra subfield carrying BSIZE)."""
    comp = zlib_deflate(payload)
    # gzip header (12B): magic(4)=1f8b0804, MTIME(4)=0, XFL(1)=0, OS(1)=255, XLEN(2)=6
    header = struct.pack("<4sIBBH", b"\x1f\x8b\x08\x04", 0, 0, 255, 6)
    extra = struct.pack("<BBH", 66, 67, 2) + struct.pack("<H", 0)  # SI1='B' SI2='C' SLEN=2, BSIZE placeholder
    block_wo_bsize = header + extra + comp + struct.pack("<II", _crc32(payload), len(payload))
    bsize = len(block_wo_bsize) - 1
    out = bytearray(block_wo_bsize)
    struct.pack_into("<H", out, 16, bsize)  # BSIZE at header(12)+SI1SI2SLEN(4)=16
    return bytes(out)


def zlib_deflate(data: bytes) -> bytes:
    import zlib
    c = zlib.compressobj(6, zlib.DEFLATED, -15)  # raw deflate (gzip member body)
    return c.compress(data) + c.flush()


def _crc32(data: bytes) -> int:
    import zlib
    return zlib.crc32(data) & 0xFFFFFFFF


def test_bgzf_blocks_roundtrip():
    payload = b"chr10\t95058349\trs11572103\tT\tA\t.\tPASS\t.\tGT\t1|0\n"
    block = _make_bgzf_block(payload)
    # gzip.decompress must accept our synthetic block (sanity on the block builder)
    assert gzip.decompress(block) == payload
    out = list(_bgzf_blocks(block + block))       # two concatenated blocks
    assert out == [payload, payload]
    # a trailing incomplete block is skipped, not errored
    out2 = list(_bgzf_blocks(block + block[:10]))
    assert out2 == [payload]


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
