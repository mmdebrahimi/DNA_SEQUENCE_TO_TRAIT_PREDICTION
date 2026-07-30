"""Spec round-trip tests for the pure-Python PLINK1 reader (dna_decode/pigment/plink_io.py).

PLINK1 .bed is a FIXED public binary spec, so we verify the decoder against a HAND-CONSTRUCTED .bed whose
bytes are computed directly from the spec (magic 0x6c1b01, SNP-major, 2 bits/sample LSB-first, codes
00=a1a1/01=missing/10=het/11=a2a2). No external data needed. Runnable via pytest OR standalone.
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.pigment.plink_io import (  # noqa: E402
    PlinkFormatError,
    find_variants,
    genotype_string,
    read_bed_variants,
    read_bim,
    read_fam,
)

# 5 samples, 2 variants. a2-dosages:
#   variant0 = [0,1,2,None,1]  -> codes [00,10,11,01,10]
#   variant1 = [2,2,0,1,None]  -> codes [11,11,00,10,01]
# Pack LSB-first, 4 samples/byte, ceil(5/4)=2 bytes/variant (hand-computed from the spec):
#   v0 byte0 = 00|(10<<2)|(11<<4)|(01<<6) = 0x78 ; byte1 = 10 = 0x02
#   v1 byte0 = 11|(11<<2)|(00<<4)|(10<<6) = 0x8f ; byte1 = 01 = 0x01
_BED = b"\x6c\x1b\x01" + bytes([0x78, 0x02, 0x8f, 0x01])
_BIM = "15\trs_A\t0\t100\tA\tG\n5\trs_B\t0\t200\tC\tT\n"
_FAM = "\n".join(f"F{i} S{i} 0 0 1 -9" for i in range(1, 6)) + "\n"


def _fixture():
    d = Path(tempfile.mkdtemp())
    (d / "x.bed").write_bytes(_BED)
    (d / "x.bim").write_text(_BIM, encoding="utf-8")
    (d / "x.fam").write_text(_FAM, encoding="utf-8")
    return d


def test_read_bim():
    bim = read_bim(_fixture() / "x.bim")
    assert len(bim) == 2
    assert bim[0].chrom == "15" and bim[0].pos == 100 and bim[0].a1 == "A" and bim[0].a2 == "G"
    assert bim[1].vid == "rs_B" and bim[1].index == 1


def test_read_fam():
    assert read_fam(_fixture() / "x.fam") == ["S1", "S2", "S3", "S4", "S5"]


def test_decode_dosages_matches_spec():
    d = _fixture()
    got = read_bed_variants(d / "x.bed", n_samples=5, indices=[0, 1])
    assert got[0] == [0, 1, 2, None, 1]
    assert got[1] == [2, 2, 0, 1, None]


def test_seek_single_variant():
    d = _fixture()
    # reading only variant 1 must give the same as reading both (direct seek, not sequential)
    assert read_bed_variants(d / "x.bed", 5, [1])[1] == [2, 2, 0, 1, None]


def test_find_variants_chrom_prefix_tolerant():
    bim = read_bim(_fixture() / "x.bim")
    assert [v.index for v in find_variants(bim, chrom="chr15", pos=100)] == [0]
    assert [v.index for v in find_variants(bim, chrom="5", pos=200)] == [1]
    assert find_variants(bim, chrom="15", pos=999) == []


def test_genotype_string():
    assert genotype_string(0, "A", "G") == "A/A"
    assert genotype_string(1, "A", "G") == "A/G"
    assert genotype_string(2, "A", "G") == "G/G"
    assert genotype_string(None, "A", "G") is None


def test_bad_magic_raises():
    d = Path(tempfile.mkdtemp())
    (d / "b.bed").write_bytes(b"\x00\x00\x00\x78")
    try:
        read_bed_variants(d / "b.bed", 5, [0])
    except PlinkFormatError:
        return
    raise AssertionError("expected PlinkFormatError on bad magic")


def test_truncated_block_raises():
    d = Path(tempfile.mkdtemp())
    (d / "t.bed").write_bytes(b"\x6c\x1b\x01" + bytes([0x78]))  # variant needs 2 bytes at n=5, only 1 present
    try:
        read_bed_variants(d / "t.bed", 5, [0])
    except PlinkFormatError:
        return
    raise AssertionError("expected PlinkFormatError on truncated block")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
