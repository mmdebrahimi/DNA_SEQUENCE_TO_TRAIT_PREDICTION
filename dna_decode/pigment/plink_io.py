"""Minimal pure-Python PLINK1 .bed/.bim/.fam reader — extract specific variants' genotypes.

Needed by the Darwin's Ark dog coat-colour scoring (scripts/dog_coat_darwins_ark_validate.py): that cohort
ships genotypes as a PLINK1 binary set, and this host has no `plink`/`bcftools` native (the Docker-free
posture, same reason `scripts/fetch_1000g_region.py` decodes BGZF by hand). We only need a HANDFUL of
variants (the ~5 coat-colour causal loci) out of millions, so we seek directly to each variant's block in
the SNP-major .bed rather than loading the whole matrix.

PLINK1 .bed spec (SNP-major, mode byte 0x01):
  - 3 magic bytes: 0x6c 0x1b 0x01.
  - Then, per VARIANT (in .bim order): ceil(n_samples/4) bytes; each byte packs 4 samples, 2 bits each,
    LEAST-significant pair first.
  - 2-bit code -> genotype (a1/a2 are .bim columns 5/6):
        00 = a1/a1   (homozygous first allele)
        01 = MISSING
        10 = a1/a2   (heterozygous)
        11 = a2/a2   (homozygous second allele)
    So a2-allele dosage = {00:0, 10:1, 11:2, 01:None}. We return a1/a2 bases from the .bim so the caller
    maps dosage onto whichever base is the coat-allele it counts (strand/allele harmonization done upstream).

Pure-python, offline, deterministic. Spec-verifiable (synthetic round-trip test), no external data needed.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_MAGIC = b"\x6c\x1b\x01"
# 2-bit code -> a2 dosage (None = missing). Index by the 2-bit value.
_CODE_TO_A2_DOSAGE = {0b00: 0, 0b01: None, 0b10: 1, 0b11: 2}


class PlinkFormatError(ValueError):
    """Malformed .bed magic / truncated block (never a silent wrong genotype)."""


@dataclass(frozen=True)
class BimVariant:
    index: int          # 0-based row in the .bim (== variant block order in the .bed)
    chrom: str
    vid: str
    pos: int
    a1: str             # PLINK allele 1 (.bim col 5) — the "counted" allele in the 2-bit encoding
    a2: str             # PLINK allele 2 (.bim col 6)


def read_bim(bim_path: str | Path) -> list[BimVariant]:
    out: list[BimVariant] = []
    for i, line in enumerate(Path(bim_path).read_text(encoding="utf-8", errors="replace").splitlines()):
        f = line.split()
        if len(f) < 6:
            continue
        try:
            pos = int(f[3])
        except ValueError:
            pos = -1
        out.append(BimVariant(index=i, chrom=str(f[0]), vid=f[1], pos=pos, a1=f[4], a2=f[5]))
    return out


def read_fam(fam_path: str | Path) -> list[str]:
    """Return sample IIDs in .fam order (column 2)."""
    out: list[str] = []
    for line in Path(fam_path).read_text(encoding="utf-8", errors="replace").splitlines():
        f = line.split()
        if len(f) >= 2:
            out.append(f[1])
    return out


def find_variants(bim: list[BimVariant], *, chrom: str, pos: int) -> list[BimVariant]:
    """All .bim variants at a (chrom, pos). Chromosome compared with/without a 'chr' prefix."""
    def norm(c: str) -> str:
        return c[3:] if c.lower().startswith("chr") else c
    nc = norm(str(chrom))
    return [v for v in bim if norm(v.chrom) == nc and v.pos == pos]


def read_bed_variants(bed_path: str | Path, n_samples: int, indices: list[int]) -> dict[int, list[int | None]]:
    """{variant_index: [a2-dosage per sample in .fam order]} for the requested variant indices.

    Seeks directly to each variant's block (SNP-major), so cost is O(len(indices)), not O(n_variants).
    a2-dosage is 0/1/2 (copies of the .bim a2 allele) or None (missing).
    """
    if n_samples <= 0:
        raise PlinkFormatError("n_samples must be > 0")
    bytes_per_variant = (n_samples + 3) // 4
    out: dict[int, list[int | None]] = {}
    with open(bed_path, "rb") as fh:
        magic = fh.read(3)
        if magic != _MAGIC:
            raise PlinkFormatError(f"bad .bed magic {magic!r} (expected SNP-major {_MAGIC!r}; "
                                   "individual-major .bed is unsupported — re-export with plink --make-bed)")
        for vi in indices:
            fh.seek(3 + vi * bytes_per_variant)
            block = fh.read(bytes_per_variant)
            if len(block) != bytes_per_variant:
                raise PlinkFormatError(f"variant {vi}: truncated block ({len(block)}/{bytes_per_variant} bytes)")
            dosages: list[int | None] = []
            for s in range(n_samples):
                byte = block[s // 4]
                code = (byte >> (2 * (s % 4))) & 0b11
                dosages.append(_CODE_TO_A2_DOSAGE[code])
            out[vi] = dosages
    return out


def genotype_string(dosage: int | None, a1: str, a2: str) -> str | None:
    """a2-dosage -> a diploid base string ('A/G') using the .bim alleles, or None if missing.

    dosage counts a2 copies: 0 -> a1/a1, 1 -> a1/a2, 2 -> a2/a2.
    """
    if dosage is None:
        return None
    if dosage == 0:
        return f"{a1}/{a1}"
    if dosage == 1:
        return f"{a1}/{a2}"
    return f"{a2}/{a2}"
