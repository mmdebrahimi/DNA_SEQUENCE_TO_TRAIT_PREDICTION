"""TBProfiler/tbdb lineage barcode (Napier-2020 lineage) loader + pin-verify (NON-FROZEN).

Ratified F: use the barcode DATA (a position -> lineage-allele table) applied directly to our VCFs,
NOT the TBProfiler caller (avoids the wrapper-trap). Pinned at tbdb commit 618cf0ff + sha256 in
`data/raw/tb_lineage_barcode/CHECKSUMS`. 1111 lineage-defining SNPs — far finer sublineage resolution
than Coll-2014's 62, so the downstream clonality collapse is more conservative (higher effective_n).

barcode.bed columns (TBProfiler BED): chrom, start(0-based), end(=VCF POS, 1-based), lineage, allele,
major-lineage, spoligotype, RD. We key on POS (single chromosome) + the derived `allele`.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

BARCODE_DIR = Path("data/raw/tb_lineage_barcode")
BARCODE_FILE = "barcode.bed"
CHECKSUMS_FILE = "CHECKSUMS"
PINNED_COMMIT = "618cf0ff5f22886971bd437929d2c49defa6c7bf"


class BarcodePinError(RuntimeError):
    pass


@dataclass(frozen=True)
class BarcodeSNP:
    pos: int        # 1-based H37Rv position (VCF POS)
    lineage: str    # e.g. "lineage4.2.2.1", "La1.8"
    allele: str     # the derived lineage-defining allele


def verify_pin(barcode_dir: Path = BARCODE_DIR) -> bool:
    """Recompute sha256 of barcode.bed vs the pinned CHECKSUMS; raise BarcodePinError on drift."""
    pins: dict[str, str] = {}
    for ln in (barcode_dir / CHECKSUMS_FILE).read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if ln and not ln.startswith("#"):
            p = ln.split()
            if len(p) >= 3:
                pins[p[-1]] = p[0]
    want = pins.get(BARCODE_FILE)
    got = hashlib.sha256((barcode_dir / BARCODE_FILE).read_bytes()).hexdigest()
    if want != got:
        raise BarcodePinError(f"barcode drift: pinned {str(want)[:12]} got {got[:12]}")
    return True


def load_barcode(barcode_dir: Path = BARCODE_DIR, verify: bool = True) -> list[BarcodeSNP]:
    if verify:
        verify_pin(barcode_dir)
    out: list[BarcodeSNP] = []
    for ln in (barcode_dir / BARCODE_FILE).read_text(encoding="utf-8").splitlines():
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 5:
            continue
        try:
            pos = int(f[2])  # BED end = 1-based VCF POS
        except ValueError:
            continue
        out.append(BarcodeSNP(pos=pos, lineage=f[3].strip(), allele=f[4].strip()))
    return out


def barcode_available(barcode_dir: Path = BARCODE_DIR) -> bool:
    return (barcode_dir / BARCODE_FILE).exists()
