"""WHO M. tuberculosis mutation catalogue v2 (2023) loader + pin-verify + RIF/INH join (NON-FROZEN).

The catalogue is pinned at commit 0bb39143 with per-file sha256 in `data/raw/who_tb_catalogue/CHECKSUMS`
(the repo has NO releases, so SHA + checksum is the only reproducible pin). This module:
  - verifies the on-disk files against the pinned checksums (refuses drift), and
  - JOINs grade-bearing master rows <-> the genomic-coordinate file on the shared `variant` key to
    build the RIF + INH determinant table.

Determinant scope (brainstorm ratification A): **all WHO grade-1/2 ("Associated with resistance") loci
per drug** — for RIF that is essentially rpoB; for INH it spans katG + inhA/fabG1 promoter + grade-1/2
inhA-coding + ahpC. We do NOT hand-narrow to a gene whitelist (a narrowed subset isn't "the catalogue").

Files are gitignored (37 MB master); only CHECKSUMS is tracked. Real-catalogue consumers/tests degrade
to skip when the files are absent (offline-safe), mirroring the project's BLAST/vf_diff pattern.
"""
from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

CAT_DIR = Path("data/raw/who_tb_catalogue")
MASTER_FILE = "WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
COORDS_FILE = "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
CHECKSUMS_FILE = "CHECKSUMS"

GRADE_COL = "INITIAL CONFIDENCE GRADING"
GRADE_1_2_PREFIXES = ("1)", "2)")  # "1) Assoc w R" / "2) Assoc w R - Interim"
CHROM = "NC_000962.3"

# Reuse-table / CRyPTIC drug code -> catalogue `drug` value.
DRUG_CATALOGUE_NAME = {
    "rifampicin": "Rifampicin", "rif": "Rifampicin",
    "isoniazid": "Isoniazid", "inh": "Isoniazid",
}


class CataloguePinError(RuntimeError):
    """Raised when an on-disk catalogue file does not match its pinned sha256."""


@dataclass(frozen=True)
class Determinant:
    drug: str
    gene: str
    variant: str
    grade: str
    tier: str
    chrom: str
    pos: int
    ref: str
    alt: str


def _parse_checksums(cat_dir: Path) -> dict[str, str]:
    """{filename: sha256} from the CHECKSUMS manifest (skips `#` comment lines)."""
    out: dict[str, str] = {}
    for ln in (cat_dir / CHECKSUMS_FILE).read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split()
        if len(parts) >= 3:
            out[parts[-1]] = parts[0]
    return out


def verify_pins(cat_dir: Path = CAT_DIR) -> dict[str, bool]:
    """Recompute sha256 of each pinned file; raise CataloguePinError on any mismatch. Returns {file: True}."""
    pins = _parse_checksums(cat_dir)
    result: dict[str, bool] = {}
    for fname, want in pins.items():
        fp = cat_dir / fname
        if not fp.exists():
            raise CataloguePinError(f"pinned file missing: {fname}")
        got = hashlib.sha256(fp.read_bytes()).hexdigest()
        if got != want:
            raise CataloguePinError(f"checksum drift on {fname}: pinned {want[:12]} got {got[:12]}")
        result[fname] = True
    return result


def _load_coords(cat_dir: Path) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    with open(cat_dir / COORDS_FILE, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            v = (row.get("variant") or "").strip()
            if not v:
                continue
            try:
                pos = int(row["position"])
            except (KeyError, ValueError):
                continue
            out.setdefault(v, []).append({
                "chrom": (row.get("chromosome") or "").strip(),
                "pos": pos,
                "ref": (row.get("reference_nucleotide") or "").strip(),
                "alt": (row.get("alternative_nucleotide") or "").strip(),
            })
    return out


def is_grade_1_2(grade: str) -> bool:
    return (grade or "").startswith(GRADE_1_2_PREFIXES)


def grade_1_2_count(cat_dir: Path = CAT_DIR) -> int:
    """Distinct grade-1/2 `variant`s across the whole master (regression pin == 438)."""
    seen: set[str] = set()
    with open(cat_dir / MASTER_FILE, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if is_grade_1_2(row.get(GRADE_COL, "")):
                v = (row.get("variant") or "").strip()
                if v:
                    seen.add(v)
    return len(seen)


def load_determinants(drug: str, cat_dir: Path = CAT_DIR, verify: bool = True) -> list[Determinant]:
    """Grade-1/2 determinants for `drug` (any reuse-code/name), joined master<->coords on `variant`.

    One Determinant per (variant, coord-row) — a variant with several nucleotide encodings yields
    several rows, each a (pos, ref, alt) to match against a VCF.
    """
    if verify:
        verify_pins(cat_dir)
    want_drug = DRUG_CATALOGUE_NAME.get(drug.strip().lower(), drug)
    coords = _load_coords(cat_dir)
    out: list[Determinant] = []
    with open(cat_dir / MASTER_FILE, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("drug") or "").strip() != want_drug:
                continue
            grade = (row.get(GRADE_COL) or "").strip()
            if not is_grade_1_2(grade):
                continue
            variant = (row.get("variant") or "").strip()
            for c in coords.get(variant, []):
                out.append(Determinant(
                    drug=want_drug, gene=(row.get("gene") or "").strip(), variant=variant,
                    grade=grade, tier=(row.get("tier") or "").strip(),
                    chrom=c["chrom"], pos=c["pos"], ref=c["ref"], alt=c["alt"],
                ))
    return out


def catalogue_available(cat_dir: Path = CAT_DIR) -> bool:
    return (cat_dir / MASTER_FILE).exists() and (cat_dir / COORDS_FILE).exists()
