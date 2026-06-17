"""TB VCF parsing + acquisition for the CRyPTIC organism-routed decoder cell (NON-FROZEN).

Two CRyPTIC per-isolate VCFs, both vs H37Rv NC_000962.3:
  - **masked** VCF: variant sites with the CRyPTIC mask applied (~217 KB/isolate). CANONICAL for
    determinant CALLS — matches the verified coordinate alignment (rpoB S450L -> 761155 C>T,
    katG S315T -> 2155168 C>G) + the cached subset.
  - **regenotyped (regeno)** VCF: per-position re-genotyping (~177 MB/isolate; explicit 0/0 ref-calls,
    ./. no-calls, alt-calls). Source of truth for CALLABILITY at a determinant window.

Call rule (brainstorm C1, verified against the cached VCFs): a NON-REFERENCE call requires
`FILTER==PASS` AND a `GT` allele index >= 1. There is **no `GCP` FORMAT field** — the quality floor is
the `MIN_DP`/`MIN_FRS`/`MIN_GCP` FILTERs, so PASS subsumes it. `DP`/`FRS`/`GT_CONF_PERCENTILE` are
exposed as provenance only, never a second hard floor.

Callability rule (brainstorm C3, ratified "full regeno"): a genomic position is CALLABLE iff the regeno
VCF carries a PASS record at that position whose `GT` is an explicit call (0/0 or alt), NOT `./.` and not
absent. A determinant window is callable iff all its queried positions are callable. An uncallable
determinant window -> the cell ABSTAINs (never susceptible-by-absence).

Acquisition caches to **D:** (the external 4.4 TB drive) by default — the full-cohort regeno fetch is
~1.6 TB and must NOT land on the chronically-full C:.
"""
from __future__ import annotations

import csv
import gzip
import urllib.request
from dataclasses import dataclass
from pathlib import Path

CHROM = "NC_000962.3"  # H37Rv — single chromosome, so VCF POS alone keys a site.

REUSE_CSV = Path("data/raw/cryptic/CRyPTIC_reuse_table_20240917.csv")
FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reproducibility/"
# D: cache (4.4 TB free) — NOT C: (the ~1.6 TB regeno cohort would exhaust C:'s ~11 GB).
DEFAULT_CACHE = Path("D:/dna_decode_cache/cryptic")

_NOCALL = {"./.", ".", "", ".|."}


@dataclass(frozen=True)
class VariantCall:
    """A single non-reference PASS call at one H37Rv position (provenance kept, not a quality gate)."""
    pos: int
    ref: str
    alt: str          # the specific ALT allele selected by the GT index (multi-allelic resolved)
    gt: str
    dp: int | None = None
    frs: float | None = None
    gt_conf_percentile: float | None = None


def _split_fmt(fmt: str, sample: str) -> dict[str, str]:
    keys = fmt.split(":")
    vals = sample.split(":")
    return {k: vals[i] for i, k in enumerate(keys) if i < len(vals)}


def _gt_alleles(gt: str) -> list[int]:
    """Allele indices from a GT field; non-numeric (./.) -> empty."""
    out = []
    for a in gt.replace("|", "/").split("/"):
        if a.isdigit():
            out.append(int(a))
    return out


def _iter_records(vcf_text: str):
    """Yield (pos, ref, alt_field, filt, fmt, sample) for non-header records on CHROM."""
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 10:
            continue
        if f[0] != CHROM:
            continue
        try:
            pos = int(f[1])
        except ValueError:
            continue
        yield pos, f[3], f[4], f[6].strip(), f[8], f[9]


def parse_masked_calls(vcf_text: str) -> dict[int, VariantCall]:
    """Parse a masked VCF -> {pos: VariantCall} for PASS, non-reference (GT allele >= 1) records.

    Multi-allelic ALT (`A,C,T`) is resolved to the single allele the GT points at. `0/0`-with-ALT,
    `./.`, FILTER-failed, and reference records are excluded (NOT calls).
    """
    calls: dict[int, VariantCall] = {}
    for pos, ref, alt_field, filt, fmt, sample in _iter_records(vcf_text):
        if filt != "PASS":
            continue
        fields = _split_fmt(fmt, sample)
        alleles = _gt_alleles(fields.get("GT", ""))
        nonref = [a for a in alleles if a >= 1]
        if not nonref:
            continue  # reference call or no-call — not a resistance call
        alts = alt_field.split(",")
        idx = nonref[0] - 1
        if idx < 0 or idx >= len(alts):
            continue
        alt = alts[idx].strip()
        if alt in (".", "", "<NON_REF>"):
            continue
        calls[pos] = VariantCall(
            pos=pos, ref=ref, alt=alt, gt=fields.get("GT", ""),
            dp=_as_int(fields.get("DP")),
            frs=_as_float(fields.get("FRS")),
            gt_conf_percentile=_as_float(fields.get("GT_CONF_PERCENTILE")),
        )
    return calls


def callable_positions(regeno_text: str, positions) -> dict[int, bool]:
    """{pos: callable?} from a regeno VCF. Callable iff a PASS record exists at pos with an explicit
    (non-`./.`) GT. A position ABSENT from the regeno VCF is uncallable (conservative)."""
    want = set(positions)
    seen: dict[int, bool] = {p: False for p in want}
    for pos, _ref, _alt, filt, fmt, sample in _iter_records(regeno_text):
        if pos not in want:
            continue
        gt = _split_fmt(fmt, sample).get("GT", "")
        seen[pos] = (filt == "PASS") and (gt not in _NOCALL)
    return seen


def is_window_callable(regeno_text: str, lo: int, hi: int) -> bool:
    """A determinant window [lo, hi] (inclusive) is callable iff EVERY position in it is callable."""
    if hi < lo:
        raise ValueError(f"window hi {hi} < lo {lo}")
    flags = callable_positions(regeno_text, range(lo, hi + 1))
    return all(flags.values())


def _as_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# --- acquisition (thin; caches to D: by default) ------------------------------------------------

def _vcf_url(rel_path: str) -> str:
    rel = rel_path.strip()
    return FTP_BASE + (rel.split("reproducibility/", 1)[1] if "reproducibility/" in rel
                       else rel.lstrip("./"))


def fetch_vcf(rel_path: str, kind: str, cache_dir: Path = DEFAULT_CACHE, timeout: int = 120) -> bytes | None:
    """Fetch + cache one decompressed VCF. `kind` in {'masked','regeno'} -> a D: subdir. Skip-existing."""
    sub = cache_dir / kind
    sub.mkdir(parents=True, exist_ok=True)
    cache = sub / (rel_path.replace("/", "_").replace("..", "").strip("_"))
    if cache.exists():
        return cache.read_bytes()
    try:
        with urllib.request.urlopen(_vcf_url(rel_path), timeout=timeout) as r:
            data = gzip.decompress(r.read())
        cache.write_bytes(data)
        return data
    except Exception as e:  # noqa: BLE001
        print(f"    [vcf fetch fail] {kind} {rel_path[:60]}: {type(e).__name__}: {e}")
        return None


def reuse_rows(csv_path: Path = REUSE_CSV) -> list[dict]:
    with open(csv_path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def vcf_paths_for(row: dict) -> tuple[str, str]:
    """(masked_rel, regeno_rel) from a reuse-table row (columns VCF + REGENOTYPED_VCF)."""
    return (row.get("VCF") or "").strip(), (row.get("REGENOTYPED_VCF") or "").strip()
