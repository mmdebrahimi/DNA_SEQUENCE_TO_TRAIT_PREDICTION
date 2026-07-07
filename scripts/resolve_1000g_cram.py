#!/usr/bin/env python
"""Resolve a 1000G sample -> its 30x high-coverage CRAM URL (the read-level/structural-surface substrate).

The CYP2D6 structural + hybrid-identity surfaces need per-sample CRAMs (read depth + PSV pileups), NOT the
phased VCF panel that `scripts/fetch_1000g_region.py` slices. This resolves an arbitrary 1000G sample name
(e.g. NA12156, HG00276) to the exact `*.final.cram` URL by parsing the canonical 1000G high-coverage
`.sequence.index` files, so `scripts/cyp2d6_structural_probe.py --compute` (and the hybrid pileup step) can
build their `sample/truth/cram_url` cohort TSVs reproducibly instead of hand-curled URLs.

Covers BOTH the 2504 panel and the additional-698-related index (family members like the CEPH trios). The
CRAM field in the index is `ftp://ftp.sra.ebi.ac.uk/...`; samtools/curl read it over http, so we normalize
`ftp://` -> `http://` (empirically the http mirror serves the same bytes + range requests, which the
Docker-free/remote-slice path relies on).

Usage:
    uv run python scripts/resolve_1000g_cram.py --sample NA12156
    uv run python scripts/resolve_1000g_cram.py --samples NA12156 HG00276 --truths '*1/*4' '*4/*5' \
        --out data/pgx_1000g/cyp2d6_extra_cohort.tsv     # emit a --compute-ready cohort TSV
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

# Canonical 1000G 30x high-coverage sequence indices (NYGC/EBI). The 2504 panel + the 698 related samples.
_BASE = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
_INDICES = (
    f"{_BASE}/1000G_2504_high_coverage.sequence.index",
    f"{_BASE}/additional_698_related/1000G_698_related_high_coverage.sequence.index",
)


def _normalize_scheme(url: str) -> str:
    """ftp://ftp.sra.ebi.ac.uk/... -> http://ftp.sra.ebi.ac.uk/... (the http mirror serves range requests
    the remote-slice path needs; ftp:// is what the index records)."""
    return "http://" + url[len("ftp://"):] if url.startswith("ftp://") else url


def parse_cram_url_from_index(index_text: str, sample: str) -> str | None:
    """PURE: given a 1000G sequence.index body + a sample name, return the sample's `*.final.cram` URL
    (http-normalized) or None. `sample` must appear as a WHOLE whitespace-delimited field (so NA12156 does
    NOT match NA121560); the CRAM is the field ending in `.final.cram` (or any `.cram`)."""
    for line in index_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        if sample not in fields:
            continue
        for f in fields:
            if f.endswith(".cram"):
                return _normalize_scheme(f)
    return None


def _fetch(url: str, timeout: int = 90) -> str:
    """Fetch an index over http via curl (urllib intermittently hangs on the big EBI files on this host —
    same rationale as fetch_1000g_region.py)."""
    r = subprocess.run(["curl", "-sSL", "--max-time", str(timeout), url],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"curl rc={r.returncode} for {url}: {r.stderr[:160]}")
    return r.stdout


def resolve_1000g_cram_url(sample: str, indices: tuple[str, ...] = _INDICES,
                           fetch=_fetch) -> str | None:
    """Resolve a 1000G sample -> its 30x CRAM URL, searching the 2504 index then the 698-related index.
    `fetch` is injectable for offline testing. Returns None if the sample is in neither index."""
    for idx in indices:
        hit = parse_cram_url_from_index(fetch(idx), sample)
        if hit:
            return hit
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Resolve 1000G sample(s) -> 30x CRAM URL(s).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--sample", help="single sample -> print its CRAM URL")
    g.add_argument("--samples", nargs="+", help="multiple samples -> emit a cohort TSV (needs --truths + --out)")
    ap.add_argument("--truths", nargs="+", help="per-sample truth diplotypes (same order/len as --samples)")
    ap.add_argument("--out", type=Path, help="write a sample/truth/cram_url cohort TSV here")
    args = ap.parse_args(argv)

    if args.sample:
        url = resolve_1000g_cram_url(args.sample)
        if not url:
            print(f"NOT_FOUND: {args.sample} in neither 1000G high-coverage index", file=sys.stderr)
            return 1
        print(url)
        return 0

    if not (args.truths and args.out) or len(args.truths) != len(args.samples):
        print("--samples requires --truths (same length) + --out", file=sys.stderr)
        return 2
    lines = ["sample\ttruth\tcram_url"]
    missing = []
    for s, t in zip(args.samples, args.truths):
        url = resolve_1000g_cram_url(s)
        if not url:
            missing.append(s)
            continue
        lines.append(f"{s}\t{t}\t{url}")
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[cohort TSV ({len(lines)-1}/{len(args.samples)} resolved) -> {args.out}]")
    if missing:
        print(f"NOT_FOUND: {', '.join(missing)}", file=sys.stderr)
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
