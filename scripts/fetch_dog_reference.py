"""Fetch a small canFam4 (UU_Cfam_GSD_1.0) reference slice for the masked-reconstruction probe.

Uses NCBI E-utilities efetch to pull a bounded region of one RefSeq chromosome accession, so the
download is kilobytes (a slice) not the whole ~2.4 Gb genome. Caches to D:/dna_decode_cache/dog_ref/.
Reversible network read; no money.

canFam4 chromosome accessions (UU_Cfam_GSD_1.0 = GCF_011100685.1): chr1 = NC_051804.1.

    uv run python scripts/fetch_dog_reference.py --accession NC_051804.1 --start 20000000 --stop 20006000
"""
from __future__ import annotations

import argparse
import sys
import urllib.parse
import urllib.request
from pathlib import Path

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_CACHE = Path("D:/dna_decode_cache/dog_ref")


def fetch_slice(accession: str, start: int, stop: int) -> str:
    params = {
        "db": "nuccore", "id": accession, "seq_start": str(start), "seq_stop": str(stop),
        "rettype": "fasta", "retmode": "text",
    }
    url = f"{_EUTILS}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"User-Agent": "dna_decode/1.0 (masked-recon probe)"})
    with urllib.request.urlopen(req, timeout=60) as r:
        text = r.read().decode("utf-8")
    seq = "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))
    seq = seq.upper()
    if not seq or any(c not in "ACGTN" for c in set(seq)):
        raise RuntimeError(f"unexpected efetch payload for {accession}:{start}-{stop} (len={len(seq)})")
    return seq


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fetch a small canFam4 reference slice via NCBI efetch.")
    ap.add_argument("--accession", default="NC_051804.1", help="RefSeq chromosome accession (canFam4)")
    ap.add_argument("--start", type=int, default=20000000)
    ap.add_argument("--stop", type=int, default=20006000)
    ap.add_argument("--out", default=None, help="output .txt (raw sequence); default = cache dir")
    args = ap.parse_args(argv)

    out = Path(args.out) if args.out else _CACHE / f"{args.accession}_{args.start}_{args.stop}.txt"
    if out.exists():
        print(f"cached: {out} ({len(out.read_text().strip())} bp)")
        return 0
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        seq = fetch_slice(args.accession, args.start, args.stop)
    except Exception as e:  # noqa: BLE001
        print(f"error: {e}", file=sys.stderr)
        return 1
    out.write_text(seq)
    n_valid = sum(c in "ACGT" for c in seq)
    print(f"wrote {out}  ({len(seq)} bp, {n_valid} ACGT, {len(seq) - n_valid} N)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
