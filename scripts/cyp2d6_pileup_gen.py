#!/usr/bin/env python
"""Generate the CYP2D6/CYP2D7 PSV pileup files (`<sample>.d6.txt` / `.d7.txt`) from a real CRAM.

The read-level hybrid-identity surface (`scripts/cyp2d6_psv_evidence.py` + `dna_decode.pgx.cyp2d6_hybrid_
identity`) CONSUMES per-sample pileups but never generated them — the generation was an ad-hoc bash loop.
This commits it: one `samtools mpileup` over the CYP2D6 PSV span + one over the CYP2D7 span (both from the
committed Cyrius `data/cyp2d6_psv/CYP2D6_SNP_38.txt`), via Docker samtools with ENA reference auto-fetch (no
full-reference download), reformatted to the `pos<TAB>ref<TAB>bases` lines `_read_pileup` expects.

Flag contract (pinned): `samtools mpileup -B -q 0 -Q 0` (permissive — matches the depth path + the evidence
script's flag_contract). NO `-f`: without a reference the pileup base column carries the LITERAL read bases
(ACGT), which is exactly what `_base_counts` counts. `REF_PATH` still supplies the reference for CRAM DECODE.

Usage:
    uv run python scripts/cyp2d6_pileup_gen.py --sample NA12156 --pileup-dir data/cyp2d6_pileups \
        --cram http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239310/NA12156.final.cram
    # resolve the CRAM URL first with scripts/resolve_1000g_cram.py --sample NA12156
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts.cyp2d6_psv_evidence import load_psvs  # noqa: E402

_SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
_ENA_REF = "https://www.ebi.ac.uk/ena/cram/md5/%s"
_PAD = 200  # flank the PSV span so edge positions get full pileup depth


def psv_spans(psvs: list[dict]) -> tuple[tuple[str, int, int], tuple[str, int, int]]:
    """(CYP2D6 region, CYP2D7 region) = (chrom, min-pad, max+pad) over the PSV coords of each paralog."""
    chrom = psvs[0]["chrom"]
    d6 = [p["pos_d6"] for p in psvs]
    d7 = [p["pos_d7"] for p in psvs]
    return (chrom, min(d6) - _PAD, max(d6) + _PAD), (chrom, min(d7) - _PAD, max(d7) + _PAD)


def mpileup_to_pos_lines(raw: str) -> str:
    """samtools mpileup `chrom pos ref depth bases quals` -> `pos<TAB>ref<TAB>bases` (what `_read_pileup`
    parses: field 0 = pos digit, field 2 = base string)."""
    out = []
    for line in raw.splitlines():
        f = line.split("\t")
        if len(f) >= 5 and f[1].isdigit():
            out.append(f"{f[1]}\t{f[2]}\t{f[4]}")
    return "\n".join(out) + ("\n" if out else "")


def generate(sample: str, cram: str, pileup_dir: Path, psvs: list[dict], docker_run=None) -> Path:
    """Write `<sample>.d6.txt` + `<sample>.d7.txt` into pileup_dir. `docker_run` injectable for tests."""
    if docker_run is None:
        from tools.docker_runner import run as docker_run
    pileup_dir.mkdir(parents=True, exist_ok=True)
    d6r, d7r = psv_spans(psvs)
    for region, suffix in ((d6r, "d6"), (d7r, "d7")):
        reg = f"{region[0]}:{region[1]}-{region[2]}"
        out = docker_run(
            _SAMTOOLS_IMAGE,
            ["samtools", "mpileup", "-B", "-q", "0", "-Q", "0", "-r", reg, cram],
            env={"REF_PATH": _ENA_REF}, capture_output=True, check=False, timeout=600)
        (pileup_dir / f"{sample}.{suffix}.txt").write_text(
            mpileup_to_pos_lines(out.stdout or ""), encoding="utf-8")
    return pileup_dir


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Generate CYP2D6/CYP2D7 PSV pileups from a CRAM (Docker samtools).")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cram", required=True, help="CRAM URL (resolve via scripts/resolve_1000g_cram.py)")
    ap.add_argument("--pileup-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    psvs = load_psvs()
    generate(args.sample, args.cram, args.pileup_dir, psvs)
    d6 = args.pileup_dir / f"{args.sample}.d6.txt"
    d7 = args.pileup_dir / f"{args.sample}.d7.txt"
    print(f"[wrote {d6.name} ({sum(1 for _ in d6.open())} pos) + {d7.name} "
          f"({sum(1 for _ in d7.open())} pos) -> {args.pileup_dir}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
