"""Read-level PGx site genotyping from a 1000G 30x CRAM — recovers sites the phased VCF drops.

The NYGC 30x PHASED VCF panel (what `scripts/fetch_1000g_region.py` slices) FILTERS OUT some star-allele-
defining variants that ARE present in the underlying reads — a documented artifact (Star Allele Search,
PMC10811916: "several star-allele-defining variants present in the Phase 3 10x dataset were absent from the
NYGC phased 30x VCF files"). Our TPMT *6/*12/*40 residual-silent samples are exactly this: the sentinels are
correct but the phased VCF never genotypes those sites, so the concordance harness can't exercise them.

This tool genotypes ARBITRARY sites directly from the CRAM READS (the full data), Docker-free-adjacent (uses
the same Docker samtools biocontainer + ENA-reference auto-fetch pattern as `cyp2d6_pileup_gen.py`, no full
reference download, range-requests over the remote CRAM). It is the general form of the CYP2D6 pileup path:
give it any gene's SENTINELS (or a raw site list) and it returns per-site read-level genotype calls.

    uv run python scripts/pgx_cram_genotype.py --sample NA18603 --gene tpmt \
        --cram http://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239377/NA18603.final.cram
    # resolve the CRAM URL first: scripts/resolve_1000g_cram.py --sample NA18603

Pure per-site allele counting (mpileup `-B -q 0 -Q 0`, no `-f` -> the base column carries LITERAL read bases
ACGT, matching `cyp2d6_pileup_gen`). NOT a full star-allele caller — it answers "is this non-core allele's
ALT present in the reads, and at what genotype?" which is exactly what proves a panel-limited site recovers.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

_SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
_ENA_REF = "https://www.ebi.ac.uk/ena/cram/md5/%s"
_PAD = 100
_MIN_DEPTH = 8          # below this the call is UNCALLABLE (not enough reads)
_HET_LO, _HET_HI = 0.15, 0.85   # ALT fraction thresholds for 0/0 | 0/1 | 1/1


@dataclass
class SiteCall:
    label: str
    chrom: str
    pos: int
    ref: str
    alt: str
    depth: int
    ref_count: int
    alt_count: int
    alt_fraction: float
    genotype: str            # "0/0" | "0/1" | "1/1" | "UNCALLABLE"
    alt_present: bool        # >=2 ALT reads -> the non-core allele IS in the reads (panel-recovered)
    note: str = ""

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in
                ("label", "chrom", "pos", "ref", "alt", "depth", "ref_count", "alt_count",
                 "alt_fraction", "genotype", "alt_present", "note")}


def _norm_chrom(c: str) -> str:
    return c if str(c).lower().startswith("chr") else f"chr{c}"


def parse_mpileup(raw: str) -> dict[int, tuple[str, str]]:
    """samtools mpileup -> {pos: (ref_col, bases_col)}. field 1=pos, 2=ref, 4=bases."""
    out: dict[int, tuple[str, str]] = {}
    for line in raw.splitlines():
        f = line.split("\t")
        if len(f) >= 5 and f[1].isdigit():
            out[int(f[1])] = (f[2], f[4])
    return out


def call_site(bases: str, ref: str, alt: str) -> tuple[int, int, int, float, str, bool]:
    """Count ALT vs REF literal bases -> (depth, ref_count, alt_count, alt_frac, genotype, alt_present)."""
    counts = Counter(c.upper() for c in bases if c.upper() in "ACGT")
    depth = sum(counts.values())
    ref_count, alt_count = counts.get(ref.upper(), 0), counts.get(alt.upper(), 0)
    denom = ref_count + alt_count
    frac = round(alt_count / denom, 3) if denom else 0.0
    if depth < _MIN_DEPTH:
        return depth, ref_count, alt_count, frac, "UNCALLABLE", False
    if frac < _HET_LO:
        gt = "0/0"
    elif frac <= _HET_HI:
        gt = "0/1"
    else:
        gt = "1/1"
    return depth, ref_count, alt_count, frac, gt, alt_count >= 2


def genotype_sites(sample: str, cram: str, sites: list[dict], docker_run=None) -> list[SiteCall]:
    """One mpileup over the sites' span; return a SiteCall per site. `docker_run` injectable for tests."""
    if not sites:
        return []
    if docker_run is None:
        from tools.docker_runner import run as docker_run
    chrom = _norm_chrom(sites[0]["chrom"])
    lo = min(s["pos"] for s in sites) - _PAD
    hi = max(s["pos"] for s in sites) + _PAD
    out = docker_run(_SAMTOOLS_IMAGE,
                     ["samtools", "mpileup", "-B", "-q", "0", "-Q", "0", "-r",
                      f"{chrom}:{lo}-{hi}", cram],
                     env={"REF_PATH": _ENA_REF}, capture_output=True, check=False, timeout=600)
    piles = parse_mpileup(out.stdout or "")
    calls: list[SiteCall] = []
    for s in sites:
        _refcol, bases = piles.get(s["pos"], ("", ""))
        depth, rc, ac, frac, gt, present = call_site(bases, s["ref"], s["alt"])
        calls.append(SiteCall(s["label"], _norm_chrom(s["chrom"]), s["pos"], s["ref"], s["alt"],
                              depth, rc, ac, frac, gt, present,
                              note=("no pileup at site" if not bases else "")))
    return calls


def sites_from_gene(gene: str) -> list[dict]:
    """Load a gene's SENTINELS as the site list (the panel-recovery target set)."""
    import importlib
    mod = importlib.import_module(f"dna_decode.pgx.{gene}_catalog")
    return [{"label": s.implies, "chrom": s.chrom, "pos": s.pos, "ref": s.ref, "alt": s.alt}
            for s in getattr(mod, "SENTINELS", [])]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="pgx_cram_genotype",
                                 description="Read-level PGx site genotyping from a 1000G CRAM (recovers "
                                             "phased-VCF-dropped sites).")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cram", required=True, help="CRAM URL (resolve via scripts/resolve_1000g_cram.py)")
    ap.add_argument("--gene", help="genotype this gene's SENTINELS (e.g. tpmt)")
    ap.add_argument("--sites", help="raw sites as chrom:pos:ref:alt:label,... (overrides --gene)")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    if args.sites:
        sites = []
        for tok in args.sites.split(","):
            c, p, r, a, *lab = tok.split(":")
            sites.append({"label": lab[0] if lab else f"{c}:{p}", "chrom": c, "pos": int(p),
                          "ref": r, "alt": a})
    elif args.gene:
        sites = sites_from_gene(args.gene)
    else:
        print("error: give --gene <g> or --sites chrom:pos:ref:alt:label,...", file=sys.stderr)
        return 2

    calls = genotype_sites(args.sample, args.cram, sites)
    out = {"sample": args.sample, "n_sites": len(calls),
           "n_alt_present": sum(c.alt_present for c in calls),
           "calls": [c.as_dict() for c in calls]}
    if args.json:
        print(json.dumps(out, indent=2))
    else:
        print(f"{args.sample}: {out['n_alt_present']}/{out['n_sites']} sites with ALT present in reads")
        for c in calls:
            print(f"  {c.label:6} {c.chrom}:{c.pos} {c.ref}>{c.alt}  depth={c.depth} "
                  f"alt={c.alt_count}/{c.depth} ({c.alt_fraction}) -> {c.genotype}"
                  f"{'  [ALT PRESENT]' if c.alt_present else ''}{'  ' + c.note if c.note else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
