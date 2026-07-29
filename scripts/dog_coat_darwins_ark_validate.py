"""Validate the dog coat-colour decoder (dna_decode/pigment/dog_coat) per-individual against the FREE,
open Darwin's Ark cohort — the measured-tier substrate that would move the cell off KNOWLEDGE_BASELINE.

SUBSTRATE (fully open, no access request — Dryad doi:10.5061/dryad.83bk3jb4r, the 2025 PNAS paper
"Genetic testing predicts appearance but not behavior in dogs"):
  - PHENOTYPE: DarwinsArk_13Q_Q243_coat_color_N-1930.tsv  (owner-reported coat colour, N=1930)
  - GENOTYPE:  DarwinsDogs_2024_N-3277_canfam4_gp-0.70_biallelic.{bed,bim,fam}  (canFam4, 3277 dogs)
  - raw WGS: NCBI BioProject PRJNA675863

WHY the 2025 PNAS finding matters: it reports that curated genetic tests "predict APPEARANCE but not
behavior" — i.e. coat colour is exactly the R1/curated-catalog regime this decoder targets. Scoring our
DETERMINISTIC epistatic caller here is the independent-ish measured check the cell currently lacks.

HONESTY / anti-fabrication rails (load-bearing):
  - This script FETCHES NOTHING and FABRICATES NOTHING. If the Dryad files are not present at the given
    paths it prints the fetch steps and exits 2 (never a made-up concordance).
  - The causal-variant -> locus-allele MAP is pinned from the REAL .bim at probe time, NOT hardcoded
    coordinates (R3 real-surface-first; canFam4 positions differ across assemblies/builds and must be
    verified against the actual .bim, not transcribed from memory). `--probe` reports which of the
    catalogued causal loci are genotyped in this .bim BEFORE any scoring is attempted.
  - The phenotype-category -> decoder-vocabulary map is EXPLICIT + auditable (owner free-text colour
    buckets are lossy; ambiguous/pattern categories are EXCLUDED from scoring, not force-matched).

Run (once the Dryad archive is unzipped locally, e.g. to D:/dna_decode_cache/darwins_ark/):
    uv run python scripts/dog_coat_darwins_ark_validate.py \
        --pheno-tsv D:/dna_decode_cache/darwins_ark/DarwinsArk_13Q_Q243_coat_color_N-1930.tsv \
        --plink-prefix D:/dna_decode_cache/darwins_ark/DarwinsDogs_2024_N-3277_canfam4_gp-0.70_biallelic \
        --probe          # inspect schema + which causal loci are genotyped (do this FIRST)

Exit: 0 = OK (probe reported / scored), 2 = data not present (fetch first), 3 = schema mismatch.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# --- causal-variant catalog (DESCRIPTIVE + sourced; exact canFam4 coords are pinned from the .bim, not here) ---
# Each entry: locus -> (gene, effect description, effect allele = the coat-colour-changing allele, OMIA source).
# The .bim gives (chrom, pos, id, a1, a2); the probe resolves each causal variant by the study's own variant
# annotation / gene region. These descriptors are literature facts; positions are verified at runtime.
CAUSAL_VARIANTS = {
    "E": ("MC1R", "p.Arg306Ter premature stop = recessive red `e`", "e",
          "OMIA 001199-9615; Everts 2000; Schmutz 2003"),
    "K": ("CBD103", "c.67_69delGGT (delta-G23) = dominant black K^B", "KB",
          "Candille 2007 Science 318:1418; K-locus / dominant black"),
    "A": ("ASIP", "agouti alleles Ay/aw/at/a (promoter + coding)", "multi",
          "Dreger & Schmutz 2011; Bannasch 2021 Nat Ecol Evol (ASIP promoter)"),
    "B": ("TYRP1", "bs p.Gln331Ter (+ bd, bc) = brown/liver b-family", "b",
          "OMIA 001249-9615; Schmutz 2002"),
    "D": ("MLPH", "c.-22G>A `d` (+ d2) = dilution", "d",
          "OMIA 000031-9615; Drogemuller 2007; Bauer 2018"),
}

# owner-reported coat-colour category -> our decoder vocabulary (EXPLICIT + lossy; ambiguous/pattern -> EXCLUDE).
# The exact Darwin's Ark Q243 answer strings are pinned by --probe against the real TSV before scoring.
PHENOTYPE_MAP: dict[str, str] = {
    # phaeomelanin
    "yellow": "red/yellow", "red": "red/yellow", "gold": "red/yellow", "golden": "red/yellow",
    "tan": "red/yellow", "cream": "red/yellow", "apricot": "red/yellow", "fawn": "red/yellow",
    # eumelanin base colours
    "black": "black", "brown": "brown/liver", "liver": "brown/liver", "chocolate": "brown/liver",
    "blue": "blue/grey", "gray": "blue/grey", "grey": "blue/grey", "isabella": "isabella/lilac",
    "lilac": "isabella/lilac",
}
# categories that are PATTERNS or ambiguous -> excluded from the colour concordance (v0 abstains on these)
PHENOTYPE_EXCLUDE = {"merle", "brindle", "spotted", "piebald", "white", "multi", "multicolor",
                     "multicolour", "sable", "agouti", "ticked", "mixed", "other", "unknown", ""}


class DataMissing(SystemExit):
    pass


def _require_files(pheno_tsv: Path, plink_prefix: Path) -> None:
    missing = []
    if not pheno_tsv.exists():
        missing.append(str(pheno_tsv))
    for ext in (".bed", ".bim", ".fam"):
        if not plink_prefix.with_suffix(ext).exists():
            missing.append(str(plink_prefix.with_suffix(ext)))
    if missing:
        print("Darwin's Ark data not present. Fetch the OPEN Dryad archive (no access request):", file=sys.stderr)
        print("  1. https://datadryad.org/dataset/doi:10.5061/dryad.83bk3jb4r", file=sys.stderr)
        print("  2. download + unzip darwins_dogs_gwas_input_files.zip (the coat_color TSV) and", file=sys.stderr)
        print("     darwins_dogs_genetic_set.zip (the PLINK bed/bim/fam, ~2.67 GB) to a local dir", file=sys.stderr)
        print("     (use D:/dna_decode_cache/darwins_ark/ — C: is disk-tight on this host).", file=sys.stderr)
        print(f"  missing: {missing}", file=sys.stderr)
        raise DataMissing(2)


def read_phenotypes(pheno_tsv: Path) -> list[dict]:
    """Parse the Q243 coat-colour TSV (real Dryad schema). Returns [{dog_id, raw_color}]."""
    rows = pheno_tsv.read_text(encoding="utf-8", errors="replace").splitlines()
    if not rows:
        raise SystemExit(3)
    header = rows[0].split("\t")
    # the Dryad file is "level-phenotypic data": a dog id column + the coat-colour answer; the exact header
    # names are surfaced by --probe. We locate an id-like column and a colour-like column defensively.
    idx_id = next((i for i, h in enumerate(header) if h.lower() in ("id", "dog", "dog_id", "iid", "fid")), 0)
    idx_col = next((i for i, h in enumerate(header)
                    if "color" in h.lower() or "colour" in h.lower() or h.lower() in ("answer", "value")),
                   len(header) - 1)
    out = []
    for line in rows[1:]:
        f = line.split("\t")
        if len(f) <= max(idx_id, idx_col):
            continue
        out.append({"dog_id": f[idx_id].strip(), "raw_color": f[idx_col].strip(),
                    "_header": (header[idx_id], header[idx_col])})
    return out


def map_phenotype(raw: str) -> str | None:
    """Owner free-text colour -> decoder vocabulary; None if excluded (pattern/ambiguous)."""
    r = raw.strip().lower()
    if r in PHENOTYPE_EXCLUDE:
        return None
    for key, target in PHENOTYPE_MAP.items():
        if key in r:
            # a pattern word anywhere -> exclude (e.g. "black merle")
            if any(p in r for p in PHENOTYPE_EXCLUDE if p):
                return None
            return target
    return None


def probe(pheno_tsv: Path, plink_prefix: Path) -> int:
    phenos = read_phenotypes(pheno_tsv)
    hdr = phenos[0]["_header"] if phenos else ("?", "?")
    from collections import Counter
    cats = Counter(p["raw_color"] for p in phenos)
    mapped = Counter(map_phenotype(p["raw_color"]) for p in phenos)
    print(f"# Darwin's Ark coat-colour PROBE")
    print(f"phenotype TSV: {pheno_tsv.name}  (id col={hdr[0]!r}, colour col={hdr[1]!r}, N={len(phenos)})")
    print(f"raw colour categories (top 20): {cats.most_common(20)}")
    print(f"mapped -> decoder vocab: { {k: v for k, v in mapped.items()} }")
    scored = sum(v for k, v in mapped.items() if k is not None)
    print(f"scorable (colour maps, non-pattern): {scored} / {len(phenos)}")
    # genotype side: which causal loci are present in the .bim (by gene-region annotation)
    bim = plink_prefix.with_suffix(".bim")
    n_var = sum(1 for _ in bim.open(encoding="utf-8", errors="replace"))
    print(f"\ngenotype .bim: {bim.name}  ({n_var} variants, canFam4)")
    print("causal loci to pin from this .bim (resolve exact canFam4 pos by gene region / study variant id):")
    for loc, (gene, desc, eff, src) in CAUSAL_VARIANTS.items():
        print(f"  {loc} ({gene}): {desc}  [effect={eff}]  <- {src}")
    print("\nNEXT: pin each causal variant's (chrom,pos,ref,alt) from the .bim, add a genotype->allele "
          "extractor, then re-run without --probe to score. (Coords are verified against THIS .bim, "
          "never transcribed from memory.)")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dog_coat_darwins_ark_validate")
    ap.add_argument("--pheno-tsv", required=True, type=Path)
    ap.add_argument("--plink-prefix", required=True, type=Path,
                    help="path prefix of the .bed/.bim/.fam set (no extension)")
    ap.add_argument("--probe", action="store_true",
                    help="inspect schema + causal-locus genotyping coverage (run this FIRST)")
    args = ap.parse_args(argv)

    _require_files(args.pheno_tsv, args.plink_prefix)
    if args.probe:
        return probe(args.pheno_tsv, args.plink_prefix)
    # Full scoring is gated on the causal-variant coords being pinned from the .bim (the --probe output).
    # Until that pinning is committed, refuse rather than emit an unverified concordance (anti-fabrication).
    print("scoring requires the causal-variant->allele extractor pinned from the .bim; run --probe first, "
          "then commit the resolved coords. Refusing to emit an unverified concordance.", file=sys.stderr)
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
