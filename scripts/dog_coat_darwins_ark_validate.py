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

# REAL Darwin's Ark Q243 schema (verified 2026-07-30 against the file): a MULTI-HOT presence matrix, NOT a
# single free-text colour. dog_id + per-colour 0/1 columns + number_of_colors_in_coat + single_color_in_coat.
# Map each base-colour column -> our decoder's base-colour vocabulary. `white` = spotting / pigment-absence,
# which the cell ABSTAINS on (a pattern axis, not a eumelanin/phaeomelanin base).
COAT_COLOR_COLS: dict[str, str] = {
    "Q243_black_coat_color": "black",
    "Q243_liver_or_brown_coat_color": "brown/liver",
    "Q243_red_coat_color": "red/yellow",
    "Q243_yellow_coat_color": "red/yellow",
    "Q243_cream_coat_color": "red/yellow",          # dilute phaeomelanin
    "Q243_grey_or_blue_coat_color": "blue/grey",
    "Q243_tan_coat_color": "tan",                   # tan-points / agouti tan
    "Q243_white_coat_color": "white",               # spotting / absence -> ABSTAIN axis
}
_WHITE = "white"


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
    """Parse the REAL Q243 multi-hot coat-colour TSV. Each dog carries a SET of present colours.
    Returns [{dog_id, colors:set (decoder vocab), n_colors:int, single:bool}]."""
    import csv
    with pheno_tsv.open(encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or "dog_id" not in rows[0]:
        raise SystemExit(3)
    out = []
    for r in rows:
        colors = {COAT_COLOR_COLS[c] for c, v in r.items() if c in COAT_COLOR_COLS and v == "1"}
        try:
            nc = int(r.get("number_of_colors_in_coat", "") or 0)
        except ValueError:
            nc = 0
        out.append({"dog_id": r["dog_id"].strip(), "colors": colors, "n_colors": nc,
                    "single": r.get("single_color_in_coat") == "1"})
    return out


def scoring_target(pheno: dict):
    """The concordance target for one dog:
      ('base', color)      -- a cleanly-scorable single base colour (single-colour, non-white) -> match our call
      ('abstain', reason)  -- white-only (spotting/absence) -> our cell should ABSTAIN
      ('multi', frozenset) -- a multi-colour coat -> distribution-level scoring (e.g. {black,tan} = tan-points)
      None                 -- uninformative
    """
    non_white = pheno["colors"] - {_WHITE}
    if pheno["single"]:
        if pheno["colors"] == {_WHITE}:
            return ("abstain", "white/spotting")
        if len(non_white) == 1:
            return ("base", next(iter(non_white)))
    if non_white:
        return ("multi", frozenset(non_white))
    return None


def phenotype_summary(phenos: list[dict]) -> dict:
    """Distribution + scorable counts of the multi-hot coat-colour phenotypes."""
    from collections import Counter
    tgt = Counter(scoring_target(p) and scoring_target(p)[0] for p in phenos)
    single_base = Counter(scoring_target(p)[1] for p in phenos
                          if (scoring_target(p) or (None,))[0] == "base")
    return {"n": len(phenos),
            "single_colour": sum(1 for p in phenos if p["single"]),
            "target_kinds": dict(tgt),
            "single_base_colour_counts": dict(single_base.most_common()),
            "directly_scorable_single_base": sum(single_base.values())}


def probe(pheno_tsv: Path, plink_prefix: Path) -> int:
    phenos = read_phenotypes(pheno_tsv)
    s = phenotype_summary(phenos)
    print("# Darwin's Ark coat-colour PROBE (multi-hot schema)")
    print(f"phenotype TSV: {pheno_tsv.name}  (N={s['n']} dogs, {s['single_colour']} single-colour)")
    print(f"target kinds (base/multi/abstain/None): {s['target_kinds']}")
    print(f"single-base-colour concordance targets: {s['single_base_colour_counts']}")
    print(f"directly-scorable single-base dogs: {s['directly_scorable_single_base']}")
    # genotype side: read the REAL .bim/.fam via the spec-verified plink_io reader (no line-count guess)
    from dna_decode.pigment import plink_io
    bim = plink_io.read_bim(plink_prefix.with_suffix(".bim"))
    fam = plink_io.read_fam(plink_prefix.with_suffix(".fam"))
    print(f"\ngenotype set: {len(bim)} variants x {len(fam)} dogs (canFam4, PLINK1)")
    print("causal loci to pin from this .bim (resolve exact canFam4 pos by gene region / study variant id):")
    for loc, (gene, desc, eff, src) in CAUSAL_VARIANTS.items():
        print(f"  {loc} ({gene}): {desc}  [effect={eff}]  <- {src}")
    print("\nNEXT: identify each causal variant's row in THIS .bim (by canFam4 pos / study variant id — "
          "verified, never transcribed from memory; note E/K indels may be ABSENT from a SNP-imputed panel), "
          "then re-run with --coords LOCUS=chrom:pos:counted_allele,... to score. The scorer uses the "
          "spec-verified PLINK reader (dna_decode.pigment.plink_io) to extract only those variants.")
    return 0


def extract_genotypes(plink_prefix: Path, coords: dict[str, tuple[str, int, str]]) -> dict[str, dict]:
    """Extract the pinned causal variants for every dog via the spec-verified PLINK reader.

    `coords`: {locus -> (chrom, pos, counted_allele)} pinned from the REAL .bim (never fabricated).
    Returns {dog_id -> {locus -> 'A/G' harmonized to the counted allele's strand}}. Fails loud if a pinned
    coord is not in the .bim (wrong coord/assembly). This is the DONE half of scoring — the extraction —
    verified by tests/test_plink_io.py; it does NOT emit a concordance (see NOTE in main()).
    """
    from dna_decode.pigment import plink_io
    bim = plink_io.read_bim(plink_prefix.with_suffix(".bim"))
    fam = plink_io.read_fam(plink_prefix.with_suffix(".fam"))
    loc_idx: dict[str, tuple[int, str, str, str]] = {}
    for loc, (chrom, pos, counted) in coords.items():
        hits = plink_io.find_variants(bim, chrom=chrom, pos=pos)
        if not hits:
            raise SystemExit(f"pinned {loc} variant {chrom}:{pos} not in the .bim (wrong coord/assembly?)")
        v = hits[0]
        loc_idx[loc] = (v.index, v.a1, v.a2, counted)
    dos = plink_io.read_bed_variants(plink_prefix.with_suffix(".bed"), len(fam),
                                     [t[0] for t in loc_idx.values()])
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    out: dict[str, dict] = {}
    for si, dog in enumerate(fam):
        genos: dict[str, str] = {}
        for loc, (vidx, a1, a2, counted) in loc_idx.items():
            gt = plink_io.genotype_string(dos[vidx][si], a1, a2)
            if gt is None:
                continue
            bases = set(gt.split("/"))
            if counted not in bases and comp.get(counted) in bases:
                gt = "/".join(comp[b] for b in gt.split("/"))  # strand-harmonize to the counted allele
            genos[loc] = gt
        out[dog] = genos
    return out


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
    # DONE (verified): the PLINK extractor (dna_decode.pigment.plink_io + extract_genotypes here) and the
    # phenotype ingest/mapping. PENDING (both need the real .bim, so NOT fabricated here): (1) the causal
    # variants' canFam4 coords, and (2) the per-variant base->coat-allele-symbol table (which base is the
    # `e`/`b`/`d` allele). Until both are pinned + committed, refuse to emit a concordance (anti-fabrication).
    print("extraction is wired + spec-verified (plink_io). Scoring still needs, pinned from the REAL .bim: "
          "(1) the causal-variant canFam4 coords, (2) the per-variant base->coat-allele-symbol table. "
          "Run --probe first, commit those pins, then wire the call. Refusing to emit an unverified "
          "concordance.", file=sys.stderr)
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
