"""Single-defining-SNP -> actionable-allele concordance vs GeT-RM, for the single-SNP PGx cells.

SLCO1B1 / CYP4F2 (and UGT1A1's promoter tag) are SINGLE-SNP cells: their whole job is "does this one
defining SNP correctly detect the actionable star-allele?" -- not a full star-allele diplotype call. The
star-allele concordance harness (pgx_getrm_concordance.py) doesn't fit them. This scores exactly that
question against the free GeT-RM Consolidated truth (tests/data/pgx_getrm/getrm_consolidated_truth.tsv):

  for each sample: GeT-RM diplotype -> actionable-allele DOSAGE (0/1/2, i.e. how many of the two haplotypes
  carry an allele in the SNP's actionable set) vs our SNP's genomic ALT dosage (0/1/2 from the 1000G VCF).

The actionable sets are grounded in PharmVar/CPIC (which star alleles the SNP tags), NOT guessed:
  SLCO1B1 rs4149056 (521T>C) -> {*5,*15,*17}   (521C decreased-function carriers)
  CYP4F2  rs2108622 (V433M)  -> {*3}
  UGT1A1  rs887829  (promoter *28-LD-tag) -> {*28,*37}   +  rs4148323 (G71R) -> {*6}

Ambiguous GeT-RM truth (parenthetical `(*37)` = uncertain) is SKIPPED (conservative, mirrors the star
harness). Reports dosage-exact concordance + carrier (any-vs-none) concordance. Honest: this validates the
CELL'S defining-SNP CALL against the accepted star-allele truth -- a single-SNP proxy, not a full caller.

    uv run python scripts/pgx_single_snp_concordance.py --gene slco1b1
"""
from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TRUTH = REPO / "tests" / "data" / "pgx_getrm" / "getrm_consolidated_truth.tsv"

GENE_SNPS: dict[str, list[dict]] = {
    "slco1b1": [{"rsid": "rs4149056", "chrom": "12", "pos": 21178615, "ref": "T", "alt": "C",
                 "actionable": {"*5", "*15", "*17"}, "label": "*5-family (521C, decreased)"}],
    "cyp4f2":  [{"rsid": "rs2108622", "chrom": "19", "pos": 15879621, "ref": "C", "alt": "T",
                 "actionable": {"*3"}, "label": "*3 (V433M)"}],
    "ugt1a1":  [{"rsid": "rs887829", "chrom": "2", "pos": 233759924, "ref": "C", "alt": "T",
                 "actionable": {"*28", "*37"}, "label": "*28/*37 (reduced; irinotecan)"},
                {"rsid": "rs4148323", "chrom": "2", "pos": 233760498, "ref": "G", "alt": "A",
                 "actionable": {"*6"}, "label": "*6 (G71R; EAS)"}],
}


def _norm_chrom(c: str) -> str:
    return c[3:] if str(c).lower().startswith("chr") else str(c)


def actionable_dosage(diplotype: str, actionable: set[str]) -> int | None:
    """0/1/2 = how many haplotypes carry an actionable star. None if truth is ambiguous (parenthetical)."""
    if "(" in diplotype or " or " in diplotype:
        return None
    haps = diplotype.split("/")
    if len(haps) != 2:
        return None
    dose = 0
    for hap in haps:
        stars = {t.strip() for t in re.split(r"[+]", hap) if t.strip().startswith("*")}
        if stars & actionable:
            dose += 1
    return dose


def vcf_alt_dosage(vcf: Path, chrom: str, pos: int, alt: str) -> dict[str, int]:
    """{sample: ALT copy count (0/1/2)} at the SNP site."""
    samples: list[str] = []
    out: dict[str, int] = {}
    for line in vcf.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 10 or not cols[1].isdigit():
            continue
        if _norm_chrom(cols[0]) != _norm_chrom(chrom) or int(cols[1]) != pos:
            continue
        alts = cols[4].split(",")
        ai = alts.index(alt) + 1 if alt in alts else -1
        if ai < 0:
            continue
        fmt = cols[8].split(":")
        gi = fmt.index("GT") if "GT" in fmt else 0
        for s, cell in zip(samples, cols[9:]):
            gt = cell.split(":")[gi]
            nums = [int(a) for a in gt.replace("|", "/").split("/") if a.isdigit()]
            out[s] = sum(1 for n in nums if n == ai)
    return out


def score_gene(gene: str) -> dict:
    vcf = REPO / "data" / "pgx_1000g" / f"{gene}_1000g.vcf"
    truth_gene = gene.upper()
    truth = {}
    for line in TRUTH.read_text(encoding="utf-8").splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 3 and f[1] == truth_gene:
            truth[f[0]] = f[2]
    snp_results = []
    for snp in GENE_SNPS[gene]:
        dose = vcf_alt_dosage(vcf, snp["chrom"], snp["pos"], snp["alt"])
        n = dosage_hit = carrier_hit = ambiguous = 0
        mism = []
        for sample, dip in truth.items():
            if sample not in dose:
                continue
            truth_dose = actionable_dosage(dip, snp["actionable"])
            if truth_dose is None:
                ambiguous += 1
                continue
            n += 1
            our = dose[sample]
            if our == truth_dose:
                dosage_hit += 1
            else:
                mism.append({"sample": sample, "truth": dip, "truth_dose": truth_dose, "our_alt": our})
            if (our > 0) == (truth_dose > 0):
                carrier_hit += 1
        snp_results.append({
            "rsid": snp["rsid"], "actionable": sorted(snp["actionable"]), "label": snp["label"],
            "n_scored": n, "ambiguous_skipped": ambiguous,
            "dosage_concordance": round(dosage_hit / n, 4) if n else None, "dosage_hits": f"{dosage_hit}/{n}",
            "carrier_concordance": round(carrier_hit / n, 4) if n else None, "carrier_hits": f"{carrier_hit}/{n}",
            "mismatches": mism[:10],
        })
    return {"gene": truth_gene, "method": "single_defining_SNP_vs_GeT-RM_actionable_allele",
            "date": datetime.date.today().isoformat(), "snps": snp_results,
            "honesty": "single-SNP CELL call vs accepted GeT-RM star-allele truth; a defining-SNP proxy, "
                       "not a full star-allele caller. Ambiguous (parenthetical) truth skipped."}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="pgx_single_snp_concordance")
    ap.add_argument("--gene", required=True, choices=list(GENE_SNPS))
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)
    res = score_gene(args.gene)
    out = REPO / "wiki" / f"pgx_single_snp_concordance_{args.gene}_2026-07-28.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print(f"# {res['gene']} single-SNP -> GeT-RM actionable-allele concordance")
        for s in res["snps"]:
            print(f"  {s['rsid']} -> {s['label']}: dosage {s['dosage_hits']} ({s['dosage_concordance']}) | "
                  f"carrier {s['carrier_hits']} ({s['carrier_concordance']}) | ambiguous-skipped {s['ambiguous_skipped']}")
    print(f"[-> {out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
