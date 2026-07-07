#!/usr/bin/env python
"""Run the deterministic PGx decoder on a REAL PGP-UK individual (GRCh37 VCF) — no whole-genome liftOver.

PGP-UK (Personal Genome Project UK) publishes FREE, open-consent, non-application-gated individual human
VCFs in ENA (`PRJEB17529`). They are **GRCh37/hg19**; the PGx catalog is GRCh38 and matches by position, so
a direct run returns all-INDETERMINATE (the documented gotcha, see wiki/pgp_uk_realization_handoff.md).

This does a TARGETED position-liftover instead of a 100 MB whole-VCF liftOver: rsIDs are build-stable, so we
know each defining/sentinel variant's GRCh37 position (resolved once via Ensembl GRCh37 REST, hardcoded +
provenance-stamped below), extract only those ~22 genotypes from the individual's VCF, relabel them to their
GRCh38 catalog positions, and run the normal `dna_decode.pgx` callers. No dep-install, no Docker, no chain
file. First real-people deterministic PGx decode on a free non-gated cohort (teed up by the mosfaer session's
handoff to DNA-11).

HONEST FRAMING: this is a DEPLOYMENT/robustness demonstration — the decoder runs end-to-end on an arbitrary
real-world individual VCF from an INDEPENDENT source (different sequencing pipeline than the 1000G/GeT-RM
benchmark) and produces sane CPIC calls with the same honest abstention. PGP-UK ships no GeT-RM truth, so
this is NOT a new accuracy-vs-truth number (that lives in the GeT-RM concordance cells).

Usage:
    uv run python scripts/pgx_decode_pgp_uk.py --vcf D:/dna_decode_cache/pgp_uk/FR07961000.vcf.gz --sample-id FR07961000
"""
from __future__ import annotations

import argparse
import gzip
import json
import tempfile
from pathlib import Path

from dna_decode.pgx import (
    cyp2b6_catalog as c6,
    cyp2c8_catalog as c8,
    cyp2c9_catalog as c9,
    cyp2c19_catalog as c19,
    cyp2d6_catalog as c2d6,
    cyp3a5_catalog as c3,
    dpyd_catalog as dp,
    slco1b1,
    tpmt_catalog as tp,
    vkorc1,
)
from dna_decode.pgx.runner import (
    call_cyp2b6,
    call_cyp2c8,
    call_cyp2c19,
    call_cyp2c9,
    call_cyp2d6,
    call_cyp3a5,
    call_dpyd,
    call_tpmt,
)

# GRCh37 positions for every defining + sentinel PGx variant, resolved via Ensembl GRCh37 REST 2026-07-05
# and cross-checked against the GRCh38 catalog coords (rsIDs are build-stable). Hardcoded for offline use.
GRCH37_POS: dict[str, int] = {
    "rs4244285": 96541616, "rs4986893": 96540410, "rs12248560": 96521657,   # CYP2C19 *2/*3/*17
    "rs28399504": 96522463, "rs12769205": 96535124,                          # CYP2C19 sentinels *4/*35
    "rs1799853": 96702047, "rs1057910": 96741053,                            # CYP2C9 *2/*3
    "rs28371686": 96741058, "rs7900194": 96702066, "rs2256871": 96708974, "rs28371685": 96740981,  # CYP2C9 sentinels
    "rs11572103": 96818106, "rs11572080": 96827030, "rs1058930": 96818119,   # CYP2C8 *2/*3/*4
    "rs776746": 99270539, "rs10264272": 99262835, "rs41303343": 99250394,    # CYP3A5 *3/*6/*7
    "rs3745274": 41512841,                                                     # CYP2B6 *6 (516G>T)
    "rs1800460": 18139228, "rs1142345": 18130918,                            # TPMT *3B/*3C
    "rs9923231": 31107689,                                                     # VKORC1
    "rs4149056": 21331549,                                                     # SLCO1B1 *5
    # CYP2D6 (chr22) defining variants — GRCh37 positions resolved via Ensembl GRCh37 REST 2026-07-07,
    # cross-checked vs the GRCh38 catalog (uniform +396002 offset for SNVs / +396003 for the 3 indels).
    "rs3892097": 42524947, "rs1065852": 42526694, "rs16947": 42523943,        # CYP2D6 *4/*10/*2
    "rs1135840": 42522613, "rs28371706": 42525772, "rs28371725": 42523805,    # CYP2D6 486/*17/*41
    "rs59421388": 42523610, "rs769258": 42526763,                             # CYP2D6 *29/*35
    "rs35742686": 42524244, "rs5030655": 42525086, "rs5030656": 42524176,     # CYP2D6 *3/*6/*9 (indels)
    # DPYD (chr1) — the four CPIC-actionable fluoropyrimidine-toxicity haplotypes; GRCh37 positions
    # resolved via Ensembl GRCh37 REST 2026-07-07, cross-checked vs the Ensembl GRCh38 catalog coords.
    "rs3918290": 97915614, "rs55886062": 97981343,                            # DPYD *2A/*13 (no function)
    "rs67376798": 97547947, "rs75017182": 98045449,                           # DPYD c.2846A>T/HapB3 (decreased)
}


def _all_variants():
    """Yield (rsid, chrom, grch38_pos, ref, alt) for every catalog site we need to lift."""
    for cat in (c19, c9, c8, c3, c6):
        for d in cat.CORE_DEFINING:
            yield d.rsid, d.chrom, d.pos, d.ref, d.alt
        for s in getattr(cat, "SENTINELS", []):
            yield s.rsid, s.chrom, s.pos, s.ref, s.alt
    for d in tp.COMPONENTS:
        yield d.rsid, d.chrom, d.pos, d.ref, d.alt
    for d in c2d6.COMPONENTS:
        yield d.rsid, d.chrom, d.pos, d.ref, d.alt
    for d in dp.CORE_DEFINING:
        yield d.rsid, d.chrom, d.pos, d.ref, d.alt
    yield vkorc1.RSID, vkorc1.CHROM, vkorc1.POS, vkorc1.REF, vkorc1.ALT
    yield slco1b1.RSID, slco1b1.CHROM, slco1b1.POS, slco1b1.REF, slco1b1.ALT


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def extract_sites(vcf_gz: Path, needed_by_pos: dict[tuple[str, int], str]) -> tuple[dict[str, tuple], set]:
    """Stream a (possibly gzipped) GRCh37 VCF. Return ({rsid: (ref, alt, gt)} for needed sites,
    set-of-chromosomes-that-had-ANY-record) — the second value lets us tell 'person is ref here' (chrom
    covered, site absent -> assumed 0/0) from 'VCF doesn't cover this chrom' (never assume ref)."""
    opener = gzip.open if str(vcf_gz).endswith(".gz") else open
    want_chroms = {c for c, _ in needed_by_pos}
    found: dict[str, tuple] = {}
    chroms_seen: set = set()
    with opener(vcf_gz, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            i = line.find("\t")
            chrom = _norm_chrom(line[:i])
            if chrom not in want_chroms:
                continue
            chroms_seen.add(chrom)
            cols = line.rstrip("\n").split("\t")
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            rs = needed_by_pos.get((chrom, pos))
            if rs is None:
                continue
            gt = "."
            if len(cols) >= 10 and "GT" in cols[8].split(":"):
                gt = cols[9].split(":")[cols[8].split(":").index("GT")]
            found[rs] = (cols[3], cols[4], gt)
    return found, chroms_seen


def build_grch38_minivcf(found: dict[str, tuple], chroms_seen: set, sample_id: str) -> tuple[str, list]:
    """Relabel extracted GRCh37 genotypes to their GRCh38 catalog positions -> a mini single-sample VCF.
    SNP sites are carried when VCF REF/ALT are allele-concordant with the catalog (same-strand); an insertion
    (*7) is carried when the VCF ALT contains the catalog ALT (alt="*" sentinels carry the VCF's own REF/ALT).
    A site ABSENT-as-variant but on a chromosome the VCF COVERS is emitted 0/0 (a called WGS VCF omits ref
    sites -> absence = confident reference; dna-pgx flags this as assumed-reference). A site on a chromosome
    the VCF never touched is OMITTED (never assume ref -> the caller reports no_input, honestly)."""
    header = ("##fileformat=VCFv4.2\n"
              f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")
    rows, audit = [], []
    for rs, chrom, pos38, ref, alt in _all_variants():
        if rs not in found:
            if chrom in chroms_seen:
                # person is reference here (called WGS VCF omits ref sites) -> encode 0/0
                out_ref = "N" if alt == "*" else ref
                out_alt = "." if alt == "*" else alt  # sentinel with no variant -> no ALT to encode
                if alt == "*":
                    audit.append({"rsid": rs, "state": "ref_no_sentinel"})
                    continue  # a no-variant sentinel simply doesn't fire; omit
                rows.append(f"{chrom}\t{pos38}\t{rs}\t{out_ref}\t{out_alt}\t.\tPASS\t.\tGT\t0/0")
                audit.append({"rsid": rs, "state": "assumed_ref_0/0_wgs"})
            else:
                audit.append({"rsid": rs, "state": "chrom_not_covered_omitted"})
            continue
        vref, valt, gt = found[rs]
        valts = valt.split(",")
        concordant = (alt == "*") or (vref == ref and alt in valts)
        if not concordant:
            audit.append({"rsid": rs, "state": "allele_mismatch", "vcf": f"{vref}>{valt}", "catalog": f"{ref}>{alt}"})
            continue
        out_ref = vref if alt == "*" else ref
        out_alt = valt if alt == "*" else alt
        rows.append(f"{chrom}\t{pos38}\t{rs}\t{out_ref}\t{out_alt}\t.\tPASS\t.\tGT\t{gt}")
        audit.append({"rsid": rs, "state": "lifted", "gt": gt})
    return header + "\n".join(rows) + "\n", audit


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Deterministic PGx decode of a real PGP-UK (GRCh37) individual.")
    ap.add_argument("--vcf", type=Path, required=True, help="PGP-UK individual VCF (.vcf or .vcf.gz, GRCh37)")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--out", type=Path, default=None, help="write the provenance JSON here")
    args = ap.parse_args(argv)

    needed_by_pos = {}
    for rs, chrom, _p38, _ref, _alt in _all_variants():
        p37 = GRCH37_POS.get(rs)
        if p37 is not None:
            needed_by_pos[(chrom, p37)] = rs
    print(f"[pgp-uk] {args.sample_id}: extracting {len(needed_by_pos)} GRCh37 PGx sites from {args.vcf.name}...",
          flush=True)
    found, chroms_seen = extract_sites(args.vcf, needed_by_pos)
    mini_vcf, audit = build_grch38_minivcf(found, chroms_seen, args.sample_id)
    tmp = Path(tempfile.mktemp(suffix=".vcf"))
    tmp.write_text(mini_vcf, encoding="utf-8")

    genes = {"cyp2c19": call_cyp2c19, "cyp2c9": call_cyp2c9, "cyp2c8": call_cyp2c8,
             "cyp3a5": call_cyp3a5, "tpmt": call_tpmt, "cyp2b6": call_cyp2b6}
    results = {}
    for g, fn in genes.items():
        rec = fn(tmp, sample_id=args.sample_id)
        results[g] = {"diplotype": rec["diplotype"], "phenotype": rec.get("phenotype"),
                      "phenotype_abbrev": rec.get("phenotype_abbrev"),
                      "phenotype_status": rec["phenotype_status"], "phasing": rec["phasing"]}
    # CYP2D6 — the flagship pharmacogene. SNP-diplotype surface (activity-score phenotype). LOAD-BEARING
    # HONESTY: a SNP VCF cannot see the CYP2D6 structural alleles (*5 del / *xN dup / *13/*36/*68 hybrids) ->
    # the caller stamps cnv_hybrid_unassessed=True; the copy-number half is resolvable only from a BAM/CRAM
    # (PGP-UK ships called VCFs, not reads), so this is a SNP-proxy diplotype, honestly flagged.
    d6 = call_cyp2d6(tmp, sample_id=args.sample_id)
    results["cyp2d6"] = {"diplotype": d6["diplotype"], "phenotype": d6.get("phenotype"),
                         "phenotype_abbrev": d6.get("phenotype_abbrev"),
                         "phenotype_status": d6["phenotype_status"], "phasing": d6.get("phasing"),
                         "activity_score": d6.get("activity_score"),
                         "cnv_hybrid_unassessed": d6.get("cnv_hybrid_unassessed", True)}
    # DPYD — the fluoropyrimidine (5-FU/capecitabine) toxicity gene. CPIC activity-score over the four
    # actionable DPD-deficiency haplotypes. All-SNP, no structural blind spot (unlike CYP2D6). A DPYD IM/PM
    # call is clinically load-bearing (dose reduction / avoidance), so this extends the real-human coverage
    # to the highest-stakes actionable pharmacogene.
    dd = call_dpyd(tmp, sample_id=args.sample_id)
    results["dpyd"] = {"diplotype": dd["diplotype"], "phenotype": dd.get("phenotype"),
                       "phenotype_abbrev": dd.get("phenotype_abbrev"),
                       "phenotype_status": dd["phenotype_status"], "phasing": dd.get("phasing"),
                       "activity_score": dd.get("activity_score")}
    vk = vkorc1.call_vkorc1(tmp); vk["sample_id"] = args.sample_id
    results["vkorc1"] = {"genotype": vk["cdna_genotype"], "sensitivity": vk["warfarin_sensitivity"]}
    sl = slco1b1.call_slco1b1(tmp); sl["sample_id"] = args.sample_id
    results["slco1b1"] = {"genotype": sl["variant_genotype"], "function": sl["function"],
                          "myopathy_risk": sl["myopathy_risk"]}

    n_real = sum(1 for g in ("cyp2c19", "cyp2c9", "cyp2c8", "cyp2d6", "dpyd", "cyp3a5", "tpmt", "cyp2b6")
                 if results[g]["phenotype"] not in (None, "Indeterminate")
                 and results[g]["phenotype_status"] == "ok")
    out = {"schema": "pgx-pgp-uk-realization-v0", "sample_id": args.sample_id,
           "cohort": "PGP-UK (Personal Genome Project UK), ENA PRJEB17529, open-consent GRCh37",
           "vcf": str(args.vcf), "assembly_source": "GRCh37 -> targeted position-liftover -> GRCh38 catalog",
           "n_genes_real_call": n_real, "results": results, "lift_audit": audit,
           "honest_tier": ("DEPLOYMENT/robustness demonstration on a real independent-cohort individual VCF; "
                           "PGP-UK ships no GeT-RM truth so this is NOT an accuracy-vs-truth number.")}
    print(json.dumps({"sample_id": args.sample_id, "n_genes_real_call": n_real,
                      "results": results}, indent=2))
    if args.out:
        args.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[provenance -> {args.out}]")
    return 0 if n_real >= 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
