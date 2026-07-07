"""THE genuine GeT-RM independent number: CYP2C19 caller vs the GeT-RM NGS consensus on real 1000G genomes.

Truth labels: the GeT-RM NGS consensus diplotypes (consensus of Astrolabe + Stargazer + Aldy; Gaedigk 2022
J Mol Diagn), as published in the ursaPGx benchmark table `star-allele-comparison_common.tsv`
(coriell-research/ursaPGx, column CYP2C19_getrm_ngs) for the ~87 samples overlapping the 1000 Genomes 30x
panel by Coriell ID. Genotypes: the public 1000G 30x phased VCF region (already fetched via Docker bcftools
-> data/pgx_1000g/cyp2c19_1000g.vcf.gz). Our caller is INDEPENDENT of the three consensus tools.

Honest tier: this is the strongest validation the project can give star-allele CALLING -- vs the field's
accepted GeT-RM consensus truth set, on real population genomes, with an independent caller. (The consensus
is itself caller-derived, so it is "vs the accepted reference," a notch below a wet-lab phenotype, but it is
the standard star-allele benchmark.) The v0 caller covers the CORE SNP set (*1/*2/*3/*17); samples whose
GeT-RM consensus uses a NON-CORE allele are scored SEPARATELY -- the v0.1 sentinel layer should WITHHOLD
them, not mis-call. The headline is core-comparable exact-diplotype concordance.
"""
from __future__ import annotations

import csv
import datetime
import json
import re
from pathlib import Path

from dna_decode.pgx.caller import call_diplotype

REPO = Path(__file__).resolve().parent.parent
# committed truth set (vendored from ursaPGx); fall back to the gitignored fetch dir if present
TRUTH = REPO / "tests" / "data" / "pgx_getrm" / "star-allele-comparison_common.tsv"
if not TRUTH.exists():
    TRUTH = REPO / "data" / "pgx_getrm" / "star-allele-comparison_common.tsv"


def _c9_caller(plain_vcf, sample):
    from dna_decode.pgx import cyp2c9_catalog as c9
    return call_diplotype(plain_vcf, sample=sample, defining=c9.CORE_DEFINING, sentinels=c9.SENTINELS,
                          reference_allele=c9.REFERENCE_ALLELE, phenotype_fn=c9.diplotype_phenotype,
                          gene=c9.GENE)


def _c8_caller(plain_vcf, sample):
    from dna_decode.pgx import cyp2c8_catalog as c8
    return call_diplotype(plain_vcf, sample=sample, defining=c8.CORE_DEFINING, sentinels=c8.SENTINELS,
                          reference_allele=c8.REFERENCE_ALLELE, phenotype_fn=c8.diplotype_phenotype,
                          gene=c8.GENE)


def _c3a5_caller(plain_vcf, sample):
    from dna_decode.pgx import cyp3a5_catalog as c3
    return call_diplotype(plain_vcf, sample=sample, defining=c3.CORE_DEFINING, sentinels=c3.SENTINELS,
                          reference_allele=c3.REFERENCE_ALLELE, phenotype_fn=c3.diplotype_phenotype,
                          gene=c3.GENE)


def _tpmt_caller(plain_vcf, sample):
    from dna_decode.pgx import tpmt_catalog as tp
    from dna_decode.pgx.compound_caller import assemble_compound_diplotype
    return assemble_compound_diplotype(plain_vcf, tp.COMPONENTS, tp.COMPOUND_RULES,
                                       reference_allele=tp.REFERENCE_ALLELE,
                                       phenotype_fn=tp.diplotype_phenotype, gene=tp.GENE, sample=sample)


def _c2b6_caller(plain_vcf, sample):
    from dna_decode.pgx import cyp2b6_catalog as c6
    return call_diplotype(plain_vcf, sample=sample, defining=c6.CORE_DEFINING, sentinels=c6.SENTINELS,
                          reference_allele=c6.REFERENCE_ALLELE, phenotype_fn=c6.diplotype_phenotype,
                          gene=c6.GENE)


def _c2d6_caller(plain_vcf, sample):
    from dna_decode.pgx import cyp2d6_catalog as c2d6
    from dna_decode.pgx.cyp2d6_caller import assemble_cyp2d6_diplotype
    return assemble_cyp2d6_diplotype(plain_vcf, c2d6.COMPONENTS, c2d6.STAR_PRIORITY,
                                     reference_allele=c2d6.REFERENCE_ALLELE,
                                     phenotype_fn=c2d6.diplotype_phenotype, gene=c2d6.GENE, sample=sample)


# Per-gene config: which 1000G region VCF, which truth column, the core SNP set, *38-equivalence, caller.
# *38 is the TRUE variant-free reference (NORMAL function, phenotype==*1; rs3758581 distinguishes, but that
# is phenotype-irrelevant) -> a GeT-RM *38 scores phenotype-equivalent to *1. (CYP2C9 has no *38 equivalent.)
GENES = {
    "cyp2c19": {"vcf": REPO / "data" / "pgx_1000g" / "cyp2c19_1000g.vcf.gz",
                "truth_col": "CYP2C19_getrm_ngs", "core": {"*1", "*2", "*3", "*17"},
                "ref_equiv": {"*38": "*1"}, "caller": lambda v, s: call_diplotype(v, sample=s),
                "out_stem": "pgx_getrm_concordance_2026-06-25", "fetch_note": "Docker bcftools"},
    "cyp2c9":  {"vcf": REPO / "data" / "pgx_1000g" / "cyp2c9_1000g.vcf.gz",
                "truth_col": "CYP2C9_getrm_ngs", "core": {"*1", "*2", "*3"},
                "ref_equiv": {}, "caller": _c9_caller,
                "out_stem": "pgx_getrm_concordance_cyp2c9_2026-06-25", "fetch_note": "Docker bcftools"},
    "cyp2c8":  {"vcf": REPO / "data" / "pgx_1000g" / "cyp2c8_1000g.vcf.gz",
                "truth_col": "CYP2C8_getrm_ngs", "core": {"*1", "*2", "*3", "*4"},
                "ref_equiv": {}, "caller": _c8_caller,
                "out_stem": "pgx_getrm_concordance_cyp2c8_2026-07-05",
                "fetch_note": "pure-Python tabix-over-HTTP (scripts/fetch_1000g_region.py; no Docker)"},
    "cyp3a5":  {"vcf": REPO / "data" / "pgx_1000g" / "cyp3a5_1000g.vcf.gz",
                "truth_col": "CYP3A5_getrm_cons", "core": {"*1", "*3", "*6", "*7"},
                "ref_equiv": {}, "caller": _c3a5_caller,
                "truth_file": REPO / "tests" / "data" / "pgx_getrm" / "getrm_cyp3a5_consensus.tsv",
                "out_stem": "pgx_getrm_concordance_cyp3a5_2026-07-05",
                "fetch_note": "pure-Python tabix-over-HTTP (scripts/fetch_1000g_region.py; no Docker); GeT-RM CDC CYP3A4/5 table"},
    "tpmt":    {"vcf": REPO / "data" / "pgx_1000g" / "tpmt_1000g.vcf.gz",
                "truth_col": "TPMT_getrm_cons", "core": {"*1", "*3A", "*3B", "*3C"},
                "ref_equiv": {}, "caller": _tpmt_caller,
                "truth_file": REPO / "tests" / "data" / "pgx_getrm" / "getrm_tpmt_consensus.tsv",
                "out_stem": "pgx_getrm_concordance_tpmt_2026-07-05",
                "fetch_note": "pure-Python tabix-over-HTTP (no Docker); GeT-RM CDC consolidated PGx table; COMPOUND *3A caller"},
    "cyp2b6":  {"vcf": REPO / "data" / "pgx_1000g" / "cyp2b6_1000g.vcf.gz",
                "truth_col": "CYP2B6_getrm_cons", "core": {"*1", "*6"},
                "ref_equiv": {}, "caller": _c2b6_caller,
                "truth_file": REPO / "tests" / "data" / "pgx_getrm" / "getrm_cyp2b6_consensus.tsv",
                "out_stem": "pgx_getrm_concordance_cyp2b6_2026-07-05",
                "fetch_note": "pure-Python tabix-over-HTTP (no Docker); GeT-RM CDC consolidated PGx table; SINGLE-SNP *6-proxy (785 absent from panel)"},
    "cyp2d6":  {"vcf": REPO / "data" / "pgx_1000g" / "cyp2d6_1000g.vcf",
                "truth_col": "CYP2D6_getrm_cons",
                "core": {"*1", "*2", "*3", "*4", "*6", "*9", "*10", "*17", "*29", "*35", "*41"},
                "ref_equiv": {}, "caller": _c2d6_caller,
                "out_stem": "pgx_getrm_concordance_cyp2d6_2026-07-06",
                "fetch_note": "pure-Python tabix-over-HTTP (scripts/fetch_1000g_region.py; no Docker); "
                              "SNP surface only — structural alleles BAM-required + EXCLUDED"},
}
REF_EQUIV = GENES["cyp2c19"]["ref_equiv"]   # back-compat alias (CYP2C19 *38==*1)

# --- CYP2D6 tiered-honest truth classification (the last-major-pharmacogene structural surface) ---
_CYP2D6_STRUCTURAL_STARS = {"*5", "*13", "*36", "*61", "*63", "*68"}
_CYP2D6_CORE = {"*1", "*2", "*3", "*4", "*6", "*9", "*10", "*17", "*29", "*35", "*41"}


def _cyp2d6_stars(text: str) -> list[str]:
    """Extract primary star tokens (strip suballele letters + parentheses) from one truth string."""
    out = []
    for tok in re.split(r"[/|+()]", text):
        m = re.match(r"\s*(\*\d+)", tok)
        if m:
            out.append(m.group(1))
    return out


def _classify_cyp2d6_truth(raw: str) -> tuple[str, tuple[str, ...]]:
    """Classify a CYP2D6 GeT-RM truth string into a tier + normalized primary-star tuple.

    Tiers (brainstorm: raw_truth kept separately; NEVER collapse ambiguity into a match):
      structural  — gene deletion *5, duplication *xN/x2, or CYP2D6-CYP2D7 hybrid *13/*36/*61/*63/*68.
                    NOT VCF-decodable -> EXCLUDED from the scored denominator (may be silently mis-called).
      ambiguous   — a parenthetical ALTERNATIVE annotation (e.g. '*2 (*35)') -> excluded (uncertain truth).
      noncore_snp — a SNP allele outside the core set (*14/*15/*21/*40/*46) -> residual (mis-called).
      core_snp    — every allele in the SNP core set -> the scored tier.
    """
    s = raw.strip()
    stars = _cyp2d6_stars(s)
    if re.search(r"[xX]\d", s) or any(st in _CYP2D6_STRUCTURAL_STARS for st in stars):
        return "structural", tuple(stars)
    if "(" in s:   # parenthetical alternative -> genuinely ambiguous truth; do NOT score
        return "ambiguous", tuple(stars)
    def key(a):
        d = "".join(ch for ch in a.lstrip("*") if ch.isdigit())
        return (int(d) if d else 999, a)
    norm = tuple(sorted(stars, key=key))
    if all(st in _CYP2D6_CORE for st in stars):
        return "core_snp", norm
    return "noncore_snp", norm


def _run_cyp2d6(cfg) -> int:
    """Dedicated tiered-honest CYP2D6 concordance (isolated from the 7 other genes' shared path)."""
    VCF = cfg["vcf"]
    if VCF.suffix == ".gz":
        VCF = _gunzip_to_plain(VCF)
    elif VCF.with_suffix(VCF.suffix + ".gz").exists() and not VCF.exists():
        VCF = _gunzip_to_plain(VCF.with_suffix(VCF.suffix + ".gz"))
    TRUTH_F = cfg.get("truth_file", TRUTH)
    if not VCF.exists() or not TRUTH_F.exists():
        print(f"ERROR: need {VCF} + {TRUTH_F} (fetch the CYP2D6 region first)")
        return 2

    samples_in_vcf = set()
    for line in VCF.read_text(encoding="utf-8").splitlines():
        if line.startswith("#CHROM"):
            samples_in_vcf = set(line.rstrip("\n").split("\t")[9:])
            break

    from dna_decode.pgx.cyp2d6_catalog import diplotype_phenotype
    caller = cfg["caller"]
    rows = []
    with open(TRUTH_F, encoding="utf-8") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            cor = rec["Coriell"].strip()
            raw = (rec.get(cfg["truth_col"], "") or "").strip()
            if not raw or raw in ("NA", ".") or cor not in samples_in_vcf:
                continue
            kind, norm = _classify_cyp2d6_truth(raw)
            r = caller(VCF, cor)
            pred = _norm(r.diplotype or "")
            match = (pred == norm) if kind == "core_snp" else None
            # phenotype-level agreement (secondary; *35==*2==*1 are all normal-function -> forgiving view)
            pheno_match = None
            if kind == "core_snp" and len(norm) == 2 and r.allele1:
                pheno_match = diplotype_phenotype(*norm) == diplotype_phenotype(r.allele1, r.allele2)
            # honest diagnosis of a core mis-call: predicted HOMOZYGOUS at a defining SNP while truth is a
            # het-with-*1 is the SIGNATURE of hidden structural confounding (a *4 tandem/dup makes the
            # defining SNP read 1|1) -> the cnv_hybrid_unassessed caveat in action, NOT a SNP-logic bug.
            diag = None
            if kind == "core_snp" and match is False and r.allele1 == r.allele2 and "*1" in norm:
                diag = "likely_structural_confound (predicted homozygous; truth het-with-*1 -> hidden CNV/hybrid; cnv_hybrid_unassessed)"
            rows.append({"sample": cor, "raw_truth": raw, "normalized_truth": "/".join(norm),
                         "predicted": "/".join(pred) if pred else None, "tier": kind,
                         "match": match, "phenotype_match": pheno_match, "miscall_diagnosis": diag,
                         "phenotype_status": r.phenotype_status, "phasing": r.phasing})

    core = [r for r in rows if r["tier"] == "core_snp"]
    noncore = [r for r in rows if r["tier"] == "noncore_snp"]
    structural = [r for r in rows if r["tier"] == "structural"]
    ambiguous = [r for r in rows if r["tier"] == "ambiguous"]
    core_hits = sum(1 for r in core if r["match"])
    pheno_hits = sum(1 for r in core if r["phenotype_match"])
    concordance = round(core_hits / len(core), 4) if core else None

    rep = {
        "schema": "pgx-cyp2d6-getrm-concordance-v0",
        "gene": "CYP2D6",
        "analysis_date": datetime.date.today().isoformat(),
        "truth_source": ("GeT-RM consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) via ursaPGx "
                         "star-allele-comparison_common.tsv, column CYP2D6_getrm_cons"),
        "genotype_source": f"1000 Genomes 30x phased panel (CYP2D6 chr22 region, {cfg['fetch_note']})",
        "caller_is_independent_of_consensus_tools": True,
        "is_snp_surface_only": True,
        "n_overlap_samples": len(rows),
        # TIERED denominators (brainstorm #2/#4 — never one inflated number):
        "core_snp_n": len(core),
        "noncore_snp_n": len(noncore),
        "structural_excluded_n": len(structural),
        "ambiguous_truth_excluded_n": len(ambiguous),
        "core_snp_diplotype_concordance": concordance,
        "core_snp_diplotype_hits": f"{core_hits}/{len(core)}",
        "core_snp_phenotype_concordance": (round(pheno_hits / len(core), 4) if core else None),
        "core_snp_phenotype_hits": f"{pheno_hits}/{len(core)}",
        "core_miscall_diagnoses": [{"sample": r["sample"], "raw_truth": r["raw_truth"],
                                    "predicted": r["predicted"], "diagnosis": r["miscall_diagnosis"]}
                                   for r in core if r["match"] is False],
        "structural_note": ("structural alleles (*5 deletion / *xN duplication / *13/*36/*61/*63/*68 "
                            "CYP2D6-CYP2D7 hybrids) are NOT VCF-decodable -> EXCLUDED from the scored "
                            "denominator; they are NOT withheld and may be SILENTLY MIS-CALLED "
                            "(cnv_hybrid_unassessed). Full typing needs a BAM/CRAM + Cyrius-class caller."),
        "honesty_tier": ("GeT-RM CONSENSUS core-diplotype concordance on the SNP-DECODABLE subset, "
                         "independent caller. Structural + ambiguous truth EXCLUDED (tiered denominators; "
                         "raw + normalized truth both retained). SNP surface only — NOT full CYP2D6 typing."),
        "core_rows": core, "noncore_rows": noncore,
        "structural_rows": structural, "ambiguous_rows": ambiguous,
    }
    out_json = REPO / "wiki" / f"{cfg['out_stem']}.json"
    out_md = REPO / "wiki" / f"{cfg['out_stem']}.md"
    out_json.write_text(json.dumps(rep, indent=2), encoding="utf-8")

    core_set = "/".join(sorted(_CYP2D6_CORE, key=lambda a: (int("".join(c for c in a if c.isdigit()) or 999), a)))
    L = [f"# CYP2D6 caller vs GeT-RM consensus on real 1000G ({rep['analysis_date']})", "",
         f"**Truth:** {rep['truth_source']}", f"**Genotypes:** {rep['genotype_source']}",
         "**Surface:** SNP-decodable star alleles ONLY (structural alleles BAM-required; see note).", "",
         f"- Overlap samples: **{rep['n_overlap_samples']}**  (tiered below — no single inflated denominator)",
         f"- **Core-SNP diplotype concordance: {rep['core_snp_diplotype_hits']} "
         f"({rep['core_snp_diplotype_concordance']})**  (truth in {core_set})",
         f"- Core-SNP PHENOTYPE concordance: {rep['core_snp_phenotype_hits']} "
         f"({rep['core_snp_phenotype_concordance']})  (*35==*2==*1 all normal-function)",
         f"- Non-core SNP alleles (residual, mis-called; *14/*15/*21/*40/*46): **{rep['noncore_snp_n']}**",
         f"- Structural EXCLUDED (BAM-required; *5/*xN/*13/*36/*68): **{rep['structural_excluded_n']}**",
         f"- Ambiguous-truth EXCLUDED (parenthetical alternative annotation): **{rep['ambiguous_truth_excluded_n']}**",
         "", f"_{rep['honesty_tier']}_", "", f"_{rep['structural_note']}_", "",
         "## Core-SNP samples (the scored tier)", "",
         "| sample | raw truth | normalized | predicted | match |", "|---|---|---|---|---|"]
    for r in sorted(core, key=lambda x: x["sample"]):
        L.append(f"| {r['sample']} | `{r['raw_truth']}` | {r['normalized_truth']} | {r['predicted']} "
                 f"| {'OK' if r['match'] else 'X'} |")
    if rep["core_miscall_diagnoses"]:
        L += ["", "### Core mis-call diagnosis"]
        for m in rep["core_miscall_diagnoses"]:
            L.append(f"- **{m['sample']}** truth `{m['raw_truth']}` -> predicted `{m['predicted']}`: "
                     f"{m['diagnosis'] or 'SNP-level disagreement'}")
    for title, bucket in (("Non-core SNP (residual — mis-called)", noncore),
                          ("Structural (EXCLUDED — BAM-required; may be silently mis-called)", structural),
                          ("Ambiguous truth (EXCLUDED — uncertain annotation)", ambiguous)):
        if bucket:
            L += ["", f"## {title}", "", "| sample | raw truth | normalized | SNP-proxy predicted |",
                  "|---|---|---|---|"]
            for r in sorted(bucket, key=lambda x: x["sample"]):
                L.append(f"| {r['sample']} | `{r['raw_truth']}` | {r['normalized_truth']} | {r['predicted']} |")
    L.append("")
    out_md.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:14]))
    print(f"[report -> {out_md} + {out_json}]")
    return 0


def _norm(diplo: str) -> tuple[str, ...]:
    """'*17|*1' / '*1/*1' -> sorted ('*1','*17'). Returns () if unparseable."""
    if not diplo or diplo.strip() in ("", "NA", "."):
        return ()
    parts = [p.strip() for p in re.split(r"[/|]", diplo.strip()) if p.strip()]
    def key(a):
        m = re.match(r"\*(\d+)", a)
        return (int(m.group(1)) if m else 999, a)
    return tuple(sorted(parts, key=key))


def _gunzip_to_plain(gz: Path) -> Path:
    """call_diplotype reads .read_text(); the VCF is bgzipped -> stage a plain-text copy once."""
    import gzip
    plain = gz.with_suffix("")  # cyp2c19_1000g.vcf
    if not plain.exists() or plain.stat().st_mtime < gz.stat().st_mtime:
        plain.write_bytes(gzip.open(gz, "rb").read())
    return plain


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="CYP2C9/CYP2C19 caller vs GeT-RM consensus on real 1000G.")
    ap.add_argument("--gene", default="cyp2c19", choices=list(GENES))
    args = ap.parse_args(argv)
    cfg = GENES[args.gene]
    if args.gene == "cyp2d6":   # dedicated tiered-honest path (structural surface); isolated from the rest
        return _run_cyp2d6(cfg)
    VCF, CORE, REF_EQUIV = cfg["vcf"], cfg["core"], cfg["ref_equiv"]
    truth_col, caller, out_stem = cfg["truth_col"], cfg["caller"], cfg["out_stem"]
    TRUTH_F = cfg.get("truth_file", TRUTH)   # per-gene truth override (CYP3A5 uses its own GeT-RM table)
    gene_label = args.gene.upper()

    if not VCF.exists() or not TRUTH_F.exists():
        print(f"ERROR: need {VCF} + {TRUTH_F} (fetch the {gene_label} region first)")
        return 2
    plain_vcf = _gunzip_to_plain(VCF)

    # 1000G sample set present in the VCF
    samples_in_vcf = set()
    for line in plain_vcf.read_text(encoding="utf-8").splitlines():
        if line.startswith("#CHROM"):
            samples_in_vcf = set(line.rstrip("\n").split("\t")[9:])
            break

    rows = []
    with open(TRUTH_F, encoding="utf-8") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            cor = rec["Coriell"].strip()
            truth = _norm(rec.get(truth_col, ""))
            if not truth or cor not in samples_in_vcf:
                continue
            r = caller(plain_vcf, cor)
            pred = _norm(r.diplotype or "")
            # *38 -> *1 phenotype-equivalent view of the truth (for the metabolizer-correct check)
            truth_pheno = tuple(sorted((REF_EQUIV.get(a, a) for a in truth),
                                       key=lambda a: (int(re.match(r"\*(\d+)", a).group(1))
                                                      if re.match(r"\*(\d+)", a) else 999, a)))
            is_core = all(a in CORE for a in truth)
            withheld = r.phenotype_status == "phenotype_withheld"
            involves_38 = any(a in REF_EQUIV for a in truth)
            # bucket: core_exact | star38_equiv | withheld | genuine_miscall
            if is_core:
                bucket = "core_exact" if pred == truth else "core_mismatch"
            elif withheld:
                bucket = "withheld"
            elif involves_38 and pred == truth_pheno:
                bucket = "star38_equiv"          # phenotype-correct; only the *38 star-label differs
            else:
                bucket = "genuine_miscall"
            rows.append({
                "sample": cor, "getrm": "/".join(truth),
                "predicted": "/".join(pred) if pred else None,
                "phenotype_status": r.phenotype_status, "withheld": withheld,
                "is_core_truth": is_core, "bucket": bucket,
                "match": (pred == truth) if is_core else None,
            })

    core = [r for r in rows if r["is_core_truth"]]
    noncore = [r for r in rows if not r["is_core_truth"]]
    core_hits = sum(1 for r in core if r["match"])
    star38_equiv = sum(1 for r in rows if r["bucket"] == "star38_equiv")
    noncore_withheld = sum(1 for r in noncore if r["withheld"])
    genuine_miscall = sum(1 for r in rows if r["bucket"] == "genuine_miscall")
    pheno_correct = core_hits + star38_equiv

    rep = {
        "schema": f"pgx-{args.gene}-getrm-concordance-v0",
        "gene": gene_label,
        "analysis_date": datetime.date.today().isoformat(),
        "truth_source": ("GeT-RM NGS consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) via the ursaPGx "
                         f"benchmark star-allele-comparison_common.tsv, column {truth_col}"),
        "genotype_source": f"1000 Genomes 30x phased panel ({gene_label} region, {cfg.get('fetch_note', 'region slice')})",
        "caller_is_independent_of_consensus_tools": True,
        "n_overlap_samples": len(rows),
        "n_core_comparable": len(core),
        "core_diplotype_concordance": round(core_hits / len(core), 4) if core else None,
        "core_diplotype_hits": f"{core_hits}/{len(core)}",
        "n_noncore_truth": len(noncore),
        "star38_phenotype_equivalent": star38_equiv,
        "noncore_correctly_withheld": noncore_withheld,
        "genuine_silent_miscall": genuine_miscall,
        "genuine_silent_miscall_pct": round(100 * genuine_miscall / len(rows), 1) if rows else None,
        "phenotype_correct_incl_star38": f"{pheno_correct}/{len(rows)}",
        "correct_or_abstains": f"{pheno_correct + noncore_withheld}/{len(rows)}",
        "honesty_tier": ("GeT-RM CONSENSUS concordance on real 1000G genomes, independent caller. The "
                         "strongest star-allele-CALLING validation tier available (vs the field's accepted "
                         "consensus truth set). v0 covers the CORE SNP set; non-core-truth samples are scored "
                         "separately (the v0.1 sentinel layer should WITHHOLD, not mis-call)."),
        "core_rows": core,
        "noncore_rows": noncore,
    }
    out_json = REPO / "wiki" / f"{out_stem}.json"
    out_md = REPO / "wiki" / f"{out_stem}.md"
    out_json.write_text(json.dumps(rep, indent=2), encoding="utf-8")
    def _star_num(a):
        d = "".join(ch for ch in a.lstrip("*") if ch.isdigit())
        return (int(d) if d else 999, a)
    core_set = "/".join(sorted(CORE, key=_star_num))

    L = [f"# {gene_label} caller vs GeT-RM consensus on real 1000G ({rep['analysis_date']})", "",
         f"**Truth:** {rep['truth_source']}", f"**Genotypes:** {rep['genotype_source']}", "",
         f"- Overlap samples scored: **{rep['n_overlap_samples']}**",
         f"- **Core-comparable diplotype concordance: {rep['core_diplotype_hits']} "
         f"({rep['core_diplotype_concordance']})**  (GeT-RM truth in {core_set})",
         f"- Phenotype-correct incl. *38==*1: **{rep['phenotype_correct_incl_star38']}** "
         f"(+{star38_equiv} *38 phenotype-equivalent samples)",
         f"- Correctly WITHHELD by sentinel: **{noncore_withheld}**",
         f"- **Genuine silent mis-call: {rep['genuine_silent_miscall']}/{rep['n_overlap_samples']} "
         f"({rep['genuine_silent_miscall_pct']}%)** -- non-core alleles beyond the v0 SNP set "
         f"(+ sentinels where present); the honest residual blind spot.",
         f"- Correct-or-abstains: **{rep['correct_or_abstains']}**", "",
         f"_{rep['honesty_tier']}_", "",
         "## Core-comparable samples (GeT-RM truth in the v0 SNP set)", "",
         "| sample | GeT-RM | predicted | match |", "|---|---|---|---|"]
    for r in sorted(core, key=lambda x: x["sample"]):
        L.append(f"| {r['sample']} | {r['getrm']} | {r['predicted']} | {'OK' if r['match'] else 'X'} |")
    if noncore:
        L += ["", "## Non-core-truth samples (v0.1 sentinel SHOULD withhold)", "",
              "| sample | GeT-RM | core-proxy | phenotype_status |", "|---|---|---|---|"]
        for r in sorted(noncore, key=lambda x: x["sample"]):
            L.append(f"| {r['sample']} | {r['getrm']} | {r['predicted']} | {r['phenotype_status']} |")
    L.append("")
    out_md.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:14]))
    print(f"[report -> {out_md} + {out_json}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
