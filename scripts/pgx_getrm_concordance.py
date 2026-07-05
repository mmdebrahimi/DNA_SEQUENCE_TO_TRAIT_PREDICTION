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
}
REF_EQUIV = GENES["cyp2c19"]["ref_equiv"]   # back-compat alias (CYP2C19 *38==*1)


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
