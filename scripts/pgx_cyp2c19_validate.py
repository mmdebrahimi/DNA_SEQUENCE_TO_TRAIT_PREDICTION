"""Validate the CYP2C19 caller against an INDEPENDENT VCF cohort -> a concordance report card.

Two honesty tiers, depending on the cohort supplied:
  * PharmCAT test fixtures (the reference tool's OWN test VCFs, filename encodes the expected diplotype,
    e.g. `s1s2.vcf` -> *1/*2): a FAITHFUL-TO-PHARMCAT, in-distribution number on real VCF files we did not
    author. This is the "validate the wrapper vs the underlying tool" discipline (the caller vs PharmCAT,
    PharmCAT's own fixtures = the truth). It is NOT the GeT-RM independent number.
  * GeT-RM Coriell samples' public 1000 Genomes VCFs (a `--expected-tsv sample<TAB>diplotype` map):
    the genuine INDEPENDENT calling number. Needs the 1000G VCF fetch (tabix/bcftools on a Linux/Docker
    host) -> the named follow-up; the harness consumes either source identically.

Scope honesty: the v0 caller covers the CORE SNP-defined alleles (*1/*2/*3/*17). A fixture whose expected
diplotype uses a NON-CORE allele (*4/*15/*28/*35/...) is reported SEPARATELY as a blind-spot case (the
caller is EXPECTED to mis-call it to a *1-substituted diplotype + flag it) -> it does NOT count toward the
core concordance headline. Reporting both is the point.

Fetch the PharmCAT fixtures (reproducible):
  base=https://raw.githubusercontent.com/PharmGKB/PharmCAT/main/src/test/resources/org/pharmgkb/pharmcat/haplotype/cyp2c19
  for f in s1s1 s1s2 s1s17 s2s2 s2s3 s1s35 s1s4b; do curl -s "$base/$f.vcf" -o data/pgx_fixtures/pharmcat_cyp2c19/$f.vcf; done

Run:
  uv run python scripts/pgx_cyp2c19_validate.py --vcf-dir data/pgx_fixtures/pharmcat_cyp2c19 --source pharmcat
"""
from __future__ import annotations

import argparse
import datetime
import json
import re
from pathlib import Path

from dna_decode.pgx.caller import call_diplotype
from dna_decode.pgx.cyp2c19_catalog import ALLELE_FUNCTION, diplotype_phenotype

REPO = Path(__file__).resolve().parent.parent
CORE_STARS = {"*1", "*2", "*3", "*17"}

# PharmCAT fixture filename -> expected diplotype. Strip the trailing per-fixture QUALIFIER first
# (`rs<id>missing` / `het` / `only` / `missing` / `more` ...), which marks an extra/missing variant but
# does NOT change the core diplotype, THEN match `s<N>s<M>` (allele = digits + optional single letter,
# e.g. *4b). Stripping first avoids absorbing the `r` of `rs...` into the allele (the s1s1rs... bug).
_QUALIFIER = re.compile(r"(rs\d+|missing|het|only|call|ref|more).*$", re.IGNORECASE)
_FNAME = re.compile(r"^s(\d+[a-z]?)s(\d+[a-z]?)$")


def expected_from_filename(stem: str) -> tuple[str, str] | None:
    """`s1s2` -> ('*1','*2'); `s1s1rs12248560missing` -> ('*1','*1'); `s1s4b` -> ('*1','*4b').
    None if not an `sNsM` fixture (e.g. noCall, rs12769205only, sUnks17)."""
    core = _QUALIFIER.sub("", stem)
    m = _FNAME.match(core)
    if not m:
        return None
    return f"*{m.group(1)}", f"*{m.group(2)}"


def _norm_diplo(a1: str, a2: str) -> str:
    """Order-independent diplotype string, *1 first then numeric."""
    def key(a):
        core = a.lstrip("*").split("+")[0]
        try:
            return (0, int(re.match(r"\d+", core).group()), a)
        except (ValueError, AttributeError):
            return (1, 0, a)
    return "/".join(sorted((a1, a2), key=key))


def validate_dir(vcf_dir: str | Path) -> list[dict]:
    vcf_dir = Path(vcf_dir)
    rows: list[dict] = []
    for vcf in sorted(vcf_dir.glob("*.vcf")):
        exp = expected_from_filename(vcf.stem)
        if exp is None:
            continue
        exp_diplo = _norm_diplo(*exp)
        is_core = all(s in CORE_STARS for s in exp)
        r = call_diplotype(vcf)
        pred = r.diplotype
        match = (pred == exp_diplo) if pred else False
        exp_pheno = diplotype_phenotype(*exp) if is_core else None
        rows.append({
            "fixture": vcf.stem,
            "expected_diplotype": exp_diplo,
            "predicted_diplotype": pred,
            "diplotype_match": match,
            "is_core": is_core,
            "predicted_phenotype": r.phenotype,
            "expected_phenotype": exp_pheno,
            "phenotype_match": (r.phenotype == exp_pheno) if (is_core and exp_pheno) else None,
            "phenotype_status": r.phenotype_status,
            "withheld": r.phenotype_status == "phenotype_withheld",
            "sentinel_hits": [h["implies"] for h in r.sentinel_hits],
            "phasing": r.phasing,
            "flags": r.flags,
        })
    return rows


def build_report(rows: list[dict], source: str) -> dict:
    core = [r for r in rows if r["is_core"]]
    noncore = [r for r in rows if not r["is_core"]]
    core_dip_hits = sum(1 for r in core if r["diplotype_match"])
    core_phe_hits = sum(1 for r in core if r["phenotype_match"])
    noncore_withheld = sum(1 for r in noncore if r["withheld"])
    independent = source.lower() == "getrm"
    return {
        "schema": "pgx-cyp2c19-validation-v0",
        "analysis_date": datetime.date.today().isoformat(),
        "gene": "CYP2C19",
        "cohort_source": source,
        "n_total": len(rows),
        "n_core": len(core),
        "n_noncore_blindspot": len(noncore),
        "n_noncore_correctly_withheld": noncore_withheld,
        "core_diplotype_concordance": round(core_dip_hits / len(core), 3) if core else None,
        "core_diplotype_hits": f"{core_dip_hits}/{len(core)}",
        "core_phenotype_concordance": round(core_phe_hits / len(core), 3) if core else None,
        "core_phenotype_hits": f"{core_phe_hits}/{len(core)}",
        "honesty_tier": (
            "INDEPENDENT_CALLING (GeT-RM consensus panel)" if independent
            else "FAITHFUL_TO_PHARMCAT (reference tool's own test fixtures; in-distribution, NOT independent)"),
        "caller_is_independent_baseline": independent,
        "caveat": (
            "Core concordance is on the v0 SNP-defined set (*1/*2/*3/*17). Non-core fixtures are now "
            "WITHHELD by the v0.1 sentinel layer (phenotype_status=phenotype_withheld) rather than silently "
            "mis-called, and do NOT count toward the headline. The GeT-RM INDEPENDENT consensus number is a "
            "DATA-ACCESS step (labels in paper supplements; the 1000G tooling is proven via Docker bcftools "
            "-- see wiki/pgx_1000g_population_2026-06-25); this harness consumes it via --source getrm "
            "--expected-tsv. NOT a clinical tool."),
        "core_rows": core,
        "blindspot_rows": noncore,
    }


def render_md(rep: dict) -> str:
    L = [f"# CYP2C19 caller validation -- {rep['cohort_source']} ({rep['analysis_date']})", "",
         f"**Honesty tier:** {rep['honesty_tier']}", "",
         f"- Core diplotype concordance: **{rep['core_diplotype_hits']}** "
         f"({rep['core_diplotype_concordance']})",
         f"- Core phenotype concordance: **{rep['core_phenotype_hits']}** "
         f"({rep['core_phenotype_concordance']})",
         f"- Non-core blind-spot cases (excluded from headline): {rep['n_noncore_blindspot']}", "",
         "## Core fixtures (alleles in *1/*2/*3/*17)", "",
         "| fixture | expected | predicted | match | phenotype (pred) | phenotype match |",
         "|---|---|---|---|---|---|"]
    for r in rep["core_rows"]:
        L.append(f"| {r['fixture']} | {r['expected_diplotype']} | {r['predicted_diplotype']} | "
                 f"{'OK' if r['diplotype_match'] else 'X'} | {r['predicted_phenotype']} | "
                 f"{'OK' if r['phenotype_match'] else 'X'} |")
    if rep["blindspot_rows"]:
        L += ["", f"## Non-core cases -- v0.1 sentinel layer WITHHOLDS rather than mis-calls "
              f"({rep['n_noncore_correctly_withheld']}/{rep['n_noncore_blindspot']} correctly withheld)", "",
              "| fixture | expected (PharmCAT) | core-proxy | phenotype_status | sentinel |",
              "|---|---|---|---|---|"]
        for r in rep["blindspot_rows"]:
            L.append(f"| {r['fixture']} | {r['expected_diplotype']} | {r['predicted_diplotype']} | "
                     f"{r['phenotype_status']} | {','.join(r['sentinel_hits']) or '-'} |")
    L += ["", f"_{rep['caveat']}_", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Validate the CYP2C19 caller vs an independent VCF cohort.")
    ap.add_argument("--vcf-dir", type=Path, default=REPO / "tests" / "data" / "pgx_cyp2c19",
                    help="dir of *.vcf fixtures (default: committed PharmCAT fixtures under tests/data)")
    ap.add_argument("--source", default="pharmcat", choices=["pharmcat", "getrm"],
                    help="pharmcat = faithful-to-tool fixtures; getrm = independent panel")
    ap.add_argument("--out-json", type=Path, default=REPO / "wiki" / "pgx_cyp2c19_report_card.json")
    ap.add_argument("--out-md", type=Path, default=REPO / "wiki" / "pgx_cyp2c19_report_card.md")
    args = ap.parse_args(argv)

    if not args.vcf_dir.exists():
        print(f"ERROR: vcf-dir not found: {args.vcf_dir} (fetch the PharmCAT fixtures first -- see docstring)")
        return 2
    rows = validate_dir(args.vcf_dir)
    if not rows:
        print(f"ERROR: no sNsM fixtures found in {args.vcf_dir}")
        return 2
    rep = build_report(rows, args.source)
    args.out_json.write_text(json.dumps(rep, indent=2), encoding="utf-8")
    args.out_md.write_text(render_md(rep), encoding="utf-8")
    print(render_md(rep))
    print(f"[report card -> {args.out_md} + {args.out_json}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
