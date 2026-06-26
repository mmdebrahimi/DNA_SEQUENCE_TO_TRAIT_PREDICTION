"""PGx trio Mendelian co-segregation QC (Unit B) -- a free internal validation of the diplotype CALLING.

The 1000 Genomes 30x panel includes 602 trios. A correctly-called diplotype must be Mendelian-consistent:
the child's two alleles must descend one-from-each-parent. This checks the PGx caller's calls across all
trios whose 3 members are in the on-disk region VCFs (CYP2C19 + CYP2C9) and counts consistent / violation /
uncallable. A clean consistency rate strengthens the CALLING claim INDEPENDENTLY of GeT-RM (different axis:
inheritance physics, not a consensus panel).

PERF: the VCF is read ONCE per gene into an in-memory {pos: per-sample-GT} map (the 1000G lines carry 3202
sample columns -- a per-sample re-read would re-split billions of fields; the wide-line-split footgun). Each
trio member's diplotype is then assembled from memory via the shared `assemble_diplotype`.

Honest scope: validates CALLING consistency (Mendelian), not phenotype. A VIOLATION definitively flags a
calling error; consistency is necessary-not-sufficient (parent+child could share a mis-call). Independent
axis from the GeT-RM panel. Needs the 1000G ped + the region VCFs (Docker-fetched).
"""
from __future__ import annotations

import datetime
import gzip
import json
from pathlib import Path

from dna_decode.pgx import cyp2c9_catalog as c9
from dna_decode.pgx import cyp2c19_catalog as c19
from dna_decode.pgx.caller import _interpret, assemble_diplotype

REPO = Path(__file__).resolve().parent.parent
PED = REPO / "data" / "pgx_1000g" / "1000G_3202_samples.ped"
GENES = {
    "CYP2C19": {"vcf": REPO / "data" / "pgx_1000g" / "cyp2c19_1000g.vcf.gz",
                "defining": c19.CORE_DEFINING, "sentinels": c19.SENTINELS,
                "ref": c19.REFERENCE_ALLELE, "pheno": c19.diplotype_phenotype, "gene": c19.GENE},
    "CYP2C9":  {"vcf": REPO / "data" / "pgx_1000g" / "cyp2c9_1000g.vcf.gz",
                "defining": c9.CORE_DEFINING, "sentinels": c9.SENTINELS,
                "ref": c9.REFERENCE_ALLELE, "pheno": c9.diplotype_phenotype, "gene": c9.GENE},
}


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def _diplotypes_for_gene(cfg) -> dict[str, str | None]:
    """ONE pass over the VCF -> {sample: diplotype string} for every sample, using the shared assembler."""
    want_pos = {(d.chrom, d.pos): ("def", d) for d in cfg["defining"]}
    for s in cfg["sentinels"]:
        want_pos[(s.chrom, s.pos)] = ("sent", s)
    samples: list[str] = []
    rows: dict[int, dict] = {}   # pos -> {ref, alt, gts:[...]}
    with gzip.open(cfg["vcf"], "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = cols[9:]
                continue
            key = (_norm_chrom(cols[0]), int(cols[1])) if cols[1].isdigit() else None
            if key not in want_pos:
                continue
            fmt = cols[8].split(":")
            gi = fmt.index("GT") if "GT" in fmt else 0
            rows[int(cols[1])] = {"ref": cols[3], "alt": cols[4],
                                  "gts": [c.split(":")[gi] for c in cols[9:]]}
    idx = {s: i for i, s in enumerate(samples)}
    out: dict[str, str | None] = {}
    for sample, si in idx.items():
        calls = {}
        for d in cfg["defining"]:
            r = rows.get(d.pos)
            calls[d.star] = (_interpret(d, r["ref"], r["alt"], r["gts"][si]) if r
                             else _interpret(d, d.ref, d.alt, None))
        sent_counts = {}
        for s in cfg["sentinels"]:
            r = rows.get(s.pos)
            if r and s.alt in r["alt"].split(","):
                ai = r["alt"].split(",").index(s.alt) + 1
                nums = [int(a) for a in r["gts"][si].replace("|", "/").split("/") if a.isdigit()]
                sent_counts[s.rsid] = sum(1 for n in nums if n == ai)
        res = assemble_diplotype(calls, sentinel_counts=sent_counts, reference_allele=cfg["ref"],
                                 phenotype_fn=cfg["pheno"], sentinels=cfg["sentinels"], gene=cfg["gene"])
        out[sample] = res.diplotype
    return out


def _mendelian_ok(child: list[str], father: list[str], mother: list[str]) -> bool:
    """child {c1,c2} consistent if one allele comes from each parent's diplotype (either assignment)."""
    c1, c2 = child
    return (c1 in father and c2 in mother) or (c2 in father and c1 in mother)


def main(argv=None) -> int:
    if not PED.exists():
        print(f"ERROR: 1000G ped not found at {PED}")
        return 2
    trios = []
    for line in PED.read_text(encoding="utf-8").splitlines():
        f = line.split()
        if len(f) < 5 or f[0] == "FamilyID":
            continue
        if f[2] != "0" and f[3] != "0":
            trios.append((f[1], f[2], f[3]))   # (child, father, mother)

    gene_results = {}
    for gene, cfg in GENES.items():
        if not cfg["vcf"].exists():
            gene_results[gene] = {"status": "vcf_absent"}
            continue
        diplo = _diplotypes_for_gene(cfg)
        usable = [t for t in trios if all(m in diplo for m in t)]
        consistent = violation = uncallable = 0
        viols = []
        for child, father, mother in usable:
            cc, cf, cm = diplo[child], diplo[father], diplo[mother]
            if not (cc and cf and cm):
                uncallable += 1
                continue
            if _mendelian_ok(cc.split("/"), cf.split("/"), cm.split("/")):
                consistent += 1
            else:
                violation += 1
                if len(viols) < 20:
                    viols.append({"child": child, "c": cc, "f": cf, "m": cm})
        n = consistent + violation
        gene_results[gene] = {
            "status": "ok", "n_trios_in_vcf": len(usable), "n_callable": n,
            "consistent": consistent, "violation": violation, "uncallable": uncallable,
            "mendelian_consistency": round(consistent / n, 4) if n else None,
            "violations_sample": viols,
        }

    rep = {
        "schema": "pgx-trio-mendelian-v0", "analysis_date": datetime.date.today().isoformat(),
        "cohort": "1000 Genomes 30x panel trios (20130606_g1k_3202_samples ped)",
        "n_trios_total": len(trios),
        "honesty": ("Validates CALLING consistency (Mendelian inheritance), NOT phenotype. A VIOLATION "
                    "definitively flags a calling error; consistency is necessary-not-sufficient. Independent "
                    "axis from the GeT-RM panel."),
        "genes": gene_results,
    }
    out_json = REPO / "wiki" / "pgx_trio_mendelian_2026-06-25.json"
    out_md = REPO / "wiki" / "pgx_trio_mendelian_2026-06-25.md"
    out_json.write_text(json.dumps(rep, indent=2), encoding="utf-8")
    L = [f"# PGx trio Mendelian co-segregation QC ({rep['analysis_date']})", "",
         f"**Cohort:** {rep['cohort']} ({rep['n_trios_total']} trios)", ""]
    for gene, r in gene_results.items():
        if r.get("status") != "ok":
            L.append(f"- {gene}: {r.get('status')}"); continue
        L.append(f"- **{gene}: Mendelian consistency {r['consistent']}/{r['n_callable']} "
                 f"({r['mendelian_consistency']})** across {r['n_trios_in_vcf']} trios "
                 f"(violations {r['violation']}, uncallable {r['uncallable']})")
    L += ["", f"_{rep['honesty']}_", ""]
    out_md.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"[report -> {out_md} + {out_json}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
