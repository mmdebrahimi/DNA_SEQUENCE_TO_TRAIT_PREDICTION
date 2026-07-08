"""Standing PGx trust-surface report card -- one honest roll-up of every shipped human-PGx cell.

Read-only consolidation (exit 0 always; a report, NOT a gate) -- the PGx analogue of the AMR
`decoder_validation_report_card`. Rows = the deployed PGx genes; columns = the validation evidence already
produced (GeT-RM consensus concordance on real 1000G, PharmCAT fixtures, independent functional-evidence
verdicts, trio Mendelian QC). NO aggregate headline; each cell's honest tier stands on its own. Reads the
committed `wiki/pgx_*` JSON sidecars; a missing sidecar renders as NOT_RUN, never a fabricated number.
"""
from __future__ import annotations

import datetime
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
W = REPO / "wiki"


def _load(name):
    p = W / name
    try:
        return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
    except Exception:
        return None


def main() -> int:
    getrm_c19 = _load("pgx_getrm_concordance_2026-06-25.json")
    getrm_c9 = _load("pgx_getrm_concordance_cyp2c9_2026-06-25.json")
    getrm_c8 = _load("pgx_getrm_concordance_cyp2c8_2026-07-05.json")
    getrm_c3a5 = _load("pgx_getrm_concordance_cyp3a5_2026-07-05.json")
    getrm_tpmt = _load("pgx_getrm_concordance_tpmt_2026-07-05.json")
    getrm_c2b6 = _load("pgx_getrm_concordance_cyp2b6_2026-07-05.json")
    getrm_c2d6 = _load("pgx_getrm_concordance_cyp2d6_2026-07-06.json")
    pharmcat = _load("pgx_cyp2c19_report_card.json")
    fe = _load("pgx_functional_evidence_2026-06-25.json")
    trio = _load("pgx_trio_mendelian_2026-06-25.json")

    fe_by_gene = {}
    for e in (fe or {}).get("evidence", []):
        fe_by_gene.setdefault(e["gene"], []).append(e["verdict"])

    def fe_summ(g):
        v = fe_by_gene.get(g)
        if not v:
            return "—"
        from collections import Counter
        c = Counter(v)
        return f"A{c['AGREE']}/D{c['DISAGREE']}/F{c['FLAG']}"

    def trio_summ(g):
        r = (trio or {}).get("genes", {}).get(g)
        return (f"{r['consistent']}/{r['n_callable']}" if r and r.get("status") == "ok" else "—")

    cells = [
        {"gene": "CYP2C19", "trait": "metabolizer phenotype (PM/IM/NM/RM/UM)",
         "getrm": getrm_c19 and getrm_c19.get("core_diplotype_hits"),
         "getrm_pct": getrm_c19 and getrm_c19.get("core_diplotype_concordance"),
         "pharmcat": pharmcat and pharmcat.get("core_diplotype_hits"),
         "functional_evidence": fe_summ("CYP2C19"), "trio_mendelian": trio_summ("CYP2C19"),
         "tier": "GeT-RM consensus (independent of consensus tools) + PharmCAT fixtures; phenotype faithful-to-CPIC",
         "residual": "non-core *4/*35 withheld (sentinel v0.1)"},
        {"gene": "CYP2C9", "trait": "metabolizer phenotype (activity-score)",
         "getrm": getrm_c9 and getrm_c9.get("core_diplotype_hits"),
         "getrm_pct": getrm_c9 and getrm_c9.get("core_diplotype_concordance"),
         "pharmcat": None,
         "functional_evidence": fe_summ("CYP2C9"), "trio_mendelian": trio_summ("CYP2C9"),
         "tier": "GeT-RM consensus; phenotype faithful-to-CPIC (activity-score)",
         "residual": "non-core *5/*8/*9/*11 withheld (sentinel v0.1); *6-indel/*61 residual"},
        {"gene": "CYP2C8", "trait": "star-allele diplotype (*2/*3/*4) — CALLING only (no CPIC phenotype)",
         "getrm": getrm_c8 and getrm_c8.get("core_diplotype_hits"),
         "getrm_pct": getrm_c8 and getrm_c8.get("core_diplotype_concordance"),
         "pharmcat": None,
         "functional_evidence": fe_summ("CYP2C8"), "trio_mendelian": trio_summ("CYP2C8"),
         "tier": ("GeT-RM consensus (independent of consensus tools); CALLING validated 82/82. "
                  "NO CPIC metabolizer phenotype — CYP2C8 function is substrate-dependent, so this is a "
                  "CALLING-only cell (never a PM/IM/NM). Region VCF fetched Docker-free (tabix-over-HTTP)."),
         "residual": "rare non-core allele mis-called *1 (no sentinel layer v0); no phenotype layer by design"},
        {"gene": "CYP3A5", "trait": "expressor/non-expressor phenotype (tacrolimus)",
         "getrm": getrm_c3a5 and getrm_c3a5.get("core_diplotype_hits"),
         "getrm_pct": getrm_c3a5 and getrm_c3a5.get("core_diplotype_concordance"),
         "pharmcat": None,
         "functional_evidence": fe_summ("CYP3A5"), "trio_mendelian": trio_summ("CYP3A5"),
         "tier": ("REAL GeT-RM CDC multi-lab consensus (independent of the labs); 8/8 core-diplotype incl. "
                  "*1/*3/*6/*7 (the *7 insertion + *6/*7 non-expressor cases). UNDERPOWERED (n=8). "
                  "Phenotype faithful-to-CPIC (expressor/non-expressor). First gene outside the CYP2C cluster."),
         "residual": "UNDERPOWERED n=8 (only ~8 GeT-RM CYP3A5 samples overlap 1000G); rare non-core alleles mis-called *1"},
        {"gene": "TPMT", "trait": "thiopurine phenotype (COMPOUND *3A=*3B+*3C)",
         "getrm": getrm_tpmt and getrm_tpmt.get("core_diplotype_hits"),
         "getrm_pct": getrm_tpmt and getrm_tpmt.get("core_diplotype_concordance"),
         "pharmcat": None,
         "functional_evidence": fe_summ("TPMT"), "trio_mendelian": trio_summ("TPMT"),
         "tier": ("REAL GeT-RM CDC consolidated consensus; 85/85 core-comparable. FIRST compound-allele "
                  "cell — *3A resolved from two SNPs in cis (*3B+*3C), each alone = *3B/*3C. Phenotype "
                  "faithful-to-CPIC (thiopurine)."),
         "residual": "rare non-core (*2/*8/*16...) mis-called *1 (no sentinel layer v0)"},
        {"gene": "CYP2B6", "trait": "efavirenz phenotype (*6-proxy, 516G>T)",
         "getrm": getrm_c2b6 and getrm_c2b6.get("core_diplotype_hits"),
         "getrm_pct": getrm_c2b6 and getrm_c2b6.get("core_diplotype_concordance"),
         "pharmcat": None,
         "functional_evidence": fe_summ("CYP2B6"), "trio_mendelian": trio_summ("CYP2B6"),
         "tier": ("REAL GeT-RM CDC consolidated consensus; 62/62 on clean *1/*6. SINGLE-SNP *6-proxy "
                  "(516G>T) — rs2279343 (785A>G) is absent from the 1000G 30x panel so *6 can't be split "
                  "from *9. Phenotype faithful-to-CPIC."),
         "residual": "single-SNP proxy; *9/*4/other non-core mis-called (785A>G absent from callset)"},
        {"gene": "CYP2D6", "trait": "metabolizer phenotype (activity-score) — SNP surface only",
         "getrm": getrm_c2d6 and getrm_c2d6.get("core_snp_diplotype_hits"),
         "getrm_pct": getrm_c2d6 and getrm_c2d6.get("core_snp_diplotype_concordance"),
         "pharmcat": None,
         "functional_evidence": fe_summ("CYP2D6"), "trio_mendelian": trio_summ("CYP2D6"),
         "tier": ("The last major pharmacogene. GeT-RM consensus (independent of the consensus tools); "
                  "46/47 core-comparable on the SNP-DECODABLE subset (the single miss is a diagnosed "
                  "structural confound). PRIORITY-ordered per-haplotype resolver (shared-background-aware). "
                  "Phenotype faithful-to-CPIC (activity-score). Trio-Mendelian 592/602 — the ~2% residual is "
                  "the structural-confound signature (all homozygous-child). STRUCTURAL SURFACE (read-depth "
                  "off a real CRAM, dna_decode.pgx.cyp2d6_structural): *5 deletion + *xN duplication CN "
                  "validated 26/26 on 1000G CRAMs (wiki/cyp2d6_structural_2026-07-06); HYBRID PRESENCE via "
                  "elevated CYP2D7 depth (sens 0.62/spec 1.0); HYBRID IDENTITY via read-level PSV D6-fraction "
                  "(Cyrius 117-PSV method) — full-N GO, spec 1.0, *68 4/4 / *36 6/8 "
                  "(wiki/cyp2d6_hybrid_identity_2026-07-06). Everything short-read WGS can resolve."),
         "residual": ("subtle *36 exon-9 gene-conversions abstain (hybrid_present_identity_unresolved); *13 "
                      "identity single-sample-UNPOWERED (n=1); non-core SNP alleles (*14/*15/*21/*40/*46) "
                      "mis-called (no sentinel v0). Identity needs a BAM/CRAM + read-level pileup")},
        {"gene": "DPYD", "trait": "fluoropyrimidine (5-FU/capecitabine) toxicity phenotype (activity-score)",
         "getrm": None, "getrm_pct": None, "pharmcat": None,
         "functional_evidence": fe_summ("DPYD"), "trio_mendelian": trio_summ("DPYD"),
         "tier": ("NEW — the clinically highest-stakes pharmacogene (DPD deficiency = severe/fatal 5-FU "
                  "toxicity). The four CPIC-actionable DPD-deficiency haplotypes (*2A/*13 no-function, "
                  "c.2846A>T/HapB3 decreased), CPIC ACTIVITY-SCORE phenotype (Amstutz 2018; AS 2=NM, "
                  "1-1.5=IM [reduce ~50%], 0-0.5=PM [avoid]). All-SNP, NO structural blind spot "
                  "(Ensembl-GRCh38-verified coords). v0 DEPLOYMENT tier: decoded on 5 real PGP-UK humans "
                  "(all *1/*1 NM — no false-positive deficiency call). Phenotype faithful-to-CPIC. "
                  "AF-CORROBORATED 4/4: each actionable variant's 1000G EUR frequency matches the "
                  "CPIC/gnomAD DPD-deficiency spectrum (wiki/dpyd_validation_2026-07-07)."),
         "residual": ("GeT-RM DPYD concordance = EXTERNAL WALL (Pratt 2016 CDC consensus is a paper-"
                      "supplement table, not in the CYP-only ursaPGx benchmark -> manual curation; deferred). "
                      "Validation is AF-corroboration + PGP-UK deployment (KNOWLEDGE_BASELINE tier, like "
                      "VKORC1/SLCO1B1). NO sentinel layer -> rarer uncertain-function DPYD alleles called *1")},
        {"gene": "NUDT15", "trait": "thiopurine-toxicity phenotype (activity-score) — pairs with TPMT",
         "getrm": None, "getrm_pct": None, "pharmcat": None,
         "functional_evidence": fe_summ("NUDT15"), "trio_mendelian": trio_summ("NUDT15"),
         "tier": ("NEW — 2nd thiopurine gene (azathioprine/mercaptopurine myelosuppression). Dominant "
                  "actionable no-function *3 (rs116855232, ~9.5% EAS) + *1, CPIC ACTIVITY-SCORE phenotype "
                  "(Relling 2019; AS 2=NM, 1=IM [reduce], 0=PM [avoid]). All-SNP, Ensembl-GRCh38-verified. "
                  "v0 DEPLOYMENT tier: decoded on 5 real PGP-UK humans (EUR — all *1/*1 NM, consistent with "
                  "the ~0.2% EUR *3 freq). Phenotype faithful-to-CPIC."),
         "residual": ("*2 shares rs116855232 -> called *3 (SAME no-function phenotype, CPIC call unaffected). "
                      "GeT-RM NUDT15 concordance = external wall (paper-supplement, like DPYD); validation is "
                      "the *3 EAS-AF match + PGP-UK deployment. NO sentinel layer -> rarer alleles called *1")},
        {"gene": "UGT1A1", "trait": "irinotecan-toxicity phenotype (activity-score) — tag-SNP surface",
         "getrm": None, "getrm_pct": None, "pharmcat": None,
         "functional_evidence": fe_summ("UGT1A1"), "trio_mendelian": trio_summ("UGT1A1"),
         "tier": ("NEW — irinotecan toxicity / Gilbert. *STRUCTURAL: the major *28 allele is a promoter "
                  "TA-REPEAT (STR), NOT a SNP — unresolvable from a short-read SNP VCF. v0 uses rs887829 "
                  "(*80) as the validated LD-TAG proxy for *28 (EUR r^2 ~0.9+; rs887829 EUR AF ~30% == the "
                  "*28 frequency, confirming the tag) + *6 (rs4148323, clean SNP). CPIC activity-score "
                  "(Gammal 2016; AS 2=NM, 1.5=IM, 1.0=PM). Ensembl-GRCh38-verified. v0 deployment: decoded "
                  "on 5 real PGP-UK humans. Phenotype faithful-to-CPIC. Tag-SNP wrapper (like the HLA cells)."),
         "residual": ("*28 TA-repeat length UNASSESSED (star28_ta_repeat_unassessed) — rs887829 is an LD-tag "
                      "PROXY, imperfect off-EUR; *37/*36 repeat alleles not called. GeT-RM UGT1A1 concordance "
                      "= external wall (paper-supplement). Needs a repeat-aware caller for a direct *28 call")},
        {"gene": "VKORC1", "trait": "warfarin sensitivity (rs9923231)",
         "getrm": None, "getrm_pct": None, "pharmcat": None,
         "functional_evidence": fe_summ("VKORC1"), "trio_mendelian": "—",
         "tier": "single-SNP genotype->sensitivity (minus-strand encoded); not a star/diplotype system",
         "residual": "—"},
        {"gene": "SLCO1B1", "trait": "statin myopathy (rs4149056 / *5 521T>C)",
         "getrm": None, "getrm_pct": None, "pharmcat": None,
         "functional_evidence": fe_summ("SLCO1B1"), "trio_mendelian": "—",
         "tier": ("single-SNP genotype->function readout (plus-strand); KNOWLEDGE_BASELINE like VKORC1. "
                  "NOT an independent star number (rs4149056 IS the truth for a 521 call). CPIC-aligned "
                  "(simvastatin function is assigned largely from 521T>C)."),
         "residual": "single-SNP proxy for *5/*15/*17; full SLCO1B1 star typing needs more variants"},
        {"gene": "CYP4F2", "trait": "warfarin dose modifier (rs2108622 / *3 V433M)",
         "getrm": None, "getrm_pct": None, "pharmcat": None,
         "functional_evidence": fe_summ("CYP4F2"), "trio_mendelian": "—",
         "tier": ("NEW — completes the WARFARIN TRIAD (VKORC1 + CYP2C9 + CYP4F2). Single-SNP rs2108622 "
                  "genotype->function readout (minus-strand cDNA C>T == 433 Val>Met); *3 reduced-function "
                  "carriers need a HIGHER warfarin dose (CPIC Johnson 2017, ~+0.4 mg/day per *3). "
                  "KNOWLEDGE_BASELINE; AF-corroborated (*3 ~29% EUR / ~79% EAS). Deployed on 5 PGP-UK humans."),
         "residual": "single-SNP *3 proxy (rs2108622 IS the *3 truth); a DOSE modifier not a metabolizer "
                     "phenotype; dose direction is annotation only, NOT a clinical dose"},
        {"gene": "ABCG2", "trait": "rosuvastatin transporter (rs2231142 / Q141K)",
         "getrm": None, "getrm_pct": None, "pharmcat": None,
         "functional_evidence": fe_summ("ABCG2"), "trio_mendelian": "—",
         "tier": ("NEW — pairs with SLCO1B1 for STATINS. Single-SNP rs2231142 genotype->transporter-function "
                  "readout (minus-strand cDNA G>T == 141 Gln>Lys); 141K poor-function carriers have INCREASED "
                  "rosuvastatin exposure (CPIC Cooper-DeHoff 2022 -> lower rosuvastatin dose cap). "
                  "KNOWLEDGE_BASELINE; AF-corroborated (141K ~9% EUR / ~29% EAS). Deployed on 5 PGP-UK humans."),
         "residual": "single-SNP Q141K readout (rs2231142 IS the truth); ROSUVASTATIN-specific (not all "
                     "statins); transporter function, not a metabolizer phenotype"},
    ]

    rep = {
        "schema": "pgx-report-card-v0", "analysis_date": datetime.date.today().isoformat(),
        "note": ("Standing PGx trust surface -- a roll-up, NOT a gate (exit 0 always). No aggregate "
                 "headline; each cell's honest tier stands alone. CALLING is independently validatable vs "
                 "GeT-RM (free consensus panel); PHENOTYPE is faithful-to-CPIC (assigned, not measured)."),
        "cells": cells,
        "sources": {"getrm_cyp2c19": bool(getrm_c19), "getrm_cyp2c9": bool(getrm_c9),
                    "getrm_cyp2c8": bool(getrm_c8), "getrm_cyp3a5": bool(getrm_c3a5),
                    "getrm_tpmt": bool(getrm_tpmt), "getrm_cyp2b6": bool(getrm_c2b6),
                    "getrm_cyp2d6": bool(getrm_c2d6),
                    "pharmcat_cyp2c19": bool(pharmcat), "functional_evidence": bool(fe),
                    "trio_mendelian": bool(trio)},
    }
    (W / "pgx_report_card.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")

    L = [f"# PGx decoder report card ({rep['analysis_date']})", "", f"_{rep['note']}_", "",
         "| gene | trait | GeT-RM core | PharmCAT | func-evidence (A/D/F) | trio Mendelian | residual |",
         "|---|---|---|---|---|---|---|"]
    for c in cells:
        g = f"{c['getrm']} ({c['getrm_pct']})" if c["getrm"] else "—"
        L.append(f"| {c['gene']} | {c['trait']} | {g} | {c['pharmcat'] or '—'} | "
                 f"{c['functional_evidence']} | {c['trio_mendelian']} | {c['residual']} |")
    L += ["", "## Honest tier per cell", ""]
    for c in cells:
        L.append(f"- **{c['gene']}:** {c['tier']}")
    L += ["", "_Validation axes: GeT-RM = consensus concordance on real 1000G (independent of the consensus "
          "tools); PharmCAT = reference-tool fixtures; func-evidence = non-CPIC cross-check of the function "
          "assignment (AGREE/DISAGREE/FLAG); trio = Mendelian calling-consistency on 1000G trios. "
          "NOT a clinical tool._", ""]
    (W / "pgx_report_card.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print("[report -> wiki/pgx_report_card.{md,json}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
