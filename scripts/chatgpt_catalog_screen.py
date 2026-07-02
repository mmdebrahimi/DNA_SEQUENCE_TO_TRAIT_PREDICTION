"""Screen the ChatGPT genome->phenotype data-source catalog through the F1 8-gate scorecard.

The catalog is overwhelmingly HUMAN-CLINICAL; its own headline warning (label leakage / evidence circularity)
IS this project's central wall. This screen routes each representative source to a verdict + a usage class.
Gate verdicts grounded in the 2026-07-02 web checks (DepMap free/measured; ClinVar type-1 circularity) +
known access facts. Writes wiki/chatgpt_catalog_screen_<date>.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts.dataset_candidate_scorecard import Candidate, rank, score  # noqa: E402

P = "pass"; F = "fail"; U = "unknown"

# Only LABEL-substrate sources are scorecard-scoreable. Feature/benchmark sources are categorized separately.
CANDIDATES = [
    Candidate(
        "DepMap / CCLE", "human cancer cell lines",
        "measured drug sensitivity (PRISM) + CRISPR dependency (quantitative)",
        {"G1_accessible": P, "G2_non_circular": P, "G3_sampling_independent": P, "G4_unit_joinable": P,
         "G5_provenance_separable": P, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        depth_estimate=1000,
        notes="WEB-VERIFIED free CSVs (depmap.org): PRISM 4,518 compounds x ~578 lines, CRISPRGeneEffect, "
              "CCLE mutations/expression for >1,000 lines, joined by depmap_id. The CANCER ANALOG of the yeast "
              "win. De-confound by lineage/tissue-of-origin (like yeast clades). EXECUTOR-PILOTABLE NOW.",
        sources=["https://depmap.org/portal/download/", "https://doi.org/10.6084/m9.figshare.9393293.v4"]),
    Candidate(
        "GDC / TCGA (open tier)", "human tumors",
        "somatic variants x survival / subtype",
        {"G1_accessible": P, "G2_non_circular": P, "G3_sampling_independent": U, "G4_unit_joinable": P,
         "G5_provenance_separable": P, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        depth_estimate=10000,
        notes="Open tier has somatic MAF + clinical/survival. G3 UNKNOWN: tumor purity/heterogeneity + "
              "site/batch confounds. Controlled tier (germline/BAM) is USER-gated.",
        sources=["https://gdc.cancer.gov/"]),
    Candidate(
        "ClinVar / ClinGen", "human germline variants",
        "curated pathogenicity assertions",
        {"G1_accessible": P, "G2_non_circular": F, "G3_sampling_independent": P, "G4_unit_joinable": P,
         "G5_provenance_separable": U, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        has_curated_catalog=True,
        notes="G2 FAIL: WEB-VERIFIED type-1 circularity -- assertions derived from the in-silico predictors a "
              "decoder competes with (ACMG PP3/BP4). REJECT as a LEARNED label. Usable only as a DETERMINISTIC "
              "catalog (WHO-TB analog) validated against a TEMPORAL-HOLDOUT / independent measured label.",
        sources=["https://www.ncbi.nlm.nih.gov/clinvar/", "https://pubmed.ncbi.nlm.nih.gov/36413997/"]),
    Candidate(
        "CIViC", "human cancer variants",
        "curated actionability/oncogenicity assertions",
        {"G1_accessible": P, "G2_non_circular": F, "G3_sampling_independent": P, "G4_unit_joinable": P,
         "G5_provenance_separable": U, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        has_curated_catalog=True,
        notes="CC0 (free). Same circularity shape as ClinVar (curated interpretations). REJECT as label; "
              "deterministic-catalog use with independent validation only.",
        sources=["https://civicdb.org/"]),
    Candidate(
        "All of Us / UK Biobank / FinnGen / dbGaP / EGA", "human cohorts",
        "measured disease/biomarker phenotypes + EHR",
        {"G1_accessible": F, "G2_non_circular": P, "G3_sampling_independent": U, "G4_unit_joinable": P,
         "G5_provenance_separable": P, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        depth_estimate=500000,
        notes="G1 FAIL by design: CONTROLLED ACCESS (researcher agreement / DAC approval / institutional "
              "affiliation, possibly cost). These are NOT executor-fetchable -- they are the USER-AUTHORITY "
              "acquisition path (reproducibility-freeze forward-path #1). Real measured phenotypes (G2 pass). "
              "G3 UNKNOWN: EHR ascertainment/site confounds. Highest supervised-training value IF acquired.",
        sources=["https://www.researchallofus.org/", "https://www.ukbiobank.ac.uk/"]),
]

# Not scorecard-scored (not phenotype-LABEL substrates) -- categorized in the artifact:
FEATURE_SOURCES = ["GTEx", "ENCODE", "SCREEN", "eQTLGen", "4D Nucleome", "FANTOM5", "Factorbook",
                   "Ensembl VEP", "CADD", "SpliceAI"]
BENCHMARK_SOURCES = ["GIAB", "CAGI", "precisionFDA"]
COMMERCIAL_GATED = ["HGMD", "COSMIC", "OncoKB", "Tempus", "23andMe"]


def main() -> int:
    from datetime import date
    ranked = rank(CANDIDATES)
    by_name = {c.name: c for c in CANDIDATES}
    L = ["# ChatGPT genome->phenotype catalog — screened through the 8-gate scorecard (2026-07-02)", "",
         "The catalog is USEFUL, and its own central warning (label leakage / evidence circularity) VALIDATES "
         "this project's hard-won discipline. But routed through the scorecard, most of it is (a) "
         "USER-authority-gated acquisition, (b) circular-label catalogs needing independent validation, or "
         "(c) human-variant FEATURE sources — NOT free executor-actionable phenotype substrates. Exactly ONE "
         "new source clears all 8 gates like the yeast substrate: **DepMap/CCLE**.", "",
         "## Label-substrate sources (scorecard-scored)", "",
         "| candidate | verdict | decoder | depth | fails | unknowns | usage class |",
         "|---|---|---|---|---|---|---|"]
    usage = {
        "DepMap / CCLE": "A — EXECUTOR-PILOTABLE NOW (yeast analog)",
        "GDC / TCGA (open tier)": "A' — pilotable, verify G3 confounds first",
        "ClinVar / ClinGen": "B — deterministic catalog + temporal-holdout validation (NOT a label)",
        "CIViC": "B — deterministic catalog + independent validation (NOT a label)",
        "All of Us / UK Biobank / FinnGen / dbGaP / EGA": "C — USER-AUTHORITY acquisition (freeze fwd-path #1)",
    }
    for s in ranked:
        L.append(f"| {s['name']} | **{s['verdict']}** | {s['decoder_type']} | "
                 f"{s['depth_estimate'] or 'catalog'} | {', '.join(s['fails']) or '—'} | "
                 f"{', '.join(s['unknowns']) or '—'} | {usage.get(s['name'],'')} |")
    L += ["", "## Non-label sources (categorized, NOT scored — they aren't genotype→phenotype label substrates)",
          "",
          f"- **D — Feature sources** (regulatory/expression features for a HUMAN variant-effect model, not "
          f"organismal-phenotype labels): {', '.join(FEATURE_SOURCES)}. Free; relevant only if the project "
          f"pivots to human noncoding variant-effect prediction.",
          f"- **E — Benchmark / anti-circularity harnesses**: {', '.join(BENCHMARK_SOURCES)}. Align with the "
          f"project's leakage discipline; adopt as an evaluation harness if going human.",
          f"- **Commercial / gated**: {', '.join(COMMERCIAL_GATED)}. License/consent-gated; USER-authority.", "",
          "## Notes per scored source", ""]
    for s in ranked:
        c = by_name[s["name"]]
        L.append(f"- **{c.name}** — {c.notes} Sources: {', '.join(c.sources)}")
    L += ["", "## Headline",
          "**DepMap/CCLE is the one new free, measured, joinable, deep, non-circular substrate** — the cancer "
          "analog of the yeast win, pilotable exactly like it. Everything else is user-gated acquisition "
          "(biobanks), circular-as-label catalogs (ClinVar/CIViC — deterministic-catalog use only, with "
          "independent validation), or human-variant feature sources. The catalog's circularity warning is the "
          "project's own wall, independently confirmed."]
    out = REPO / "wiki" / f"chatgpt_catalog_screen_{date.today().isoformat()}.md"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    for s in ranked:
        print(f"{s['verdict']:7} {s['name']}  [{usage.get(s['name'],'')}]")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
