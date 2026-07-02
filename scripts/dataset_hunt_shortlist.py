"""F3 of the dataset hunt: instantiate the F2-researched candidates, score them with the F1 scorecard,
and emit the ranked shortlist `wiki/dataset_hunt_shortlist_<date>.md`.

Gate verdicts are from the F2 web sweep (2026-07-02) where marked WEB-VERIFIED; others are FROM-KNOWLEDGE
first-pass reads flagged UNKNOWN on the gates that a follow-up sweep must confirm (honest — never asserted).
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts.dataset_candidate_scorecard import Candidate, GATES, rank, score  # noqa: E402

P = "pass"; F = "fail"; U = "unknown"

CANDIDATES = [
    # --- WEB-VERIFIED this run (F2 2026-07-02) ---
    Candidate(
        "1002 Yeast Genomes (Peter 2018)", "S. cerevisiae (yeast)",
        "growth/fitness across ~35 lab conditions (quantitative)",
        {"G1_accessible": P, "G2_non_circular": P, "G3_sampling_independent": P, "G4_unit_joinable": P,
         "G5_provenance_separable": P, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        depth_estimate=1011,
        notes="VCF matrix + phenotype matrix both free (1002genomes.u-strasbg.fr; ENA PRJEB13017). "
              "Phenotype = controlled lab growth across conditions (common-garden -> sampling-independent). "
              "Deepest single-species substrate found.",
        sources=["http://1002genomes.u-strasbg.fr/", "https://www.nature.com/articles/s41586-018-0030-5",
                 "https://www.ebi.ac.uk/ena/browser/view/PRJEB13017"]),
    Candidate(
        "DGRP2 (Drosophila Genetic Reference Panel)", "D. melanogaster (fruit fly)",
        "31 harmonized quantitative organismal phenotypes (12 studies; DGRPool aggregates more)",
        {"G1_accessible": P, "G2_non_circular": P, "G3_sampling_independent": P, "G4_unit_joinable": P,
         "G5_provenance_separable": P, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        depth_estimate=205,
        notes="205 inbred lines, full genotypes public (dgrp2.gnets.ncsu.edu). Single source population "
              "(Raleigh NC) -> G5 met by line-level holdout, not geography. DGRPool = harmonized phenotypes.",
        sources=["http://dgrp2.gnets.ncsu.edu/", "https://www.nature.com/articles/nature10811",
                 "https://elifesciences.org/reviewed-preprints/88981"]),
    Candidate(
        "CaeNDR (Caenorhabditis Natural Diversity Resource)", "C. elegans (+briggsae/tropicalis) worm",
        "quantitative wild-isolate traits (drug/toxin response etc.)",
        {"G1_accessible": P, "G2_non_circular": P, "G3_sampling_independent": P, "G4_unit_joinable": U,
         "G5_provenance_separable": P, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        depth_estimate=300,
        notes="MIT-licensed, AWS Open Data, per-strain VCF. G4 UNKNOWN: confirm a HOSTED per-strain trait "
              "corpus exists (the GWAS tool is BYO-phenotype; Andersen-lab published trait sets likely qualify).",
        sources=["https://caendr.org/", "https://registry.opendata.aws/caendr/",
                 "https://pmc.ncbi.nlm.nih.gov/articles/PMC10767927/"]),
    # --- FROM-KNOWLEDGE first-pass (NOT web-verified this run; UNKNOWN on the risk gate) ---
    Candidate(
        "Rice 3000 Genomes + IRRI phenotypes", "Oryza sativa (rice)",
        "agronomic quantitative traits",
        {"G1_accessible": P, "G2_non_circular": P, "G3_sampling_independent": U, "G4_unit_joinable": P,
         "G5_provenance_separable": P, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        depth_estimate=3000,
        notes="G3 UNKNOWN: field/environment confound is the classic plant-GWAS risk (multi-site trials); "
              "need controlled/BLUP phenotypes. FROM-KNOWLEDGE, verify.",
        sources=["https://snp-seek.irri.org/"]),
    Candidate(
        "ClinVar", "Homo sapiens",
        "variant -> disease/pathogenicity (curated)",
        {"G1_accessible": P, "G2_non_circular": U, "G3_sampling_independent": P, "G4_unit_joinable": P,
         "G5_provenance_separable": U, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        has_curated_catalog=True,
        notes="Deterministic-path (curated catalog). G2 UNKNOWN: clinical assertions can be predictor-derived "
              "(ACMG uses in-silico tools) -> circularity risk, the exact wall the HIV cell had to dodge.",
        sources=["https://www.ncbi.nlm.nih.gov/clinvar/"]),
    Candidate(
        "Mouse Phenome Database / Collaborative Cross", "Mus musculus (mouse)",
        "large measured phenome across strains",
        {"G1_accessible": P, "G2_non_circular": P, "G3_sampling_independent": P, "G4_unit_joinable": U,
         "G5_provenance_separable": P, "G6_depth_or_catalog": U, "G7_genotype_fetchable": U,
         "G8_label_not_censored": P},
        depth_estimate=100,
        notes="G4/G7 UNKNOWN: per-strain genotype<->phenotype join granularity + fetchable per-strain "
              "genotypes need confirming. FROM-KNOWLEDGE, verify.",
        sources=["https://phenome.jax.org/"]),
    Candidate(
        "Arabidopsis flowering-time (1001G + AraPheno)", "A. thaliana (plant)",
        "flowering time (quantitative)",
        {"G1_accessible": P, "G2_non_circular": P, "G3_sampling_independent": P, "G4_unit_joinable": P,
         "G5_provenance_separable": F, "G6_depth_or_catalog": P, "G7_genotype_fetchable": P,
         "G8_label_not_censored": P},
        depth_estimate=1003,
        notes="CLOSED NEGATIVE for embeddings (G2 2026-06-12: embedding learned population structure, not "
              "the causal signal). G5 marked FAIL for the LEARNED path (structure not separable from signal). "
              "Deterministic Mendelian sub-traits only; do NOT re-run embeddings here.",
        sources=["https://arapheno.1001genomes.org/"]),
]


def main() -> int:
    from datetime import date
    ranked = rank(CANDIDATES)
    by_name = {c.name: c for c in CANDIDATES}
    lines = ["# Dataset hunt — ranked shortlist (F3, 2026-07-02)", "",
             "Scored by `scripts/dataset_candidate_scorecard.py` (the F1 8-gate rubric). "
             "PASS = all gates pass + a viable decoder paradigm; VERIFY = no fail but gaps to confirm; "
             "REJECT = a gate fails or no paradigm. Gate legend: " +
             "; ".join(f"**{k}** {v}" for k, v in GATES.items()) + ".", "",
             "| rank | candidate | creature | verdict | decoder | depth | fails | unknowns |",
             "|---|---|---|---|---|---|---|---|"]
    for i, s in enumerate(ranked, 1):
        lines.append(f"| {i} | {s['name']} | {s['creature']} | **{s['verdict']}** | {s['decoder_type']} | "
                     f"{s['depth_estimate'] or ('catalog' if s['has_curated_catalog'] else '?')} | "
                     f"{', '.join(s['fails']) or '—'} | {', '.join(s['unknowns']) or '—'} |")
    lines += ["", "## Notes + sources per candidate", ""]
    for s in ranked:
        c = by_name[s["name"]]
        lines.append(f"- **{c.name}** ({c.creature}) — {c.phenotype}. {c.notes} "
                     f"Sources: {', '.join(c.sources)}")
    lines += ["", "## Headline",
              f"Top substrate: **{ranked[0]['name']}** ({ranked[0]['verdict']}, {ranked[0]['decoder_type']}, "
              f"depth {ranked[0]['depth_estimate']}). It is the F4 pilot-fetch target.",
              "", "Honest scope: the top 3 are WEB-VERIFIED (F2 2026-07-02); the rest are FROM-KNOWLEDGE "
              "first-pass reads with UNKNOWN on their risk gate — a follow-up sweep must confirm before they "
              "rank as PASS."]
    out = REPO / "wiki" / f"dataset_hunt_shortlist_{date.today().isoformat()}.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    for s in ranked:
        print(f"{s['verdict']:7} {s['decoder_type']:14} {s['name']}")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
