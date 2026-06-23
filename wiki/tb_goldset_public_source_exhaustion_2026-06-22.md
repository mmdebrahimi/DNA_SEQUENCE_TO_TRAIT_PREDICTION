# TB independent gold set — public-source exhaustion (verified 5×, 2026-06-22)

**Conclusion: a free, public, per-isolate, MEASURED-phenotype, non-CRyPTIC TB RIF/INH set does not exist on
the open surface. Verified directly across 5 sources (incl. both of ChatGPT deep-search's top picks). The
independent TB number is an EXTERNAL wall (author contact / DUA / data acquisition), NOT a code or
search-effort problem.** This is the project's banked finding — "the binding constraint is LABELS, not
models" (reproducibility freeze 2026-06-13) — now empirically proven for TB.

## The structural pattern (the same wall every time)
Recent public TB WGS papers publish, per isolate: **the GENOTYPIC resistance call (TB-Profiler/MTBseq) +
genomes/accessions**, and the **measured phenotype only in AGGREGATE**. The per-isolate MEASURED pDST — the
only non-circular label for scoring our WHO-catalogue rule — is consistently behind **author-request or a
DUA**. Using the public per-isolate (genotypic) label = rule-vs-rule = circular (gate G1).

## The 5 verified sources
| Source | Per-isolate label that IS public | Measured pDST status | Verdict |
|---|---|---|---|
| **Thorpe 2024** (Sci Rep, Thailand, ~59) | genotypic (TB-Profiler category, Table S1) | aggregate only | NO-GO (circular label) |
| **India `PRJNA1155695`** (Front Microbiol 2025, 2207) | accessions public (Suppl S3); label not | not per-isolate public | NO-GO (author contact) |
| **TB Portals** (~3,305) | genotypic predictions (GAP tool) public | phenotypic DST in CLINICAL channel | NO-GO (DUA/access-gated; clinical split out from the genomic DUA) |
| **Ethiopia childhood 2026** (`PRJNA1104194`/`1204469`) | genotypic mutations (MOESM1 docx); BioSample carries NO phenotype | aggregate in paper + "from corresponding author on request" | NO-GO (author contact) |
| **Thailand 2025** (Thawong, `PMC12344131`, 2005) | genotypic `DR` (TB-Profiler v6.2.2, Table S1) | **Table S3 = aggregate counts only** (INH 1732/1811 R; RIF 1720/1813 R) | NO-GO (per-isolate phenotype not public) |

Each was checked at the source (EuropePMC full text + supplements; NCBI SRA/BioSample for Ethiopia). The
adapter + leakage check + ingester are all ready the instant ANY per-isolate measured-label export lands —
the blocker is exclusively the label, never the pipeline.

## Why this is not "try the next source"
ChatGPT's deep search ranked Ethiopia childhood + Thailand as the two best public candidates; both verified
NO-GO above for the SAME reason. The remaining candidates (Ethiopia EPTB, Latvia 2023/2024, India) carry the
same or worse blocker per ChatGPT's own notes (author contact / all-MDR-no-balance / leakage). The wall is
structural to how the field publishes, not specific to any one dataset — CRyPTIC + the WHO 2023 catalogue
already swept the public phenotypic TB WGS, so what remains public is the genotypic re-call.

## Recommendation: BANK TB at the in-distribution baseline
- **Shippable honest TB product:** the in-distribution baseline computed 2026-06-22 from the CRyPTIC parquet
  (RIF raw 0.916/0.974, INH 0.889/0.989; lineage-collapsed RIF 0.41 / INH 0.349; labelled in-distribution,
  NOT independent). See `wiki/tb_cryptic_parquet_baseline_2026-06-22.md`.
- **The independent number is a USER-DECISION external wall.** Two ways it ever gets cleared, both yours:
  1. **Acquire a non-public label** — an author reply (Thorpe/India/Ethiopia all said "on request"), a TB
     Portals clinical-DUA grant, or a collaborator's lab pDST. Any one of these → one command to the number
     (the whole ingestion/leakage/scoring path is built + verified).
  2. **Prospective-lock** (already shipped) — the free path that never hits this wall; accrues over time.
- **Do not** re-survey public sources or treat a genotypic per-isolate label as independent — that's the
  circular trap this document exists to prevent.

## Provenance
Verified 2026-06-22 (Soraya). Tools: EuropePMC REST (full text + supplementaryFiles), NCBI eutils
(SRA runinfo + BioSample). Adapter `scripts/build_tbportals_candidates.py`; ingester
`scripts/build_tb_goldset_manifest.py`; leakage `dna_decode/organism_rules/tb_goldset.py`. Companion:
`wiki/tb_goldset_thorpe2024_assessment_2026-06-22.md`, `wiki/tb_portals_goldset_runbook_2026-06-22.md`,
`research_outputs/tb_goldset_source_shortlist_2026-06-22.md`.
