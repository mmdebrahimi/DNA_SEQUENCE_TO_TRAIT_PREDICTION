# Candidate sources for the TB independent gold set (golden baseline) — 2026-06-22

> Source survey for a NON-CRyPTIC TB genotype↔phenotype set (measured RIF/INH DST + downloadable genomes) to
> score the frozen TB rule out-of-distribution. Web-search-derived (MEDIUM confidence — verify the exact
> accessions + DST availability + CRyPTIC overlap before committing labor). Designed to be MERGED with a
> parallel ChatGPT shortlist. The CRyPTIC-overlap check is code: `tb_goldset.assert_independent`
> (`scripts/build_tb_goldset_manifest.py`); run it on any candidate's ENA accessions.

## The 4 criteria each candidate must clear
1. **Measured DST** (lab RIF/INH R/S, not a software prediction).
2. **Downloadable genomes** (ENA/SRA reads or assemblies → masked VCF vs H37Rv NC_000962.3).
3. **Post-2022 / non-CRyPTIC** (CRyPTIC's June-2022 release + WHO-2023 catalogue swept most earlier public TB).
4. **Provenance-checkable** (per-isolate ENA accession so the leakage check can prove non-overlap).

## Ranked shortlist (MERGED — Soraya survey + ChatGPT survey, 2026-06-22)

| # | Source | DST? | Genomes? | Access | ~N | post-2022 / non-CRyPTIC? | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | **India PZA/TB WGS study, BioProject `PRJNA1155695`** (Frontiers Microbiol 2024; ChatGPT find) | yes — **MGIT 960 pDST**, RIF 1.0 / INH 0.1 µg/mL | yes — SRA reads; Suppl. **S3** = per-isolate Project/Sample/Biosample/SRA accessions | **free public** (NCBI/SRA) | **2,207** | 2018–2020 → leakage check MANDATORY | **best powered+public** ⚠ per-isolate R/S label NOT confirmed published (S1/S2 aggregate) → likely needs AUTHOR CONTACT for the phenotype table |
| 2 | **Thorpe et al. 2024, Sci Rep 14:5201** (Thailand) | resistance profile in paper | yes — ENA, Table S1 | **free public** | ~59 | **2024 deposit** — cleanest temporal independence | **best low-effort first cut** (verify Table S1 has per-isolate RIF/INH R/S) |
| 3 | **NIAID TB Portals** ✅ **DUA SIGNED 2026-06-22 — ACTIVE PATH** | yes — "Sequenced DST" per specimen (phenotypic + molecular; we keep phenotypic only) | ~3,305 genomes + public SRA | free, **DUA signed** | ~3,305 | non-CRyPTIC consortium; overlap-check via SRA | **chosen powered option** — adapter + runbook shipped (`scripts/build_tbportals_candidates.py`, `wiki/tb_portals_goldset_runbook_2026-06-22.md`); awaiting the user's DST+accession export download |
| 4 | **Latvia MDR-TB, ENA `PRJEB59824`** (BMC ID 2023; ChatGPT find) | yes — phenotype+genotype in suppl. | yes — ENA | free public | ~63 (all MDR) | older → leakage risk | **R-only SENSITIVITY SMOKE TEST** (no specificity — all RIF/INH R) |
| 5 | **ReSeqTB** (platform.reseqtb.org) | curated WGS+DST | ~4,636 | free, registration + DUA | ~4,636 | aggregates public incl. likely CRyPTIC + pre-2022 | overlap-heavy; later |
| 6 | ENA **PRJEB7727** (Sierra Leone, 91) / **PRJEB37609** (Lancet Microbe 2020) | yes | free public | ~91 / TBD | older, overlap risk | fallbacks |

## Recommendation (merged, in order)
1. **India `PRJNA1155695` first** (ChatGPT's "start here") — IF the per-isolate RIF/INH phenotype is obtainable
   (Suppl. S3 has the accessions but the per-isolate R/S label is NOT confirmed published — **verified
   2026-06-22**; likely an author-contact item). Powered + fully public + good R/S mix if labels land.
2. **Thorpe 2024** as the **low-effort first cut** if you want something fully public + post-2023 NOW (small,
   ~59; verify per-isolate DST in Table S1).
3. **Latvia `PRJEB59824`** as a **sensitivity-only smoke test** (all MDR → tests the R path end-to-end; do
   NOT report specificity from it).
4. **TB Portals** for the powered independent number once the pipeline works (sign the free DUA).
- MVP balance target (ChatGPT): ~10 each of RIF-S/INH-S, RIF-S/INH-R, RIF-R/INH-R (+ optional RIF-R/INH-S).
  If a source can't balance, record the limitation — never fake balance.

## The make-or-break (both surveys agree)
Per isolate you must prove: **measured DST label + genome/read accession + not in CRyPTIC.** Genomes +
accessions are abundant; the scarce/binding item is the **per-isolate measured R/S label** (often aggregate
in papers, or DUA-gated, or needs author contact). That is the real wall — not variant-calling, not modeling.

## Machinery shipped to consume any winning source (2026-06-22)
- **Alias-aware leakage check** — `tb_goldset.assert_independent_aliased` (compares run+sample+biosample, SRA
  `SRR/SAMN` AND ENA `ERR/SAMEA`; ChatGPT's accession-namespace catch). Verified on the real CRyPTIC table.
- **Pre-flight validator** — `scripts/validate_tb_goldset_candidates.py` (schema + labels + class balance;
  VCF-existence opt-in so you validate labels/accessions before variant-calling).
- **Ingester + scorer** — `build_tb_goldset_manifest` (now reads `run_/sample_/biosample_accession` aliases)
  → `score_tb_independent_goldset`.

## Open items (still worth a targeted lookup)
- India `PRJNA1155695`: obtain the per-isolate RIF/INH DST table (author contact / a companion metadata file?).
- Thorpe 2024 Table S1: confirm per-isolate R/S (not just a profile).
- Any 2024–2025 single-study TB WGS+DST ENA deposit with ≥150 isolates **and per-isolate labels** = would beat all.

## Provenance
Survey 2026-06-22 (Soraya). Sources: tbportals.niaid.nih.gov/access-data + datasharing.tbportals.niaid.nih.gov;
platform.reseqtb.org; Thorpe et al. 2024 Sci Rep 14:5201 (nature.com/articles/s41598-024-55865-1); ENA
PRJEB7727; ENA PRJEB37609 (Lancet Microbe 2020). MEDIUM confidence — confirm specifics before acquisition.
