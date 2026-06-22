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

## Ranked shortlist

| # | Source | DST? | Genomes? | Access | ~N | post-2022 / non-CRyPTIC? | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | **Thorpe et al. 2024, Sci Rep 14:5201** (Thailand clinical isolates) | yes (lineage + resistance in paper) | yes — ENA accessions in Table S1 | **free public** (ENA) | ~59 (+1456 ref) | **2024 deposit** — strongest temporal independence; check Table S1 vs CRyPTIC | **BEST low-effort first cut** |
| 2 | **NIAID TB Portals** (tbportals.niaid.nih.gov) | yes — "Sequenced DST" (measured, per specimen) | yes — ~3,305 genomes + SRA metadata | free but **DAR + DUA** (you sign) | ~3,305 | different consortium (non-CRyPTIC); overlap-check via SRA acc | **BEST powered option** |
| 3 | **ReSeqTB** (platform.reseqtb.org) | yes — curated WGS+DST knowledgebase | yes — ~4,636 public MTBC isolates + UVP pipeline | free, **registration + DUA** | ~4,636 | aggregates public sets incl. likely CRyPTIC + pre-2022 → leakage check ESSENTIAL; temporal independence weaker | useful but overlap-heavy |
| 4 | **Sierra Leone cohort PRJEB7727** (91 strains, phenotypic DST + WGS) | yes | yes — ENA PRJEB7727 | free public | ~91 | older; overlap/pre-2022 risk — check | fallback |
| 5 | **Lancet Microbe 2020 PRJEB37609** (direct-from-sample WGS + DST) | yes | yes — ENA PRJEB37609 | free public | (cohort-size TBD) | 2020 → high CRyPTIC-overlap risk | low priority |

## Recommendation (what I'd do)
1. **Lowest-effort first cut → candidate #1 (Thorpe 2024).** Fully public ENA reads + a 2024 deposit (cleanest
   "not in CRyPTIC"). Pull Table S1 accessions + the per-isolate RIF/INH DST from the paper → run the
   leakage check → if ≥~20/class survive, variant-call → score. Underpowered but a genuinely independent first signal.
2. **Powered number → candidate #2 (TB Portals).** Sign the genomic DAR/DUA; pull "Sequenced DST" + genomes +
   SRA accessions; leakage-check against CRyPTIC; score. Largest clean option.
3. Use #3–#5 to top up only what survives the leakage check.

## Open items to verify (good targets for the ChatGPT merge)
- Exact ENA accession list + per-isolate RIF/INH DST table for **Thorpe 2024 Table S1** (is DST published per isolate, or only "resistance profile"?).
- TB Portals "Sequenced DST" granularity (does it give RIF/INH R/S per isolate with the SRA accession?) + DUA turnaround.
- Any **2024–2025 single-study TB WGS+DST ENA deposit** with ≥150 isolates (a powered, fully-public, clearly-post-CRyPTIC set would beat all of the above) — my searches kept tripping a usage filter on TB+resistance phrasing, so this is the most likely gap ChatGPT fills.

## Provenance
Survey 2026-06-22 (Soraya). Sources: tbportals.niaid.nih.gov/access-data + datasharing.tbportals.niaid.nih.gov;
platform.reseqtb.org; Thorpe et al. 2024 Sci Rep 14:5201 (nature.com/articles/s41598-024-55865-1); ENA
PRJEB7727; ENA PRJEB37609 (Lancet Microbe 2020). MEDIUM confidence — confirm specifics before acquisition.
