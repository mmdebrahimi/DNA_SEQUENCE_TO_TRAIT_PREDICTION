# ETEC Reference Substrate — von Mentzer 2021 — 2026-05-30

> Soraya `--advance` run `2026-05-30-1200-ep4-pathotype`, step 1. Resolves the ETEC gap from the bounded-slice accession analysis (`horesh_bounded_slice_accession_resolution_2026-05-30.md`): all 181 clean Horesh ETEC rows are Sanger lane IDs (need re-assembly). This finds the accession-bearing ETEC alternative.

## Finding

ETEC for the bounded slice should come from the **von Mentzer 2021 7-lineage reference genomes**, NOT Horesh's Sanger-lane-ID ETEC rows.

- Source: von Mentzer et al. 2021, *Sci Rep* 11:9256 (doi 10.1038/s41598-021-88316-2; PMC8085198). Complete PacBio long-read reference genomes — 8 strains, 7 lineages (L1–L7), curated chromosomes + plasmids.
- BioProject: **PRJEB33365** (ENA secondary ERP116152). NOTE: PRJEB33365 is the *umbrella* ETEC collection (1108 SRA experiments); the 8 complete reference genomes are the published subset.
- 8 strains: E925, E1649, E36, E2980, E1441, E1779, E562, E1373 (full table in `etec_reference_vonmentzer_strains_2026-05-30.csv`).

## ⚠️ Label-provenance nuance (provenance-split discipline)

ETEC reference labels are **toxin-typed** (LT / ST), and the resolver keys on those same toxin genes (`eltA/eltB`, `estA`). So ETEC validation is **closer to resolver-conformance than external validity** — weaker than the other two arms:

| Slice arm | Label basis | Independence from resolver markers |
|---|---|---|
| ExPEC | isolation site (blood/urine) | **STRONG** — fully independent → genuine external validity |
| EPEC | Hazen-DECA curated (clinical/epi) | **MEDIUM-STRONG** — curated, largely independent |
| ETEC | toxin typing (LT/ST) | **WEAK** — same gene family the resolver uses → near-conformance |

Implication: report ETEC under the resolver-conformance column primarily; do NOT claim strong external-validity prediction skill for ETEC. ExPEC remains the cleanest external-validity arm. This keeps the bounded slice honest.

## RESOLVED 2026-05-30 (Soraya run 1807-ep4-etec-gca2) — accessions enumerated from supplementary

Downloaded the paper's supplementary XLSX (Springer CDN; PMC `/bin/` 403s curl) and parsed (openpyxl):

**(a) 8 reference chromosome accessions** (MOESM5 "Additional file 3") — now in `etec_reference_vonmentzer_strains_2026-05-30.csv` `chromosome_accession` column:

| strain | lineage | chromosome (ENA) | strain | lineage | chromosome |
|---|---|---|---|---|---|
| E925 | L1 | LR883050 | E1779 | L5 | LR883006 |
| E1649 | L2 | LR882973 | E562 | L6 | LR883000 |
| E36 | L3 | LR882997 | E1373 | L7 | LR882990 |
| E2980 | L3 | LR882978 | E1441 | L4 | LR883012 |

FASTA fetch is direct: `https://www.ebi.ac.uk/ena/browser/api/fasta/<LRxxxxxx>`.

**(b) BONUS — a 558-genome by-accession ETEC collection** (MOESM4 "Additional_file_2_new_phylo") → `etec_vonmentzer_collection_gca_2026-05-30.csv` (558 rows: `strain_id, gca_accession, lineage, pathotype, toxin_profile, country`). All carry **GCA_ accessions** (NCBI), pathotype=ETEC, lineage-typed (L1–L7), well-distributed. This is far more than the 8 references — it makes ETEC the **best-resourced arm** of the slice, fully by-accession, no re-assembly.

**Net:** the bounded slice is now ExPEC 135 + EPEC 125 + ETEC (8 complete refs OR up to 558 GCA) — all by-accession. The toxin-label-circularity caveat (ETEC = near-conformance, not strong external validity) still applies and governs how ETEC results are reported.

## Net effect on the bounded slice

The slice can be **3-class** after all — ExPEC + EPEC (260 Horesh WGS-accession rows, `horesh_bounded_slice_wgs_accession_candidates_2026-05-30.csv`) + ETEC (8 von Mentzer references, this memo) — all by-accession, no re-assembly. But with honest per-arm validity weighting (ExPEC strongest, ETEC weakest).

## Sources
- von Mentzer 2021: https://www.nature.com/articles/s41598-021-88316-2 / PMC8085198
- BioProject: https://www.ebi.ac.uk/ena/browser/view/PRJEB33365
