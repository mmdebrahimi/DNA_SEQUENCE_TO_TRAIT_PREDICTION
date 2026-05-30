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

## Remaining lookup (1 step; workhorse or a later web pass)

Per-strain GCA assembly accessions did not enumerate cleanly via ENA portal (`result=assembly` empty) or the NCBI BioProject page (shows SRA reads, not assemblies) — the complete genomes are likely registered as individual NCBI GCA records. Resolution recipe:
- `datasets summary genome taxon "Escherichia coli" --assembly-source genbank` filtered to the 8 strain names, OR
- search NCBI Assembly for each strain alias (E925 etc.) + "von Mentzer", OR
- ENA `result=sample&query=study_accession="PRJEB33365"` → map the 8 sample aliases → linked assembly.

This is the only un-resolved piece; the substrate decision (use von Mentzer for ETEC) stands regardless.

## Net effect on the bounded slice

The slice can be **3-class** after all — ExPEC + EPEC (260 Horesh WGS-accession rows, `horesh_bounded_slice_wgs_accession_candidates_2026-05-30.csv`) + ETEC (8 von Mentzer references, this memo) — all by-accession, no re-assembly. But with honest per-arm validity weighting (ExPEC strongest, ETEC weakest).

## Sources
- von Mentzer 2021: https://www.nature.com/articles/s41598-021-88316-2 / PMC8085198
- BioProject: https://www.ebi.ac.uk/ena/browser/view/PRJEB33365
