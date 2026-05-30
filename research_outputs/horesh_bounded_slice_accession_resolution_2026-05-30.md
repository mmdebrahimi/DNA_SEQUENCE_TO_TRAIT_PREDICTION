# Horesh Bounded-Slice — Accession Resolution — 2026-05-30

> Solo laptop work executing move #1 of the bounded-slice decision (`pathotype_horesh_bounded_slice_decision_2026-05-30.md`): identify which clean H1-passing rows carry **direct WGS accessions** (runnable by accession, no re-assembly) vs Sanger lane IDs. Pure analysis on the local F1 CSV — no workhorse runtime needed.

## Candidate manifest

`research_outputs/horesh_bounded_slice_wgs_accession_candidates_2026-05-30.csv` — **260 rows**, columns: `sample_id, wgs_accession, pathotype_class, source_study, ST, isolation, country`. The workhorse materializer can fetch every row by `wgs_accession` from NCBI/ENA directly; no read re-assembly. (Tracked in `research_outputs/` — NOT `data/`, which is gitignored — so it syncs to the workhorse via `git pull`, not a relay bundle.)

## Findings (clean, supported-class rows by locator type)

| Class | clean total | **WGS accession** | Sanger lane ID | other |
|---|---|---|---|---|
| ExPEC | 1,574 | **135** | 1,319 | 120 |
| EPEC | 269 | **125** | 141 | 3 |
| ETEC | 183 | **0** | 181 | 2 |

Direct-WGS-accession sources: Salipante 2014 (129), Hazen 2016 (69), Hazen 2013 (57), Subashchandrabose 2013 (4), Chen 2013 (1).

## ⚠️ Refinement to the bounded-slice decision: the accession-only slice is 2-class, not 3

**ETEC has ZERO direct-WGS-accession rows** — all 181 clean ETEC are Sanger lane IDs (von Mentzer / Kallonen-style) that would need re-assembly from reads. So an **accession-only** bounded slice (the cheap path, no re-assembly) covers **ExPEC + EPEC only**.

Implications:
- The fast slice = **ExPEC (135) + EPEC (125) = 260 genomes**, all fetchable by accession. More than enough for a smoke + a real 2-class evaluation.
- **ETEC requires a decision:** either (a) accept ETEC needs a read→assembly step (von Mentzer reference genomes are the cleanest ETEC source anyway — they have their own assemblies; check those accessions separately), or (b) ship the first slice as ExPEC+EPEC and add ETEC in a second pass. Recommend (b): don't block the slice on ETEC re-assembly.
- ExPEC labels are isolation-source-derived (independent of markers) → genuine genotype→phenotype test. EPEC clean labels are dedicated-source-study curated (Hazen DECA) → also independent. So a 2-class ExPEC+EPEC slice is a legitimate external-validity evaluation, not a circularity artifact.

## Suggested 5-row smoke (all accession-bearing, ExPEC+EPEC)

| sample_id | wgs_accession | class | source |
|---|---|---|---|
| horesh_JSIS00000000 | JSIS00000000 | ExPEC | Salipante 2014 |
| horesh_JSMY00000000 | JSMY00000000 | ExPEC | Salipante 2014 |
| horesh_AIEY00000000 | AIEY00000000 | EPEC | Hazen 2013 |
| horesh_AIEX00000000 | AIEX00000000 | EPEC | Hazen 2013 |
| horesh_JSLK00000000 | JSLK00000000 | ExPEC | Salipante 2014 |

(2 EPEC + 3 ExPEC — exercises both classes + the by-accession fetch path. `JSIS00000000` = E. coli upec-276, already confirmed resolvable on NCBI.)

## Workhorse consumption (when synced)
1. Point the materializer's first-N selection at this CSV's `wgs_accession` column instead of the F1 `Assembly_name`.
2. Run the 5-row smoke above; then scale to the full 260 if exit criteria pass.
3. Report resolver-conformance vs external-validity separately (per the decision memo).
4. Leave ETEC out of slice 1; revisit via von Mentzer reference assemblies.

## Provenance
- Source: `data/external/horesh2021_F1_genome_metadata.csv` (local).
- WGS-accession test: `Assembly_name` matches `^[A-Z]{4,6}\d{8,}` (INSDC WGS master prefix). Sanger lane test: `^\d+_\d+#\d+`.
- Solo laptop analysis; no VirulenceFinder/blastn needed (those stay workhorse-only).
