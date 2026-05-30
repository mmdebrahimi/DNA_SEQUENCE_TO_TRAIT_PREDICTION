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

## Remaining lookup — NARROWED 2026-05-30 (Soraya run 1135-ep4-etec-gca, PARTIAL)

Three web passes (Nature [auth-gated], PMC8085198, ENA/NCBI portal) did NOT enumerate the 8 per-strain chromosome accessions — they live in the paper's **supplementary Additional File 4/5** (binary/Excel; not HTML-indexed). What IS now pinned:

- Deposition: **EMBL/ENA, study ERP116152** (= BioProject PRJEB33365), submission ERA2030910.
- Records are **ENA sequence accessions in the `LR88xxxx` range**, NOT GCA assemblies. Plasmids are `LR883051`+ (paper Table 2); the 8 chromosomes are adjacent LR accessions in the same submission.
- Once the 8 chromosome LR accessions are known, FASTA fetch is direct: `https://www.ebi.ac.uk/ena/browser/api/fasta/<LRxxxxxx>`.

**Exact remaining step (money-free; pick one):**
1. Download Additional File 4 or 5 from PMC8085198 (supplementary) → read the per-strain chromosome accession column.
2. OR `datasets summary genome taxon "Escherichia coli"` + grep the 8 strain aliases (E925 …) — if NCBI mirrored them as GCA.
3. OR ENA: enumerate all sequence records under study ERP116152 and match the 8 strain aliases.

The substrate DECISION (use von Mentzer for ETEC) stands regardless; only the literal accession list is pending. Do NOT guess the LR numbers — confirm from a source.

## Net effect on the bounded slice

The slice can be **3-class** after all — ExPEC + EPEC (260 Horesh WGS-accession rows, `horesh_bounded_slice_wgs_accession_candidates_2026-05-30.csv`) + ETEC (8 von Mentzer references, this memo) — all by-accession, no re-assembly. But with honest per-arm validity weighting (ExPEC strongest, ETEC weakest).

## Sources
- von Mentzer 2021: https://www.nature.com/articles/s41598-021-88316-2 / PMC8085198
- BioProject: https://www.ebi.ac.uk/ena/browser/view/PRJEB33365
