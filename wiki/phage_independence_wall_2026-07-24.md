# Phage receptor cell — independence is `blocked:external` (2026-07-24)

**Goal:** turn the phage receptor cell from IN_DISTRIBUTION (KNOWLEDGE_BASELINE) → a real INDEPENDENT tier by
scoring the BASEL-2021-built caller on a HELD-OUT set with MEASURED receptors + fetchable genomes.

**Verdict: `blocked:external`.** The code wall is closed (a ready scorer, below); the remaining gap is data
reachability via the available fetch tools. R2/R3 pre-bar scan (no fabrication):

## The two candidate independent substrates

| Substrate | Reachable? | Fair independence test? | Why |
|---|---|---|---|
| **2026 Morinière benchmark** (Arkin/Mutalik, LBNL; 255 E. coli phages, 193 receptor-assigned, **19 receptor classes**, consistent host, **different lab**, + 49 independently-validated phages) | **NO via my tools** | **YES** (spans my receptor classes on one host) | biorxiv article + PDF return **HTTP 403** to WebFetch; the v2 GitHub repo / BioProject / supp Datasets S6+S7 URLs live inside the 403'd paper and are not surfaced by 5 web searches or `gh search`; a guessed supp media URL 404s. Open data in principle — just not fetchable by me. |
| **2025 completed BASEL** (Humolli/Maffei; 37 new phages; **BioProject PRJNA1207239**; Tables S1/S2) | **YES** | **NO** | 34/37 use **O16 O-antigen** (deliberately isolated on an O-antigen-**restored** K-12), a receptor class my 2021-trained catalogue structurally **excludes** (2021 K-12 lacks O-antigen). Worse, the receptor is **host-context-dependent** — the same genus (e.g. Teseptimavirus) is LPS_core in 2021 but O-antigen in 2025 by host strain. An OOD-boundary probe, not a fair independence test. |

## The code wall is CLOSED — only data remains

`dna_decode/phage/receptor_caller.py::independent_validate(ref_manifest, ref_dir, test_manifest, test_dir,
receptor_map=...)` scores any held-out `(accession, receptor)` manifest against the BASEL-2021 reference by
genome-homology transfer (reuses `call_receptor`; the test set's receptor is a MEASURED column, not
catalogue-inferred → genuinely independent). `_load_labeled_manifest` reads the measured receptor column +
an optional `receptor_map` to rename the study's receptor vocabulary onto our classes. Proven on synthetic
cross-set data (`tests/test_phage_receptor_caller.py`). **The moment a fair `(accession, receptor_class)`
table + genomes land, the independent number is one call away.**

## Unblock (user action — like the TB Portals case)

1. Open the 2026 benchmark **v2**: `https://www.biorxiv.org/content/10.64898/2026.04.02.716166v2` → the
   **Data / Code Availability** section. Grab: the **GitHub repo** URL, any **BioProject / GenBank** accession
   for the 255 genomes, and **supplementary Datasets S6 + S7 / Tables S1 + S3** (the per-phage receptor
   assignments). Drop the receptor table (+ the genome accessions, or the genomes) on disk.
2. Build a `receptor_map` from their 19 receptor classes → our classes (`RECEPTOR_CLASSES` in
   `dna_decode/data/phage_receptor.py`). Most map 1:1 (FhuA/BtuB/LamB/OmpC/OmpF/OmpA/Tsx/FadL/LPS_core/ECA/
   NfrA/LptD/FepA); their finer LPS/O-antigen sub-classes fold to `LPS_core`/`O_antigen`.
3. Fetch the genomes to `data/phage_ref/<indep>/` + write `<indep>_manifest.tsv` (`accession\treceptor`),
   then run `independent_validate("data/phage_ref/basel_manifest.tsv", "data/phage_ref/basel",
   "<indep>_manifest.tsv", "data/phage_ref/<indep>", receptor_map=MAP)` → the independent number.
4. The tier moves IN_DISTRIBUTION → INDEPENDENT (LBNL benchmark = different lab, different host prep, not in
   the BASEL catalogue provenance) once that number lands.

## What was NOT done (honesty)

- No fabricated independent number. The 2025 OOD-boundary run was deliberately NOT reported as an
  "independent PASS" — its O-antigen/host-context confound makes the receptor labels non-comparable to the
  catalogue's claim.
- The scorer is validated on synthetic data only (real independent data is the external gap).
