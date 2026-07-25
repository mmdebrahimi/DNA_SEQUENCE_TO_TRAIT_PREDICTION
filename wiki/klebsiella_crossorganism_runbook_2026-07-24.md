# Klebsiella cross-organism generalization — decision-ready runbook (2026-07-24, overnight)

**The axis:** does the deterministic RBP/depolymerase → phenotype paradigm (validated for E. coli phage →
outer-membrane receptor) transfer CROSS-ORGANISM to *Klebsiella pneumoniae* phage → capsule (K/KL-type)?
A new host + new phenotype test of the paradigm's breadth. Free, reversible, mine (phage/host-tropism lane),
non-colliding (no sibling session owns Klebsiella).

**Why it's a runbook, not a result:** the data is reachable IN PRINCIPLE but browser-gated for my headless
tools (Mendeley public-API returns empty file lists; Nature redirects to an auth cookie page). Same class as
the LBNL/biorxiv cases — one browser session unblocks it. So this is a morning drop-and-run.

## Data sources (located; grab any ONE)

| Source | What it has | URL |
|---|---|---|
| **Beamud 2023** "Genetic determinants of host tropism in Klebsiella phages" (Cell Reports) | Mendeley package: **Depolymerase_DB**, a **capsular-tropism** table, candidate receptors, predicted RBP PDBs | https://data.mendeley.com/datasets/c696dvvynf/2 |
| **Ferriol-González 2024** (Microbiol Spectrum) | **71 phages × 77 KL reference serotypes** infection matrix + **depolymerase amino-acid sequences** (supp) | https://pmc.ncbi.nlm.nih.gov/articles/PMC11448410/ |
| **Nat Comms 2025** s41467-025-63861-w "Unlocking data in Klebsiella lysogens to predict capsular type-specificity of phage depolymerases" | a curated **depolymerase → KL-type** dataset + a prediction model (a natural baseline) | https://www.nature.com/articles/s41467-025-63861-w |

**What to drop on disk:** a table/FASTA mapping each phage **depolymerase protein sequence → the KL(K)-type it
targets** (e.g. `depolymerase.faa` + a `depolymerase\tKL_type` TSV). Save to `data/phage_ref/klebsiella/`.

## The build (ready to execute the moment the data lands)

1. Ingest the depolymerase→KL-type table + sequences.
2. Build a **depolymerase k-mer caller** — the exact `rbp_caller` architecture (`protein_kmers` +
   `nearest_rbp_receptor`), retargeted: nearest-depolymerase transfer → predicted KL-type. Wheel-only.
3. **Leave-one-out** over the depolymerase set → cross-ORGANISM number: does depolymerase homology predict
   KL-type? Report per-KL-type + overall.
4. Compare to the Nat Comms 2025 model (if its predictions are in the dataset) — an independent baseline.

## Honest hypothesis (UNFALSIFIED — the cheap test decides it)

Two-sided prior:
- **Against:** the E. coli cross-lab RBP number came in at **0.364** (`wiki/phage_rbp_crosslab_result`; best-case
  = best-match, a genuine k-mer-transfer generalization limit). Cross-ORGANISM is strictly harder than cross-lab,
  so naive k-mer transfer may fail here too.
- **For:** depolymerases are **more modular** than tail fibers — the capsule-degrading enzymatic domain is a
  cleaner sequence→function unit, and the field (Beamud/Ferriol/NatComms) reports depolymerase sequence
  clustering DOES track KL-type. So a homology caller may transfer BETTER on capsule than on OMP receptors.

The leave-one-out number decides it. Either outcome is a real finding: TRANSFERS = the paradigm generalizes
cross-organism on modular determinants; FAILS = the paradigm is host-specific even for determinant-scan.

## Overnight status
Data-fetch attempted (2 searches + 3 fetch attempts); reachable-but-browser-gated → parked as this runbook
rather than rabbit-holed. The overnight window was instead spent on the CONTAINED, in-hand diagnostic that
closed the E. coli cross-lab extraction-noise caveat (see `phage_rbp_crosslab_result_2026-07-24.md`).
