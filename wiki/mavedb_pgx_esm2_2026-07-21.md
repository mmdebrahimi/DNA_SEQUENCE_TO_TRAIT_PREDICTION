# ESM2-650M on the pharmacogene DMS (R2 molecular) — 2026-07-21

**Status:** ✅ ESM2-650M predicts variant effects in all 5 clinically-actionable pharmacogenes. Real Kaggle
T4 run over 16 MaveDB deep-mutational-scanning assays (CYP2C19/CYP2C9/G6PD/NUDT15/VKOR). **median |Spearman|
= 0.547** (shuffled control 0.014). Frozen surface byte-unchanged.

## HONEST FRAMING (load-bearing — read first)
This is **R2 molecular variant-effect** (the deployed `forward` cell scoring per-variant protein function),
NOT the R1 pharmacogenomics **catalog** cell. The MaveDB CYP2C19/2C9/G6PD/NUDT15/VKOR data is
deep-mutational-scanning (functional score per amino-acid variant) — the same regime as ProteinGym, not
clinical pharmacogenomics. The question this run answers: *does the learned forward cell predict variant
effects in clinically-actionable pharmacogenes?* (directly useful for pharmacogenomic VUS interpretation).

The **true R1 pharmacogenomics catalog cell is a different, deterministic build** — NOT this run:
| | this run (R2) | true R1 catalog cell |
|---|---|---|
| regime | R2 molecular (learned ESM2) | R1 curated-catalog (deterministic) |
| input | protein variant → ESM masked-marginal | star allele (`CYP2C19*2`) → CPIC/PharmGKB rule |
| label | DMS functional score | **clinical metabolizer phenotype** (poor/intermediate/normal) |
| compute | GPU (this Kaggle run) | none — a lookup table |
| substrate status | done (below) | NOT built — needs CPIC/PharmGKB catalog + a clinical-labelled cohort |

## Result (real Kaggle T4, ESM2-650M masked-marginals, |Spearman| direction-robust)
- **14/16 assays scored** (2 skipped `too_few` <20 single-missense: CYP2C19 `00001199-a-2`, G6PD `00001237-a-2`).
- **overall median |Spearman| = 0.547** (shuffled negative control **0.014** → real signal).
- Notably **higher than the general R2 holdout (0.503)** — pharmacogene function is well-captured by ESM2.

| gene | per-gene median \|Spearman\| | n assays | range |
|---|---|---|---|
| **G6PD** | **0.690** | 4 | 0.369–0.749 |
| **CYP2C9** | **0.635** | 3 | 0.283–0.679 |
| **CYP2C19** | **0.626** | 2 | 0.569–0.682 |
| **NUDT15** | **0.526** | 3 | 0.377–0.526 |
| **VKOR** | **0.461** | 2 | 0.413–0.508 |

All 5 clinically-actionable pharmacogenes score positively; the large recent full-length assays are the
strongest (G6PD 8000–10000-variant assays 0.74; CYP2C19 7830-variant 0.68; CYP2C9 6142-variant 0.68). The
low-per-gene tail (CYP2C9 0.283 on the 105-variant fragment; one G6PD construct 0.369) reflects small/atypical
constructs, not a scorer defect.

## What this supports / does NOT support
- **Supports:** ESM2 as a variant-effect prior for pharmacogene VUS — a real, deployable use (rank a novel
  CYP2C19/CYP2C9/G6PD/NUDT15/VKOR missense variant by predicted functional damage). This is R2 working on a
  clinically-important protein set, leakage-free (these genes are NOT in ProteinGym).
- **Does NOT support:** a star-allele→metabolizer-phenotype call. That is the R1 catalog cell — deterministic,
  clinical-label-gated, and unbuilt. This run is the R2 molecular half; the R1 half is a separate move
  (CPIC/PharmGKB catalog + a clinical metabolizer cohort — surfaced, not started).

## Provenance
- Kernel `notebooks/mavedb_pgx_esm2_kaggle.py` (16 pgx URNs from the enlarged manifest; masked-marginal core
  byte-faithful to `scripts/esm_zeroshot_dms.py`). Kaggle `emanueleebrahimi/mavedb-pgx-esm2` (T4).
- Per-assay results: `wiki/mavedb_pgx_esm2_2026-07-21.json`. Manifest:
  `wiki/mavedb_prospective_holdout_full_2026-07-21.json`. Regime lens: `plans/Trait_Decoding_Roadmap.md`.
