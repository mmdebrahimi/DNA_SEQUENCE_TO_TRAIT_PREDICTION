# Real blaTEM CDS end-to-end forward demo — genome edit → codon → ESM2 → measured DMS (R3, 2026-07-14)

Closes the forward cell's real-surface gap: the genome-level path (`predict_genome_edit`) was unit-validated
on the genetic code; this drives it on a **REAL blaTEM coding sequence** against the **measured** wet-lab DMS.

## Substrate (real, verified)

- **CDS: PZ538321.1** (E. coli class A β-lactamase, 861 nt, fetched from NCBI, committed at
  `data/forward_ref/blatem_3349172526.fna`). **Its translation is byte-identical (0 diffs, 286 aa) to the
  ProteinGym `BLAT_ECOLX_Stiffler_2015` target protein** — so the coordinate frame is real, not synthetic.
- **Predictor: ESM2-650M** zero-shot, from the cached masked-marginal table (no model re-run).
- **Label: measured DMS** ampicillin-fitness (Stiffler 2015).

## What was run

Every single-nucleotide substitution of the real CDS sense region (858 positions × 3 alt bases = **2,574
edits**) was driven through `predict_genome_edit` → codon → AA change → ESM2, classified silent/missense/
nonsense; every MISSENSE edit that realizes a DMS-measured variant was cross-checked against its measured
fitness.

| quantity | value |
|---|---|
| single-nt edits enumerated | 2,574 |
| consequence split | 614 silent / 1,874 missense / 86 nonsense |
| missense edits realizing a DMS variant | **1,715** |
| **Spearman(ESM2, measured DMS)** on the single-nt-accessible set | **0.761** |

Consistent with the full-DMS ESM2 number (0.732) — the single-nt-accessible subset is the honest natural
mutation set a real point edit can reach (a codon-distant AA needs ≥2 nt changes and is excluded by
construction, not by cherry-picking).

## Concrete real edits (direction agrees end-to-end)

| nucleotide edit | codon | AA change | ESM2 score | measured DMS | read |
|---|---|---|---:|---:|---|
| **C542G** | CCT→CGT | **P181R** | −13.78 | −3.56 | both strongly damaging (core Pro→Arg) |
| **G628A** | GAG→AAG | **E210K** | +2.69 | +0.23 | both benign |

## Honesty rails

- **Real coordinate frame** — the CDS translation is asserted == the DMS protein (`translation_matches_dms_
  target=True`); a frame/allele error would fail loudly. Pinned offline by
  `tests/test_forward_genome_edit.py::test_committed_real_blatem_cds_translates_to_286aa`.
- **Single-nt-accessible subset only** — this is the biologically real point-edit set, not the full DMS
  (which includes multi-nt-change AAs). Reported as such, never as "the whole DMS."
- **ESM2 is the learned Regime-B predictor** (fitness-aligned molecular direction); the antagonistic
  clinical-resistance direction stays with the Regime-A determinant catalogue.

## Status

`GENOME_PATH_REAL_SURFACE_VALIDATED`. The forward "edit E. coli genotype → predict phenotype" capability is
now validated **end-to-end from a real nucleotide edit** — the north-star loop closed on the molecular
regime: real CDS → real point edit → codon → learned predictor → measured wet-lab fitness, ρ=0.76 over 1,715
real variants. Run: `uv run python scripts/blatem_genome_demo.py`. Frozen decoder surface byte-unchanged
(`verify_lock OK`); `dna_decode/forward` NON-frozen.
