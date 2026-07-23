# GEMME toolchain built + the 3-way ESM2+GEMME+ProSST hybrid on held-out MaveDB (2026-07-23)

**Status:** ✅ the GEMME evolution modality is fully built (toolchain + validation + held-out integration) and
the 3-way hybrid answers the open question: **adding evolution LIFTS the ESM2+ProSST hybrid on a leakage-free
held-out set.** Frozen AMR surface byte-unchanged (READ-only).

## C1 — the GEMME toolchain (built + validated)

`run_gemme` was a stub ("subprocess wiring deferred to the first real Linux run"). Finalized against the
self-contained **GEMME Docker image** (`elodielaine/gemme:gemme` — JET2 + R + MUSCLE + python2.7 bundled;
www.lcqb.upmc.fr/GEMME), the same Docker route as Mash/AMRFinder/Bakta. Pipeline: ColabFold MSA (cached) →
`a3m_to_aligned_fasta` (drop lowercase inserts → query-length match-column FASTA) → `docker run … python2.7
gemme.py ali.fasta -r input -f ali.fasta` → parse `*_normPred_evolCombi.txt` → `{mutation: score}`.

**Validated end-to-end on TEM-1** (GEMME's own aliBLAT tutorial protein): MSA depth 3831 → 5434 variants →
**Spearman 0.7191 vs the wet-lab DMS** (literature GEMME-TEM1 ≈ 0.72; project ESM2 0.7315). Real, strong, sane.

## C2 — the 3-way on held-out MaveDB

GEMME needs Docker (Kaggle can't run it), so GEMME per-variant tables were computed LOCALLY for 17 held-out
proteins and shipped as a **Kaggle dataset** (mirroring the ProSST pre-quantized-tokens pattern). The Kaggle
hybrid notebook loads them and computes the **3-way ESM2+GEMME+ProSST rank-average** where GEMME covers every
scored variant, comparing 2-way vs 3-way PAIRED on that subset.

**Result — N GROWN 13 → 25 (29 GEMME tables shipped; 25 with full variant coverage). The larger run
CONFIRMS + TIGHTENS the finding ~100×.**

| Decoder | median \|Spearman\| (N=25 subset) |
|---|---|
| ESM2 | (subset) |
| GEMME | 0.556 |
| **2-way (ESM2+ProSST)** | 0.575 |
| **3-way (ESM2+GEMME+ProSST)** | **0.593** |

**PAIRED (the correct statistic):**

| Comparison | N=13 (first run) | **N=25 (larger run)** |
|---|---|---|
| **3-way > 2-way** | 10/13, p=0.046 | **21/25, p=0.00046** |
| **3-way > GEMME-alone** | 11/13, p=0.011 | **22/25, p=8e-05** |
| GEMME-alone > 2-way | 3/13 | 6/25 |

Doubling the GEMME set took the 3-way-beats-2-way result from *marginal* (p=0.046) to *highly significant*
(p=0.00046) — the direction held (84% win, 21/25) and the p dropped ~100×. The median lift is +0.018.
(The N=25 covered subset is less-hard than the N=13 one — 2-way baseline 0.575 vs 0.423 — because more
higher-signal proteins are included; the within-subset comparison is what matters.)

## Findings

1. **Adding evolution (GEMME) genuinely lifts the hybrid** — the 3-way beats the 2-way ESM2+ProSST on 10/13
   held-out proteins (sign-test p = 0.046, median delta +0.015). This CONFIRMS the ProteinGym modality-hybrid
   finding (ESM2+GEMME+ProSST best) on a leakage-free held-out set.
2. **GEMME contributes ORTHOGONAL signal, it is not simply "better."** GEMME alone beats the 2-way on only
   3/13 (its higher median is inflated by a few big wins, e.g. MLH1 2-way 0.375 → 3-way 0.523). Yet the 3-way
   beats BOTH the 2-way AND GEMME-alone significantly — evolution adds information the sequence+structure pair
   misses. The three modalities (sequence ⊕ evolution ⊕ structure) are complementary.
3. **The lift is modest in magnitude but now HIGHLY significant.** The median delta is +0.018 (small — a few
   low-signal proteins like CXCR4 0.05, RAS 0.15 are where GEMME slightly hurts), but at N=25 the direction is
   decisive (21/25, p=0.00046). The larger run (user-directed "toward the full 90 to tighten the p-value")
   did exactly that: the first N=13 run was marginal (p=0.046); doubling the GEMME set confirmed it at p<0.001.

## Honest scope

- **Within-subset comparison ONLY.** The GEMME-covered subset is HARDER than average (2-way baseline 0.423
  here vs 0.602 over all 90) — do NOT compare the 3-way's 0.498 to the overall 2-way 0.602. The valid
  comparison is 2-way 0.423 → 3-way 0.498 ON THE SAME 13 proteins (the notebook computes both paired).
- N=13 (17 GEMME tables built; 4 excluded for partial variant coverage). More proteins would tighten it;
  the builder is restartable (`scripts/build_gemme_holdout_tables.py`) toward the full manifest.
- GEMME is deterministic (0 learned params); MSAs from the free ColabFold MMseqs2 API (cache-first,
  single-query etiquette). Ran locally via Docker.

## Reproduce

```
# C1: GEMME on any protein (needs Docker + the pulled elodielaine/gemme:gemme image)
python -c "from dna_decode.forward.gemme_scorer import run_gemme; from dna_decode.forward.msa_fetch import fetch_msa; ..."
# held-out GEMME tables (restartable):
MSYS_NO_PATHCONV=1 uv run python scripts/build_gemme_holdout_tables.py --limit 30
uv run python scripts/kaggle_upload_dataset.py wiki/gemme_holdout_tables.json gemme-holdout-tables
uv run python scripts/kaggle_push_poll.py push notebooks/mavedb_holdout_hybrid_kaggle_filled.py \
  mavedb-holdout-hybrid --gpu --dataset gemme-holdout-tables
```

Artifacts: `wiki/mavedb_holdout_hybrid_threeway_2026-07-23.json` (with the `threeway` block). Tests:
`tests/test_gemme_runner.py` (4, offline). GEMME tables are gitignored (2.4 MB, on the Kaggle dataset).
Builds on `wiki/mavedb_holdout_hybrid_2026-07-23.md` (the N=76 2-way).
