# J2 CPU smoke — the FIRST real ESM-2 run on real ProteinGym DMS on this machine (2026-07-08)

Answering "why can't `--advance` capture new data?" — it can. This captured real public data + ran the real
model, with **no GPU, no D:, no Kaggle, no money**: fetched the ProteinGym v1.0 reference + the smallest
assays from the public HF mirror (`ICML2022/ProteinGym`) and ran the **canonical** scorer
(`scripts/esm_zeroshot_dms.py::score_assay`) on CPU. Reusable: `scripts/j2_cpu_smoke.py`.

## Result — ESM-2 **35M** (smallest model), CPU, 6 smallest assays

| assay | seq len | rho | shuffled |
|---|---|---|---|
| IF1_ECOLI_Kelsic_2016 | 72 | **+0.477** | +0.009 |
| SUMO1_HUMAN_Weile_2017 | 101 | **+0.496** | +0.005 |
| TAT_HV1BR_Fernandes_2016 | 86 | −0.008 | −0.016 |
| F7YBW8_MESOW_Aakre_2015 | 93 | −0.007 | −0.041 |
| CCDB_ECOLI_Adkar_2012 | 101 | −0.015 | +0.039 |
| CCDB_ECOLI_Tripathi_2016 | 101 | +0.036 | +0.059 |

**Mini-median |Spearman| = 0.025** (shuffled 0.027).

## What this proves (and does NOT)

- **The J2 pipeline is CORRECT on real data.** IF1 (0.477) + SUMO1 (0.496) with shuffled ≈ 0 are
  unambiguous real signal — the real ESM-2 loads via transformers, the tokenizer/mask/`aa_ids` wiring is
  right on the real model, and `score_assay`'s masked-marginals correlate with real wet-lab DMS. This is the
  first real-model evidence beyond the synthetic unit tests (R3 real-surface).
- **The mini-median (0.025) is NOT a benchmark number** — and claiming "0.48 validated" here would be
  dishonest. Two honest reasons the median is near-chance:
  1. **Smallest-assay bias:** the CPU-cheap subset is exactly the atypical short constructs (viral TAT,
     toxin-antitoxin CCDB/F7YBW8) whose DMS phenotype (binding/toxicity) a conservation-style masked-marginal
     captures poorly — even for big models. The smallest-6 is a biased sample, not the 217-assay benchmark.
  2. **35M is the weakest ESM-2** (published 0.48 is 650M). The 650M-on-CPU test over the same 6 assays is the
     controlled check (see below).

## ESM-2 650M on the same 6 assays (CPU) — model-size control

Same 6 assays, **ESM-2 650M** (CPU, ~100–190 s/assay):

| assay | 35M rho | **650M rho** | 650M shuffled |
|---|---|---|---|
| IF1_ECOLI_Kelsic_2016 | +0.477 | **+0.599** | +0.010 |
| SUMO1_HUMAN_Weile_2017 | +0.496 | +0.509 | −0.009 |
| CCDB_ECOLI_Adkar_2012 | −0.015 | **+0.472** | +0.018 |
| CCDB_ECOLI_Tripathi_2016 | +0.036 | **+0.340** | −0.043 |
| F7YBW8_MESOW_Aakre_2015 | −0.007 | **+0.335** | −0.139 |
| TAT_HV1BR_Fernandes_2016 | −0.008 | +0.017 | +0.003 |

**650M mini-median = 0.406** (vs 35M's 0.025).

**Finding: model size rescues the "chance" assays.** The three assays 35M near-chanced on (CCDB×2, F7YBW8)
jump to 0.34–0.47 with 650M — so the near-chance 35M median was mostly **model capacity**, not the pipeline
or the assays. Only TAT_HV1BR stays low on both (+0.017) — a genuinely hard short HIV-Tat construct.
**Honest caveats:** (1) these are still the 6 *smallest* assays (biased short/hard), so 0.406 is not the
full-set number; (2) two assays have elevated shuffled |rho| (F7YBW8 0.139, CCDB_Tripathi 0.043) — a
few-effective-variant noise artifact on the smallest assays, so *those* per-assay rhos are noisier (the
median shuffled 0.014 is clean). The real 0.48-comparable number is 650M over the FULL 217-set (the GPU run).

## Honest bottom line

The real number that supports the 0.48 claim still needs the **650M+ model over the FULL assay set** — i.e.
the free-GPU Kaggle run (the smallest-6/35M CPU proxy can't stand in for it). But the pipeline is now
**validated on real data end-to-end on this machine**, and `scripts/j2_cpu_smoke.py` lets any laptop
reproduce a real per-assay number for free. Reproduce: `python scripts/j2_cpu_smoke.py --k 6 --data-dir <dir>`.

## ProteinGym v1.1 benchmark catalog (wiki/proteingym_v1.1_substitutions_catalog.tsv, 2026-07-08)

Structured map of the full current benchmark (217 substitution assays), built from the v1.1 reference — for
planning the GPU-run sharding + the free-CPU-scorable subset. Zero-thrash (metadata only).

- **Taxon:** Human 96 · Prokaryote 50 · Eukaryote 40 · Virus 31.
- **Size buckets:** ≤200 aa 96 · 201–400 aa 44 · 401–1022 aa 61 · >1022 aa 16 (windowing-only).
- **Free-CPU-tractable (≤400 aa & ≥20 single mutants): 140 / 217** — most of the benchmark can be
  characterized on a laptop with a small/mid ESM-2, no GPU.
- **Best free-CPU substrate = the v1.1 Tsuboyama domains** (TCRG1 37 aa/621 mut, PIN1 39/686, YNZC 39/714 …)
  — tiny + thousands of mutants each. NOT on the HF v1.0 mirror (404); they need the v1.1 archive (~1 GB),
  which is Kaggle-bound anyway. Staging them is the highest-value next free capture once the archive host is
  reachable.

NOTE: this catalog is a benchmark MAP, not new scored numbers — the scored data comes from the 35M-full run
(in progress) + the 650M control above.
