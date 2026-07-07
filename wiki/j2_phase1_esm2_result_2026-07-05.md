# J2 Phase 1 — the REAL ESM2-650M number landed (2026-07-05, free GPU)

**The learned-representation thesis, confirmed on the real model.** Ran ESM2-650M (zero-shot
masked-marginals) over the ProteinGym deep-mutational-scan benchmark on a **free Kaggle GPU** and correlated
its per-variant scores against the wet-lab measured effects. Result **matches the published field number**:

| metric | value |
|---|---|
| assays scored | **201** of 217 (16 skipped: too-long >1022 aa or <20 joinable variants) |
| **median \|Spearman\| (ESM2-650M vs wet-lab DMS)** | **0.4914** |
| shuffled negative control (median) | **0.0136** (~0 → signal real, not artifact) |
| PASS bar (pre-registered ≥ 0.45) | **PASS** |
| field-match (≥ 0.48, the published ESM2-650M number) | **YES — matches/slightly exceeds** |
| device | CUDA (Tesla P100, via the self-heal below) |

Strongest joins (signed Spearman): `A4GRB6_PSEAI_Chen_2020` +0.738 · `BLAT_ECOLX_Firnberg_2014` +0.737 ·
`BLAT_ECOLX_Stiffler_2015` +0.731 · `GRB2_HUMAN_Faure_2021` +0.726 · `TCRG1_MOUSE_Tsuboyama_2023` +0.719.

Full per-assay data: `wiki/j2_phase1_esm2_result_2026-07-05.json`.

## What it confirms
- The DMS falsifier (2026-07-04) proved the substrate with an **AlphaMissense proxy** at median |ρ| **0.417**.
  This run confirms it on the **REAL learned protein model** (ESM2-650M itself) at **0.491** — J2 Phase 1
  done, on our own model, matching the field. The learned/JEPA protein-representation direction is
  empirically green-lit at the molecular-variant-effect layer.
- Scope unchanged: this is the **molecular** layer (protein function is strongly sequence-determined). It
  does NOT rescue the complex-organismal-phenotype direction (the 0-for-5 de-confounded negative stands).
- **Phase 2** ("beat 0.48" via fine-tune / a trained JEPA) is the genuine research bet — deliberately NOT
  run here; gated on whether beating the field powers a concrete decoder feature.

## How it ran (free, zero local GPU)
Executed entirely from the laptop via the **Kaggle API** on Kaggle's **free GPU** — no browser, no money,
kernel kept private:
- Data: the committed ProteinGym subset uploaded as the Kaggle dataset `emanueleebrahimi/proteingym-j2`
  (a single forward-slash zip after a Windows `Compress-Archive` backslash bug was fixed with a Python
  `zipfile` rebuild).
- Kernel: `scripts`-style headless script pushed via `kaggle kernels push`; polled to completion; result
  pulled via `kaggle kernels output`.
- **Self-heal (load-bearing):** Kaggle's free single-GPU was a **Tesla P100 (sm_60)**, and Kaggle's stock
  PyTorch dropped Pascal support (`no kernel image available for execution on the device`). The kernel now
  reinstalls an official sm_60-capable `torch==2.5.1` (cu121) **and removes the mismatched
  torchvision/torchaudio** (stale ABI would crash the transformers import; ESM needs neither), then runs a
  fast CUDA smoke-test before the real work. This makes the kernel run on ANY free Kaggle GPU (P100 or T4).

The runner + runbook + offline tests are the committed deliverable: `notebooks/j2_phase1_esm2_proteingym.py`
+ `notebooks/J2_PHASE1_RUNBOOK.md` + `tests/test_j2_phase1_notebook.py`. Frozen AMR surface byte-unchanged.
