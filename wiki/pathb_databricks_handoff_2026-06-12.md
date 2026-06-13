# Path B Databricks Handoff (2026-06-12)

## Purpose

Freeze the current Arabidopsis Path B Databricks state after the first real primary-model packet, two successful seed-stability repeats, and the failed `window_budget=128` scale-up attempt on the T4 cluster.

## Executive summary

- Primary-model environment is **real and runnable** on Databricks GPU.
- `PlantCaduceus_l32` probe is **green** on `Eng-GPU-Cluster-01`.
- Full FT10 inputs are staged under the requested `RCA Analysis` volume root.
- Three real primary-model packets at `window_budget=64` are complete.
- `window_budget=128` was attempted and failed on T4 GPU memory.
- Scientific status is now:
  - **not environment-blocked**
  - **not a clean win**
  - **stable mixed / ambiguous, with repeatable global `r2` advantage at `window_budget=64`**

## Active Databricks target

- Workspace:
  - `eng-601`
- Cluster:
  - `Eng-GPU-Cluster-01`
  - cluster id `0610-172603-xlr5p5ei`
  - single-node `Standard_NC4as_T4_v3 [T4]`
- Active workspace bundle root:
  - `/Workspace/Users/zsvc-athena-eng@aero.bombardier.com/plant_trait_databricks_bundle`
- Active volume root:
  - `/Volumes/dev_insightzone/iz-engineering-data-analytics/stress/RCA Analysis`

## Databricks notebooks

- Install/runtime:
  - `https://adb-4018984958712813.13.azuredatabricks.net/?o=4018984958712813#notebook/441461760253312`
- Probe:
  - `https://adb-4018984958712813.13.azuredatabricks.net/?o=4018984958712813#notebook/791402853440459`
- Bounded run:
  - `https://adb-4018984958712813.13.azuredatabricks.net/?o=4018984958712813#notebook/791402853440460`

## Confirmed staged inputs

Under:

- `/Volumes/dev_insightzone/iz-engineering-data-analytics/stress/RCA Analysis/`

Confirmed present:

- `phenotype/FT10_pheno_261.csv`
- `groups/groups_ft10_pc9.csv`
- `annotation/GCF_000001735.4_TAIR10.1_genomic.gff`
- `kinship/kinship_ibs_mac5.hdf5`
- `pseudogenomes/` with `1003` files
- manifest:
  - `inputs_manifest/g2_dry_manifest_2026-06-09_ft10_full1003.json`

## Confirmed probe success

Probe succeeded on GPU with:

- `mamba_ssm` importable
- CUDA visible
- `Tesla T4` visible
- `PlantCaduceus_l32` config load OK
- tokenizer load OK
- model load OK
- forward pass OK

Practical verdict:

- `ready_for_primary_pathb_probe`

## Load-bearing code fixes already applied

### Notebook/runtime bootstrap

`01_probe_plantcaduceus.py`

- default `bundle_root` fixed to `zsvc-athena-eng`
- bootstrap install path if `mamba_ssm` absent
- helper imports switched to file-path loading
- `sys.modules` registration added
- Hugging Face caches redirected to `/local_disk0/rca_analysis/hf_cache`

`02_run_pathb_primary_bounded.py`

- default `bundle_root` fixed
- bootstrap install path added
- stale cached modules purged before runner load
- post-`restartPython()` path hardened repeatedly:
  - rehydrate widget/runtime state directly from widgets
  - re-import `Path`, `os`, `sys`, `yaml`, `importlib`, `importlib.util`
  - avoid dependence on first-cell-only helper state

### Primary model wrapper

`plant_trait_runtime/models/foundation.py`

- filter tokenizer outputs to model-supported keys
- drop unsupported fields like `token_type_ids`
- synthesize all-ones mask when `attention_mask` is absent

### Config

`config/datasources.yaml`

- `plantcaduceus.embedding_dim` fixed:
  - old: `768`
  - new: `2048`

## Real bugs encountered and already resolved

1. Wrong workspace user path
2. `mamba_ssm` install / visibility confusion after notebook restart
3. `ModuleNotFoundError: scripts`
4. Dataclass loader failure from file-path imports without `sys.modules`
5. `token_type_ids` unsupported by `Caduceus.forward`
6. `KeyError: 'attention_mask'`
7. Stale runtime modules surviving reruns
8. Stale `embedding_dim = 768`
9. Multiple post-restart notebook state-loss bugs in `02_run_pathb_primary_bounded.py`

## Completed scientific results

### Packet 1: canonical first primary-model packet

Output files:

- `dbfs:/Volumes/dev_insightzone/iz-engineering-data-analytics/stress/RCA Analysis/outputs/runs/plant_trait_primary_bounded_full1003.json`
- `dbfs:/Volumes/dev_insightzone/iz-engineering-data-analytics/stress/RCA Analysis/outputs/runs/plant_trait_primary_bounded_full1003.md`

Metrics:

- embedding:
  - `r2 = -0.0455`
  - `spearman = 0.2129`
  - `within_group_r2 = -0.1725`
- structure-only:
  - `r2 = -0.4490`
  - `spearman = 0.4843`
  - `within_group_r2 = 0.0375`

Interpretation:

- embedding beats structure-only on global `r2`
- structure-only beats embedding on rank correlation and within-group fit
- verdict: **mixed / ambiguous first packet**

### Packet 2: seed-stability repeat

Run:

- `window_budget = 64`
- `seed = 7`

Output files:

- `dbfs:/Volumes/dev_insightzone/iz-engineering-data-analytics/stress/RCA Analysis/outputs/runs/plant_trait_primary_bounded_full1003_seed7.json`
- `dbfs:/Volumes/dev_insightzone/iz-engineering-data-analytics/stress/RCA Analysis/outputs/runs/plant_trait_primary_bounded_full1003_seed7.md`

Metrics:

- embedding:
  - `r2 = -0.0355`
  - `spearman = 0.2211`
  - `within_group_r2 = -0.1097`
- structure-only:
  - `r2 = -0.4490`
  - `spearman = 0.4843`
  - `within_group_r2 = 0.0375`

Interpretation:

- same directional pattern as seed `42`
- embedding global `r2` remains better than structure-only
- structure-only still better on `spearman` and `within_group_r2`
- verdict: **mixed, slightly more stable than a single one-off**

### Packet 3: second seed-stability repeat

Run:

- `window_budget = 64`
- `seed = 99`

Output files:

- `dbfs:/Volumes/dev_insightzone/iz-engineering-data-analytics/stress/RCA Analysis/outputs/runs/plant_trait_primary_bounded_full1003_seed99.json`
- `dbfs:/Volumes/dev_insightzone/iz-engineering-data-analytics/stress/RCA Analysis/outputs/runs/plant_trait_primary_bounded_full1003_seed99.md`

Metrics:

- embedding:
  - `r2 = -0.0305`
  - `spearman = 0.2258`
  - `within_group_r2 = -0.1281`
- structure-only:
  - `r2 = -0.4490`
  - `spearman = 0.4843`
  - `within_group_r2 = 0.0375`

Interpretation:

- same directional pattern as seeds `42` and `7`
- embedding global `r2` remains better than structure-only
- structure-only still better on `spearman` and `within_group_r2`
- verdict: **stable mixed / ambiguous across three real packets**

### Packet 4: higher-window follow-up

Run:

- `window_budget = 128`
- `seed = 42`

Result:

- failed with **CUDA OOM** on the single-node `Tesla T4` cluster
- practical failure mode: attempted extra allocation of `512 MiB` after the process had effectively exhausted available device memory

Interpretation:

- this is a **hardware/memory ceiling**, not a notebook-control-flow or model-load failure
- `window_budget=128` is not currently runnable unchanged on this T4 environment

## Current follow-up state at handoff

Resolved bounded packet set:

1. `window_budget=64`, `seed=42`, `run_suffix=full1003`
2. `window_budget=64`, `seed=7`, `run_suffix=full1003_seed7`
3. `window_budget=64`, `seed=99`, `run_suffix=full1003_seed99`

Attempted but unresolved scale-up:

4. `window_budget=128`, `seed=42`, `run_suffix=full1003_wb128`
   - failed on T4 memory (`CUDA OOM`)

## Resume instructions

### First check

Open these existing outputs under:

- `/Volumes/dev_insightzone/iz-engineering-data-analytics/stress/RCA Analysis/outputs/runs/`

Files that should now exist:

- `plant_trait_primary_bounded_full1003.json`
- `plant_trait_primary_bounded_full1003.md`
- `plant_trait_primary_bounded_full1003_seed7.json`
- `plant_trait_primary_bounded_full1003_seed7.md`
- `plant_trait_primary_bounded_full1003_seed99.json`
- `plant_trait_primary_bounded_full1003_seed99.md`

### Resume decision

1. Treat the `window_budget=64` result as **stable mixed / ambiguous** across three real seeds.
2. Do **not** retry `window_budget=128` unchanged on the same T4 cluster.
3. If a higher-window follow-up is still needed, only do it with:
   - a larger GPU, or
   - a more memory-efficient inference path, or
   - a reduced per-pass window load
4. Otherwise freeze G2 as:
   - primary model runnable
   - bounded evidence collected
   - still ambiguous, not confirmed

## What not to do

- Do not reopen Windows fallback work
- Do not go back to BIS
- Do not use the old 3-pack follow-up job shape again
- Do not revert the neutral Databricks-facing naming

## Load-bearing local files

- `C:\Users\b0652085\PycharmProjects\PythonProject\DNA_AI_Decoder\plant_trait_databricks_bundle\notebooks\00_install_primary_runtime_combined.py`
- `C:\Users\b0652085\PycharmProjects\PythonProject\DNA_AI_Decoder\plant_trait_databricks_bundle\notebooks\01_probe_plantcaduceus.py`
- `C:\Users\b0652085\PycharmProjects\PythonProject\DNA_AI_Decoder\plant_trait_databricks_bundle\notebooks\02_run_pathb_primary_bounded.py`
- `C:\Users\b0652085\PycharmProjects\PythonProject\DNA_AI_Decoder\plant_trait_databricks_bundle\plant_trait_runtime\models\foundation.py`
- `C:\Users\b0652085\PycharmProjects\PythonProject\DNA_AI_Decoder\plant_trait_databricks_bundle\config\datasources.yaml`

## Short port sentence

> Path B Databricks is past environment bring-up and now has three real primary-model packets on `PlantCaduceus_l32` (`full1003`, `seed=7`, `seed=99`), all showing the same stable mixed pattern: embedding better than structure-only on global `r2`, but worse on `spearman` and `within_group_r2`; the attempted `window_budget=128` follow-up failed on T4 CUDA memory, so the current honest verdict is stable mixed / ambiguous, not cleanly positive or blocked.
