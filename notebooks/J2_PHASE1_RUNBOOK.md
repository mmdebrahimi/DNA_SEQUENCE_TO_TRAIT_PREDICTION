# J2 Phase 1 — get the real ESM-2 number on a FREE GPU (Kaggle / Colab)

**Goal:** run the REAL learned protein model (ESM2-650M, zero-shot masked-marginals) over the joinable
ProteinGym DMS assays → the real median |Spearman| (target ~0.48, the published field number). This is
**inference only** — free MIT weights, ~3–4 GB VRAM, fits a free T4/P100 (16 GB) easily. No training, no
money, no cross-machine relay.

**Honest scope (unchanged):** Phase 1 mostly *re-confirms* a known field number (the DMS falsifier already
got 0.417 via an AlphaMissense proxy — the substrate is proven). The real bet is **Phase 2** (fine-tune /
JEPA to *beat* 0.48), which is genuine training and worth gating on "does beating it power a decoder
feature." Run Phase 1 because it's free + nails the real substrate number end-to-end on our own model.

Runner: `notebooks/j2_phase1_esm2_proteingym.py` (self-contained — upload just this one file; scoring core
is drift-guarded against `scripts/esm_zeroshot_dms.py` by `tests/test_j2_phase1_notebook.py`).

---

## Step 0 — get the ProteinGym data into the kernel

The runner needs `DMS_substitutions.csv` (or `pg_reference.csv`) + the per-assay DMS CSVs. **Best path: you
already have the exact validated data on `D:/dna_decode_cache/proteingym`** (217 assays, all present —
verified). Upload it once as a **private Kaggle Dataset**:

1. Kaggle → **Datasets → New Dataset** → drag the `D:/dna_decode_cache/proteingym` folder (it has
   `pg_reference.csv` + `pg_dms/DMS_ProteinGym_substitutions/`). ~1 GB, well under the 20 GB limit.
2. In your notebook → **Add Data** → attach it. It mounts at `/kaggle/input/<your-slug>`.
3. Pass `--data-dir /kaggle/input/<your-slug>` (the runner tolerates both the `pg_reference.csv`+`pg_dms/`
   layout and the official `DMS_substitutions.csv`+`DMS_ProteinGym_substitutions/` layout).

Alternatives: **Colab** → upload the same folder to Google Drive, `drive.mount`, point `--data-dir` at it.
**Fetch fallback** → `--fetch` downloads the reference file from the official ProteinGym GitHub, but the
large per-assay CSVs still need attaching (reference alone can't score) — see
[ProteinGym repo](https://github.com/OATML-Markslab/ProteinGym) / its Resources page, or the HF mirror
[ICML2022/ProteinGym](https://huggingface.co/datasets/ICML2022/ProteinGym).

---

## Path A — single free kernel (simplest)

**Kaggle** (Notebook → Settings → Accelerator = **GPU T4 ×2** or **P100**; Internet ON for the weight
download). One cell:

```python
!pip -q install "transformers>=4.40" "torch" 2>/dev/null   # usually pre-installed on Kaggle
!python /kaggle/working/j2_phase1_esm2_proteingym.py \
    --data-dir /kaggle/input/<your-proteingym-slug> \
    --model facebook/esm2_t33_650M_UR50D \
    --out /kaggle/working/j2_full.json
```
(Upload `j2_phase1_esm2_proteingym.py` via **Add Data → Upload** or paste it into a cell as a file.)

**Colab** (Runtime → Change runtime type → **T4 GPU**):
```python
from google.colab import drive; drive.mount('/content/drive')
!python j2_phase1_esm2_proteingym.py --data-dir /content/drive/MyDrive/proteingym --out j2_full.json
```

**Expected:** a per-assay table then `median |Spearman| = 0.4x (shuffled ~0.0x)`. **PASS ≥ 0.45**;
**≥ 0.48 = matches the field.** Wall-clock ~30–90 min on a free T4 (proteins > 1022 aa are skipped, not
windowed — a known follow-up that lifts the score slightly).

---

## Path B — DISTRIBUTE across two free kernels (≈ half the wall-clock)

Run **shard 0/2 on Kaggle** and **shard 1/2 on Colab AT THE SAME TIME**, then merge (merge needs no GPU):

```python
# Kaggle kernel:
!python j2_phase1_esm2_proteingym.py --data-dir /kaggle/input/<slug> --shard 0/2 --out shard0.json
# Colab kernel (simultaneously):
!python j2_phase1_esm2_proteingym.py --data-dir /content/drive/MyDrive/proteingym --shard 1/2 --out shard1.json

# then, anywhere (download both JSONs to one place; no GPU needed):
!python j2_phase1_esm2_proteingym.py --merge shard0.json shard1.json
# -> MERGED: N assays, median |Spearman| = 0.4x (shuffled 0.0x) -> PASS
```

Sharding is deterministic + strided (`sorted(ids)[i::n]`) → the two halves are disjoint and cover all
assays, so the merged median is the true full-cohort number. Stack more kernels with `--shard i/3`, `i/4`…
(Kaggle gives 30 GPU-hr/week; Colab ~15–30; genuinely free, no card.)

---

## What you paste back
The final line: `median |Spearman| = 0.4x (shuffled 0.0x)` (single kernel) or the `MERGED: …` line
(distributed). That's the real J2 Phase 1 number. If it lands ≥ 0.48 we've matched the field on our own
model; then the Phase 2 "beat it" decision (GPU-heavier, possibly money) is the next fork.

---

## Phase 2 — BEAT 0.48 on free compute (added 2026-07-08)

Phase 1 re-confirms the ~0.48 baseline. **Phase 2 beats it**, still free (inference + ensembling, no
training). Three levers now in the runner (all CPU-tested; see `wiki/j2_phase2_beat_0.48_plan_2026-07-08.md`
for the pre-registered plan + bar):

1. **ESM2-3B in fp16** (`--model facebook/esm2_t36_3B_UR50D --dtype float16`) — the certain lift
   (650M ≈0.48 → 3B ≈0.51 published); fits a free T4 16 GB.
2. **Long-protein windowing** (`--long-mode window`) — scores proteins >1022 aa (previously dropped) in a
   centered window. Completeness, not a guaranteed median lift.
3. **Ensemble** (650M + 3B) — run each with `--keep-scores`, then combine on CPU with `--ensemble-merge`.

```python
# The certain beat (one GPU cell):
!python j2_phase1_esm2_proteingym.py --data-dir /kaggle/input/<slug> \
    --model facebook/esm2_t36_3B_UR50D --dtype float16 --long-mode window --out j2_3b.json

# Ensemble (two GPU cells + one CPU merge):
!python j2_phase1_esm2_proteingym.py --data-dir /kaggle/input/<slug> --model facebook/esm2_t33_650M_UR50D \
    --dtype float16 --long-mode window --keep-scores --out j2_650m.json
!python j2_phase1_esm2_proteingym.py --data-dir /kaggle/input/<slug> --model facebook/esm2_t36_3B_UR50D \
    --dtype float16 --long-mode window --keep-scores --out j2_3b.json
!python j2_phase1_esm2_proteingym.py --ensemble-merge j2_650m.json j2_3b.json   # no GPU
```

**PASS (pre-registered):** median |Spearman| **> 0.48** on the full joinable cohort, shuffled < 0.05.
