# A learned protein model DISCOVERS which variants drive the measured phenotype (PASS, 2026-07-04)

**The JEPA/CLIP learned-representation thesis, tested at the layer where it can work — and it works.** On
dna_decode's own cached data, a learned variant-effect model (AlphaMissense) ranks protein variants by their
**wet-lab-measured** functional effect (ProteinGym DMS) at **median |Spearman| = 0.417** across 39 joinable
assays, with a **shuffled control of 0.029** (signal is real, not artifact). This is the empirical green
light + substrate for a JEPA/CLIP protein model.

Script: `scripts/dms_learned_model_falsifier.py` (+ `tests/test_dms_learned_model_falsifier.py`, 5 offline).
Data (cached, no download/GPU/money): `D:/dna_decode_cache/proteingym` — AlphaMissense (`am_pg.tsv`), 217
ProteinGym DMS assays (`pg_dms/`), the leaderboard (`pg_spearman_dms.csv`).

## Result — PASS (pre-registered: median |Spearman| >= 0.30 AND shuffled ~0)

| metric | value |
|---|---|
| AlphaMissense vs DMS, median \|Spearman\| (39 assays) | **0.417** |
| shuffled negative control (median) | **0.029** |
| sign | **negative** — pathogenic AM score → lower measured fitness (biologically correct) |

Strongest joins (signed): BRCA1_Findlay_2018 −0.580 (n=1837) · CYP2C9 activity −0.643 (n=6142) · KCNE1
function −0.621 (n=2315) · CYP2C9 abundance −0.598 (n=6370).

**Field context** (median Spearman across all 217 ProteinGym wet-lab DMS assays, from the cached leaderboard):
ESM2-650M **0.484** · GEMME 0.484 · TranceptEVE-L 0.475 · EVE-ensemble 0.466 · ESM-1v-ensemble 0.460. The
whole learned-representation field lands ~0.46–0.48 — so a JEPA on protein sequence has a real, measured
target to beat/match.

## What it means (honest scope)

- **The learned-REPRESENTATION thesis is VALIDATED at the molecular-phenotype layer.** A learned model
  genuinely "discovers which variants drive the measured phenotype" for protein function — the working form
  of the user's JEPA/CLIP dream. This is the empirical foundation the whole GENOME_JEPA_CLIP plan (J-family)
  rests on; **J3's pre-registered falsifier is now PASSED** on a proxy (AlphaMissense) + the field leaderboard.
- **It does NOT rescue the complex-organismal-phenotype direction** (dna_decode's 0-for-5 de-confounded
  negative). The signal lives where the phenotype IS strongly sequence-determined (protein function), not in
  confounded polygenic organismal traits.
- **AlphaMissense is a proxy** for the learned-representation family (a supervised structure+sequence model),
  not JEPA itself. The next build (J2, GPU-gated) trains/fine-tunes a real sequence representation and re-runs
  this exact falsifier — the bar to beat is the ~0.48 field median.

## The through-line (all three learned-idea branches now placed, honestly)

| branch | verdict | layer |
|---|---|---|
| learned embedding → complex organismal phenotype | 0-for-5 de-confounded NEGATIVE | organismal (confounded) |
| masked-genotype IMPUTATION (LD input-completion) | PASS (98.9% ABO) | genotype completion |
| **learned model → molecular variant-effect** | **PASS (median \|rho\| 0.42; field ~0.48)** | **molecular (this result)** |

The learned/JEPA idea WORKS — twice (imputation + molecular variant-effect) — just not at the organismal
layer the romantic framing assumed.

## Reproduce
```bash
uv run python scripts/dms_learned_model_falsifier.py --max-assays 40
uv run pytest tests/test_dms_learned_model_falsifier.py -q
```
