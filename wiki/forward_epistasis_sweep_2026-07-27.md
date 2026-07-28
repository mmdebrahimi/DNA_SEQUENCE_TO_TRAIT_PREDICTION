# Epistasis sweep: joint-ESM2 vs additive-null across proteins + mutation orders (2026-07-27, overnight)

**Extends the GB1-doubles epistasis test (row 579) to 5 proteins × mutation orders 2–6.** Question: does the
joint-ESM2 advantage over the additive null (a) vary by protein, (b) GROW with mutation order (the literature's
"denser = more epistasis" prediction)? ESM2-650M on a free Kaggle T4 (4500 multi-mutants, ProteinGym assays).

## Per-protein (overall)

| protein | seq | n | additive ρ | joint ρ | Δ (joint−additive) |
|---|---:|---:|---:|---:|---:|
| GFP_AEQVI_Sarkisyan_2016 (mutation-dense) | 238 | 1500 | 0.101 | 0.087 | −0.014 |
| HIS7_YEAST_Pokusaeva_2019 | 220 | 1500 | 0.342 | 0.347 | +0.005 |
| GRB2_HUMAN_Faure_2021 | 217 | 300 | 0.656 | 0.661 | +0.005 |
| SPG1_STRSG_Olson_2014 (GB1) | 448 | 300 | 0.228 | 0.232 | +0.004 |
| F7YBW8_MESOW_Aakre_2015 (ParD) | 93 | 900 | 0.543 | **0.260** | **−0.283** |

## Δ (joint−additive) by mutation order

| protein | k2 | k3 | k4 | k5 | k6 |
|---|---|---|---|---|---|
| GFP | +0.003 | −0.002 | +0.028 | +0.012 | +0.009 |
| HIS7 | −0.004 | +0.026 | +0.009 | +0.001 | +0.013 |
| F7YBW8 (ParD) | +0.054 | −0.086 | −0.126 | | |
| GRB2 / GB1 | +0.005 / +0.004 | | | | |

## Findings

1. **The "joint advantage grows with mutation order" hypothesis is FALSIFIED.** Per-order Δ is small
   (±0.03) and non-monotonic — no consistent growth with order on GFP or HIS7. ESM2's conditional
   (mutated-background) scoring does not systematically extract more higher-order epistasis at higher orders.

2. **The additive null is validated CROSS-PROTEIN.** For 4/5 proteins Δ ≈ ±0.005 — joint scoring adds
   essentially nothing on average. Extends the GB1-doubles result (row 579, +0.0096) to a 5-protein panel.
   **`predict_multi_effect`'s additive null is the robust default.**

3. **NEW: joint can be MUCH WORSE than additive, and degrade WITH order (F7YBW8/ParD: Δ = −0.283; k2 +0.054
   → k4 −0.126).** On this small 93-aa protein, masking a position in a heavily-mutated background produces
   worse predictions than the WT-anchored additive. **Hypothesis (unfalsified mechanism):** the multi-mutant
   background is far out-of-distribution for ESM2 (many substitutions in a short sequence), so its conditional
   log-probs become unreliable, whereas the additive keeps each single anchored to the real WT. This is a
   *robustness* argument FOR the additive null — joint scoring has a real downside risk OOD.

4. **Mutation-dense GFP is hard for ALL predictors** (additive 0.101, joint 0.087) — consistent with
   "simple baselines rival PLMs in mutation-dense design" (bioRxiv 2026); here even the additive barely works.

## Verdict (strengthens rows 578/579)

The additive-null multi-mutant is the **robust default across proteins and mutation orders**: joint scoring
adds ~0 on average, never reliably grows with order, and can hurt badly out-of-distribution. A joint/
epistasis-aware mode is not worth deploying for a zero-shot predictor.

## Scope / honesty (H8)
- Small per-order n (300); one protein (F7YBW8) drives the striking negative — a strong signal but a single
  case. The OOD *mechanism* is an unfalsified hypothesis, not established fact; the observed Δ's are data.
- ProteinGym assays, ESM2-650M zero-shot, $0 Kaggle T4. Reproduce: `scripts/epistasis_sweep_prep.py` +
  the sweep kernel (`emanueleebrahimi/epistasis-sweep-bundle`).
