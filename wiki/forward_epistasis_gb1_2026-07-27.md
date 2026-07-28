# Epistasis validation: joint-ESM2 vs the additive-null multi-mutant (GB1, 2026-07-27)

**Does scoring the fully-mutated sequence JOINTLY beat the additive null the forward cell's
`predict_multi_effect` uses?** Tested on the canonical GB1 double-mutant DMS (Olson 2014 / ProteinGym
`SPG1_STRSG_Olson_2014`; 2500 sampled doubles, all with both measured singles), ESM2-650M on a free Kaggle T4.

## Result (n=2500 GB1 doubles, Spearman vs measured double-mutant fitness)

| predictor | ρ |
|---|---|
| ESM2 **additive** (WT-background single deltas summed — the deployed `predict_multi_effect` scoring) | 0.298 |
| ESM2 **joint** (each mutation masked in the OTHER mutation's background) | 0.308 |
| **joint − additive** | **+0.0096** |
| **measured-additive** (sum of the two real single-mutant fitnesses s1+s2) | **0.958** |

## Verdict

1. **ESM2 encodes only MINIMAL usable higher-order epistasis here.** Joint beats additive by +0.0096 — real
   and directionally consistent with Tsui et al. 2024 (ESM2 encodes *modest* higher-order epistasis) but
   marginal. **The additive-null `predict_multi_effect` is the right deployed choice**: joint scoring costs a
   forward pass per multi-mutant at inference and returns ~+0.01 Spearman for a zero-shot predictor — not
   worth it. The additive null is validated as the sensible default.

2. **The dominant signal is additivity itself.** Summing the two *measured* single-mutant fitnesses predicts
   the double at **ρ=0.958** — GB1 double fitness is ~96% additive, with epistasis as a small residual (the
   classic Olson 2014 finding). This is why the additive null is a strong baseline.

3. **Zero-shot ESM2 (~0.30) is far below the supervised measured-additive (0.958)** — expected: ESM2 predicts
   from sequence with no access to the measured single effects; the measured-additive is a supervised oracle.

## Scope / honesty (H8)

- **Confirmatory of the literature**, not a discovery — the R2 pre-bar scan established this up front
  (Tsui 2024: ESM2 encodes modest epistasis; additive is the standard ESM2 multi-mutant scoring; "simple
  baselines rival PLMs in mutation-dense design", bioRxiv 2026). The value is a **project-owned number that
  validates the additive-null scope** of the row-578 multi-mutant capability on real data.
- **+0.0096 is small** (n=2500); the honest claim is "joint captures modest epistasis, additive is the strong
  cheap baseline", NOT "joint clearly wins".
- One protein (GB1), doubles only. Higher-order / mutation-dense regimes (htFuncLib) may differ (the 2026
  bioRxiv shows simple baselines rival PLMs there too).

## Reproduce
Bundle `scripts`-style prep from ProteinGym `SPG1_STRSG_Olson_2014.csv`; ESM2-650M kernel on Kaggle T4
(`emanueleebrahimi/gb1-epistasis-bundle` + kernel run), $0.
