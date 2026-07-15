# Forward-cell dosage head — cross-organism generalization sweep + a rank-vs-dose finding (2026-07-15)

The dosage head was validated on one protein (PTEN). This sweeps it across the forward cell's tree-of-life
panel (E. coli / human / yeast / Arabidopsis), each scored by its best available method (cached ESM2 >
AlphaMissense-if-human > BLOSUM), to establish it as a real capability and to ask whether calibrated coverage
AND magnitude-informativeness generalize. All local/instant (no GPU). `scripts/forward_dosage_sweep.py`,
target 80% coverage, 15 shuffled 50/25/25 splits.

## Result — 10 proteins, 4 organisms, 2 kingdoms

| protein | organism | method | n | coverage | narrowing | pt Spearman | verdict |
|---|---|---|---:|---:|---:|---:|---|
| TEM-1 β-lactamase | E. coli | ESM2 | 4996 | 0.801 | **0.301** | 0.726 | CALIBRATED_DOSAGE |
| CcdB | E. coli | ESM2 | 1663 | 0.817 | **−0.313** | 0.492 | calibrated-UNINFORMATIVE |
| IF1 | E. coli | BLOSUM | 1367 | 0.798 | −0.165 | 0.176 | calibrated-UNINFORMATIVE |
| PTEN | human | ESM2 | 7260 | 0.802 | 0.130 | 0.517 | CALIBRATED_DOSAGE |
| TPMT | human | AlphaMissense | 3648 | 0.809 | 0.175 | 0.565 | CALIBRATED_DOSAGE |
| CYP2C9 | human | AlphaMissense | 6370 | 0.804 | 0.190 | 0.603 | CALIBRATED_DOSAGE |
| MSH2 | human | AlphaMissense | 16749 | 0.801 | 0.100 | 0.415 | CALIBRATED_DOSAGE |
| RL40A (ubiquitin) | yeast | ESM2 | 1253 | 0.826 | 0.088 | 0.530 | CALIBRATED_DOSAGE |
| SR43C | Arabidopsis | ESM2 | 889 | 0.791 | 0.236 | 0.586 | CALIBRATED_DOSAGE |
| MBD11 | Arabidopsis | BLOSUM | 1155 | 0.794 | −0.000 | 0.168 | calibrated-UNINFORMATIVE |

**Calibrated 10/10 · informative 7/10.**

## Two findings

1. **Calibration generalizes universally.** Held-out coverage sits at the 80% target on every protein across
   E. coli / human / yeast / Arabidopsis (range 0.791–0.826) — the split-conformal guarantee is organism-
   agnostic, so the forward cell can emit a *calibrated* magnitude interval for any organism.

2. **Dosage-informativeness is DISTINCT from rank quality (the non-obvious finding).** Interval narrowing —
   how much the score pins the magnitude vs the predict-the-mean marginal — does NOT simply track Spearman:
   - Strong methods on most proteins narrow the interval (TEM-1 ESM2 +0.30, SR43C ESM2 +0.24, human AM +0.10
     to +0.19).
   - **A moderate rank score can FAIL to narrow the dose interval**: CcdB-ESM2 ranks at Spearman 0.49 yet has
     *negative* narrowing (−0.31) — the conditioned interval is WIDER than the marginal at 80%. Weak BLOSUM
     scores (IF1/MBD11, Spearman ~0.17) narrow ~0 as expected.
   - **Reading:** "ranks variants well" (Spearman) ≠ "pins the magnitude" (narrowing). A decoder that orders
     edits correctly may still not tighten the calibrated dose interval — the DMS distribution shape governs
     whether conditioning helps at a given coverage. The informativeness gate surfaces this honestly (a
     negative narrowing is the check WORKING, not a bug).

## Status

The dosage head is now a validated, generalized capability: calibrated magnitude + prediction intervals
across 4 organisms / 2 kingdoms, with an informativeness gate that honestly separates "calibrated" from
"informative." 6 offline dosage tests. Frozen decoder surface byte-unchanged (`verify_lock OK`);
`dna_decode/forward` NON-frozen. Run: `uv run python scripts/forward_dosage_sweep.py`.
