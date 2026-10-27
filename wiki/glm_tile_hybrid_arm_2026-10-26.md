# Tile-rung hybrid arm (conv ⊕ max_tx) — 2026-10-26

**VERDICT: `HYBRID_ADDS`**

hybrid +0.0339 over max_tx, 10/10 seeds, CI95 [0.0238, 0.0409]

## The precondition, measured first

Spearman between the conv and `max_tx` prediction vectors, per seed: median **0.3862** (min 0.3590, max 0.4790).

Near 1 would mean the two scorers are redundant and a rank-average cannot add; near 0 means they are orthogonal, which is the condition under which the combine is measured to pay.

## Arms

| arm | median spearman | per-seed |
|---|---|---|
| `control_shuffled_conv` | +0.2703 | +0.279, +0.270, +0.268, +0.256, +0.280, +0.271, +0.282, +0.253, +0.290, +0.264 |
| `conv` | +0.3234 | +0.324, +0.353, +0.329, +0.323, +0.287, +0.313, +0.290, +0.299, +0.331, +0.360 |
| `hybrid_conv_maxtx` | +0.4254 | +0.427, +0.449, +0.432, +0.419, +0.416, +0.425, +0.410, +0.414, +0.449, +0.426 |
| `hybrid_conv_maxtx_decomp` | +0.4404 | +0.444, +0.459, +0.447, +0.426, +0.439, +0.441, +0.429, +0.428, +0.460, +0.433 |
| `hybrid_conv_maxtx_gc` | +0.4217 | +0.423, +0.445, +0.428, +0.409, +0.421, +0.432, +0.411, +0.420, +0.444, +0.399 |
| `hybrid_conv_maxtx_gc_decomp` | +0.4383 | +0.438, +0.454, +0.442, +0.420, +0.439, +0.443, +0.427, +0.430, +0.455, +0.416 |
| `max_tx` | +0.3968 | +0.397, +0.400, +0.392, +0.379, +0.405, +0.396, +0.399, +0.390, +0.411, +0.376 |

## Paired margins (same splits, per seed)

- **hybrid vs max_tx**: median +0.0339, 10/10 seeds, CI95 [0.0238, 0.0409]
- **hybrid+gc vs max_tx**: median +0.0302, 10/10 seeds, CI95 [0.0226, 0.0351]
- **hybrid+decomposition vs max_tx**: median +0.0472, 10/10 seeds, CI95 [0.0404, 0.0522]
- **hybrid+gc+decomposition vs max_tx**: median +0.0408, 10/10 seeds, CI95 [0.0374, 0.0467]
- **CONTROL (shuffled conv) vs max_tx**: median -0.1232, 0/10 seeds, CI95 [-0.1275, -0.1187]

## Does a THIRD channel stack on top of `hybrid_conv_maxtx`?

**STACKING VERDICT: `WEAK_STACK`**

hybrid_conv_maxtx_decomp [WEAK_STACK]: +0.0145 over hybrid_conv_maxtx, 10/10 seeds, CI95 [0.0108, 0.0172] -- REPRODUCIBLE (10/10 clears the 8-seed bar) but short of the +0.02 margin; hybrid_conv_maxtx_gc [NO_STACK]: -0.0039 over hybrid_conv_maxtx, 4/10 seeds, CI95 [-0.0095, 0.0026] -- this third channel carries nothing the measured 2-channel lacks; hybrid_conv_maxtx_gc_decomp [WEAK_STACK]: +0.0103 over hybrid_conv_maxtx, 9/10 seeds, CI95 [0.0039, 0.0156] -- REPRODUCIBLE (9/10 clears the 8-seed bar) but short of the +0.02 margin

| combo | verdict | margin over the 2-channel |
|---|---|---|
| `hybrid_conv_maxtx_decomp` | `WEAK_STACK` | +0.0145, 10/10 seeds, CI95 [0.0108, 0.0172] |
| `hybrid_conv_maxtx_gc` | `NO_STACK` | -0.0039, 4/10 seeds, CI95 [-0.0095, 0.0026] |
| `hybrid_conv_maxtx_gc_decomp` | `WEAK_STACK` | +0.0103, 9/10 seeds, CI95 [0.0039, 0.0156] |

**Predicted BEFORE the run:** decomposition-as-third is predicted to LOSE, and to lose by MORE than gc did. Measured pairwise orthogonality (10-seed median) puts decomposition <-> max_tx at +0.7812, against gc <-> max_tx's +0.5508, so decomposition is the MORE redundant of the two available third channels -- and gc-as-third already lost (0.4217 against the 2-channel's 0.4254). The expected verdict is therefore NO_STACK for both new combos.

**What a negative does and does not close:** a negative closes THESE FOUR SCORERS on THIS rung, not stacking in general. conv and gc are both composition-family; max_tx and decomposition are both readings of the same Promoter Calculator output, max_tx being a scalar reduction of it. So the four span only TWO modality families, and a channel from a NEW family is untested here. A negative is evidence about the available channels, never about the mechanism.


## Pairwise orthogonality — every channel pair any combo uses

Near 1 = redundant, so a rank-average cannot add. Near 0 = orthogonal, the condition under which the combine is measured to pay. Ordered most-redundant first.

| pair | median | min | max |
|---|---|---|---|
| `max_tx` ↔ `decomposition` | +0.7812 | +0.7451 | +0.8146 |
| `gc` ↔ `decomposition` | +0.6041 | +0.5828 | +0.6278 |
| `max_tx` ↔ `gc` | +0.5508 | +0.5478 | +0.5754 |
| `conv` ↔ `gc` | +0.4593 | +0.4100 | +0.6557 |
| `conv` ↔ `decomposition` | +0.4213 | +0.3890 | +0.5440 |
| `conv` ↔ `max_tx` | +0.3862 | +0.3590 | +0.4790 |


## Honest limits

- The combine is label-free and parameter-free (a rank-average), so it RANKS and does not dose -- the same deployability class as the inverse cell.
- Both inputs are out-of-sample on each seed's held-out peaks, so no second-level leakage is introduced; but the pool is MEASURED peak tiles, not generator output.
- One rung, one experiment, one combine rule. A negative does not show that no learned representation can add here -- only that THIS encoder, rank-averaged with THIS scalar, does not.
- `max_tx` and `conv` are both fitted by ridge on the same split, so their errors may be correlated through the split rather than through the biology.
- The STACKING question puts more than one superset combo against ONE pre-registered bar and the aggregate stacking verdict is the best of them. That is a multiple-comparison exposure: with k combos the chance one clears the margin by noise exceeds the bar's nominal rate, so read the per-combo block, not only the aggregate.
- Orthogonality here is Spearman between PREDICTION vectors. Two scorers can be rank-correlated and still carry complementary error structure; a rank-average is only sensitive to the rank relationship, which is why this is the measure the arm uses -- it is not a general claim that the channels share information.
- SCOPE on the direction finding, and it bounds what this arm showed: `max_tx` and `decomposition` are NOT two modalities, they are two readings of ONE Promoter Calculator output -- scripts/glm_tile_predictions.py passes the same decomposition dict to both `tile_features()` and `max_tx_scalar()`, the latter being a scalar REDUCTION of it (`log1p(max Tx_rate)`). So their +0.7812 redundancy is mechanically expected, and the pair that is most redundant while still stacking is a REDUCTION/SUPERSET pair. Pairwise prediction-vector Spearman is refuted as a screen ACROSS THAT PAIR; it is NOT shown invalid between genuinely different modality families, which remain untested here. Read the result as 'finer-grained access to a channel already present pays', not as 'a new kind of information pays'.
- The CI is a normal approximation on the MEAN with z=1.96, while every point estimate reported beside it is the MEDIAN -- two different estimands -- and at n=10 the correct critical value is t(9)=2.262, so these intervals are ~15% narrower than their nominal level. RE-DERIVED under t and the conclusion SURVIVES on all three combos (decomposition [0.0103, 0.0177] and gc+decomposition [0.0030, 0.0165] still exclude zero; gc still includes it), with sign tests corroborating independently (p=0.00098 / 0.0107 / 0.828). Nothing published here is wrong. It is NOT fixed because `paired()` lives in promcalc_decomposition_arm and is SHARED with the published parent arm, so correcting it moves parent numbers and needs its own pre-registration.
