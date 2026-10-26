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
| `hybrid_conv_maxtx_gc` | +0.4217 | +0.423, +0.445, +0.428, +0.409, +0.421, +0.432, +0.411, +0.420, +0.444, +0.399 |
| `max_tx` | +0.3968 | +0.397, +0.400, +0.392, +0.379, +0.405, +0.396, +0.399, +0.390, +0.411, +0.376 |

## Paired margins (same splits, per seed)

- **hybrid vs max_tx**: median +0.0339, 10/10 seeds, CI95 [0.0238, 0.0409]
- **hybrid+gc vs max_tx**: median +0.0302, 10/10 seeds, CI95 [0.0226, 0.0351]
- **CONTROL (shuffled conv) vs max_tx**: median -0.1232, 0/10 seeds, CI95 [-0.1275, -0.1187]

## Honest limits

- The combine is label-free and parameter-free (a rank-average), so it RANKS and does not dose -- the same deployability class as the inverse cell.
- Both inputs are out-of-sample on each seed's held-out peaks, so no second-level leakage is introduced; but the pool is MEASURED peak tiles, not generator output.
- One rung, one experiment, one combine rule. A negative does not show that no learned representation can add here -- only that THIS encoder, rank-averaged with THIS scalar, does not.
- `max_tx` and `conv` are both fitted by ridge on the same split, so their errors may be correlated through the split rather than through the biology.
