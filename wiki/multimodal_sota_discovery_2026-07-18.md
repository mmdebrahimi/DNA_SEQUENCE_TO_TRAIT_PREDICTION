# Ready-made multimodal SOTA (SaProt / VenusREM / …) vs our ESM2+ProSST — don't adopt

**Date:** 2026-07-18 · **Method:** paired per-protein abs-Spearman from ProteinGym's PRECOMPUTED columns
(N=95 assays with all columns present) · **Script:** `scripts/multimodal_sota_discovery.py` · **Data:**
`wiki/multimodal_sota_discovery_2026-07-18.json`

## Question + reframe (R2)

The directive: "discover ready-made multimodal SOTA (SaProt/VenusREM)". The decision this informs is whether
to ADOPT one of those off-the-shelf models over our validated 2-way `ESM2+ProSST` (+0.067 vs ESM2, 93%,
`wiki/prosst_lift_2026-07-18.md`). That is answerable **paired, per-protein, with ZERO install** — SaProt,
VenusREM, ESCOTT, S3F, ProtSSN are all on the ProteinGym leaderboard with precomputed columns. The heavy
install only matters if the answer is "yes, adopt one".

## Result — NONE of the ready-made SOTA beats our ESM2+ProSST (paired)

Paired vs our 2-way `ESM2+ProSST` (rank-average of the precomputed ESM2-650M + ProSST-2048 columns —
apples-to-apples; our own ProSST reproduces its column at 1.0), median |Spearman| over N=95:

| model | modality | median | Δ vs our 2-way | win | sign-p | beats us? |
|---|---|---:|---:|---:|---:|:--:|
| **our 3-way** ESM2+ProSST+GEMME | seq+struct+evo | 0.5411 | **+0.008** | 59/95 | 0.024 | ✅ (marginal) |
| **our 2-way** ESM2+ProSST | seq+struct | **0.5304** | baseline | — | — | — |
| VenusREM | struct+MSA-retrieval (leaderboard #1) | 0.5309 | −0.011 | 39/95 | **0.10 (n.s.)** | ❌ |
| SaProt-650M-AF2 | seq+Foldseek-3Di structure | 0.4988 | −0.030 | 32/95 | 0.002 | ❌ |
| ESCOTT | evolution+structure | 0.5049 | −0.029 | 26/95 | ~0 | ❌ |
| S3F-MSA | structure+MSA | 0.5027 | −0.016 | 33/95 | 0.004 | ❌ |
| ProtSSN-ensemble | structure GNN | 0.4749 | −0.045 | 13/95 | ~0 | ❌ |
| MSA-Transformer | evolution | 0.4456 | −0.069 | 19/95 | ~0 | ❌ |

**Our simple ESM2+ProSST rank-average matches or beats every ready-made multimodal SOTA on this subset,
paired.** VenusREM (the strongest) is statistically indistinguishable-to-slightly-worse (−0.011, p=0.10,
loses 56/95). The only thing that beats our 2-way is our own 3-way, and only marginally (+0.008 — consistent
with `wiki/three_way_lift_2026-07-18.md`).

## What they are + deployment cost (the reason not to adopt)

Even setting performance aside, the ready-made SOTA are **heavier to deploy** than our ESM2+ProSST (which
runs locally on `transformers`+torch with ProteinGym/ColabFold structure tokens):

- **SaProt-650M** — structure-aware PLM fusing residue tokens with **Foldseek 3Di** structural-alphabet
  tokens. Needs the Foldseek binary + a 3D structure (AF2/ColabFold) + pLDDT masking. *And it scores −0.030
  below our 2-way* — strictly worse, more infra.
- **VenusREM** — retrieval-enhanced structure+MSA hybrid; ranks **1st on the full ProteinGym leaderboard**
  (by the published aggregate metric over all 217 assays). Deployment: a conda env + **HMMER + EVCouplings**
  for MSA + downloading the a2m/a3m homology alignments **and the Uniref100 database**. For that cost it
  does **not** beat our 2-way paired on the shared subset.

## Verdict — DON'T adopt; our validated hybrid is at the frontier

- **Performance:** no ready-made multimodal model beats our ESM2+ProSST paired on the shared subset; VenusREM
  ties. Our 2-way is already competitive with the leaderboard.
- **Cost:** the ready-made models add Foldseek / HMMER+EVCouplings / Uniref100-DB infra for zero paired gain.
- **So:** keep the validated, locally-runnable **ESM2+ProSST** 2-way as the deployable hybrid. Adopting
  SaProt/VenusREM is not worth the deployment weight.

## Honest scope

- **VenusREM tops the FULL published leaderboard** (217 assays, aggregate 5-category-mean metric). This
  comparison is the **N=95 struct+MSA-available subset, PAIRED** — a different (stricter, per-protein) lens.
  So the honest claim is "on the shared subset, paired, none beats our 2-way; VenusREM ties" — NOT "we beat
  the #1 model everywhere". VenusREM may edge us on the full 217 by the aggregate metric; it does not justify
  its deployment cost for our use.
- Our-2-way / our-3-way here use ProteinGym's precomputed ESM2/ProSST columns (apples-to-apples with the
  ready-made columns); our own ProSST reproduces its column at 1.0, so our-own-2-way ≈ this.
- All from precomputed columns — no model was run, no install. `scripts/multimodal_sota_discovery.py`.

Frozen decoder surface byte-unchanged (`verify_lock OK`).
