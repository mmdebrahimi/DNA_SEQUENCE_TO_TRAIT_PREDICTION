# Multimodal Family C + epoch synthesis (2026-07-27)

**`--until-mvp` Family C ("extend the forward cell with a lifting modality") verdict: `blocked` — on TWO
independent walls (R4 owned-elsewhere + R2 closed-negative). Do NOT build.**

## The two walls (v1.12 classification)

1. **R4 — owned-elsewhere.** `dna_decode/forward/` is **Soraya-fwd's owned frontier** (session-board:
   "Soraya-fwd owns `forward/`"). Any modality build inside it collides with a live sibling session. →
   hand-off/verify lane only, never duplicate-build.

2. **R2 — closed negative (empirically exhausted).** Soraya-fwd already ran the full modality-hybrid
   exploration and closed it WITH DATA on ProteinGym paired comparisons:
   - **ESM2 (sequence) + ProSST (structure) = the validated frontier**: +0.067 median, win 52/56=93%,
     sign-p 1e-11, positive on every phenotype.
   - **+GEMME (evolution) 3-way**: does NOT beat the 2-way (+0.0035, win 31/56=55%, sign-p 0.50 n.s.;
     structure already captures the lift, evolution on top is redundant / cancels on Stability).
   - **MSA-Transformer (evolution)**: reproduces PG's column but does NOT lift the hybrid (n.s.).
   - **SaProt / VenusREM / ESCOTT / S3F / ProtSSN (ready-made SOTA)**: NONE beats ESM2+ProSST paired
     (VenusREM #1 ties -0.011 p=0.10; all others worse) — and cost far more to deploy.

   Adversarial check (per the "brainstorm a fresh negative before committing" rail): sequence / structure /
   evolution are the three orthogonal protein modalities for variant-effect, ALL covered; cross-modality
   options (SaProt=seq+struct, VenusREM=struct+MSA) are tested-and-tied. **No untried orthogonal modality
   with a plausible lift remains.** A 4th modality (MD dynamics) is neither free nor evidenced to lift DMS.

## Multimodal epoch synthesis (the honest terminal)

The multimodal ambition splits exactly along the project's R2/R3 regime boundary — and both halves are now
RESOLVED:

| Multimodal half | Regime | State |
|---|---|---|
| **Molecular** (variant-effect: DNA/protein sequence ⊕ structure ⊕ evolution → molecular phenotype) | R2 (learned wins) | **DONE + CLOSED** — the forward cell's ESM2+ProSST hybrid is built/validated/deployed (Soraya-fwd); the modality space is exhausted with paired negatives |
| **Organism** (DNA-encoder + 2nd modality → organism/individual phenotype) | R3 (closed negative) | **CLOSED** — Family A falsifier (2026-07-27): the DNA arm has variant-effect signal (eQTL sign auROC 0.80) but ties a linear baseline cross-individual; inherits the R3 population-structure wall |

**Therefore the multimodal epoch's reversible/free executor frontier is genuinely exhausted.** The molecular
multimodal that works is shipped; the organism multimodal is closed. The ONLY unrealized multimodal is
organism-level, and it is gated on the SAME acquisition fork as unlocking R3: a dbGaP / UK Biobank (or
equivalent) paired controlled-access dataset + a de-confounded within-ancestry design — an AUTHORITY + MONEY
fork, not executor work.

## Recommendation

- **Bank Family C** (already done + closed + owned-elsewhere — no reversible build remains).
- The multimodal north-star is realized in the regime where it can be (molecular) and closed in the regime
  where free data can't support it (organism).
- The single move that would extend multimodal past this is the **organism-paired-data acquisition** — a user
  authority/money decision (draft anchor already exists: `wiki/data_acquisition_voi_memo_2026-07-18.md` +
  the label-acquisition memo series).
