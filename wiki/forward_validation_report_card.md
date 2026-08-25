# Forward-cell validation report card (standing trust surface)

Read-only roll-up of the forward variant-effect cell's DMS-validated numbers (the molecular analogue of
the AMR `decoder_validation_report_card`). DMS is the one place the project's label wall does not bind
(free, independent, per-variant magnitude labels). Zero-shot unless stated. **No aggregate headline** —
each capability carries its own honest tier.

| capability | regime | tier | metric | scope |
|---|---|---|---|---|
| blosum62 (deterministic, no deps) | B_molecular | DMS_VALIDATED_BENCHMARK_WIDE | ProteinGym median |Spearman| 0.2012 | n=209 assays; 61 below 0.15 — modest, the honest floor |
| esm2-650M (learned, universal) | B_molecular | DMS_VALIDATED_BENCHMARK_WIDE | ProteinGym median |Spearman| 0.49275 | n=194 assays; 153 above 0.3 — the sequence baseline at scale |
| modality hybrid (ESM2+ProSST[+GEMME]) | B_molecular | DMS_VALIDATED_CEILING | ESM2 baseline median 0.4926; hybrid ceiling ~0.547 (ESM2+GEMME+ProSST) | n=95 struct+MSA assays; +0.056 paired, win 90.5% (prereg bar met); RANKS not doses |
| alphamissense / esm-if / prosst (structure/pathogenicity) | B_molecular | DMS_VALIDATED_PER_PROTEIN | per-protein DMS Spearman (AM human-only; ESM-IF/ProSST structure tier) | 18 proteins in the method leaderboard; ESM-IF does NOT beat ESM2 on PTEN |
| multi-mutant additive-null (predict_multi_effect) | B_molecular | DMS_VALIDATED | GB1 doubles: joint-ESM2 beats additive by only +0.0096 | additive null validated as the robust deployed default; GB1 doubles ~96% additive (measured s1+s2 -> 0.958) |
| epistasis characterization (joint vs additive) | B_molecular | CHARACTERIZED | 5 proteins x orders 2-6: joint ~= additive (delta ~+-0.005); 'grows with order' FALSIFIED | novel: joint can be worse OOD (ParD WITHIN-ORDER delta -0.053, degrades with order); additive is robust. The pooled -0.283 was a mutation-order POOLING artifact -- corrected 2026-08-25, wiki/forward_epistasis_pooling_correction_2026-08-25.md |
| inverse (effect->edit proposal) | B_molecular | DEPLOYABLE_RANK_ONLY | rank/percentile inverse works 4/4 assays (beats null); learned oracle earns keep 3/4; magnitude NOT deployable | RANKS not doses (percentile pts); magnitude needs the target's own DMS (circular by construction) |
| regime router (predict_edit) | router | DEPLOYED | A determinant -> AMR catalogue / B molecular -> DMS methods / C organismal -> ABSTAIN | resistance NEVER routed to a likelihood model; organism-polygenic ABSTAINS (closed negative) |

## Honest scope
- **DMS-validated, zero-shot** — per-variant Spearman vs measured DMS; not a calibrated dose (the inverse is RANK-only; the hybrid RANKS, does not dose).
- **Regime B only** — resistance (A) routes to the frozen AMR catalogue (never a likelihood model); organism-polygenic (C) ABSTAINS (closed negative).
- **Multi-mutant = additive null** (validated robust across proteins/orders; joint epistasis adds ~0 and can hurt out-of-distribution).
- Sources are the committed `wiki/forward_*` artifacts; regenerate with `scripts/build_forward_report_card.py`.