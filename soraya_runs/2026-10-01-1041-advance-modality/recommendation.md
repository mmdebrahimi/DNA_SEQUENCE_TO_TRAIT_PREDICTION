# Recommendation — 2026-10-01-1041-advance-modality

## The gLM question is now closed on the technical side

Every lever measured, none pays:

| lever | verdict | artifact |
|---|---|---|
| zero-shot likelihood | fails (0.4485; BLOSUM62 matches it) | `hiv_esm_vs_catalog_2026-07-09` |
| scale | regresses 650M→3B→15B | `proteingym_esm2_650m_full_2026-07-09` |
| capacity over the genotype | loses 4/4, 2 significant | `hiv_context_vs_linear_floor_2026-09-30` |
| pretrained embeddings + head | PARTIAL, loses to catalog | `esm_supervised_head_result_2026-07-03` |
| **orthogonal modality** | **no room** | `hiv_onehot_unseen_headroom_2026-10-01` |

**Do not re-propose an estimator swap or a feature channel.** What ships is what the evidence supports: the
supervised one-hot blind-spot complement (3 genes, in the L2 doubt layer) + the LD imputer, with the
deterministic catalog primary.

## The only live lever is LABELS — and it is yours

The binding constraint was never model capacity; it is **catalog incompleteness per arm**, and all three
validated genes are **one virus**. A fourth arm needs a new *label*, not a new model.

Two things to know before deciding:

- Public-label expansion is a **recorded closure** with 10 rejection gates (`negative_results_map`), and
  the gates are runnable — `scripts/screen_candidate_gates.py` screens a candidate in seconds.
- The other two target-site cells with free labels **cannot** supply a fourth arm: sarscov2-mpro is
  TN-starved (37R/5S) and fungal has no free phenotype source.

So the realistic options are (a) accept the HIV-only scope and bank the gLM arc, or (b) acquire a
non-public wet-lab/clinical label source, which clears the gates by construction. (b) is an acquisition
decision, not an executor task.

## Still open from the previous run

**AUTHORITY FORK — restrict the gentamicin `rmt` rescue to E. coli?** Edits the frozen surface and
invalidates the v2 lock. Drafted position unchanged: **do not restrict yet** — the Klebsiella over-call
evidence (PPV 0.475) fails this repo's own source-diversity bar at 98.4% one-study, while NCBI-PD says
53R/0S across 12 BioProjects and clears it. The L2 `organism_scope` warning already ships the caveat at
call time.
