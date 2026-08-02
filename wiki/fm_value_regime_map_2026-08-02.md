# When does a foundation model add value for genotype→phenotype? — the regime map (2026-08-02)

The definitive, grounded answer to the question this project (and this session) kept circling. Every
regime below is backed by a REAL, verified result — not speculation. The FM's value is not "yes" or "no";
it depends entirely on the **regime**, and the three regimes have three different answers.

## Regime A — self-supervised RECONSTRUCTION proxy (mask the genome, guess it back)

**Question:** does a frozen DNA foundation model's fill-in-the-blank skill capture genome *function*?
**Answer: NO — it's calibration, not a decoding signal. CLOSED negative.**
Evidence (this session): NT-v2 masked reconstruction on dog canFam4. After fixing a leaky baseline + using
per-base marginals, NT is at ~parity with a low-order Markov chain; scale helps only against a weak
baseline; and the skill is ANTI-structured — better on repetitive intergenic DNA, *worse* than Markov in
coding regions (`wiki/dog_nt_f2_strata_2026-07-31.md`). The FM has the genome's "spelling," not its
"meaning." **FM does not add usable value here.**

## Regime B — IN-DISTRIBUTION prediction (predict a trait where training data spans the variants)

**Question:** in a panel where you have labeled examples across the genetic variation, does an FM beat a
simple model?
**Answer: NO — the genotype is a SUFFICIENT STATISTIC; a simple model is near-optimal.**
Evidence (this session): the confound-free yeast (Bloom, n=1008) + mouse (BXD) crosses. Genomic prediction
DECODES quantitative traits (yeast 12/12; mouse brain weight r=0.57 — and this *generalizes* across
kingdoms). A nonlinear model captures epistasis and beats linear ridge *where epistasis is strong* (yeast
26/46, +0.16 on Maltose) but not on additive traits (mouse). Critically, an FM sequence-embedding of a
variant's context can only RE-ENCODE the allele identity the marker already carries — it adds no genetic
information for these individuals (`wiki/yeast_bloom_layer2_verdict_2026-07-31.md`). This mirrors the AMR
track: the deterministic curated catalog beats learned scorers in-distribution. **FM does not add value
here; optimize the model (nonlinearity) instead — and only where the architecture is epistatic.**

## Regime C — TRANSFER / zero-shot to UNSEEN variants (no training label for the variant)

**Question:** predict the effect of a variant you have NO labeled example of — where the marker/panel
gives you nothing.
**Answer: YES — this is the one regime where the FM's pretrained knowledge is the only source of signal,
and it genuinely wins.**
Evidence (the `forward`/ProteinGym cell, re-verified 2026-08-02): ESM2-650M zero-shot scores **median
Spearman 0.490 over 217/217 ProteinGym deep-mutational-scanning assays** (`our_median_spearman = 0.49`,
reproducing the published number), far above naive baselines. Combining ORTHOGONAL modalities
(ESM+GEMME+ProSST) lifts past the single-model ceiling (+0.056, wins 90.5% of proteins paired;
`wiki/forward_modality_hybrid_2026-07-17.md`). Plus clinical-variant-effect + HIV/Mpro extensions. **FM
ADDS REAL VALUE here — it is the correct tool precisely when there is no in-panel label to learn from.**

## The one-line map

| regime | is the FM the right tool? | why | project evidence |
|---|---|---|---|
| A: reconstruction proxy | **No** | measures spelling, not function | dog NT recon (closed neg) |
| B: in-distribution prediction | **No** | genotype is a sufficient statistic; simple model near-optimal | yeast+mouse crosses; AMR catalog |
| C: transfer / unseen-variant | **YES** | pretrained knowledge is the only signal source | ProteinGym ESM 0.49; modality-hybrid |

## Why this resolves the whole arc
The recurring "does the AI add value / why give up on the FM" tension dissolves once you split by regime:
the FM was being pointed at Regimes A and B (where it *can't* win by construction) and judged a failure;
its real home is Regime C, where the project's `forward` cell already shows it winning. **The honest
design rule: use deterministic catalogs / simple genomic prediction in-distribution (Regimes A/B); reach
for the foundation model only in the transfer regime (C), where nothing else has signal.**

## The one genuine OPEN frontier
Regime C is validated at the PROTEIN level (ESM on amino-acid variants). The untested extension is a
**DNA-level transfer test** — a DNA foundation model (NT/Caduceus) zero-shot on unseen *non-coding /
regulatory* DNA variant effects, where protein-level ESM doesn't apply. Prior is guarded (NT was
calibration-not-signal in Regime A), it needs a free DNA variant-effect benchmark + GPU (Kaggle), and it
is a separate substantial build — the natural next mission, not part of this map.

Verified: ProteinGym re-derivation (`wiki/proteingym_esm2_650m_full_2026-07-09.json`). Frozen AMR/forward
surfaces byte-unchanged (read-only synthesis).
