# A 346-parameter biophysical model beats every learned arm on the 150 bp tiles — and subsumes GC

**2026-10-12.** Machine-readable twin: `wiki/promcalc_vs_learned_2026-10-12.json`. Reproduce with the
isolated env at `D:/tools/promcalc_env` (see *Isolation* below).

## Headline

| arm | median Spearman (10 seeds, leave-peak-out) | params | trained on this substrate |
|---|---|---|---|
| `gc` (the gate's own `ridge_gc`) | 0.3112 | 1 feature | yes |
| `encoder-w20` | 0.3163 | 2,657 | yes |
| `encoder-multi` (previous best) | 0.3234 | 13,217 | yes |
| `encoder-pyramid` | 0.3167 | 27,649 | yes |
| **Promoter Calculator, `fwd_max_tx`** | **0.3968** | **346** | **NO** |
| Promoter Calculator, `fwd_sum_tx` | 0.3953 | 346 | NO |
| GC + biophysical feature | 0.3950 | 2 features | the ridge is |
| replicate ceiling | 0.941 | — | — |

Paired on the same per-seed splits: **+0.0686 over GC (5/5 seeds with a logged GC)** and **+0.0401 over
`encoder-multi` (4/4)**. The learned arms' paired CIs against GC all include zero. **This one does not.**
It is the first thing in this project to separate from the GC baseline on this substrate — and it is a
mechanistic model that never saw the data.

## The pre-committed retraction check — PASSED, result KEPT

The plan fixed the rule before the data landed: *if the full-set median moves more than 0.02 from the
provisional figure, the non-random subset was carrying the result and the provisional is retracted.*

| statistic | provisional (14,885 tiles, file order) | full set (44,106) | move |
|---|---|---|---|
| `fwd_max_tx` | 0.3804 | **0.3968** | +0.0164 |
| `fwd_sum_tx` | 0.3856 | **0.3953** | +0.0097 |

Both under 0.02, and both moved **up**. Per-seed test sets went from ~1,600 rows to the protocol's
~11,000. **Result kept, not footnoted.** Sweep completed 44,046/44,046 with `err=0` in 74 min.

## The question that actually mattered: does mechanism ADD to composition?

**It does, and then some — it SUBSUMES composition.**

- GC + biophysical beats GC alone by **+0.0856, winning 10/10 seeds** (0.3820-0.4103 vs 0.2747-0.3347).
- **The feature-permutation control COLLAPSED: -0.0002, 2/10 seeds.** Destroying the sequence-to-mechanism
  pairing destroys the entire gain, so it is attributable to that pairing and not to adding any column.
- Spearman(GC, log1p mech) over all 44,106 tiles = **-0.5563**: correlated, nowhere near redundant.
- **But GC + mech (0.3950) is NOT better than mech alone (0.3968).** The direction of containment is
  one-way: the mechanistic feature carries GC's information and more, while GC adds nothing on top of it.
  Report this as **subsumption, not orthogonality** — an earlier draft of this analysis said "orthogonal
  axis", which the GC-adds-nothing row refutes.

**Baseline faithfulness is proven, not assumed.** The GC arm here is the gate's `ridge_gc` IMPORTED, not
reimplemented, and the 2-feature model mirrors its exact pipeline (`StandardScaler` -> `Ridge(alpha=1.0)`
on `log1p(expression)`). A replication of the 1-feature path reproduces `ridge_gc` to within 1e-9 on
every seed tested.

## Two methodological errors found by running this, both of which produced a spectacular wrong answer

1. **A raw-GC Spearman is the SIGN-FLIP of the gate's GC arm.** Raw GC correlates with expression at
   **-0.2994** on seed 0; the gate reports **+0.3000**. Same quantity to 0.0006, opposite sign, because
   `ridge_gc` FITS the direction. A first pass compared a 2-feature fit against raw-GC and reported a
   **+0.5681 gain, 10/10 seeds** — which was arithmetic on a sign flip, not an information gain.
2. **A shuffled-LABEL control cannot validate a rank metric here, and specifying it was a planning
   error.** Shuffling training labels randomises the fitted coefficient's SIGN but not its ranking
   power: the prediction stays a monotone function of GC, and GC genuinely correlates with the test
   target at ~0.30, so the control returns **+-0.41 by construction** rather than zero. The plan demanded
   it collapse; it cannot. The correct control permutes the **FEATURE** (destroying the
   sequence-to-mechanism pairing), and that one collapses to -0.0002 as required.

Both errors pointed the same way — toward a bigger, more exciting number — which is why the control
mattered more than the estimate.

## What the model is and what it is not

`Promoter_Calculator` from the Salis lab (LaFleur / Hossain / Salis). Scans both strands for every
candidate transcription start site, decomposes each into UP element, -35, spacer, -10, extended -10,
discriminator and initial transcribed region, and returns a Gibbs binding energy plus a calibrated
initiation rate. ~73 forward candidates per 150 bp tile, ~0.26 s/tile single-core.

**Reductions.** One measured value per tile versus ~146 candidate sites means the per-site predictions
must be reduced. `max_tx` (strongest single promoter) and `sum_tx` (total initiation) land within noise
of each other, and **`min_dG` is the EXACT negation of `max_tx` at all 10 seeds** — so they are ONE
feature, reported once, not three correlated results.

**Orientation is respected by the model, and this is a real check.** Forward-strand-only (0.3968) beats
best-of-both-strands (0.3800) at full n, as it did on the subset. The assay has a defined orientation, so
reverse-strand promoters should not drive the reporter — and measurably they add noise.

**This is a NEW measurement, not a reproduction.** Verified against the primary source: the preprint
states "346 interaction energy parameters" verbatim and reports R2 0.45-0.65 on external in-vivo sets,
but contains **NO Urtecho/Kosuri comparison and NO CNN baseline whatsoever**. A research pass had
surfaced a claim of "R2=0.60 on the Urtecho data, beating the CNNs" and correctly flagged it unverified;
reading the paper shows that 0.60 belongs to a genome-integrated modular-promoter set. Nobody has run
this model on endogenous genome-wide tiles before.

## Isolation (load-bearing)

Installed at `D:/tools/promcalc_env`, NOT in the project venv: the package pulls **sklearn 1.9.1** while
this project is pinned at **1.8.0**, so an in-project install would have moved sklearn under all 6,409
tests. Verified after install that the project still reports sklearn 1.8.0 / numpy 2.4.4.

## Honest limits

- **One substrate.** Endogenous 150 bp E. coli peak tiles with a 0.941 replicate ceiling. Says nothing
  about the fragment rung or any designed library.
  > **CORRECTED 2026-10-13.** This line read *"(ceiling ~0.46-0.47)"* for the fragment rung. **0.4583
  > is that rung's SINGLE-REPLICATE agreement, not its ceiling** — the same confusion
  > `wiki/mpra_substrate_census_2026-10-12.json`'s `ceiling_field_note` records itself making and
  > fixing. The published readout is the MEAN of two replicates, so the ceiling is the
  > Spearman-Brown step-up square-rooted: `sqrt(2r/(1+r))`. **Re-derived here rather than cited, and
  > the method validates against a known answer:** tiles `r11=0.7945 -> 0.9410`, reproducing the
  > published 0.941 exactly; fragments `r11=0.4583 -> 0.7928`. So the stale figure **understated the
  > fragment rung's headroom by +0.3345**, which matters because it is the rung this memo's
  > size-lever note points away from.
- **Still 42% of the ceiling.** 0.3968 against 0.941. Beating the learned arms is not solving the task;
  the gap remains the dominant fact.
- **No claim that learning cannot win here** — only that a 346-parameter mechanistic model beats the
  learned arms this project has actually built, at 44k training rows. Training-set size is unavailable
  **on this rung**: 46,713 rows is the entire file (the census records the same verdict verbatim,
  `size_lever: "EXHAUSTED"`).
  > **TWO CORRECTIONS 2026-10-13, both to this bullet, because it was quoted onward and did damage.**
  > (1) It called training-set size **"the literature's first-ranked lever (log-linear and
  > unsaturated)"**. That ranking is **UNSOURCED** — a grep finds it in exactly two places, this memo
  > and a script written a day later quoting this memo, which makes one unattributed claim look
  > corroborated. The claim is plausible and may well be right; it carries no citation, so it is
  > stated here as a hypothesis and the superlative is struck.
  > (2) **"Unavailable on this rung" is true and was then quoted WITHOUT the scoping**, as "rows are
  > not available". They are: **321,123 LB + 318,457 M9 sheared fragments** sit in the same cache
  > directory with their own owner loader (`genomewide.load_fragments`) — **7.3x** the peak tiles, at
  > a 0.7928 ceiling. A rung-scoped exhaustion note reads as a project-level wall the moment the rung
  > is dropped, and this is the exact error class the `blocked:substrate` rail was written for: *a
  > substrate blocker is a claim about what you looked for.*
- **Not a deployed-decoder change.** Whether this feature belongs in what the project ships is a
  surface-change question and is deliberately out of scope here.
- The sweep is single-condition (the peak-tile readout); the two fragment growth conditions are untouched.
