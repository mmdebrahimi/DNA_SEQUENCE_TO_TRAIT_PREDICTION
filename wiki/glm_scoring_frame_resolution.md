# Resolving the scoring-frame confound behind the bacterial-versus-human gap

**2026-11-16.** Machine-readable companion: [`wiki/glm_scoring_frame_resolution.json`](glm_scoring_frame_resolution.json).

**PARTIALLY RESOLVED, and the decisive gate did not grade.** Read the decision first, then the correction, then the findings.

## Decisions needed

1. **The null ceiling in the Step 3 bar was mis-specified by me, and re-deriving it is a user authority call.** The gate returned `INDETERMINATE_NULL_NOT_CLEAN` on a shuffled-label arm reaching |0.0674| against a 0.05 ceiling. **CORRECTED 2026-11-17 — my first diagnosis of this was wrong on the mechanism.** I recorded it as `null_max_abs=0.05` inherited from a bar sized for ~11,200 pooled rows and applied to 3,822 active-frame rows. The row counts are not the defect. Measured on the exact active frame the gate used (14,913 tiles, reproduced to the row): tiles within a peak share **ICC 0.6931** of the label's variance, so the design effect is **5.900** and the 3,822-tile test set is only **648 independent units**. SE(rho|H0) = 1/sqrt(n_eff-1) = **0.0393**, so the frozen 0.05 bar sits at **1.27 SE** — and the statistic it bounds is a MAX of |rho| over 10 seeds, judged against a bar set as if it were one draw, with `abs()` making it two-sided against a hazard (leakage) that can only push a null UPWARD. **Measured consequence: a perfectly clean pipeline fails this criterion 89.6% of the time.** The observed 0.0674 is at the **41st percentile** of a clean null's max — 59% of clean runs are worse. All three breaches are NEGATIVE; the largest positive null is +0.0287, the 7th percentile in the leakage direction. So the INDETERMINATE verdict carries no information about leakage. **The margin bar is mis-sized the same way, which the first diagnosis missed:** +0.15 is 3.81 SE and was calibrated on the pooled 44,106-row run, so fixing only the null would have converted an honestly-ungraded result into a FALSELY-graded one. A design-derived bar (99th pct of a clean null's max-of-10) is **0.1291**; the median margin +0.2190 is 1.70x it and **9 of 10 seeds clear it**, the weakest (+0.1075) falling below. Re-deriving BOTH arms, stated before re-running, remains a USER AUTHORITY call — it edits a pre-registration.
2. **Whether the two eukaryotic artifacts should gain a sixth `honest_limits` entry.** This plan's position is a forward pointer only; both files are byte-identical. See *The published verdict is untouched* below.

## Correction carried forward

The **0.8194** single-threshold statistic reported earlier on 2026-11-16 — "87% of the achievable ranking signal", and the inference that the pooled Spearman is *overwhelmingly a detection metric* — **is retracted and remains retracted**. `rho = sqrt(3k(n-k)/(n^2-1))` is a closed form in the split proportion alone; at `n=44,106`, `k=14,913` it equals 0.8194 exactly, and uniform, Gaussian, Pareto and bimodal samples all return it to four decimals. The maximum is 0.866 at *any* median split of *any* continuous label.

The retracted value is **not deleted** from `wiki/glm_bacterial_scoring_frame_2026-11-16.json`. Deleting the evidence of an error is what lets it be re-derived. That file must **never be cited as provenance for a number**: the contract-number audit verifies by whole-number presence, so a contract citing it would find 0.8194 and certify it. No `cell_registry` contract cites it today and a test holds that line.

## Finding 1 — the label's reproducibility rises across its range, and this replaces the retracted claim

Bin on `RNA_exp_sum_1` (group 1's published summary, verified as the exact mean of its own two sub-columns on 44,106 of 44,106 rows). Score the Spearman between `RNA_exp_sum_2_1` and `RNA_exp_sum_2_2` **inside** each decile. The binning variable and the scored pair share no measurement, so bucketing cannot mechanically force the pair to agree.

| decile | n | within-group-2 agreement | implied ceiling |
|---|---|---|---|
| 0 | 4,411 | +0.1437 | 0.5013 |
| 1 | 4,411 | +0.1105 | 0.4462 |
| 2 | 4,411 | +0.1472 | 0.5066 |
| 3 | 4,411 | +0.1804 | 0.5528 |
| 4 | 4,411 | +0.2175 | 0.5978 |
| 5 | 4,411 | +0.3357 | 0.7090 |
| 6 | 4,410 | +0.3857 | 0.7461 |
| 7 | 4,410 | +0.4187 | 0.7683 |
| 8 | 4,410 | +0.5755 | 0.8547 |
| 9 | 4,410 | **+0.9323** | **0.9823** |

Trend **+0.9879**. **It must be read against the attenuation null, which reaches +0.4303 on its own** — binning on group 1 restricts group 1, which correlates with group 2, so range restriction alone produces substantial apparent slope. The observed trend is more than double that, and the mirror direction (bin group 2, score group 1) gives **+0.9515 against a null max of +0.4499**. That replication is free: it reuses the same four columns.

**The consequence reaches every bacterial number in this project.** The published **0.941 is not a single ceiling**. It is a pooled average over deciles whose own ceilings run **0.5013 to 0.9823**. An arm reported "against the replicate-derived 0.941 ceiling" is being judged against a figure that describes no part of the range.

## Finding 2 — the decisive gate is ungraded

Phase 1 reconciled exactly before the frame changed: `gc` seed 0 pooled **+0.3000** against the committed anchor, `encoder-multi` **+0.3236** against the published 0.3235, with **13,217 parameters at 11 epochs** — matching the published `epochs_ran_per_seed[0]`.

Phase 2 then ran 10 seeds on the active frame and the **null gate fired first by the bar's own ordering**, so nothing is graded. The numbers it produced are recorded in the JSON under `ungraded_numbers_DO_NOT_QUOTE_AS_A_RESULT` and are **not** reported here as a result. A bar overridden once it returns an inconvenient answer is not a bar.

What the diagnosis does and does not support: all three breaching null values are **negative** (−0.0599, −0.0622, −0.0674), 4 of 10 nulls are positive and 6 negative, and a shuffled-label model that anti-correlates has learned noise pointing the wrong way rather than fabricated signal. A conv is also a smooth function of sequence rather than an independent random vector, so its spread is legitimately wider than `1/sqrt(3822) = 0.0162`. None of that licenses changing the gate — it explains why a properly-derived ceiling for this frame would likely differ, which is the decision above.

## Finding 3 — the margin grows on the active frame

From `wiki/active_frame_rescore_2026-11-16.json`. No model trained; arm identity assigned by reconcile at **max|delta| 0.0000**, each file discriminating against the wrong series at 0.0126, split gate passing on all 10 seeds for both `n_test` and the exact held-out index set.

| arm | pooled margin | active margin | active GC | seeds |
|---|---|---|---|---|
| `asgiven-position-resolved-promoteratlas-150bp` | +0.1973 | **+0.3158** | +0.0389 | 10/10 |
| `promoteratlas-base-embeddings-frozen-150bp` | +0.1948 | **+0.3016** | +0.0389 | 10/10 |

Both arms beat GC by **more** on the active frame, because GC's pooled advantage was almost entirely detection. The readable quantity is the **margin under the same restriction**; the active-frame absolutes are range-restricted on the outcome and carry `not_a_ledger_candidate`.

**Do not merge these with the scratch numbers.** The 2026-11-16 scratch used a 4-mer ridge, not the published conv, and fit GC *within* the active subset (**−0.0464**) rather than on the pooled train set (**+0.0389**). Step 3 independently reproduced the scratch's −0.0464 with GC fit inside the active frame.

## Finding 4 — the human label is not similarly concentrated

`modal_bin_share` uses equal-**width** bins and is not scale-invariant, so the raw comparison is invalid: `expn_med` is a raw ratio on [0, 158.3] while human activity is already `log2(RNA/DNA)`. On the scale the models actually fit (`log1p`, which is `ridge_gc`'s target):

| label | modal-bin share | distinct-value share |
|---|---|---|
| bacterial, log1p (model scale) | **0.5704** | 0.9057 |
| bacterial, raw — *not comparable* | 0.9790 | 0.9057 |
| human HepG2 | 0.2203 | 0.1170 |
| human K562 | 0.2319 | 0.0559 |

The direction survives the scale correction at 2.5× (down from 4.4× raw), **but the scale-invariant measure runs the other way**: the bacterial label has concentrated mass with nearly all-distinct values, the human labels have spread mass but are heavily tied. The two labels are **differently pathological**, not one a milder version of the other.

**No human ceiling exists and none was invented.** The deposit carries **zero per-replicate columns across 6 files**, so no Spearman-Brown step-up is derivable; `human_ceiling` is `null` with a reason rather than filled with the bacterial 0.941, which both eukaryotic artifacts currently carry as context. This resolves the plan's R2 `[unverified]` flag in the negative. WTC11 refused at join rate 0.0000 (0 of 23,064 rows) and the refusal is recorded, not worked around.

## Which of the eukaryotic arm's four differences this resolves

**None of them.** The arm's first honest limit says it cannot say whether assay, element length (230 bp vs 150 bp), kingdom or library construction carries the effect. This plan adds a **fifth** candidate — the scoring frame — and shows that candidate is real on the bacterial side and **unmeasurable** on the human side. Assay, element length and kingdom are left open.

## The published verdict is untouched, and needs no reinterpretation

`verdict: SUBSTRATE_LIMITED` is unchanged in both `wiki/eukaryotic_arm_result_{HepG2,K562}_2026-11-01.json`. Its `verdict_detail` — *"the CI for the human gain lies entirely ABOVE the bacterial gain 0.0122"* — is a comparison of two measured numbers and stays true whatever this plan finds. The frame result lands **inside** `honest_limits[0]` (which difference carries the effect — now five candidates) and **strengthens** `honest_limits[1]` (the 0.0122 anchor treated as a fixed point, now frame-dependent). So no authority call is needed on the verdict field; the only one left is the sixth-`honest_limits` question above.

## Honest limits

- The decisive gate is **ungraded**. Findings 1, 3 and 4 do not depend on it, but the question the plan was built to answer has no graded answer.
- The active frame conditions on the outcome; only differentials under the same restriction are readable.
- `modal_bin_share` is not scale-invariant even on `log1p` — that transform was chosen because it is what the models fit, not because it is canonical.
- The curve's null reaches +0.4303 on its own, so the trend is interpretable **only** against it. Three per-decile ceilings below 0.55 rest on 4,411 rows each.
- One substrate, one organism, one assay series (GSE144621). Nothing here generalises to bacterial regulation at large.
- Three null replicates in the gate is not a null distribution, as the gate's own scope language says.

## Ledger

**Nothing was written to `soraya_runs/arms-ledger.json`, and that file is byte-identical to `HEAD`.** The plan asked for one criterion-level annotation key *and* an unchanged arm count of 16, and those two requirements cannot both hold: `len(criteria[<crit>])` **is** the arm count and every value is iterated as an arm dict (`scripts/declare_successor_arms.py:172-173`), so a criterion-level key reads as a 17th arm. The top-level alternative was implemented, and it broke a pinned key set — `tests/test_record_successor_arm_results.py:367` asserts `set(data) == {"schema_version", "criteria"}` — so it was reverted.

**There is no augment-only location in that file.** Recording the disclosure there requires editing a guard, which is a user-authority call and was not taken. The disclosure therefore lives in this artifact and its JSON companion only, which means **a reader who consults the ledger alone will still see the `ceiling:` criterion named after 0.941 with no indication that the figure is a pooled average over deciles whose own ceilings run 0.5013 to 0.9823.** That gap is stated here rather than closed.
