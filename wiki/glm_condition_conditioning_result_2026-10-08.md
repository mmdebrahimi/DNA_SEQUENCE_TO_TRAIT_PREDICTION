# G-C: conditioning on growth medium adds nothing measurable — a VALID null

**Verdict `NO_GAIN`** (`wiki/glm_condition_conditioning_2026-10-08.json`, schema
`glm-condition-conditioning-v1`). 297,868 coordinate-keyed fragments shared between LB and M9, 10 seeds,
position-blocked splits, median of the per-seed per-medium-mean Spearman.

| arm | per-medium mean | pooled | role |
|---|---|---|---|
| `two_model` (comparator) | **+0.2032** | +0.3086 | two independent per-medium models |
| `partial_pooling` (**PRIMARY**) | **+0.2033** | +0.3085 | gain **+0.0001** against a bar of **0.02** |
| `interaction` | +0.2032 | +0.3086 | DEGENERATE — equals the comparator by construction |
| `indicator` | +0.2029 | +0.2932 | too weak by construction; the failure is the point |
| `shuffled` | +0.2026 | +0.2027 | validity control |
| `partial_pooling_shuffled` | +0.2026 | +0.2027 | **null is CLEAN** (gain −0.0006) |
| `data_matched` | +0.2035 | +0.3077 | row-count control |

The gain is **two orders of magnitude short** of the bar — not marginal. `null_clean=true` and
`registered_protocol_used=true` are stamped by the script into its own artifact, so decision D1's two
validity predicates are satisfied by a machine-readable field rather than a human claim.

**This is a null about CONDITIONING, not about condition.** The substrate genuinely carries
condition-specific signal and that was measured first: noise-matched cross-condition agreement 0.4062 vs
within-condition 0.4676 (gap −0.0614), disattenuated true cross-condition 0.7887, so roughly 38% of the
real variance is condition-specific. Something is there; a shared linear model over length-invariant 3-mers
does not capture it.

## Three defects found in the arms themselves, two of them mine

**1. The primary arm was DEGENERATE by construction.** A full interaction `[seq ‖ m ‖ seq×m]` with a binary
`m` spans exactly the same hypothesis space as two independent per-medium models (LB uses `w`, M9 uses
`w+v`). Verified on synthetic data whose truth genuinely has medium-specific coefficients: corr 0.9999998,
max prediction difference 0.0118 — and it reproduces on the real fragments, where `interaction` and
`two_model` agree to four decimals at the median (+0.2032 both). The first design made the primary arm a
medium INDICATOR, which can only learn a global offset and is *strictly less* expressive than the
comparator; correcting to a full interaction over-shot to exactly EQUIVALENT. Neither arm can show a gain.
Only **partial pooling** — penalising the condition-specific block harder than the shared block — sits
strictly between the two endpoints and can beat both.

**2. The pooling scale was a silent NO-OP, so all three conditioned arms were the same model.** The scale
was applied in `build_design`, i.e. *before* `StandardScaler`. Standardisation divides out any constant
column factor, so the standardised matrix was unchanged and `PARTIAL_POOLING_SCALE = 0.3` did nothing
whatsoever: measured r = 0.9999970 against the unscaled model, identical to seven decimals. The fix is
placement, not magnitude — the penalty is now applied *after* standardisation
(`fit_ridge(col_scale=…)`), where scaling a standardised column by `c` means the block is penalised by
`alpha/c²`. That makes `c=1` the no-pooling endpoint and `c→0` full pooling, verified monotone:

| scale | r vs `two_model` |
|---|---|
| 1.0 | 0.999997 |
| 0.3 | 0.999736 |
| 0.05 | 0.952 |
| 0.001 | 0.892 |

A test now pins that a pre-standardisation scale is inert and a post-standardisation one is not, so the
defect cannot return silently. **Both defects were caught by writing the tests, not by the run** — the
script would otherwise have produced a confident `NO_GAIN` from three arms that were secretly one arm.

**3. An MVP predicate was unsatisfiable.** Criterion 2 named `tests/test_glm_condition.py`, a file that was
never going to exist (convention here is `test_<script>`). It failed *closed* — pytest exits 4 on a missing
file, verified — so it would have blocked the MVP rather than passing vacuously, which is the safe
direction. This is the **second** mis-specified path in this family's bar, after criterion 3's
`_2026-10-07.json`. Corrected to the real path; the criterion itself is unchanged.

## The pooling strength is SELECTED, and the curve is FLAT

The pooling scale is the arm's one real hyperparameter, and asserting it would make the headline a
statement about an arbitrary constant — especially since the asserted 0.3 had never actually been
exercised. It is now chosen on an **inner position-blocked split of the TRAINING rows only**, from the
frozen grid `(1.0, 0.5, 0.3, 0.1, 0.03)`; the test split never informs the choice. The inner split is
blocked for the same measured reason the outer one is (91.5% of consecutive fragments overlap), and a test
pins that it is a partition of the training indices with no test index on either side.

**The whole curve is reported because the whole curve is flat.** Seed 0: `{1.0: 0.2122, 0.5: 0.2122,
0.3: 0.2122, 0.1: 0.2122, 0.03: 0.2123}` — a 0.0001 spread. The argmax moves between 0.03, 0.1 and 1.0
across the 10 seeds (`[0.03, 0.03, 0.03, 0.03, 0.03, 1.0, 0.1, 0.1, 0.03, 0.03]`), i.e. **the selection is
picking noise**. Read that as *pooling strength is not a lever on this substrate*, **never** as *0.03 is
the right amount of pooling*. Reporting the curve is what makes this visible instead of hidden behind an
argmax — and 1.0 is deliberately IN the grid so that "no pooling is best" is a sayable outcome.

### One seed's primary arm silently BECAME the comparator, and that is a consequence of the flat curve

Found by an independent coverage audit of this code, then verified against the committed numbers.
**Scale 1.0 means no pooling is applied at all**, so on a seed whose argmax lands there, `partial_pooling`
is byte-identical to the degenerate `interaction` arm — the one documented above as equal to the two-model
comparator by construction. That happened on **seed 5**, where `partial_pooling` and `interaction` are
exactly equal at **+0.2161**.

It does **not** overturn the verdict: 9 of 10 seeds selected 0.03 or 0.1 (genuine partial pooling), and a
collapsed seed pulls the primary arm *toward* the comparator, which is the direction of the `NO_GAIN`
result already reported. But it has to be visible, because it is the sharp form of "the argmax is noise":
when that noise lands on the no-pooling endpoint, the primary arm stops being partial pooling at all.

**Disclosed rather than tie-broken.** Re-pointing the argmax after seeing which scales it chose would be
editing the method post-hoc. `selection_collapsed_to_comparator` + a `n_seeds_collapsed_to_comparator`
count now ship in the script so a future run states it directly; the committed artifact predates those
fields, and its per-seed `selected_scale` list is where a reader can count the collapses today.

## A partial mechanistic account, and one caveat upgraded from asserted to demonstrated

**The design is rank-deficient by construction, and sklearn says so** (`LinAlgWarning`, rcond ~1e-8). The
64 3-mer *frequencies* sum to 1.0 for every row. That is an **affine** constraint, so the raw block is full
rank 64 — but `StandardScaler` **centres** each column, which converts it into an exact linear dependency:
measured rank **63 of 64**, condition number **9.3e14**. Ridge's alpha is therefore doing load-bearing work
on every arm, which partly explains the flat curve: an extra per-block penalty has little room to move a
fit already heavily regularised on a near-collinear basis. Not a defect — ridge solves it, and the
`two_model` baseline reproduces the published numbers exactly — but it bounds how much the pooling lever
could ever have done here.

**The pooled metric's excess over the per-medium mean IS the between-media offset — now measured, not
assumed.** The honest-limits list previously asserted this from a level comparison (pooled 0.310 vs
per-medium mean 0.197). The shuffled control demonstrates it: destroying the medium label collapses pooled
to **0.2027**, essentially equal to its own per-medium mean of **0.2026**. The entire pooled-vs-per-medium
gap is the medium offset, which is exactly why `per_medium_mean` is the pre-registered primary — a gain
visible only in `pooled` would be an offset effect, not the conditioning this family is about.

## The engineering finding that made 10 seeds affordable

Measured before optimising: 3-mer feature extraction cost **~181 s** per arm pass over the full fragment
set × 2 media, against **~4 s** for the ridge fit it feeds. So **~97% of the sweep's cost was recomputing a
pure function of the sequence** — the features are identical across arms (which differ only in how the
medium flag enters) and across seeds (which differ only in the split). One pass per medium, cached: the
one-time cost is **385 s**, and a 10-seed × 7-arm sweep went from ~3.5 h to ~35 min. That is what made
decision D2 (full data, multi-seed median) affordable rather than a compute trade, and the refactor is
**bit-faithful** — `two_model` returns +0.2013 / +0.1916 on seeds 0 and 1, identical to the pre-refactor
smoke, reproduced on three separate runs.

## Honest limits

- **TWO conditions only.** A null here rules out conditioning *as tested, on two growth media*; it is weak
  evidence against conditioning in general. Whether that closes the family is an acceptance-bar call and
  therefore the user's.
- Conditioning is tested as sequence × medium **interactions over 3-mer features**, which is not the same
  as a model whose **representation** is modulated by condition. A null here does not foreclose that.
- The substrate ceiling is low (LB **0.7925** / M9 **0.8040**) because 97% of fragments carry a single
  barcode, so absolute numbers are small by construction and must be read against those per-medium
  ceilings, never against 1.0.
- Only length-invariant features apply (48–499 bp, 359 distinct lengths), so this is **not** comparable to
  the designed grid's positional-one-hot number.
- `min_gain = 0.02` is **ASSERTED**, not derived. The per-medium ceilings that contextualise it are derived.
- The flat curve has a *partial* account (rank deficiency + ridge). It is not established that collinearity
  is the whole reason pooling does not help.

Reproduce: `uv run python scripts/glm_condition_conditioning.py --seeds 10 --no-download`
(~35 min, needs `D:/dna_decode_cache/mpra`). Tests: `uv run pytest tests/test_glm_condition_conditioning.py -q`
(21 tests, offline).
