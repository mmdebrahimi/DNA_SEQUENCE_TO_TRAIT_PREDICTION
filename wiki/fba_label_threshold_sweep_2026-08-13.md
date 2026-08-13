# The carbon-panel headline survives every label cutoff — and a confidence filter I wrote destroyed it (2026-08-13)

> Closes the open gap named in `wiki/fba_conditional_carbon_2026-08-12.md`: *"the 217-gene set rests on
> one inherited cutoff (`fit < -2`) … no threshold sweep was run."* One has now been run.

## The question

Every number in the 25-carbon-source conditional result rests on one inherited label rule: a gene is
essential in a condition iff its mean RB-TnSeq fitness is `< -2` (Bernstein 2023, reused for
comparability with the shipped Keio validation). Nobody had ever checked whether the headline moves under
a different bar. It also rests on a threshold applied to *noisy* measurements — RB-TnSeq ships a
per-measurement t-statistic that the loader never read until yesterday.

**Efficiency note that made this cheap:** the FBA deletion calls do not depend on the label rule at all —
only *which genes get scored* does. So the deletions run once over the union of every gene conditionally
essential at any setting (391 genes × 25 conditions), and all 35 settings re-score from that cache.

## Result 1 — the headline is not a lucky cutoff

Sweeping the fitness bar with no confidence filter:

| fit bar | n cond.-ess. | exact-set | per-cell | null | **lift** | constant |
|---|---|---|---|---|---|---|
| < −1.0 | 318 | 27 | 0.7197 | 0.6801 | **+0.0396** | 87.1% |
| < −1.5 | 261 | 24 | 0.7146 | 0.6690 | **+0.0456** | 86.6% |
| **< −2.0 (shipped)** | **217** | **23** | **0.7368** | **0.6623** | **+0.0745** | **84.8%** |
| < −2.5 | 218 | 22 | 0.7561 | 0.6582 | **+0.0979** | 85.3% |
| < −3.0 | 195 | 17 | 0.7774 | 0.6361 | **+0.1413** | 86.1% |

**Lift over the constant null is positive at every setting, and rises monotonically as the label bar
tightens** — the shipped `< −2` sits in the middle, not at a maximum. A cleaner label set makes the model
look *better*, which is the direction a real signal should move. The constant fraction is stable at
85–87% throughout, so "the model mostly does not switch" is not a cutoff artifact either.

Adding the per-cell confidence bar changes almost nothing (20 settings, all positive lift):

| | n | exact | lift |
|---|---|---|---|
| fit < −2.0, no t bar | 217 | 23 | +0.0745 |
| fit < −2.0, \|t\| ≥ 2 | 219 | 23 | +0.0689 |
| fit < −2.0, \|t\| ≥ 3 | 214 | 23 | +0.0793 |
| fit < −2.0, \|t\| ≥ 4 | 209 | 21 | +0.0955 |

**Verdict: the headline survives all 20 sound settings.** This was the single largest open question about
the substrate, and it closes cleanly.

## Result 2 — the confidence filter I shipped yesterday was anti-selective, and the sweep caught it

The first version of `min_abs_t` required `|mean t| ≥ bar` in **every** condition (the complete-row rule,
carried over from how fitness values are handled). Under it, **15 of 15 grid settings collapsed to 100%
constant predictions and ZERO commitments**, with lift going sharply negative (−0.08 to −0.46).

That is not the labels being noisy. It is the filter removing the phenomenon:

> A conditionally essential gene is confidently essential in **one** condition (large |t|) and
> confidently **neutral** in the other 24 (t ≈ 0). Requiring high |t| *everywhere* therefore selects
> against exactly the switchers, and keeps genes with strong signal in every condition — which are, by
> construction, the ones that never switch.

The tell was visible in the numbers before the reasoning: as the fitness bar tightened, the surviving
count under `|t| ≥ 2` went *up* (24 → 42 → 58), not down. A filter whose yield grows when you make the
other criterion stricter is selecting on something other than what you think.

`min_abs_t_mode="per_cell"` is now the default and the correct semantics: a low-confidence measurement
cannot **support** an essential call, but it is no reason to discard a gene whose other conditions are
cleanly measured. `all_conditions` is kept, documented as anti-selective, and pinned by a test.

**The reporting trap this exposed:** the first verdict function pooled all axes and returned
`HEADLINE_IS_CUTOFF_DEPENDENT (5/10 powered settings beat their null)`. Every one of the 5 failures sat on
the broken axis. Pooling a degenerate instrument with a sound one manufactures a false negative — the
exact mirror of the retracted-headline error from yesterday, in the opposite direction. The verdict is now
axis-aware and names a degenerate instrument as such.

## Result 3 — what actually carries the conditional signal

Splitting every essential call behind a commitment by mechanism:

- **`infeasible`** — no growth value at all; the ATPM maintenance floor cannot be met without the gene.
  A boolean signal with no threshold involved (see `wiki/fba_infeasibility_finding_2026-08-13.md`).
- **`sub_threshold`** — a real, finite growth ratio that falls below 1% of wild type. **The only kind of
  call a threshold retune could ever move.**

Across every powered setting, the infeasible share of commitment calls is **53.9–66.2%**.

This **refines** yesterday's finding rather than repeating it. At the *gene* level 32 of 33 commitments
touched an infeasible cell; at the *call* level **34–46% are genuine sub-threshold crossings**. So
the model is *mostly* — not purely — a boolean can-grow/cannot-grow predictor for conditional purposes.
The 1% cutoff is doing real work for a substantial minority of calls, and the honest verdict is
`CONDITIONAL_SIGNAL_IS_MOSTLY_FEASIBILITY`, not "binary feasibility".

## Honest limits

- The FBA calls are **identical across all settings by construction** — only which genes get scored
  changes. This measures the label-cutoff sensitivity of the *metric*, not of the model.
- The null moves with the cutoff (0.636–0.680), so `lift` is the comparable column; raw per-cell
  agreement is not.
- `per_cell` mode withholds essential calls but does not withhold *dispensable* ones — an unconfident
  near-zero fitness is still read as "not essential". A symmetric abstain-both-ways variant would shrink
  the two-sided set further and was not run.
- All 25 conditions are aerobic carbon sources. No oxygen axis.
- The union deletion pass reuses one FBA result across settings, so a solver quirk would propagate to
  every row identically rather than averaging out. The solver audit rides in the artifact.

## Reproduce

```bash
uv run python scripts/fba_label_threshold_sweep.py
```

Needs `feba.db` (7.4 GB, not committed — figshare `10.6084/m9.figshare.25236931`, CC BY 4.0).
Sidecar: `wiki/fba_label_threshold_sweep_2026-08-13.json`.
Tests: `tests/test_fba_label_threshold_sweep.py` (9) + `tests/test_fba_fitness_browser_t.py` (7).
