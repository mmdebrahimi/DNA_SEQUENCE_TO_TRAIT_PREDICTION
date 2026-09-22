# The genotype's resistance boundary sits one dilution below the clinical cut-off — SUPPORTED, powered

**Verdict: `SUPPORTED`** on both drugs, against a bar frozen after the mechanism analysis and before any
test statistic existed (`wiki/tb_implied_boundary_acceptance_bar.json`). Run:
`uv run python scripts/tb_implied_boundary.py` (offline, read-only, no network, no Docker, seconds).

| drug | recovered clinical cut-off | top susceptible rung | carriage there | carriage below | gap | permutation max |
|---|---|---|---|---|---|---|
| rifampicin | 0.5 mg/L | 0.5 (n=188) | **0.378** | 0.029 | 0.3492 | 0.0435 |
| isoniazid | 0.1 mg/L | 0.1 (n=275) | **0.160** | 0.011 | 0.1491 | 0.0273 |

Both gaps exceed the **maximum** of 1000 permutations — RIF by 8×, INH by 5.5× — not merely the 95th
percentile.

> ## Powering caveat added 2026-09-22 — read this before quoting `SUPPORTED`
>
> The frozen bar is **silent on `PHENOTYPE_QUALITY`**, but this repo's two other CRyPTIC scorers
> (`score_tb_cryptic.py:189`, `score_tb_cryptic_parquet.py:150`) both hard-filter to `HIGH`. This run
> did not. So it is powered partly by label tiers the repo elsewhere discards:
>
> | drug | boundary rung n (as run) | HIGH-quality only | frozen floor |
> |---|---|---|---|
> | rifampicin | 188 | **84** | 100 |
> | isoniazid | 275 | 145 | 100 |
>
> The bar requires **both** drugs, so **under the repo's own TB quality standard this run reads
> `INDETERMINATE`, not `SUPPORTED`.** The two compared strata also differ: the boundary rung is
> **2.5× (RIF) / 4.4× (INH)** enriched for LOW-quality labels relative to the lower-S stratum.
>
> **The effect itself is not a label-quality artifact — it replicates within HIGH alone:** RIF gap
> **0.3363** (vs 0.3492 pooled), INH **0.1446** (vs 0.1491). So the finding holds and the *verdict
> label* does not. **A powering failure is not a refutation**, and the two must not be collapsed.
>
> The verdict field is deliberately **left as the mechanical output of the frozen bar**; the HIGH-only
> reading ships beside it under `quality_gate`. Re-sizing a frozen bar after seeing the result is an
> authority call, and this session has already had to override two mis-specified bars — a third
> silent edit is exactly what that record argues against. Found by the `/innovate` sweep
> (`wiki/innovate_sweep_2026-09-21.md`, survivor `QS-TBQUAL`).

## What this is, and why it moved here

The 2026-09-12 ECOFF arm asked whether a wild-type anchor agrees with the genotype better than a
clinical breakpoint. On the Oxford cohort it returned `WEAK_DIRECTIONAL`: an 18× relative enrichment in
the right direction that sat **inside** its own null, because the discriminating stratum was 25 isolates
carrying 2 determinants. That memo named its own fix — *"a cohort with a full dilution series at the low
end … not a bigger cohort"*.

**CRyPTIC is that cohort and it was already on disk.** 12,287 *M. tuberculosis* isolates with measured
broth-microdilution MIC on a genuine two-fold ladder (RIF 0.03→8 over 12 levels, INH 0.025→12.8 over 12)
against Oxford's collapsed gentamicin panel of `{1, 2, 4, 32}`.

**The inventory rule paid for itself before any code was written.** CRyPTIC-as-MIC-substrate is *not*
unexplored: `scripts/tb_mic_calibration.py` shipped a calibrated-MIC-interval arm on the same two drugs
on 2026-07-12 (`CALIBRATED_MIC_INTERVALS`), and it already records the exact validity argument I was
about to reconstruct — that this is wet-lab BMD MIC, not the closed BV-BRC negative (G1: 91% of BV-BRC
MIC is XGBoost-from-genome). It also left behind the per-isolate determinant cache this run reuses, so
no 2.9 GB parquet stream was needed. Its question was different (interval calibration, not anchor
choice), so this is new work — but only checking established that.

## The framing change that matters more than the sample size

Oxford's version needed **two sourced numbers** — an ECOFF and a clinical breakpoint — which makes the
result hostage to both and carries a fabrication hazard on a biological reference value.

Here neither is needed. CRyPTIC ships its own `BINARY_PHENOTYPE`, and that column turns out to be a
**pure step on the MIC ladder** — measured, not assumed: R-fraction is exactly **0.000** at every rung at
or below one value and exactly **1.000** above it. So the effective clinical cut-off is *recovered* from
the shipped data and is re-derivable by anyone. `recover_cutoff` **refuses** when the step is impure,
because an impure step would mean the phenotype is not a MIC threshold and the whole recovery argument
is void.

**No breakpoint or ECOFF is recalled from memory anywhere in this arm.** The measured quantity is *where*
the genotype's boundary falls on the ladder, not whether one asserted cut-off beats another.

## The mechanism is visible in the composition, and needs two statistics to see

The determinants sitting at the boundary rung are not a thinner version of the high-MIC ones — but the
shift shows up at a **different resolution in each drug**, and either statistic alone would have missed
one of them:

| drug | allele-level (top-1 differs) | gene-level (share shift) |
|---|---|---|
| RIF | **yes** — boundary is `Asp435Tyr` / `Leu452Pro` / `Ile491Phe`; high-MIC is dominated by `Ser450Leu` (2,898) | **vacuous** — every RIF determinant is in `rpoB`, so this statistic *cannot* move |
| INH | no — `katG_p.Ser315Thr` is commonest at both ends | **yes** — `inhA` share 0.479 at the boundary vs 0.185 high-MIC (shift 0.294) |

I flagged the top-1 comparison as deliberately weak in the code rather than quietly reporting only the
statistic that happened to fire: it reports "no shift" for isoniazid while the per-gene mix moves
substantially. Both ship.

**The literature reading of *which* mechanisms these are is deliberately not asserted.** The gene is taken
lexically from the determinant string's prefix — a pure string operation on data already in the cache —
so the claim is only that the mix differs, which is checkable from the counts.

## Three defects found in this script, all mine

**A third, found 2026-09-22: an all-susceptible ladder returned a vacuous cut-off instead of refusing.**
`recover_cutoff([-2,-1,0], ["S","S","S"]) -> 0.0`. With no resistance anywhere there is no threshold,
but the loop set `out` to the top rung, the `above` slice was then empty so its check was *skipped*, and
an information-free column yielded a confident-looking number. Harmless here (both classes present in
the real data) — but it would have poisoned any generic scanner calling this directly, which is exactly
what a sibling `/innovate` candidate proposed building before it was killed on other grounds. Now
refuses when either class is absent; proven non-vacuous by re-injection.

## Two defects found while building this, both mine

**A real bug in the cut-off recovery.** It verified that every rung *above* the candidate was resistant
but never checked that everything *below* was susceptible — so a non-monotone phenotype (an `R` rung
beneath an `S` one) would pass and return a bogus cut-off. The headline is **unchanged** by the fix,
because the real data is a clean step in both directions; what changed is that a future cohort which
isn't can no longer slip through. Proven non-vacuous by re-injecting the defect.

**An exact-tie trap in my own statistic.** `largest_shift_gene` used `max()`, but with exactly two genes
the shares sum to 1, so both shift by *identically* 0.2939 and the "winner" was whichever sorted first —
an arbitrary choice dressed as a finding. It now reports the magnitude plus **every** gene attaining it,
with a `largest_shift_is_tied` flag. Same class as the FBA lesson: every exact tie hides an arbitrary
choice.

## Honest limits

- **POWERING DEPENDS ON LABEL TIERS THIS REPO ELSEWHERE DISCARDS** — see the caveat box above. Under
  HIGH-only labels the RIF boundary rung is 84 against a frozen floor of 100, so the run reads
  `INDETERMINATE`. The effect replicates within HIGH; the powering does not.
- **IN-DISTRIBUTION.** The WHO catalogue was built partly from CRyPTIC, so the determinant calls and this
  cohort are **not independent**. This locates a boundary on the ladder; it is **not** an independent
  validation of the catalogue.
- **One organism, two drugs, one compendium** — and both drugs share the same isolates, so they are
  **not two independent replications**.
- **A determinant at a clinically-susceptible rung reads two ways and this arm cannot separate them:** a
  genuine low-level-resistance mechanism sitting under the cut-off, or a catalogue false positive. The
  rung-position dependence (a 6–18× jump at exactly the top rung rather than a flat rate across
  susceptible rungs) argues against *random* false positives, but does not exclude the second reading.
- **This does not name a numeric ECOFF** and is not a claim that any ECOFF value is correct.
- Frozen surfaces are read-only throughout; TB rules live in the non-frozen `organism_rules` package.
- The boundary rungs are **uncensored exact values** — asserted by test, not assumed. Censoring sits at
  the ladder ends (`<=0.03`, `>4`), so no boundary isolate is a censored bound.

## Reusable

**When a comparison needs an external constant, check whether the data already determines it.** The
Oxford arm spent a whole run blocked on sourcing two reference values, and its successor needed neither —
because the cohort shipped a phenotype column that *is* a threshold, making the cut-off a recoverable
property rather than an input. The refusal path is what makes that legitimate: recovery is only sound
while the step is pure, so the code must be able to say "this cannot be recovered here".

**A composition statistic has a resolution, and the wrong one reports no change.** Per-gene shares are
vacuous for a single-gene drug; top-1 allele comparison is blind to a mix shift under a constant leader.
Neither drug here would have been fully characterised by one of them.
