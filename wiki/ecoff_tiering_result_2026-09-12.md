# ECOFF-anchored tiering: the direction is right, the evidence does not clear its own control

**Verdict: `WEAK_DIRECTIONAL`** — one of four outcomes named in
`wiki/ecoff_tiering_acceptance_bar.json`, frozen before any carriage number existed. Per that bar this
**must be reported as NOT established, never as a partial win.**

Run: `uv run python scripts/ecoff_tiering_feasibility.py` then
`uv run python scripts/ecoff_tiering_evaluate.py` (offline, read-only, seconds).

## The question

A determinant-based decoder predicts **mechanism presence** — wild-type vs non-wild-type — while it is
scored against a **clinical breakpoint**, which asks whether treatment will succeed and depends on dose
and site as well as the organism. The EUCAST WGS subcommittee's position, surfaced by the 2026-09-03
prior-art scan, is that a genotype-based decoder should be anchored to the ECOFF instead.

That makes a testable prediction. The two anchors disagree only about isolates that are **clinically
susceptible but above the ECOFF**. If the ECOFF is the better anchor for this decoder, those isolates
should carry the decoder's own determinants more often than the clinically-susceptible isolates below
the ECOFF.

## Most of the question could not be asked at all

The feasibility gate (Step 2) ran first, with its outcomes named in advance. **Three of four drugs were
unscoreable before any comparison was possible:**

| drug | ECOFF | Oxford panel floor | verdict |
|---|---|---|---|
| gentamicin | 2 | 1 | **SCOREABLE** |
| ciprofloxacin | 0.06 | 0.125 | DEGENERATE_ECOFF_BELOW_PANEL |
| ceftriaxone | (0.125) tentative | 0.5 | DEGENERATE_ECOFF_BELOW_PANEL |
| tetracycline | 8 | — | NOT_IN_COHORT (no MIC column) |

For cipro and ceftriaxone the ECOFF sits **below the lowest dilution the panel reports**, so 100% of
isolates are non-wild-type and the anchor separates nothing. That is a property of this panel, not of
the ECOFF: a cohort reporting finer low-end dilutions would make both scoreable. The plan named
`DEGENERATE_ECOFF_BELOW_PANEL` before any ECOFF was known, which is what keeps it a result rather than
a discovery made halfway through.

## The one drug that could be asked

| stratum | n | determinant carriage |
|---|---|---|
| clinically-S **and** ECOFF-NWT | 25 | **0.080** (2 isolates) |
| clinically-S **and** ECOFF-WT | 2,652 | **0.0045** (12 isolates) |
| clinically-R *(control)* | 192 | **0.922** |

The direction is right, and not marginally so — **an 18-fold relative enrichment**. It still fails the
bar: the observed gap is **0.0755 against a permutation maximum of 0.1159**, so it sits inside its own
null.

**That is not a quibble, it is the correct call.** With 25 isolates against a base rate of 0.45%, a
single extra carrier moves the rate by 0.04; over 1000 shuffles chance alone lands 3 carriers in the
small stratum routinely. Two carriers cannot distinguish an 18-fold effect from noise, and a bar that
passed this would pass noise.

**The control stratum is what makes the null result interpretable.** Clinically-resistant isolates carry
a determinant 92.2% of the time, so the deployed rule plainly works on this cohort. Had that been low,
the flat susceptible strata would have said nothing about the anchor and everything about the rule.

## Why the bar was a permutation null and not an effect size

Two acceptance bars earlier in this project were mis-specified because a threshold was frozen before the
analysis that says what a correct result can look like. Here the substrate, the determinant definition
and the strata were all measured **first**; only the verdict rule was frozen, and it was frozen on a
**control that self-calibrates to the stratum sizes actually present** rather than on a number I would
have had no basis to choose. At n=25 there is no effect size I could have justified in advance.

## The determinant definition decided this, and it was nearly wrong

The rule counts a genome if the **deployed** gentamicin rule counts it: AMRFinder `Subclass` containing
`GENTAMICIN` as a **substring** (so `GENTAMICIN/KANAMYCIN/TOBRAMYCIN` counts), threshold 1.

Keying instead on `Class=AMINOGLYCOSIDE` would have pulled in **3,298 STREPTOMYCIN rows**
(`aph(3'')-Ib`, `aph(6)-Id`, `aadA`) against 343 gentamicin ones — producing a result about streptomycin
co-carriage that still read as a gentamicin finding. The cohort's aminoglycoside genes are mostly not
gentamicin genes.

Also measured: **zero `rmt`/`npmA` rows exist in this cohort**, independently reconfirming the
2026-09-03 Oxford rmt probe. The v2 rescue is inert here and cannot have moved the result either way.

## Honest limits

- **One organism, one UK region, 2008–2018, one collapsed MIC panel.** Oxford is single-source — it
  trips G7 and is `n/a` on G2 for exactly that reason (`wiki/oxford_gate_screen_2026-09-11.md`).
- **n=25 and 2 carriers.** A `WEAK_DIRECTIONAL` verdict at this size is **weak evidence against
  nothing** — it neither supports nor refutes the claim. The honest summary is that this cohort cannot
  answer the question, not that the answer is no.
- **This compares carriage, not accuracy.** It asks which anchor the genotype agrees with better; it
  does not establish that either anchor predicts clinical outcome.
- **Ceftriaxone's ECOFF is a TECOFF** (parenthesised, 908 observations from 4 sources). Its degeneracy
  here does not depend on that — 0.125 is below the 0.5 floor either way — but any ceftriaxone
  conclusion would.
- **No frozen file was touched.** `mic_tiers.py` is sha256-pinned in the v2 prospective lock; adopting
  ECOFF anchoring would retire that lock and restart the prospective clock, and remains a
  user-authority decision. A `SUPPORTED` verdict would not have changed that, and this is not one.

## What would actually answer it

A cohort with a **full dilution series at the low end**. Oxford's panel reports four gentamicin values
(1, 2, 4, 32) and nothing between 4 and 32; cipro jumps 0.5 → 8. The degeneracy and the thin stratum
are both consequences of that grid, not of the biology or the anchors. Measured MICs on a standard
two-fold series would make three of the four drugs scoreable and put hundreds of isolates in the
discriminating stratum instead of 25.

## Reusable

**A self-calibrating control beats a guessed threshold when the stratum size is not known in advance.**
Both bars this project mis-specified were fixed numbers chosen before the mechanism was understood. A
permutation null on the label being tested adapts to whatever n turns out to be, and it refuses
underpowered results without needing anyone to have predicted the power.

**Check that your predictor works on the stratum where it should, before reading a null where it
shouldn't.** The 92.2% carriage among clinically-resistant isolates is what licenses reading the
susceptible strata as a statement about the anchor rather than about the rule.
