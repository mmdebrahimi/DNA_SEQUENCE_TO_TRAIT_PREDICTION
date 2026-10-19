# Adoption of `learned-model-is-the-product` — and the evidence deliberately NOT cited

**Date:** 2026-10-19
**Register:** `soraya_runs/standing-positions.json`
**Status set:** `adopted-user-was-right` (the `ADOPTED` constant, `standing_positions.py:73`)
**Resolved by:** `project_state/trait-to-edit/big-idea.md`

## What this is, and what it is not

The position — *the GLM (a learned, generative model) is the product; the deterministic catalog was a
bridge built only so there was a number to compare against* — was raised by the user three times and
stood `open`. Three raises while open is the register's own LOUD threshold: the signal that the project
is talking past the user.

It is now **ADOPTED**, which in this register means *investigated, the USER was right, adopt it*. That
status exists precisely because the other two closing statuses (`resolved-by-evidence`,
`withdrawn-by-user`) both smuggle "the model was right" into the word *resolved*.

**ADOPTION IS NOT EVIDENCE.** No experiment tested this claim. The basis is the user's authority plus
the fact that the mission document now encodes the position as its opening premise. Nothing in this
memo or in the register should be read as a measurement.

## The evidence that looks like it settles this, and why it does not

`wiki/glm_headroom_allrt_2026-09-30.md` reports that a supervised model over position-resolved
substitution tokens reaches a **median 0.8141 leave-one-STUDY-out against 0.4750 zero-shot**, with
**7 of 8 RT drugs passing versus 0 of 8**, on a subset that is catalog-negative by construction.

(The near-identical **0.8142** belongs to a *different* artifact,
`wiki/glm_alphabet_headroom_2026-09-30.md`, and is the single-drug EFV figure rather than the all-RT
median. The two are one digit apart and easy to swap; quote each with its own artifact.)

That is a real result and it supports the learned **direction**. It does **not** resolve this position,
because:

- it is **DISCRIMINATIVE** (scoring which variant is resistant), while the position is about
  **GENERATION** (producing the edit to make);
- those are different quantities, and conflating them is this project's recorded **orientation-error 0**
  — the catalog-vs-GLM category error, corrected across four separate sessions;
- citing it as the resolving measurement is exactly the substitution `proxy_warning()` exists to catch.

So it is named here and kept **out of** the register's `tested_exactly` field on purpose. A reader who
later finds this adoption and wants the "proof" should find this paragraph rather than a number that
measures something else.

## Consistency with the mission's modality correction

Adopting this does not contradict the modality correction made in the same pass. On the promoter-tile
substrate a **346-parameter biophysical model** that was never trained on it scores rho **0.3968**,
above every purely-learned arm; the only arm whose paired CI excludes zero is the **orthogonal hybrid at
0.4254**, and its shuffled-conv control **hurts** (**-0.1232**, 0/10 seeds).

The hybrid **contains a learned channel.** So the measurement says a learned channel pays *when combined
with orthogonal non-learned ones* — not that learning is the wrong bet. It bounds which **lever** to
pull next; it does not demote the generative product.

## What stays open

`substrate-is-not-the-subject` remains **`open` at `times_raised: 2`**. The mission's north-star step-4
rewording (step 4 is CARRIED by `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN`, not
"deliberately absent") *addresses* it; it does not resolve it, and this pass makes no claim that it does.

## Honest limits

- The adoption rests on **authority, not measurement**. If a later session wants evidence for the
  generative claim, none exists yet and this memo is not it.
- The register's `times_raised: 3` carries `increments_without_id` history and a recorded count
  correction; the count is the LOUD signal, so it is worth re-reading the register rather than trusting
  the number quoted here.
- Nothing in the frozen AMR surface, the prospective lock, or any report card is touched by this.
