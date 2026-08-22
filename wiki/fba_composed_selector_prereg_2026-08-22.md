# PRE-REGISTRATION — composed expression selector (wildtype-safe by construction)

**Written:** 2026-08-22, **before the scoring run**. **Frozen on commit.** Amendments must be dated
sections appended here, never silent edits.

**What has been looked at, stated plainly:** the preflight in `wiki/fba_selector_preflight_2026-08-22.md`
measured **TRN availability and wildtype behaviour only**. No recovery count, no false-positive count, no
endpoint of any kind has been computed. That boundary is the same one used for the two prior
pre-registrations in this arc.

## 1 · Hypothesis

Two expression-gating attempts have now failed, and the preflight shows **both failed for the same
under-diagnosed reason**: the gate removed capability the model needs, collapsing the wildtype, which then
manufactured recoveries via the `wt <= 0 -> all genes essential` path.

> **H1.** A selector that is **wildtype-safe by construction** — per-gene thresholding, restricted to
> knockouts that disable nothing, and verified **jointly** so no isozyme pair is fully gated — recovers
> isozyme-masked genes without the collapse artifact.

This is not "a better threshold." It adds the joint-composition step that both prior attempts lacked.

## 2 · Intervention (frozen)

Per condition *c*:

1. **Per-gene criterion.** Gene *g* is a candidate if its expression in *c* is below the **10th percentile
   of its OWN distribution across the full 1,035-sample PRECISE-1K compendium**. (Not a genome-wide
   percentile — that is the falsified selector.)
2. **Single-gene eligibility.** Keep *g* only if `gpr_disabled_reactions(g)` is empty.
3. **JOINT verification (the new step).** Compute the reactions disabled by the whole candidate SET. While
   any reaction is disabled, remove from the set the member with the **highest** expression in *c*
   (deterministic tie-break: lexicographic gene id) and recompute. Terminate when the set disables nothing.
4. Apply the resulting set as knockouts, then run the deletion scan as usual.

Step 3 makes wildtype invariance a **constructed** property, so it must be asserted at run time, not hoped
for (see §5).

**Sensitivity range, declared in advance:** per-gene percentile ∈ {5, 10, 20}. **Primary is 10.**
Reporting the best of the three as the headline is **forbidden**.

## 3 · Target set (frozen) — unchanged from the earlier run

The 8 isozyme-masked genes from `wiki/fba_orphan_protection_2026-08-21.json`:
`ilvI` `ilvH` `ilvB` `ilvN` `aroF` `tktA` `trxA` `ompC`

Mechanistic prediction recorded now: a target can only be recovered in a condition where its **partner
isozyme** is in the gated set. A recovery in a condition where the partner was *not* gated is evidence
**against** the stated mechanism, even if the count improves.

## 4 · Endpoints (frozen, ordered)

| | endpoint | success |
|---|---|---|
| **Gate-0 (hard)** | wildtype growth unchanged (≤1e-6 abs) in **every** scored condition | any deviation ⇒ `INVALID_WILDTYPE_PERTURBED` |
| **Primary** | of the 8, how many become predicted-essential in ≥1 condition where they are experimentally essential | **≥ 4 of 8** |
| **Mechanism check** | fraction of recoveries where the partner isozyme was in the gated set | **1.0 expected**; <1.0 ⇒ mechanism disconfirmed |
| Secondary | net recall change over the 131 | reported, not a bar |
| **Guardrail** | false positives over the full gene × condition grid | **must not rise more than 20 % relative** |

**A primary hit with the guardrail breached is a FAILURE** — the rule that decided both prior runs.
**Gate-0 failing invalidates the run entirely**, regardless of every other number; that is the explicit
fix for the artifact that invalidated the first attempt.

## 5 · Determinism + validity gates (frozen)

- `processes=1`; the whole pipeline runs **twice**; identical essentiality calls required, else
  `INDETERMINATE`.
- Gate-0 is asserted **per condition, at run time**. A constructed invariant that is not checked is an
  assumption.
- Conditions with no matched expression are excluded from scoring and reported, not silently dropped.

## 6 · Pre-committed verdicts

| outcome | verdict |
|---|---|
| Gate-0 fails anywhere | **INVALID_WILDTYPE_PERTURBED** (no other number is reported as a result) |
| runs disagree | **INDETERMINATE** |
| guardrail breached | **FAILURE** regardless of the primary |
| ≥4/8 recovered, mechanism 1.0, guardrail held | **H1 SUPPORTED** |
| ≥4/8 recovered, mechanism <1.0 | **H1 SUPPORTED, MECHANISM DISCONFIRMED** |
| 1–3/8 recovered, guardrail held | **H1 WEAKLY SUPPORTED** — do not tune the percentile to reach 4 |
| 0/8 recovered | **H1 FALSIFIED** |

## 7 · Known limits, recorded in advance

1. **The TRN is not used by this design.** It was fetched and verified to cover 8/8 targets, but the
   composed selector needs only expression. A TRN-gated variant is a *separate* pre-registration; folding
   it in here would confound two changes.
2. **11 of 25 conditions** have matched expression. A target cannot be recovered where there is no data.
3. **Strain mismatch** — PRECISE-1K is K-12 MG1655; labels are `Keio` (BW25113 parent). Unchanged.
4. **The compendium reference mixes stresses and strains.** A condition-matched reference would be
   tighter; using the full compendium is the conservative, non-tuned choice.
5. **Low mRNA does not prove absent protein.** The claim is about a gating rule, never about the cell.
6. **Step 3's resolution rule (drop the highest-expressed member) is a choice, not a derivation.** It is
   fixed in advance and not tuned; an alternative rule is a separate experiment.
