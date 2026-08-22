# PRE-REGISTRATION — expression-gated GPR

**Written:** 2026-08-22, **before any run**, and while `D:` (which holds PRECISE-1K) is **disconnected**.
That is deliberate: the data is physically unreachable as this is written, so no result can have informed
any choice below. **Frozen on commit.** Any later change must appear as a dated AMENDMENT in this file,
never as a silent edit.

This project has already **retracted one positive result** (`wiki/fba_eflux_graded_RETRACTION_2026-08-17.md`)
for exactly the failure this document exists to prevent: a new endpoint chosen after the primary one
disappointed, on a run whose variance exceeded its effect.

## 1 · Hypothesis

Four constraint-based levers (gap-fill, threshold retune, pFBA, E-Flux) failed identically on conditional
gene essentiality. Measured cause, in part: **31 of 131 experimentally-essential genes (23.7 %) are
provably uncallable** by any of them (`wiki/fba_replaceability_lp_2026-08-22.md`). Of those, **8 are
isozyme-masked** — the function IS essential and only GPR redundancy hides it.

> **H1.** Deleting an isozyme is a **boolean** operation on the GPR. E-Flux scales a reaction's *bounds*,
> but an `or` of two isozymes keeps the reaction functional however the bounds move — so E-Flux could
> never have recovered a masked gene. **Gating the boolean** (treating an isozyme as *absent* when its
> expression is below threshold, not merely capacity-limited) is a structurally different intervention
> and can recover them.

This is why H1 is not "the same lever with a new knob": it changes the *operator*, not the parameter.

## 2 · Pre-declared target set (frozen)

The 8 genes in `wiki/fba_orphan_protection_2026-08-21.json` → `impact_on_experimental_deficit`
→ `genes_isozyme_masked`:

`ilvI` `ilvH` `ilvB` `ilvN` `aroF` `tktA` `trxA` `ompC`

Naming the target **in advance** is the entire point. A post-hoc count of "genes that improved" is not a
result.

Mechanistic note recorded now, so it cannot be back-fitted: `ilvG`/`ilvM` are **absent from iML1515**, and
`ACHBS`/`ACLS` carry `(ilvN and ilvB) or (ilvI and ilvH)` — two genuinely functional isozymes. So H1
predicts these specific recoveries occur **only in conditions where the partner isozyme is not expressed**.
If a recovery appears in a condition where BOTH pairs are well expressed, that is evidence **against** the
stated mechanism even if the number improves.

## 3 · Intervention (frozen)

For each condition with PRECISE-1K expression, and each gene *g*: mark *g* **absent** when its expression
is below the gate, then evaluate every GPR with the absent set removed, then run the deletion as usual.

- **Primary gate: expression below the 20th percentile of that condition's own gene distribution.**
  Per-condition percentile, not an absolute TPM — absolute cutoffs are not comparable across samples.
- **Pre-declared sensitivity range: {10, 20, 30}th percentile**, reported in full. The primary is the
  20th. Reporting the best of the three as the headline is **forbidden**.
- Expression normalisation, join key and condition set are inherited **unchanged** from
  `scripts/fba_eflux_bridge.py` (including its join-key fix). No re-derivation.

## 4 · Endpoints (frozen, ordered)

| | endpoint | success |
|---|---|---|
| **Primary** | of the 8 pre-declared genes, how many become predicted-essential in ≥1 condition where they are experimentally essential | **≥ 4 of 8** |
| Secondary | net change in overall recall over the 131 | reported, not a bar |
| **Guardrail** | false positives across the full gene × condition grid | **must not rise by more than 20 % relative to baseline** |

**Binary essentiality (ratio ≤ 0.01) is the primary readout, unchanged.** No fitness-scale transform is in
scope here; that is a separate pre-registration.

A primary hit with the guardrail breached is a **FAILURE**, not a partial success. Gating expression off
trivially makes more things look essential; the guardrail is what makes the primary meaningful.

## 5 · Determinism gate (frozen — this is the retraction's lesson)

Before any comparison is interpreted:

- all deletions run with `processes=1`;
- the **entire pipeline runs twice**, and the two runs must produce **identical essentiality calls**;
- the primary metric's difference between arms must exceed the between-run spread by **≥ 1000×**
  (`dna_decode/fba/nitrogen.py::determinism_verdict`).

**A bar cleared by a run whose variance exceeds its effect is not a result.** If the runs disagree, the
outcome is `INDETERMINATE` — never the better of the two.

## 6 · Pre-committed verdicts

| outcome | verdict | what happens next |
|---|---|---|
| ≥4/8 recovered, guardrail held, determinism passed | **H1 SUPPORTED** | report; the mechanism-class is real |
| 1–3/8 recovered, guardrail held | **H1 WEAKLY SUPPORTED** | report as underpowered; do **not** tune the gate to reach 4 |
| 0/8 recovered | **H1 FALSIFIED** | report the negative; the isozyme-masked class is not expression-explained |
| guardrail breached at the primary gate | **FAILURE** regardless of the primary | report as failure |
| runs disagree | **INDETERMINATE** | report; fix determinism before re-running |

## 7 · Known confounds, recorded in advance

1. **Strain mismatch.** PRECISE-1K is K-12 MG1655; the Fitness Browser labels are `Keio` (BW25113
   parent). Inherited from the E-Flux bridge and unchanged.
2. **Condition mismatch.** PRECISE-1K conditions are matched to the carbon panel by name; the overlap was
   ~11 of 25 in the bridge. A gene can only be recovered in a condition that has expression.
3. **A 20th-percentile gate marks ~20 % of genes absent by construction** — hence the guardrail.
4. **Absence ≠ non-function.** Low mRNA does not prove no protein; the gate is a proxy and the claim must
   stay "expression-gated GPR recovers *n* of 8", never "the isozyme is absent in the cell".

## 8 · Blocked on

`D:/dna_decode_cache/precise1k` (drive disconnected). Everything except the run is complete. **This is an
external wall, not a code wall** — no further building closes it.
