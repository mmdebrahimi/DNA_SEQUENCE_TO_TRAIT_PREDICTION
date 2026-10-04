# Klebsiella × meropenem: the first prospective Klebsiella evidence, and the rule calls every isolate R (2026-10-03)

**Second cell ever scored against the live 2026-08-31 v2 surface** (`verify_lock` OK, 5 files pinned),
and the first prospective evidence of any kind for Klebsiella.

**N=16 eligible / 0 excluded (7R / 9S): tp 7, fp 9, tn 0, fn 0 — acc 0.438, sens 1.000, spec 0.000.
`status = UNDERPOWERED`.**

**Read the specificity, not the sensitivity.** sens 1.000 here is uninformative — the rule called
*every* isolate R, so perfect sensitivity is arithmetic, not skill.

## Where it came from

The 2026-09-28 sweep wrote 442 eligible rows and only 426 were scored (Campylobacter × cipro). These
16 were the unscored remainder, located by reading that sweep's own output off D: rather than its
commit message. A **fresh sweep on 2026-10-03 reproduced the eligible set exactly** — 442 rows,
identical key set, 0 added, 0 dropped — so the rows are real and still eligible, and the 5-day re-run
was a genuine **accruing zero** (PD grew +1,958 rows; nothing newly eligible).

## Diagnosed before being reported

The deployed rule is `threshold 1, subclass_any={CARBAPENEM}` — one AMRFinder `CARBAPENEM`-subclass
determinant is enough to call R. Recomputed independently from the cached AMRFinder runs through the
deployed `call_resistance`: **tp=7 fp=9 tn=0 fn=0, identical to the artifact.**

**All 9 false positives carry exactly one carbapenemase:**

| determinant | FP isolates |
|---|---|
| `blaOXA-48` | 5 |
| `blaNDM-1` | 2 |
| `blaOXA-181` | 1 |
| `blaKPC-3` | 1 |

So this is **not** a detection failure, and **not** the broad-drug-class over-call this repo has on
record for intrinsic genes — these are genuine acquired carbapenemases that AMRFinder itself files
under `CARBAPENEM`. The rule did exactly what it says it does.

The disagreement is between **determinant presence** and **measured MIC**. **6 of the 9 are
OXA-48-family** (`blaOXA-48` ×5 + `blaOXA-181`), which is where the FPs concentrate.

> **NAMED HYPOTHESIS, not established here:** that OXA-48-family producers commonly test
> meropenem-susceptible at the clinical breakpoint absent porin loss. That is a clinical-microbiology
> claim this run did **not** test — it only measured the concentration (6/9). Do not promote it to a
> mechanism without evidence.

**What is NOT claimed:** that the labels are wrong. Every FP is determinant-positive *and*
MIC-susceptible; both readings can be simultaneously true, and calling the label an artefact would be
the lazy move. Nor is any change to the frozen surface proposed — that is an authority fork, and n=9
single-batch is nowhere near enough to motivate one.

## Against the rule's own frozen-era number

The rule's `validated` string reads: *"Klebsiella N=30 acc 0.867/sens 1.0/**spec 0.733** (acquired
carbapenemase, CARBAPENEM-subclass: blaKPC/NDM/OXA-48 …)"*.

So specificity went **0.733 → 0.000**. Under independence, P(0 of 9 susceptible isolates called
correctly | spec = 0.733) = **6.9e-06**.

**That p-value is an upper bound on surprise, not a result.** These 16 sit on consecutive `SAMN5593xxxx`
slots with a single release date — one submission batch — so they are epidemiologically correlated and
the effective N is well below 9. A single OXA-48-carrying cluster could account for the whole
discordance as one or two independent observations.

## A gap in the rule's self-description

The `validated` string declares an **FN** mode — *"blind to porin-loss-mediated R (expected FN
mode)"* — and declares **no FP mode at all**. The failure measured here is an FP mode: carbapenemase
present, phenotype susceptible. It is undeclared, which is exactly the shape the L2 doubt layer exists
to carry (qualify the call without changing it).

## Honest limits, in order of severity

1. **UNDERPOWERED by this repo's own bar** — 7R/9S against a ≥10-per-class floor. The scorer said so
   itself rather than being told.
2. **SINGLE-SOURCE** — one submission batch, so largest-source share is effectively 1.000 and it fails
   the project's own ≥5-BioProject / ≤0.60 diversity bar by the maximum margin, exactly as the
   Campylobacter cell does.
3. **Essentially one mechanism class** — carbapenemase presence. This says nothing about the rule's
   behaviour on determinant-negative isolates; there are none here (tn=0, fn=0).
4. **Not clonality-corrected.** Read 0.000 as isolate-level.
5. The lock guarantees only that these isolates became public after 2026-08-31 — leakage-free **by
   construction**, which is the one strong property this tier has and it is unaffected by the above.

## What this does and does not support

It **does** support: the frozen meropenem rule, applied prospectively to Klebsiella, called 9 of 9
susceptible isolates resistant, and every one of those 9 carried a real carbapenemase.

It does **not** support a specificity estimate for the cell. One underpowered, single-source batch
cannot replace the frozen-era N=30. The honest next step is a **source-diverse** Klebsiella meropenem
cohort (≥5 BioProjects, ≤0.60 share) — which is the same instrument
`scripts/source_diverse_validate.py` already applies elsewhere, and which refused 6 of 10 cells when
it was last run.
