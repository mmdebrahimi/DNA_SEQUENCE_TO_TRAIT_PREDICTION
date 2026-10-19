# CABBAGE arm — four scored cells, 27,879 isolates, and what they jointly establish

**2026-10-10 → 10-11.** Read this before re-deriving anything about the CABBAGE arm. Per-cell detail:
`cabbage_scored_cell_2026-10-10.md` (Campylobacter cipro), `cabbage_ecoli_gentamicin_2026-10-11.md`,
`cabbage_salmonella_cipro_2026-10-11.md`, `cabbage_salmonella_ceftriaxone_2026-10-11.md`; mechanism
probe `cabbage_qrdr_multihit_2026-10-11.json`.

## What the arm is

The EMBL-EBI AMR Portal ships **pre-computed AMRFinderPlus** per assembly (~1 KB each) alongside
measured MIC. That removes ~95 s/strain of local Docker, which is the only reason cells of this size
are affordable here — the Campylobacter cohort alone would have been ~133 hours of local AMRFinder.
Every cell runs the **frozen** `eval/amr_rules.call_resistance` unchanged, labels go through this
repo's own `external_mic_labels.tier_for_isolate` on STRICT operator-aware MIC tiers, and a
fail-closed `cohort_manifest` leakage gate excludes prior accessions.

## The four cells

| cell | N scored | acc | sens | spec | FP | source-diversity verdict |
|---|---|---|---|---|---|---|
| *Campylobacter* × ciprofloxacin | 4,265 | 0.9925 | 0.9794 | 0.9966 | 11 | `source_concentrated` (96.4% USA) |
| *E. coli* × gentamicin | 7,948 | 0.9943 | 0.9016 | 0.9951 | 39 | **`source_diverse`** (21 countries, largest 0.5844) |
| *Salmonella* × ciprofloxacin | 7,834 | 0.9543 | 0.8209 | **0.9555** | **346** | `source_concentrated` (96.0%) |
| *Salmonella* × ceftriaxone | 7,832 | 0.9935 | 0.9426 | 0.9979 | 15 | `source_concentrated` (97.6%) |

**Only one of four clears this repo's own diversity bar**, and it clears it by 0.016. Quote the
*E. coli* gentamicin cell as source-diverse; quote none of the others as such.

## Three things established

**1. The v2 gentamicin `rmt` rescue is vindicated, and the caveat it shipped with is discharged.** The
v2 lock recorded: *"zero S-labelled `rmt` carriers exist in any dataset checked, so the unchanged
specificity is an ABSENCE of counter-examples, not a bound on over-calling."* Counter-examples now
exist — 2 S-labelled `rmtB1` carriers — and the measured PPV is **0.909 [0.722, 0.975]**, the highest
of any multi-carrier determinant in that cell. Sensitivity replicates: **0.892** on the N=131 cohort
that justified v2, **0.9016** here on ~60× more isolates, source-diverse. **And the over-call is not
the v2 widening**: 37 of 39 false positives come from v1-counted determinants (`aac(3)-IId` PPV 0.500
on 40 carriers, `aac(3)-VIa` 0.000 on 8), only 2 from the rescue.

**2. `threshold=1` is a per-organism-drug-pair sufficiency claim, not a property of the rule or of a
mechanism class.** Three cells make this measurable. R-rate at **exactly one** counted determinant:

| cell | mechanism | R-rate @ 1 determinant | n |
|---|---|---|---|
| *Campylobacter* × cipro | target-site QRDR | 0.9891 | 1,006 |
| *Salmonella* × ceftriaxone | acquired β-lactamase | 0.9756 | 820 |
| *Salmonella* × cipro | target-site QRDR | **0.0183** | 493 |

The *Salmonella* cipro/ceftriaxone pair is a **within-organism, within-cohort control**: 53× apart on
shared labels, a shared portal AMRFinder run and a shared cohort, which exonerates the data and puts
the fault in the rule. The *Campylobacter* row then **refutes the obvious explanation** — it is also
target-site and single-hit predicts at 0.9891, so mechanism *class* does not decide it. What decides it
is whether one determinant clears that pair's clinical breakpoint.

**3. Every false negative in all four cells carries ZERO counted determinants.** 21 + 6 + 12 + 36 = 75
false negatives, none with a determinant the rule counts. This is the rule's own declared blind spot,
stated at call time as `undetectable_mechanisms: [efflux, porin_loss, regulatory]`, not degradation.

## Three located over-calls, all parked on user authority

Each would require editing the sha256-pinned frozen surface, which retires the active 2026-08-31
prospective lock. **None was made.** The frozen five are byte-unchanged and `verify_lock` returned
ok=True after every commit.

| determinant | cell | carriers | true R | PPV | note |
|---|---|---|---|---|---|
| `gyrA_P104S` | Campylobacter cipro | 6 | 0 | 0.000 [0.0, 0.39] | AMRFinder itself files it `Class=QUINOLONE`, so the claim is **upstream** |
| `aac(3)-IId` | E. coli gentamicin | 40 | 20 | 0.500 [0.352, 0.648] | v1-counted; largest single FP source in that cell |
| `aac(3)-VIa` | E. coli gentamicin | 8 | 0 | 0.000 [0.0, 0.324] | v1-counted |
| `gyrA_D87Y` | Salmonella cipro | 169 | 2 | 0.012 [0.003, 0.042] | ~half that cell's 346 FP on its own |

An **L2 doubt flag** is the non-lock-breaking alternative for all of them.

## Two of my own calls were wrong, and both are recorded beside the correction

- The anatomy probe's **first** run read a key the record does not have, so every bucket came back
  determinant-free and it printed a confident `ALL_FN_CARRY_NO_DETERMINANT_DECLARED_BLIND_SPOT` — a
  plausible verdict matching the prospective cell's story, from an extraction that had parsed nothing.
  The tell was 997 true positives with no determinant, impossible by definition. It now **refuses
  (exit 3)** when the TP bucket is determinant-free.
- A pre-registered **strict monotonicity** test for the determinant-count effect **fails** on Salmonella
  cipro (4-determinant stratum 0.857 below the 3-determinant 1.000, overlapping CIs). The premise holds
  decisively as a single-vs-multi contrast (0.0183 vs 0.8246, disjoint CIs), but the failed test is kept
  next to it rather than swapped for the version that passes.

A third, smaller one: an accession-ordered prefix is **measurably easier** (Campylobacter first-1,444
sens 0.9971 vs full-4,265 0.9794), so `--max-isolates` now takes a **seeded random** sample.

## Standing limits for the whole arm

- **NOT clonality-corrected** — every number is isolate-level. Constraint
  `cabbage-clonality-g8-blocked-no-mash`: Docker down, no native `mash`, and the portal parquet carries
  **no MLST/ST column**, so the Mash-free lineage unit used for TB (Napier barcode) and *E. coli* (ST)
  does not exist here.
- **Near-perfect numbers on a single-mechanism cell are mechanistically expected.** Campylobacter cipro
  R is overwhelmingly one substitution (`gyrA_T86I`, 1,000 of ~1,006 counted determinants, PPV 0.995),
  so 0.9966 tests one well-catalogued point mutation, not broad catalog coverage.
- **Per-determinant PPV is MARGINAL, not causal** — carriers may hold several determinants.
- **Portal fetch coverage varies by organism and degrades on newer accessions**: Campylobacter 8.3%
  failures, Salmonella 2.1% (zero through the first ~5,800, then rising), *E. coli* gentamicin 0.6%.
  Datapoints toward `cabbage-genomes-tree-coverage-uncensused`, not a census. Failures are excluded and
  counted, never scored wrong.
- **These are not prospective** — collection years run 2000–2024, all pre-dating the active lock. And
  the arm is **namespace-separate** from the provenance-disjoint report card, so a powered Salmonella
  number here does **not** flip that arm's UNDERPOWERED verdict.
- **The portal's AMRFinder run is EBI's, not ours.** It is evidently organism-aware (point mutations
  appear), but its settings are not under our control.

## Operational note worth keeping

Plain `nohup` background jobs do **not** survive past their launching tool call on this host — two died
15–20 minutes in with no error and never appeared in the process table. **Harness-tracked** background
tasks do (~2,400 s, with a completion notification), and that is the only reason ~28,000 isolates were
affordable. Launch long work as a foreground call with an internal `timeout 2400`; do not pipe it
through `tail`, which buffers all output until exit.
