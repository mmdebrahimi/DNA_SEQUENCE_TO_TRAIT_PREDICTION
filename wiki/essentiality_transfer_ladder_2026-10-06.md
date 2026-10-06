# Cross-organism essentiality transfer ladder — 2026-10-06

Applies the conserved-core essentiality decoder **UNCHANGED** at increasing phylogenetic distance from
the organism it was tuned on, to separate **(A)** continuous decay with distance from **(B)** a cliff at
the eukaryote boundary.

**Primary metric: `coverage_lift` = coverage(essential) − coverage(non-essential).** Both terms are
conditioned on the class, so it is base-rate robust and therefore genuinely cross-rung comparable. AUROC
is reported as a **secondary** and is deliberately **not** compared across rungs.

## Verdict

**PRIMARY — bacterial_technology_matched_subladder:** `INDETERMINATE_INSUFFICIENT_RUNGS`

> 1 scored rung(s) < min_rungs=3

**SECONDARY — full ladder (technology-confounded AND sampling-frame-confounded):** `INDETERMINATE_INSUFFICIENT_RUNGS`

Scored rungs: **2** · walled rungs: **3** · frozen bar: `{"plateau_tol": 0.06, "cliff_drop": 0.12, "min_rungs": 3}`

## Per-rung table

| rung | organism | depth | tech | cov(ess) | cov(non) | **coverage_lift** | adjusted | null p95 / MAX | state |
|---|---|---|---|---|---|---|---|---|---|
| `ecoli` | *Escherichia coli* | 8 | transposon_insertion | 0.3789 | 0.0664 | **0.3125** | 0.3213 | 0.0298 / 0.0518 | scored |
| `paeruginosa` | *Pseudomonas aeruginosa* | 5 | transposon_insertion | — | — | — | — | — | **WALL_SCHEMA_UNVERIFIED** |
| `saureus` | *Staphylococcus aureus* | 2 | transposon_insertion | — | — | — | — | — | **WALL_SCHEMA_UNVERIFIED** |
| `scerevisiae` | *Saccharomyces cerevisiae* | 1 | deletion_collection | — | — | — | — | — | **WALL_SCHEMA_UNVERIFIED** |
| `human` | *Homo sapiens* | 1 | crispr_ko | 0.1689 | 0.0056 | **0.1633** | 0.1787 | 0.0214 / 0.0446 | scored |

## Distance classes (derived, not asserted)

Ordering comes from NCBI Taxonomy lineages that were fetched and name-verified, never remembered.

- depth **8**: `ecoli`
- depth **5**: `paeruginosa`
- depth **2**: `saureus`
- depth **1**: `human`, `scerevisiae`  ← TIED: one distance class, so their spread is a consistency check, not a rung ordering

## Reconciliation

Before any rung was scored, the committed numbers were recomputed and matched:

- `ecoli`: coverage_lift **0.3125**, AUROC 0.6952 (both reproduced)
- `human`: coverage_lift **0.1633**, AUROC 0.5805 (both reproduced)

## Threshold derivation

- **plateau_tol:** ceil(max null_MAX across the two existing rungs) = ceil(0.0581) -> 0.06; the null MAX rather than p95, matching this repo's practice of clearing the maximum
- **cliff_drop:** 2 x plateau_tol -> 0.12; a drop must be twice the noise floor to be a cliff
- **min_rungs:** 3; two rungs are the motivating observation, not a shape
- **observed_two_rung_spread:** 0.1492
- **non_vacuity:** 0.1492 > cliff_drop 0.12 (a cliff is reachable) and > plateau_tol 0.06 (a plateau is not trivially satisfied)
- **HONEST_LABEL:** DESCRIPTIVE-CONSISTENCY LOCK, not a frozen-before-the-numbers endpoint test: these are derived from the two rungs that motivated the question, so a bar built from them cannot also test them.

## Honest limits

These are part of this template and are emitted on **every** run regardless of verdict.

### Label technology is confounded with distance across the full ladder
The ladder spans transposon-insertion → deletion-collection → CRISPR knockout screens, and that is
irreducible with these sources. It is the reason the **technology-matched bacterial sub-ladder is the
PRIMARY test** and the full five-rung ladder is only secondary. The bacterial sub-ladder has at most
three rungs, so its power is low — disclosed, not solved.

### Cross-rung AUROC is a DIFFERENT SAMPLING FRAME and is refused here
BAGEL CEGv2/NEGv1 are two curated extremes at base rate 0.431; the E. coli set is genome-wide at 0.0928.
Those are not samples from equivalent task populations, so comparing AUROC — or sensitivity, specificity
or precision — across rungs would compare unlike quantities. This memo leads with `coverage_lift`, which
is class-conditioned and therefore immune to that difference. **No cross-rung AUROC or sens/spec
comparison appears in this document.**

### Roughly 10% of the human miss is a catalogue PHRASING gap, not phylogeny
Measured in `wiki/essentiality_missed_vocabulary_2026-10-05`: host-specific absence is 141/566 = 0.2491
while the phrasing floor is 57/566 = 0.1007, so the ladder is not invalidated — but a conserved function
under a bacteria-specific name is a catalogue bug and must not be read as distance or as a domain
boundary. `coverage_lift_adjusted` ships beside the raw lift for exactly that reason, and it is a
**FLOOR**: the true share sits between ~0.10 and ~0.25 and vocabulary overlap cannot narrow it, so the
adjustment deliberately under-corrects.

### A rung is ONE organism
A between-rung difference is not separable from that organism's annotation quality, curation depth or
genome size. `coverage_lift` removes the base-rate confound; it does not remove this one. Note also that
the rungs differ in how their classes are sourced — a rung in `essential_plus_complement` mode has an
**annotation-defined** negative class whose size depends on the annotation rather than on a screen.

### The bar is an asserted, derived consistency lock — not an endpoint test
`MIN_JOIN_RATE`, `cliff_drop` and `plateau_tol` are asserted bars. The thresholds are *derived* from
measured permutation noise rather than guessed, but they are derived from the two rungs that **motivated**
the question, so a bar built from them cannot also test them. Per-rung numbers ship so no reader is
limited to the verdict.

## Reproduce

```bash
uv run python scripts/essentiality_transfer_ladder.py --self-check     # offline
uv run python scripts/essentiality_transfer_ladder.py --emit-memo      # artifact + this memo
```

Exit 0 = a verdict · 3 = INDETERMINATE · 2 = the reconcile gate failed, meaning nothing downstream is
trustworthy.
