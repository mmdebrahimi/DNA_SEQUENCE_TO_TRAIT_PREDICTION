# Off-pathogen cell admission gate — label-first, not trait-first (2026-07-01)

Policy for adding ANY new off-pathogen decoder cell (human or non-human trait). Distilled from the
off-pathogen brainstorm (2026-07-01): the binding constraint across this whole project is **LABELS**, so a
new cell must clear a **dataset** gate BEFORE the trait is chosen. Hunting an exciting trait and then looking
for labels recreates the OpenSNP/PGP weak-label problem.

## The gate — a candidate cell is admitted ONLY if ALL hold

1. **Measured (not self-reported) free label** — a lab/instrument-measured phenotype, publicly downloadable
   with **no DUA / no application**. (Self-report is acceptable for a PILOT/DEMO tier only, never for a
   high-confidence claim.)
2. **Stable identifiers** — a durable key joining phenotype ↔ genotype (accession / sample ID), not a
   volatile survey row.
3. **Genotype access** — the genotype (SNP/VCF/DTC) for the labelled samples is fetchable at the same bar.
4. **A deterministic rule KNOWN IN ADVANCE** — a curated, sourced genotype→phenotype rule exists before
   validation (we validate a rule; we don't fit one — the embedding/ML arm is a closed negative).
5. **Provenance + privacy posture** — rule coefficients cross-checkable to a primary source; any
   human-genotype cache stays local-only, outside the repo, never committed.

If a candidate fails ANY item → **do not add it.** No trait is worth recreating the label wall.

## Wording discipline (applies to every cell)

- Reserve **"generalizes"** for cells validated on MEASURED labels. For self-report pilots, the honest claim
  is the ENGINEERING one: *"the decoder can host deterministic, provenance-bearing non-microbial rules with
  abstention + honest tiers"* — not broad biological generalization.
- A validated-published-tool cell (IrisPlex, ABO) is an **integration/tool-validity** result, not a novel
  scientific finding — label it as such.
- Report the **deployed decision rule** (e.g. IrisPlex's 0.7 threshold), not just a permissive argmax.

## Current off-pathogen cells vs the gate (2026-07-01)

| Cell | Label | Tier | Gate status |
|---|---|---|---|
| Eye colour v0 / v0.1 | OpenSNP + PGP self-report | **PILOT/DEMO** | fails #1 (self-report) → demo only |
| ABO blood type | PGP self-report | **PILOT/DEMO** (n=6) | fails #1 → demo only |
| PGx (existing track) | GeT-RM consensus truth set | truth-set caller-concordance | passes as CALLER validation, NOT measured phenotype — frame accordingly |

## Consequence

Human self-report expansion (more OpenSNP/PGP sweeping) is **deprioritized** — it adds N inside a weak-label
regime + privacy burden without clearing the gate. The next off-pathogen cell should come through this gate
(a measured free label found FIRST), or the off-pathogen demonstration is banked at pilot strength as-is.
Non-human measured systems (e.g. yeast/plant lab-measured traits) may fit — but only after the dataset
clears items 1-5.
