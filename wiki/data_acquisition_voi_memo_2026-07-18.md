# Data-acquisition VOI memo — which label source maximally unlocks the decoder

**Date:** 2026-07-18 · **Status:** DECISION AID (the acquisition itself is a user authority/money call) ·
**Grounding:** `wiki/negative_results_map_2026-06-13.md` (the 8 rejection gates) + a 2026-07-18 web-research
pass on public genome↔measured-MIC datasets + the project's regime boundary.

## The binding constraint (recap)

The models are good enough; the **labels** are the wall. Three regimes:
- **Molecular (mutation → protein works?)** — SOLVED with free DMS; acquisition adds little.
- **Clinical resistance / MIC (mutation → drug-resistant?)** — blocked by *label* problems (circular /
  censored / surveillance-dominated). **A clean label fixes exactly the gate that blocks it.**
- **Organism traits (genome → whole-organism trait)** — blocked by a *model* ceiling (learns ancestry).
  More data does NOT fix this. Acquisition here is low-VOI.

So maximum VOI = a label source that unblocks the **middle** regime, where every ingredient except the label
is already built and validated.

## Screening rubric (a candidate must clear all four)

1. **Not circular (G1):** the phenotype is a wet-lab/clinical *measurement*, not a genome-derived call
   (AMRFinder/BlastFrost/ML-from-sequence). *This is the single most-violated gate.*
2. **Measurement, not sampling (G3) + separable provenance (G7):** an assay reading with per-record
   submitter/center/date, not a description of where the isolate came from.
3. **Uncensored where it matters (G6):** exact MICs (not `>X`/`≤X` piled at the breakpoint).
4. **Survives dedup (G4/G8):** a usable cohort (≥20/class, ≥~3 effective lineages) AFTER removing the
   surveillance ecosystem and Mash-lineage clonality.
Plus two VOI multipliers: **which regime it unlocks** + a **proven exploitation mechanism** (the supervised
complement already rescued the catalogue blind spot for the CONVERGENT-evolution pathogen HIV — a clean
measured-AST set for another convergent pathogen is the concrete way acquisition creates value).

## The free sources are exhausted — confirmed, not reopened

The 2026-07-18 research named the best free candidates; all are ones we already mined:
- **BV-BRC `PATRIC_genomes_AMR.txt`** — DOES separate measured "AMR Panel" records from "Computational
  Prediction", BUT the measured fraction is small + method-unlabelled + breakpoint-censored (our census:
  91% ML-derived, ~70% censored, within-class exact N single digits). Trips G1/G6/G8 after filtering.
- **NCBI Pathogen Detection structured antibiogram** — measured MICs linked to genomes across 100+ taxa; this
  is the source the 10 SCORED cells + the prospective-lock already use. Surveillance-dominated for the
  non-TB drugs we tested (G4). It is the deposition STANDARD a partner would use — useful to know.
- **KlebNET / Pathogenwatch** — AST is submitted-metadata-dependent + the core output is Kleborate
  (genotype-derived) → circular where it's dense (G1).
These CONFIRM the prior verdict: **no free public source clears the gates for a new drug/organism the
deterministic decoder doesn't already cover.** (The MIC-continuous negative stays closed — not reopened.)

## Ranked acquisition candidates

**🥇 Tier 1 — a CRyPTIC-style clean measured-MIC + WGS set for a CONVERGENT-evolution pathogen (non-TB).**
CRyPTIC is the existence proof: ~12k *M. tuberculosis* isolates with BMD-MIC + WGS, and it is exactly the
source that let the TB cell score. The VOI move is the bacterial/HIV analogue — a coordinated set of
isolates with **exact BMD-MIC** + genomes + per-isolate provenance, for an organism where resistance is
convergent (so the supervised complement, proven on HIV, transfers). Clears G1 (real BMD), G6 (exact MIC),
G7 (consortium metadata). **Sources to pursue:** (a) an academic clinical-microbiology consortium doing
coordinated MIC+WGS (the CRyPTIC model, replicated for e.g. *E. coli*/*Klebsiella* bloodstream isolates);
(b) a hospital clinical-micro lab willing to share deidentified isolate-level AST + genomes under a DUA;
(c) pharma AMR surveillance programs (SENTRY / ATLAS / INFORM) that hold measured MIC and sometimes genomes —
partnership/licence, likely money-gated. **Cost:** a relationship + a data-use agreement (Tier-1a/b) or a
licence (Tier-1c). **This is the recommended target.**

**🥈 Tier 2 — a wet-lab collaborator for a NEW molecular property (activity/stability on an unseen protein
family).** Extends the regime that already WORKS (molecular fitness), so it is low-risk and the tooling
(forward hybrid) is done. Free-ish (academic collaboration). Lower ceiling than Tier 1 but a safe bet.

**🥉 Tier 3 — a human biobank (UK Biobank / All of Us).** Headline VOI is large but it lands in the
organism-trait regime that fails for *model* reasons (ancestry confounding, our 0-for-5 wall), and access is
gated on institutional affiliation. **High cost, wall not actually broken by the data.** Deprioritize.

**Explicitly NOT worth acquiring:** more free public AMR/MIC (gates already tripped); more DMS (molecular
cell saturated); embedding-scale compute (closed negative).

## What converts each external wall into a code wall

A Tier-1 acquisition converts the label wall into pure engineering we've already validated: **measured AST →
`call_resistance` + the supervised complement (HIV-proven) → the report card + prospective-lock harness**, all
built. The moment a clean isolate-level measured-AST set for a convergent pathogen exists, the decoder extends
by construction — no new science, just the label.

## The concrete next step (user authority)

Not code — a *relationship*: identify a specific clinical-micro lab / AMR consortium / pharma program and open
a data-use conversation for isolate-level **exact-MIC + WGS + provenance**. The NCBI structured-antibiogram
template is the deposition format to request. Everything downstream is built and validated.

## Honest scope

This memo ranks options against the project's own gates; it does not perform an acquisition (money/partner =
user authority). The free-source verdict is a re-confirmation of the closed negative, sharpened with named
targets + the CRyPTIC precedent + the NCBI deposition standard — not a reopening.
