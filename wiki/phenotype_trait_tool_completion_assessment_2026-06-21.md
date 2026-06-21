# Bacteria/Virus phenotype→trait tool — completion assessment + wave map (2026-06-21)

Soraya conversational-executor assessment. Grounded in `plans/Trait_Decoding_Roadmap.md`,
`wiki/project_frontier_map_2026-06-19.md`, `wiki/decoder_validation_report_card.md`,
`dna_decode/data/shipped_decoder_surface.py`, and the shipped code. North star (unchanged):
**an AI DNA decoder TOOL — "DNA → what its parts do" — not papers.**

---

## TL;DR — how close are we?

**The deterministic determinant-based tool is ~75-85% of the way to a usable bacteria/virus
phenotype→trait tool, and the *architecture already spans the kingdom boundary.*** The binding
constraint is NOT model capability and NOT the decoder — it is **validation LABELS** (free,
isolate-level, de-confounded genotype↔phenotype pairs). This is the project's central, repeatedly-
verified scientific finding.

- **Bacteria: essentially DONE + validated.** 6 drugs × ~7 organisms, 10 provenance-disjoint SCORED
  cells (lineage-disclosed), reproducibility-frozen. The genome-map "honesty report" (point at ONE
  genome → evidence-tiered per-feature map) is the closest thing to the literal north star, and it just
  gained the **virulence/pathotype overlay** (this session). Bacterial AMR + virulence are the tool's
  mature core.
- **Virus: the decoder EXISTS in code, but is UNVALIDATED.** `dna_decode/data/antiviral_amr.py` +
  `scripts/flu_na_caller.py` implement an Influenza A neuraminidase (oseltamivir/zanamivir/peramivir)
  target-site resistance decoder — same deterministic pattern as the fungal ERG11 and Plasmodium
  kelch13 cells. It renders `NO_FREE_PHENOTYPE_SOURCE` on the report card because influenza NA
  resistance has no free isolate-level AST label. **The virus side is a LABEL problem, not a build
  problem.**
- **The learned/embedding arm is a CLOSED NEGATIVE** (4 de-confounded failures across the kingdom
  boundary: cipro within-lineage, pathotype, Arabidopsis flowering-time). Do NOT reopen it on free
  public labels. The deterministic decoder + the "labels not models" boundary statement ARE the
  contribution.

**So "completing the tool" does NOT mean building a learned predictor.** It means (a) UNIFYING the
shipped determinant decoders under the genome-map surface, and (b) VALIDATING the wired-but-unvalidated
cells — above all the **virus** cell, which needs a free viral genotype↔phenotype label source.

---

## What's the honest "completion" definition?

A usable v1 "bacteria/virus phenotype→trait tool" =

> Point the tool at ONE microbial genome (bacterial OR viral) → an honest, evidence-tiered per-feature
> map of what its parts do (function tiers + AMR-resistance determinants + virulence/pathotype +
> drug-resistance calls for the organism's drug panel), with every phenotype claim behind a
> validated-determinant wall, a DB-labelled unknown rate, and an explicit scope-limit when out of
> distribution.

Against that bar: bacteria ✅ (shipped + validated), virus 🟡 (decoder exists, unvalidated, not yet
wired into the genome-map surface), cross-kingdom grid 🟡 (fungal/TB/Plasmodium wired, validation
partial/blocked).

---

## Waves COMPLETED (the spine, grounded)

| Wave | What | Status |
|---|---|---|
| Phase 0-1 | cipro AMR (cached-strain + novel-genome input), cross-path concordance | ✅ shipped |
| Phase 2 | multi-drug E. coli AMR (cef/tet/gent) + the **decisive de-confound finding** (NT embedding LOSES to QRDR-POINT mechanism features; embeddings learn lineage, not mechanism) | ✅ + closed the embedding arm for concentrated mechanisms |
| Phase 3 | multi-organism AMR: E. coli + Klebsiella + Pseudomonas + S. aureus (across the gram divide); 6 drugs × 4 organisms; the "count the target determinant, not the broad class bag" principle confirmed 3× | ✅ |
| Phase 4 | non-AMR bacterial phenotype = **pathotype** (deterministic VF resolver, abstention + provenance); embedding/classifier track CLOSED (label is sampling-defined); carbon-utilization infeasible (E. coli slice too small + phylogeny trap) | ✅ resolver shipped / learned-track closed |
| Cross-kingdom | **Fungal** C. auris ERG11/FKS1 (G1 reached, LABEL_LIMITED); **TB** RIF+INH on CRyPTIC (knowledge-baseline, data-runs BLOCKED); **Virus** influenza NA + **Protozoa** Plasmodium kelch13/pfcrt (both wired in code, unvalidated) | ✅ code / 🟡 validation |
| Genome-map | v1 "Bakta honesty report" (GO) + **virulence/pathotype overlay** (this session, 1500 tests, live-verified on E. coli ST131) | ✅ shipped |
| Rigor | reproducibility freeze (10 SCORED provdisjoint cells) + lineage disclosure (clonality-corrected) + external re-validation arm (Oxford measured-MIC) + TMP-SMX experimental overlay + negative-results map (8 reusable rejection gates) | ✅ shipped |

## Waves NEXT (ranked by tool-value × not-label-blocked)

| # | Wave | Why now | Label-blocked? | Gate |
|---|---|---|---|---|
| **A** | **Unify the tool** — surface the fungal/viral/TB/Plasmodium determinant decoders as genome-map overlays (genome-map is currently bacterial-AMR + virulence only) | turns N separate CLIs into ONE "point at a genome" tool = the literal north star; pure tool-capability, no new labels needed | NO | planning chain (VCF-vs-GFF contract for TB/viral; the integrity-crux discipline) |
| **B** | **Validate the VIRUS cell** — adapt/validate a viral determinant decoder on a FREE genotype↔phenotype source. **Anchor: HIV (Stanford HIVDB)** = the canonical, historically-public viral genotype-phenotype system (RT/protease/integrase mutations → susceptibility); influenza NA stays unvalidated for lack of free labels | makes "bacteria AND virus" TRUE in the validated sense — the first validated viral cell | **MAYBE NOT** (HIVDB is the un-mined candidate; license/access verification is step 1) | full planning chain (the genotype-phenotype confound gates are exactly what the chain protects) |
| **C** | **Complete the cross-kingdom validated grid** — run the TB CRyPTIC data-runs (BLOCKED on D:/regeno cohort + hand-curated independent gold set); decide fungal C. auris label path; Plasmodium validation | closes the wired-but-unvalidated cells; TB is the highest-population unfinished cell | TB: partial (regeno compute, not labels); fungal/Plasmodium: YES (no free AST) | TB plan already saved (data-run-gated); fungal/Plasmodium need a label-source decision |
| **D** | **Strategic fork (USER decision)** — (1) label ACQUISITION (wet-lab/clinical/partnership) reopens the learned arm by construction; (2) prospective-lock validation (rigor-flavored) | the only moves that clear the label gates structurally; both are AUTHORITY decisions | n/a | user sourcing/strategy decision |

---

## Sub-plan per wave (drafted; planning-chain commands EMITTED below)

### Wave A — Unify the determinant decoders under the genome-map (strongest tool move)
**Terminal claim:** `genome_map` on a fungal/viral/TB genome surfaces its drug-resistance determinants
as overlay tiers (the fungal/viral analog of the AMR `determinant-phenotype` + the new
`virulence-determinant` tiers), behind the same coordinate-join integrity gate + presence-only wall.
**Sub-steps:**
1. `/idea-anchor` — "unify the shipped determinant decoders (fungal ERG11, influenza NA, TB RIF/INH,
   Plasmodium kelch13) as genome-map overlay tiers."
2. `/probe` — the integrity crux: fungal/viral callers are BLAST/VCF-based (ERG11 BLAST, flu NA caller,
   TB masked-VCF) — do they return coordinates the genome-map coord-join can consume? (TB/viral is the
   VCF-vs-GFF contract the v1 spike explicitly deferred.)
3. `/technical-plan` → pre-exec `/brainstorm` → `/save-plan` → `/execute-plan`.
**Falsifier:** a decoder returns only symbol-level calls (no coords) → it would symbol-fallback-NO-GO by
the genome-map's own gate (the exact trap the AMR + virulence overlays were built to avoid). Resolve in
`/probe` before building.
**Why first:** highest tool-value, NOT label-blocked, builds directly on this session's virulence-overlay
pattern (the coord-join machinery is now proven twice).

### Wave B — Validate the VIRUS cell (HIV anchor)
**Terminal claim:** a deterministic HIV RT/protease/integrase resistance decoder scored on a free,
held-out, de-confounded genotype↔phenotype set (Stanford HIVDB genotype-phenotype) → a SCORED viral
row on the report card (acc/sens/spec), OR a documented `NO_FREE_PHENOTYPE_SOURCE`/license scope-limit
if the data isn't freely usable.
**Sub-steps:**
1. **De-risk FIRST (cheap, model-invocable):** `/research` — "free, redistributable, isolate-level viral
   genotype↔phenotype drug-resistance datasets (HIV Stanford HIVDB, HCV, SARS-CoV-2 nirmatrelvir),
   screened against the project's 8 rejection gates (circular-label, study==class, sampling-defined,
   MIC-censoring, …) + license terms." (This is the GO/NO-GO; a focused web check from Soraya tripped a
   usage-policy filter this session — re-run via `/research` or an interactive WebSearch.)
2. If GO: `/idea-anchor` → `/probe` → `/technical-plan` → pre-exec `/brainstorm` → `/save-plan` →
   `/execute-plan`. Reuse the `organism_rules/` + `data/antiviral_amr.py` target-site pattern.
3. Validate held-out (HIVDB genotype-phenotype pairs are lab fold-change → de-confounded by construction,
   UNLIKE the influenza/sampling-defined blocks) + add the cell to `shipped_decoder_surface.py`.
**Falsifier:** the only free viral labels are surveillance-derived (sampling-confounded) or non-
redistributable → virus stays `NO_FREE_PHENOTYPE_SOURCE`; the influenza NA decoder remains the
honest-but-unvalidated viral artifact and the gap is a DATA gap, named not hidden.
**Why second:** it's the literal "virus" in the user's ask; HIVDB is the single best shot at the first
validated viral cell; de-risk is cheap and decisive.

### Wave C — Complete the cross-kingdom validated grid
**Terminal claim:** the TB RIF+INH v1b SCORED number (lineage-collapsed sens/spec on the prevalence-
preserving CRyPTIC cohort) + the independent gold-set number land; fungal/Plasmodium cells get an
explicit label-path decision.
**Sub-steps (TB is already planned — this is execution, not design):**
1. **TB data-runs (BLOCKED-gated by design):** stage the full per-drug prevalence-preserving CRyPTIC
   cohort (~1.6 TB regeno → D: cache) via `scripts/stage_tb_vcf_subset.py` + a regeno fetch; run
   `scripts/score_tb_cryptic.py`. **Needs the D: drive + a regeno fetch (compute, not labels).**
2. Hand-curate the post-2023 independent gold set (`wiki/tb_independent_goldset_acquisition_2026-06-17.md`)
   → `scripts/score_tb_independent_goldset.py`.
3. Fungal/Plasmodium: an AUTHORITY decision on whether to source a free label or accept the
   `LABEL_LIMITED`/`NO_FREE_PHENOTYPE_SOURCE` scope-limit (no cold-build).
**Falsifier:** TB v1b lineage-collapse blocks for lack of a lineage call on the subset
(`LINEAGE_COLLAPSE_BLOCKED_*`) → report BLOCKED honestly, never a fake metric.
**Why third:** highest-population unfinished cell (TB), but compute/label-gated; sequence after the two
unblocked waves.

### Wave D — Strategic fork (USER authority decision; not an executor task)
Per `wiki/project_frontier_map_2026-06-19.md`: (1) **label ACQUISITION** (wet-lab/clinical/partnership/
un-mined supplement) is the ONLY move that clears the label gates *by construction* and honestly reopens
the learned-decoder arm — a USER sourcing decision (draft anchor
`wiki/next_epoch_idea_anchor_prompt_2026-06-13.md`); (2) **prospective-lock** pre-registers the frozen
decoder on a not-yet-labelled cohort (leakage-impossible) — buildable but rigor-flavored, weigh against
"tool not papers." **Both are AUTHORITY decisions Soraya surfaces, not executes.**

---

## Recommendation (the autonomous call)

**Do Wave A next** (unify the determinant decoders under the genome-map) — it is the strongest pure
tool-capability move, it is NOT label-blocked, and it builds directly on this session's just-shipped
virulence-overlay coord-join machinery (now proven twice). **In parallel, de-risk Wave B** with a focused
`/research` on free viral genotype↔phenotype labels (HIV/HCV) — that single result decides whether the
"virus" half of the tool becomes a VALIDATED cell or stays an honest data-gap. Defer Wave C-TB until the
D: drive + regeno compute is available; defer Wave D to a user authority decision.

The honest meta-point (from the frontier map, reaffirmed): the project is at a **terminal honest state on
free public data** for the LEARNED arm. The remaining real growth is **tool unification (Wave A) +
new-substrate validation gated on label acquisition (Waves B/D)** — not more embedding increments.
