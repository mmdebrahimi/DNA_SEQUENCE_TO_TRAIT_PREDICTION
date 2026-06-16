# `/brainstorm` prompt — new-drug-coverage expansion (2026-06-16)

**Why this exists:** Adversarial pre-scrutiny of the "expand decoder drug coverage from in-hand measured
MIC" approach (see `wiki/new_drug_coverage_idea_anchor_prompt_2026-06-16.md`) BEFORE a technical plan is
drafted. `/brainstorm` is user-only → paste-ready below. Best run AFTER `/probe` lands the per-drug
R-count census (so Codex critiques real powering numbers, not assumptions), but usable now on the
approach. Targets the frozen-surface decision + the powering/intractability traps.

---

## Paste-ready command

```
/brainstorm New-drug-coverage expansion for the deterministic AMR decoder — adding tractable drugs from in-hand measured MIC without regressing the frozen cells.

## Problem Statement
The frozen v0.5.0 deterministic AMR decoder (dna_decode/eval/amr_rules.py::call_resistance — counts
curated AMRFinder determinants → R/S, NOT embeddings) covers 6 drugs: ciprofloxacin, ceftriaxone,
tetracycline, gentamicin, meropenem, oxacillin. Two fully-independent measured-MIC cohorts on disk
(Oxford ecoli_mic_arg, ~2900 isolates; Sci234 Sci Rep 2023, 234 isolates) carry measured MIC + complete
genotype for ~2-3× more drugs. The just-shipped fam.tsv Subclass resolver (scripts/fam_subclass_resolver.py)
reproduces deployed AMRFinder v4.2.7 family-level Class/Subclass per acquired allele, so a gene-presence
table scores by the frozen rule with no assembly (it just scored gent+cef on Sci234: cipro 0.987, cef
0.991/sens0.833/spec1.0 SCORED; gent UNDERPOWERED at 0 R isolates).
Success = add the tractable + adequately-powered subset of new drugs as honestly-validated cells, abstain
explicitly on the rest, regress NONE of the 6 frozen cells.

## Proposed Plan / Idea
1. Per-drug R/S count census over Oxford + Sci234 (CLSI/EUCAST breakpoints) → kill underpowered candidates
   (require ≥10 R AND ≥10 S pooled; gent already hit 0 R).
2. Tractability tiering: CLEAN = trimethoprim-sulfa (sul/dfr, no resolver needed), cefotaxime/ceftazidime/
   cefepime/aztreonam (ESBL — the ceftriaxone rule PATTERN via the resolver), amikacin (aminoglycoside
   subclass); CARBAPENEM-EXTEND = ertapenem/imipenem (meropenem rule; likely underpowered); ABSTAIN =
   ampicillin/amoxicillin (intrinsic blaEC → ~all R), co-amox/pip-tazo (enzyme+inhibitor), colistin
   (chromosomal mgrB/pmrB half-blind), ceftazidime-avibactam/temocillin (novel).
3. For each surviving drug: add a CLSI/EUCAST breakpoint + AMRFinder class/subclass filter + a
   determinant-count rule + binary R/S validation against in-hand MIC → a SCORED / UNDERPOWERED /
   ABSTAIN cell.
4. Frozen-surface handling (the key decision): EITHER (a) extend frozen amr_rules.DRUG_RULE +
   mic_tiers.DRUG_BREAKPOINTS/DRUG_AMRFINDER_CLASSES additively + re-validate the 6 existing cells don't
   regress, OR (b) build a SEPARATE non-frozen drug-rule overlay catalog the scorer reads alongside the
   frozen one (keeps the 2026-06-13 reproducibility freeze byte-intact). Leaning (b).
5. Start with trimethoprim-sulfa (cleanest, no resolver, MIC in both cohorts) as the proof-of-path.

## Constraints & Context
- Current behavior: 6 frozen drug rules; report card renders SCORED/UNDERPOWERED/ABSTAINS cells unioned
  from a shipped-surface registry; lineage-disclosure layer corrects clonality. amr_rules.py + mic_tiers.py
  + build_validation_report_card.py + compute_lineage_metrics.py + cohort_manifest.py are REPRODUCIBILITY-
  FROZEN (wiki/reproducibility_freeze_2026-06-13.md).
- Tech stack: Python; uv; pytest (1000+ tests); the decoder is a pure-function determinant counter; the
  resolver reads data/amrfinder_db/2026-03-24.1/fam.tsv (col3 gene_symbol / col16 class / col17 subclass,
  parent-node inheritance).
- Non-goals: NO embeddings (closed 0-for-4), NO new data acquisition (labels are in hand), NO money, NO
  cell-count padding (an uninformative or underpowered cell is worse than none).
- Key files to review: dna_decode/eval/amr_rules.py (DRUG_RULE, cipro_determinants_from_main,
  call_resistance), dna_decode/data/mic_tiers.py (DRUG_BREAKPOINTS, DRUG_AMRFINDER_CLASSES, classify_tier),
  scripts/fam_subclass_resolver.py (the new resolver + its documented family-vs-node scope-limit),
  scripts/sci234_score.py + scripts/oxford_score.py (the in-hand scorers), scripts/build_validation_report_card.py
  + dna_decode/data/shipped_decoder_surface.py (where a new cell would surface),
  wiki/external_validation_sci234_result_2026-06-16.md + wiki/oxford_external_validation_result_2026-06-15.md
  (the gent-UNDERPOWERED + cef-FN-attribution precedents), wiki/negative_results_map_2026-06-13.md.
- You are free to read any files in the repo you need.

## Specific hotspots to attack (don't validate — challenge)
1. POWERING: are R-counts in the in-hand cohorts actually ≥10/class for the CLEAN-tier drugs, or are most
   born UNDERPOWERED like gent? (E. coli amikacin-R, aztreonam-R, cefotaxime-R prevalence in a bloodstream
   cohort — likely thin for several.) Is the whole expansion mostly underpowered cells?
2. FROZEN-SURFACE: is the separate-overlay (b) actually cleaner, or does it fork the rule logic into two
   places that drift (the shared-key silent-overwrite + bundle-contract-drift lessons)? Does the report
   card's shipped-surface union handle an overlay cell without editing the frozen builder?
3. INTRINSIC-GENE / WRAPPER traps: is the trimethoprim-sulfa sul/dfr rule a real decoder or just
   re-reading AMRFinder's call (validate-wrapper-vs-underlying-tool)? Does any new β-lactam class filter
   over-match a shared AMRFinder Class token and perturb the FROZEN ceftriaxone/meropenem cells?
4. RESOLVER FIDELITY: the family-vs-node scope-limit already cost cef one FN (blaTEM-52). For aztreonam/
   ceftazidime (where TEM/SHV-ESBL variants matter more than for ceftriaxone), does family-level
   resolution under-call badly enough to make the cell misleading?
5. BREAKPOINT CURATION: are the CLSI 2024 / EUCAST 14.0 E. coli breakpoints for each new drug correct, and
   is CEFOTA(cefotaxime)-as-ceftriaxone-proxy / any cross-drug proxy defensible?
6. WORTH-IT: which of these drugs, even if SCORED, tells the user something the existing 6 don't — vs
   padding the count? Is there a drug NOT on the list that's more worth adding?
```

---

## What to expect / how it chains
- Likely Codex critiques: most CLEAN-tier drugs underpowered in these cohorts (powering kills the slate);
  overlay-vs-extend has a drift failure mode; trimethoprim-sulfa is the only safe first cell.
- After brainstorm: fold accepted findings → `/technical-plan` (the expansion), then pre-exec `/brainstorm`
  (frozen-surface = shared-state trigger), `/save-plan`, `/execute-plan`. Per planning-pipeline STOP, each
  step waits for you.
- Honest flag: if the R-count census (step 1 / `/probe`) shows the CLEAN tier is mostly underpowered, the
  right outcome may be "trimethoprim-sulfa only, bank it" — not a multi-drug expansion. Let the data decide.
