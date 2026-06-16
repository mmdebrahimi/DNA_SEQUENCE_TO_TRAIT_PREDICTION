# Milestone synthesis — external validation + first new-mechanism cell; AMR public-label track terminal (2026-06-16)

Banks the 2026-06-14 → 2026-06-16 epoch and records the strategic fork chosen for the next epoch.

## What landed this epoch
1. **First fully-INDEPENDENT external validation of the frozen deterministic decoder**, on TWO
   measured-MIC cohorts with their own genotype callers (different country / lab / AST method than the
   US-NCBI-PD tuning provenance):
   - **Oxford (UK, ~2900 isolates):** cipro acc 0.960 / sens 0.935 / spec 0.963; gent 0.990; (cef genotype-rule
     clean — the apparent gap was an OLD-AMRFinder version confound, corrected). `wiki/oxford_external_validation_result_2026-06-15.md`.
   - **Sci234 (Spain, 234 isolates):** cipro acc 0.987 / sens 1.0 / spec 0.984; cef acc 0.991 / sens 0.833 /
     spec 1.0; gent UNDERPOWERED (0 R isolates — reported honestly, not a pass). `wiki/external_validation_sci234_result_2026-06-16.md`.
   - Leakage DISJOINT by construction (both cohorts deposit 0 NCBI assemblies; tuning is 100% GCA/GCF).
2. **fam.tsv per-gene Subclass resolver** (`scripts/fam_subclass_resolver.py`) — the keystone that lets a
   gene-presence table be scored by the frozen rule without assembly; 20 tests.
3. **First NEW-MECHANISM decoder cell beyond the frozen 6: trimethoprim-sulfamethoxazole** (folate
   pathway). Scorer-local `(>=1 sul) AND (>=1 dfr)` rule in a NON-FROZEN overlay
   (`dna_decode/data/experimental_drug_rules.py`), EXPERIMENTAL_SCORED, frozen surface byte-unchanged.
   Cross-cohort: Sci234 binary acc 0.987, Oxford 0.962, and the genotype strata REPRODUCE across both labs
   (sul+dfr dominant-R 0.986/0.939; sul-only near-zero 0.0/0.024) — the AND rule generalizes.
   `wiki/external_validation_tmpsmx_result_2026-06-16.md`.

## State of the project (true current state — supersedes the stale 2026-05-26 ledger frame)
- **The deterministic AMR decoder is the terminal honest product of the public-label AMR track**, now
  externally validated on two independent measured-MIC cohorts + extended to a 4th mechanism class.
  Reproducibility-frozen (2026-06-13, commit b3761c8); 10 provenance-disjoint SCORED cells, lineage-disclosed.
- **The public-label AMR track is DONE (plateau, named).** The negative-results map closed every expansion
  (embeddings 0-for-4, MIC-continuous infeasible, pathotype label-blocked, provdisjoint grid saturated),
  and the 2026-06-16 in-hand-MIC probe census confirmed the slate is exhausted of clean new-mechanism wins —
  TMP-SMX was the last (levofloxacin = cipro-relabeled; cefepime = ceftriaxone-relabeled; rest underpowered).
  The binding constraint is **labels, not models.** Further AMR drug/cohort additions are motion, not signal.

## Strategic fork chosen for the next epoch (user decision 2026-06-16)
Of the three non-foreclosed forward paths, the user chose **ACQUISITION of a new (non-public) label source** —
the biggest lever, because it clears the label gates by construction and reopens the learned-decoder niche.
- Not chosen (still valid later): eukaryotic Path B G2 embedding test (the open learned-decoder experiment,
  pre-staged on the Precision 7780); prospective-lock validation of the frozen decoder.
- Executor-doable part of acquisition: a `/research`-backed RANKED SHORTLIST of concrete acquirable sources
  scored against the 4 gates (sampling-independent / non-circular / organism-depth / provenance-disjoint-feasible)
  + acquisition path + rough N-after-filters. The actual acquisition (MTA / contact / clinical export) is the
  user's real-world action. Framing prompt: `wiki/next_epoch_idea_anchor_prompt_2026-06-13.md`.

## Pointers
- Reproducibility freeze: `wiki/reproducibility_freeze_2026-06-13.md`. Negative-results map: `wiki/negative_results_map_2026-06-13.md`.
- New-drug-coverage probe + brainstorm (this epoch): `wiki/new_drug_coverage_{idea_anchor,brainstorm}_prompt_2026-06-16.md`.
- Executed plan: `executed_plans/TMP_SMX_External_Validation_Cell_Plan/`.
