# Execution Log — Constraint_Hierarchy_Type_A_Laws_Clade_Specialization_And_Kill_Count_Reporting_Plan
Date: 2026-10-01
Waves: 4 (max parallelism 4) — executed SEQUENTIALLY; see Mode note
Files changed: .gitignore data/mycoplasma_ref/Mgenitalium_G37_MG_382_udk_cds.fna dna_decode/constraints/__init__.py dna_decode/constraints/codon_tables.py dna_decode/constraints/registry.py dna_decode/constraints/report.py dna_decode/constraints/specializations.py dna_decode/constraints/universal.py dna_decode/forward/genome_edit.py dna_decode/forward/inverse.py dna_decode/typing/codon_map.py scripts/build_holdout_hybrid_manifest.py scripts/constraint_kill_count_probe.py scripts/forward_inverse_roundtrip.py scripts/fungal_erg11_caller.py scripts/hiv_esm_vs_catalog.py scripts/mavedb_cpu_smoke.py scripts/pear_genotype_alphabet.py scripts/pre_resistance_base_rate_census.py scripts/pre_resistance_predictive_backtest.py tests/test_codon_consumers_behavioural.py tests/test_codon_consumers_canonical.py tests/test_codon_table_single_source.py tests/test_constraint_kill_count.py tests/test_constraint_kill_count_probe.py tests/test_constraint_report_and_probe_gaps.py tests/test_constraints_codon_tables.py tests/test_constraints_codon_tables_accessors.py tests/test_constraints_registry.py tests/test_constraints_registry_resolution.py tests/test_constraints_specializations.py tests/test_constraints_universal.py tests/test_constraints_universal_and_regime_gaps.py tests/test_script_codon_consumers_canonical.py wiki/constraint_kill_count_2026-10-01.json 
Sentrux verdict: n/a — sentrux not installed
Commit: e870910 (committed directly to main; no PR — standing single-owner directive for this repo)

## Mode note
The toolkit dispatch indicated PARALLEL (waves of 3 and 4, remote present, 0 errors). SEQUENTIAL was
chosen deliberately on the standing user directive recorded in ~/.claude/session-board.md (2026-08-17
ownership consolidation: dna_decode is single-owner, commit+push directly to main), reinforced by the
board-recorded prior collision of parallel lanes in this repo's FBA/nitrogen lane. The three worktrees
on this host (lane/refseq-cache-integrity, lane/fba-structure, lane/nt-transformers-fix) are 5-6 weeks
stale, belong to other lanes, and were NOT touched.

## Tests
Baseline 5157 passed / 11 skipped / 0 failed  ->  final 5424 passed / 12 skipped / 0 failed.
+267 tests, zero regressions. The +1 skip is a pre-existing in-process import limitation of
scripts/forward_inverse_roundtrip.py, covered statically instead.

## Corrections made DURING execution (each measured, not argued)
1. NCBI tables 2/4/6 were gated [unverified] in the plan; the fetch succeeded and CORRECTED the plan --
   table 2 differs from standard in FOUR codons (TGA/ATA/AGA/AGG), not the two assumed.
2. start_codon_present was predicted inert and refuted 5 of 8 reference CDS -- all five CORRECTLY, as
   they are polyprotein-derived extracts (HIV RT/PR/IN/CA begin CCC/CCT, SARS-CoV-2 Mpro begins AGT).
   The law now abstains unless the caller asserts cds_is_complete_gene; a false refutation is the one
   failure mode that would make a constraint filter delete valid predictions.
3. The discovery guard found FIVE more codon-table copies the plan's grep missed (true count 11, not
   6) -- the hand-enumerated-list failure hitting the plan that built the guard against it.
4. dna_decode/data/amr_rules.py DOES NOT EXIST (real path dna_decode/eval/amr_rules.py), so ad-hoc
   `git diff --quiet` frozen-surface checks against it passed VACUOUSLY. The surface itself was never
   unverified (the prospective lock hashes the correct file), but the manual check was not doing what
   it claimed. The guard now asserts each frozen path exists before diffing.
5. Post-epilogue, three defects in this plan's OWN output were fixed: the probe's exit-1 branch was
   structurally unreachable (the exact shape of two previously-retracted guards), STOPS was the
   clade-blind trap surviving in the module's one derived constant, and the registration helpers
   failed open on a provenance refusal.

## Measured impact
341 of 476 (71.6%) clean real Mycoplasma genitalium G37 CDS are false-nonsense under NCBI table 1 and
clean under table 4. Not a live defect for the shipped cells (all use standard assignments) -- a
latent assumption that bites the moment the organism set widens, which is what Phase 5 does.

## Deliberately NOT shipped
Four clade properties (splicing, operons, HGT, ploidy) -- the provenance gate refused them because
citing them would mean writing references from recall. Declared in DEFERRED_CLADE_PROPERTIES with the
evidence each would need.
