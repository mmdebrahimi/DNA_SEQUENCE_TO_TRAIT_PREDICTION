# Execution Log — Meropenem_FN16_Species_Composition_Audit_Plan
Date: 2026-10-05
Waves: 6 (max parallelism 1)
Files changed: dna_decode/data/genbank_organism.py, dna_decode/eval/cohort_species.py, dna_decode/data/trust_surface.py, dna_decode/amr/cli.py, scripts/provdisjoint_species_audit.py, scripts/build_validation_report_card.py, scripts/build_identify_reference.py, tests/test_genbank_organism.py, tests/test_cohort_species.py, tests/test_provdisjoint_species_audit.py, tests/test_build_validation_report_card.py, tests/test_amr_cli_disclosure_layers.py, tests/test_evidence_surface_reachable.py, wiki/meropenem_fn16_diagnosis_2026-10-05.md, wiki/provdisjoint_species_audit_2026-10-05.{md,json}, wiki/decoder_validation_report_card.{md,json}, wiki/decisions-log.md, CLAUDE.md, LESSONS_LEARNED.md, NEXT.md, README.md, docs/ARCHITECTURE.md, .claude/testable-modules.md
Sentrux verdict: n/a — sentrux not installed
Commit: 98bcffe (plan steps 909c1f6, cb20c11, cc7da78, a555972, ee76f96, 804522e; regression fix 15fa1c4; epilogues f3f442f, bd10dfe, 98bcffe)

## Result

All 6 steps completed. Tests 5694 → 5760 (+66), 12 skipped, **0 failed**.

**The finding.** All 16 false negatives in the `Klebsiella × meropenem` SCORED cell are
*Klebsiella aerogenes*; not one is *K. pneumoniae*. The cohort is 38% off-species (37 *pneumoniae* +
23 *aerogenes* of 60) scored throughout with `-O Klebsiella_pneumoniae`. Per species: pneumoniae
`tp 12 / fp 3 / tn 22 / fn 0`, aerogenes `tp 2 / fp 0 / tn 5 / fn 16` — so on the species the rule was
validated on there are **zero misses, sens 12/12**, and the published sens 0.467 is substantially a
measurement of cohort composition. The only prior explanation on record — the frozen rule being "blind to
porin-loss-mediated R" — is false twice over: the rule counts porin truncations, and these isolates have
none.

**Not systemic.** All 10 provenance-disjoint cohorts audited: Campylobacter and all four E. coli cohorts
match what they were scored as; every Klebsiella cohort is mixed but the other four carry 1–2 *aerogenes*
(2–3%) against meropenem's 23 (38%), and meropenem is the only one where off-species isolates account for
the misses. Three Klebsiella cohorts have their cross-tab **withheld** (cached runs 3/60, 33/60, 54/60)
rather than attributed on a partial denominator.

**Disclosed, not re-scored.** The 10 `provenance_disjoint_validation_*.json` artifacts are frozen units of
the reproducibility freeze and are byte-unchanged. The pneumoniae-only matrix (sens 1.000 / spec 0.880 /
acc 0.919, N=37) ships as a *measurement*; publishing it as the cell's number is a **user authority call**.

## Deviations from the plan (as-built)

1. **Step 1 placement.** The plan specified `dna_decode/data/annotations.py` on cohesion grounds. Measured
   import cost 1.399s vs 0.241s (annotations pulls pandas) and both existing consumers are pandas-free, so
   it went into a new stdlib-only `dna_decode/data/genbank_organism.py` — which also matches the plan's own
   test filename. Verified after: `build_identify_reference` still does not pull pandas.
2. **Open Question 2 was not open.** The plan deferred "should this layer reach a CALL too?" as a design
   question. An enforced guard had already answered it (*"a layer nobody can see from the tool is not a
   disclosure"*), and deferring it produced a mid-run regression. Wired, not argued.
3. **Sidecar key shape.** The plan implied a dict keyed by `canonical_cell_key`; the sibling
   `provdisjoint_source_concentration.json` uses a LIST of `{organism, drug, …}` whose loader re-derives
   the key. Matched the sibling rather than inventing a second convention.
4. **More reuse than the survey found.** `read_selected`, `parse_cohort_dir` and `find_artifact` were
   reused alongside `reconcile_raw_metrics`, all from `compute_lineage_metrics`.
5. **A verdict constant was renamed after real data falsified it.** `SINGLE_SPECIES_AS_EXPECTED` →
   `MATCHES_SCORED_ORGANISM`: the Campylobacter cohort is legitimately 31 *C. jejuni* + 9 *C. coli* under a
   genus-level `-O` flag, so the first name asserted something false about a cohort it applied to.

## Two regressions introduced mid-run, both fixed (15fa1c4)

- The card layer had to reach a CALL (see deviation 2).
- I then wired **one of two** CLI render loops. Every reachability guard still passed because `trust_block`
  carried the block; only invoking the real bacterial CLI showed nothing printed. Same class as the
  documented salmserovar one-of-two-selection-points defect. A new guard counts the loops and is proven to
  catch exactly that bug, naming the offending loop.

## Verification

Frozen five sha256 byte-unchanged · prospective lock `2026-08-31` re-verified `ok=True` · all 10
provdisjoint artifacts byte-unchanged · report card augment-only by diff (27 cells before and after,
identical key set, every differing field under `species_composition`, state counts identical) ·
contract-number audit clean (0/215 across 25 cells) · project scope unchanged (50 entrypoints / 46 traits /
128 cells — this plan adds no decoder).

**Open, and deliberately not taken:** re-scoring the cell, and dropping lone-porin counting from the
meropenem rule. Both edit frozen state; both are recorded in `NEXT.md` as user authority calls.
