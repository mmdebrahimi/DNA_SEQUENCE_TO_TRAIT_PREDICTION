# The TB decoder was built, validated, and unreachable — `dna-tb` (2026-10-04)

## What was wrong

`dna_decode/organism_rules/` has held a working M. tuberculosis AMR decoder since 2026-06-17:
`tb_amr.score_drug`, a sha256-pinned WHO mutation catalogue v2 loader, masked-VCF determinant calling,
Napier-barcode lineage collapse, and **validated artifacts for twelve drugs**. `DrugCall` already
carries `rule_status` / `rule_scope` / `catalogue_commit` — it was *built* to be a contracted cell.

It had **no console script and zero `cell_registry` cells**. So:

- a user could not run it at all (`dna-decode tb` did not exist);
- the live evidence surface did not contain it, meaning the project's own trust surface showed no TB
  decoder while TB numbers were being cited in `CLAUDE.md`.

This is the **HCMV registry gap one step worse**. HCMV (found 2026-09-01) was at least CLI-routable and
merely uncontracted. TB was neither.

## Why the existing guards did not catch it

`test_every_trait_has_an_evidence_contract` and `test_traits_registry_matches_console_entries` both key
off `dna_decode.cli.TRAITS` and `pyproject.toml [project.scripts]`. TB appeared in **neither**, so there
was nothing for either guard to compare. A guard over a registry cannot see a decoder that never entered
the registry — the same shape as the hand-enumerated drug union that hid HCMV
(`feedback_hardcoded_exclusion_list_undercovers`, now its sixth instance).

## What shipped

`dna-tb --vcf X.vcf --drug rifampicin [--regeno Y.vcf] [--json] [--list-drugs]`, registered in all four
places (console entry, `TRAITS`, dispatch, `cell_registry` + `cell_regime`).

### Two evidence tiers, deliberately not collapsed

| drugs | tier | number |
|---|---|---|
| rifampicin, isoniazid | `NEAR_INDEPENDENT` | lineage-collapsed RIF sens **0.444** / spec **0.979**; INH **0.321** / **0.972** |
| the other ten | `KNOWLEDGE_BASELINE` | **none quoted** — CRyPTIC in-distribution only |

RIF/INH were scored on the EBI AMR-Portal cohort (N=2,845), provenance-disjoint from the CRyPTIC set the
WHO catalogue was partly built from, with measured wet-lab DST — the same evidence shape as the ten
NCBI-PD bacterial cells, hence `NEAR_INDEPENDENT` rather than `INDEPENDENT_MEASURED`. The other ten have
only a CRyPTIC-scored number, and the catalogue was built partly from CRyPTIC, so that is a knowledge
baseline. **One tier for all twelve would over-claim ten or under-claim two**, and under-claiming is as
much a trust-surface falsehood as over-claiming (the `forward` cell lesson).

**The quoted number is the lineage-collapsed one.** TB is monomorphic and heavily clonal: 2,845 isolates
collapse to ~67 barcode lineages, so raw RIF sens 0.920 is clonality-inflated against 0.444. A test pins
that the raw figure never appears unlabelled.

### A live defect found on real data, not by reading

The first version returned **S for every isolate**. `parse_masked_calls` requires `FILTER==PASS`, which
is correct for the CRyPTIC **minos** VCFs it was written against — but the minimap2/paftools.js VCFs the
AMR-Portal arm actually produced write FILTER `.`. Measured on `SAMEA1015921.vcf`: **1,382 genuine
`GT=1/1` variant rows, zero parsed calls, a confident S.**

All-S across 25 isolates is the signature of a broken match, so it was checked rather than attributed to
the cohort. The validated runner already had the fix (`run_tb_independent_amr_portal.py::_pass_mask`);
the CLI did not. After normalising, the same isolate calls **R on `rpoB_p.Ser450Leu`** — the canonical
rifampicin determinant, and the exact variant the 2026-06-17 coordinate-alignment probe verified at
761155 C>T. Validated against known biology, not "the number moved".

Across all 63 local VCFs the determinant distribution is textbook: `katG_p.Ser315Thr` 44×,
`rpoB_p.Ser450Leu` 35×, `inhA` promoter mutations as the secondary INH mechanism, rare rpoB RRDR
alternatives in single digits. Independent corroboration that the catalogue join is right.

**Promoting `.` is safe and that is argued, not assumed:** `.` means *no filters applied*, and a caller
that does filter marks failures with a filter NAME (`MIN_DP`), never `.`. Only `.` is rewritten; a real
PASS/filter vocabulary passes through untouched (pinned by test). **The consequence is disclosed**: `.`
also means there is no caller-applied quality floor, so `parse_masked_calls`'s "PASS subsumes the
MIN_DP/MIN_FRS floor" reasoning does not hold for such a VCF. Both output paths say so.

### The generalizable guard

A VCF with data rows that parses to **zero** calls now **refuses (exit 4)** rather than reporting S.
Both states would otherwise print `PREDICTION: S` — the indistinguishable-failure shape, and exactly how
this bug presented. A future caller format this normaliser does not understand now fails loudly instead
of calling an entire cohort susceptible. An *empty* VCF (zero data rows) deliberately does **not** trip
it: that is a legitimately empty file, not a parse failure.

Catalogue refusal is the sibling rail: absent or pin-mismatched → **exit 3 before any scoring**, with a
test that `score_drug` is never even reached, because "no determinant found" and "no catalogue loaded"
must not both print S.

## Two structural decisions

**`tb` is its own `track`.** Not taxonomy — the `amr` track is projected *verbatim* from the frozen
`shipped_decoder_surface` and `test_amr_projection_equals_frozen_surface_via_canonical_key` asserts set
equality, so an `amr`-track cell absent from the frozen surface breaks that guard. TB must stay out of
the frozen surface: adding it would edit a sha256-pinned file and **retire the active 2026-08-31
prospective lock**. Same resolution HIV/SARS-CoV-2 already use with their `viral` track.

**The routable drug set is read from the CATALOGUE, never from the contracts.** `all_routable_tb_drugs()`
in `routable_drugs.py` reads `DRUG_CATALOGUE_NAME` — the same source `--drug` validates against — so the
coverage test compares two *independent* sets. Deriving it from `_tb_contracts()` would compare a set to
itself and pass however many contracts were missing, which is precisely how HCMV shipped five decoders
with none. **Verified non-vacuous**: deleting one drug's contract fails the test and names it
(`missing={'moxifloxacin'}`); restoring it goes green.

## An output ambiguity fixed

`rule_status` is the library's label for how the **rule** was derived (a curated catalogue →
`KNOWLEDGE_BASELINE`); the evidence line is about the **cohort** it was scored on (provenance-disjoint
for RIF/INH). Printed bare and adjacent they read as a contradiction — `KNOWLEDGE_BASELINE` directly
above `INDEPENDENT`. Relabelled `rule deriv.:` / `cohort ev.:`; `tb_amr.py` untouched, since the
validated runners share it.

## Honest limits

- **G5 assembly attrition, 89%.** The 2,845 are the assembly-available subset of 26,941 not-leaked
  phenotyped isolates, so the cohort is **not prevalence-preserving** — status stays
  `TB_SUBSET_PLUMBING`.
- **Conservative lower bound.** The asm5 VCF route misses some determinants, and 43–44 mixed-label
  lineage clusters are *excluded* as discordant rather than majority-voted.
- **Sensitivity is low at the lineage level** — this rule misses more than half of resistant lineages.
  Specificity is high and robust (~0.97) in both the independent and in-distribution arms.
- **Ten of twelve drugs have no independent number at all.** They are a catalogue lookup.
- Not a clinical decision tool; an S call does not rule out an uncatalogued mechanism.
- A finer Mash-based lineage collapse is a **closed negative**
  (`wiki/tb_independent_mash_lineage_2026-07-09.json`): M. tuberculosis has no lineage-scale gap at Mash
  resolution, so the pinned Napier barcode is the correct partition and these numbers stand.

## Verification

Frozen five byte-unchanged (sha256 re-checked), prospective lock re-verified, AMR report card
byte-identical (TB is non-frozen and correctly invisible to it). Registry 116 → **128 cells**; 21 tests
in `tests/test_tb_cli.py`.
