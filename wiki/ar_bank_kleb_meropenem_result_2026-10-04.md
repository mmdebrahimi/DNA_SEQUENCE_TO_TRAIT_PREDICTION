# AR Bank Klebsiella × meropenem: the rule is NOT porin-blind, and that refutes its own `validated` string (2026-10-04)

**Headline: the frozen meropenem rule describes itself as *"blind to porin-loss-mediated R (expected FN
mode)"* and it is not.** AMRFinder reports `ompK35`/`ompK36` truncations with `Subclass=CARBAPENEM`,
`Type=AMR`, in `main.tsv`; the rule is `subclass_any={CARBAPENEM}, threshold=1`; so it counts them.
**Measured through the deployed `call_resistance`: 10 of 10 isolates whose ONLY carbapenem determinant is
a porin truncation are called R.**

Nothing was edited. The rule's BEHAVIOUR is unchanged and every previously-scored number already used it.
What is wrong is the rule's self-description — and that string is read by humans, including me: I used it
one message earlier to predict that sensitivity was this cell's weak axis *because of porin blindness*.
That chain is broken at the first link.

## The run

`scripts/external_cohort_revalidate.py`, frozen rule, organism triple read verbatim from the frozen cell
(`-O Klebsiella_pneumoniae` / registry `Klebsiella`). Artifact:
`wiki/external_validation_ar_bank_kleb_extval_meropenem_meropenem_kp129_2026-10-04_2026-10-04.json`.

```
strict: n_scored 48  tp 43  fp 3  tn 2  fn 0  acc 0.938  sens 1.000  spec 0.400
        n_excluded_no_assembly 81
POWERING HARD FAIL: strict scored S 5 < min_per_class 10 (underpowered)
```

**The scorer refused the specificity and that is correct** — 5 susceptible isolates cannot carry a
specificity. The G5 assembly-attrition prediction from the pre-run `/brainstorm` was exact: 129 labelled
K. pneumoniae, **81 excluded for having no assembly**, and the susceptible class collapsed 31 → 5.

## What the sensitivity does and does not mean

**sens 1.000 (43/43, fn=0) is NOT a general sensitivity estimate.** The CDC/FDA AR Isolate Bank is a
curated panel of *characterized* resistant organisms, so it is enriched for detectable mechanism by
construction: **46 of 49 run isolates (93.9%) carry an acquired `CARBAPENEM`-subclass determinant.** A
panel assembled from characterized carbapenem-resistant isolates cannot reveal a mechanism the rule
cannot see — the same structural blindness as the source-concentration finding
([[feedback_source_concentration_hides_whole_blind_spots]]) and the enrichment trap
([[feedback_enriched_archive_ascertainment_bias]]).

**But "enriched for carbapenemase" understates it**, which is how the porin finding surfaced: 10 of 49
carry **no** `bla` carbapenemase and are still called R, off porin truncations alone.

## The specificity liability, mechanistically located — and NOT a specificity claim

Of the 5 scorable susceptible isolates, split by mechanism:

| accession | label | call | carbapenem mechanism | |
|---|---|---|---|---|
| GCA_001874865.1 | S | **R** | porin-only | FALSE POSITIVE |
| GCF_001874725.1 | S | **R** | porin-only | FALSE POSITIVE |
| GCF_003204175.1 | S | **R** | porin-only | FALSE POSITIVE |
| GCA_001875745.1 | S | S | none | correct |
| GCA_003074015.1 | S | S | none | correct |

**3 of 3 false positives are porin-only; 2 of 2 correct calls carry no determinant at all.** Perfect
separation on mechanism. Biologically coherent: reduced permeability from porin loss generally confers
*reduced susceptibility*, not full meropenem resistance — the clinically resistant phenotype usually
needs porin loss **plus** a carbapenemase. Counting a porin truncation as sufficient at `threshold=1` is
therefore an over-call.

**This is a located mechanism, NOT a measured specificity.** n=5. The scorer hard-failed the powering gate
and this table does not overturn that. Read it as *"the over-call has an identified mechanism worth
testing on a powered susceptible class"*, never as *"specificity is 0.4 because of porins"*.

It is the **mirror image of the gentamicin `rmt` case**: there a determinant the rule could NOT see caused
false negatives; here a determinant the rule DOES see, but should not count alone, causes false positives.

## What this reopens

**The provenance-disjoint `fn=16` now has no explanation.** That cell (N=60, sens 0.467) was understood as
porin blindness — on this rule's own say-so. If the rule counts porin loss, those 16 missed resistant
isolates must carry *neither* an acquired carbapenemase *nor* a detectable porin truncation, and what they
do carry is unknown. That is a genuine open question, cheap to answer (their AMRFinder runs are cached),
and it was closed only by a false premise.

## Not acted on — two deliberate refusals

1. **The `validated` string is in `dna_decode/eval/amr_rules.py`, one of the five sha256-pinned files.**
   Correcting prose there would invalidate the active 2026-08-31 prospective lock and retire every
   prospective number scored against it. The correction is therefore recorded HERE and in CLAUDE.md, not
   patched. Frozen surface byte-unchanged; lock re-verified.
2. **Whether to stop counting lone porin truncations is a USER AUTHORITY call** — it edits the frozen rule,
   invalidates the lock, and is motivated by 3 isolates. Do not act on n=5.

## Honest limits

- Specificity WITHHELD (powering hard-fail, 5 susceptibles).
- Sensitivity is on an enrichment-selected reference panel; not a population estimate.
- Leakage posture DEGRADED: no BioSample preflight artifact, and 134 of 143 BioSamples appear in sibling
  AR Bank cohorts for *other* drugs. Per-drug independence is defensible for a per-drug rule but is NOT
  the drug-blind accession standard the provenance-disjoint cells meet. Do not present this beside them.
- "Panel" is a proxy for an independence unit — AR Bank cohorts carry no BioProject field.
- The porin claim rests on AMRFinder's subclass assignment; it was verified on real cached runs
  (`Type=AMR`, 114 CARBAPENEM rows) but not against an independent porin caller.
