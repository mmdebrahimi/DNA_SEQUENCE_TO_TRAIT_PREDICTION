# Full ABO O/A/B/AB decoder — free-label validated (2026-07-07)

Completes the named-deferred A-vs-B extension of the O-vs-non-O ABO cell
(`wiki/abo_pgp_result_2026-06-30.md` flagged "A-vs-B ... NOT fabricated → deferred"). The A/B tag mapping is
**sourced, not fabricated** (ClinVar RCV000019310 + the ABO 3-variant literature) and cross-checked against
the real data.

## Rule (3-variant deterministic, UK-Biobank-style)

| variant | role | A-allele (plus-strand) | B-allele |
|---|---|---|---|
| rs8176719 | 261delG O-deletion (D/I) | — | — |
| rs8176746 | c.796C>A Leu266Met | **G** | **T** |
| rs8176747 | c.803G>C Gly268Ala | **C** | **G** |

O-background subtlety (why the deletion count is needed): the common O01 allele sits on an A-type tag
background, so a DI + het-tag genotype is phenotype **B**, not AB. `call_abo_full` uses the deletion count to
disambiguate (`dna_decode/data/abo_full.py`).

## Result — 386 OpenSNP samples, self-report labels

- **4-way accuracy (O/A/B/AB): 0.902** (348/386)
- **O-vs-non-O: 0.925**
- Indeterminate: 9 (uncallable / inconsistent tags — abstains, never guesses)

Confusion highlights: `O→A` 15 + `A→O` 10 dominate the errors — the documented **non-deletional O alleles**
(the deletion SNP misses O02/weak variants) + **~15% self-report error**. `AB` recovered 22/29.

## Honest tier

- **Label tier:** self-reported ABO (near-independent, non-circular, ~15% noisy). NOT a lab assay.
- Near-Mendelian at the common alleles; A2 / cis-AB / non-deletional O / weak subgroups NOT captured.
- **NOT a clinical tool.**

Reproduce: `uv run python scripts/abo_full_opensnp_validate.py` (reuses `data/j3_abo/j3_abo_substrate.json`;
no zip re-scan). Tests: `tests/test_abo_full.py` (6).
