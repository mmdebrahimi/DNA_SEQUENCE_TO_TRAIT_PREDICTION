# TMP-SMX external validation — first NEW-MECHANISM decoder cell (EXPERIMENTAL), cross-cohort (2026-06-16)

First drug-coverage expansion beyond the frozen 6-drug surface. Trimethoprim-sulfamethoxazole
(co-trimoxazole) is a genuinely-NEW resistance mechanism — acquired folate-pathway genes (`sul*` +
`dfrA/B*`) — distinct from the fluoroquinolone / β-lactam / aminoglycoside / tetracycline cells already
covered. Scored by a **scorer-local `(>=1 sul) AND (>=1 dfr)` rule** in the NON-FROZEN overlay
(`dna_decode/data/experimental_drug_rules.py`), validated on BOTH in-hand measured-MIC cohorts. **No
frozen file touched.** This is an EXPERIMENTAL cell (not a deployed claim) on the external-validation arm.

## Result (binary R/S = primary for clean measured MIC; CLSI/EUCAST co-trimoxazole TMP-component S<=2 / R>=4)

| cohort | n | accuracy | sens (R) | spec (S) | R / S | strict tier |
|---|---|---|---|---|---|---|
| **Sci234** (Spain) | 233 | **0.987** | **0.972** | **0.994** | 71 / 162 | 0.986 / 0.993 (n=223) |
| **Oxford** (UK) | 2866 | **0.962** | **0.922** | **0.977** | 787 / 2079 | S-degenerate (no MIC<=0.5) → binary primary |

## The de-risk that matters: mechanistic strata REPRODUCE across two independent labs
The headline sens/spec could be carried by prevalence; the real test is whether the AND rule's genotype
strata separate the SAME way on both cohorts. They do:

| stratum | Sci234 n / R-rate | Oxford n / R-rate |
|---|---|---|
| **sul + dfr** (rule → R) | 70 / **0.986** | 773 / **0.939** |
| sul-only (rule → S) | 40 / **0.000** | 338 / **0.024** |
| dfr-only (rule → S) | 6 / 0.167 | 106 / 0.179 |
| neither (rule → S) | 117 / 0.009 | 1649 / 0.021 |

`sul+dfr` is the dominant-R stratum and `sul-only` is near-zero on BOTH — confirming the **AND rule (not
OR)** generalizes. An `OR` rule would mis-call the 40 Sci234 + 338 Oxford sul-only isolates (all/near-all
S). `strata_reproduced=True` on both → cells emitted SCORED (the scorer stamps INDETERMINATE if the
strata fail to reproduce).

## Verdict
- **TMP-SMX validates as a new decoder cell on two fully-independent measured-MIC cohorts** (Spain + UK),
  with the mechanism strata reproducing — the strongest form of cross-cohort evidence the project uses.
- **Binary is the honest metric** (acc 0.987 / 0.962). The 4×-margin strict tier (HIGH_S needs MIC≤0.5)
  over-excludes co-trimoxazole and is S-degenerate on Oxford (min MIC 1) — strict reported for Sci234
  (0.986/0.993), binary headlined for Oxford. Matches the Oxford/Sci234 binary-primary precedent.
- **Bounded blind-spot (documented):** the acquired-gene AND rule cannot see `folP/folA` target
  point-mutation TMP-R (the sul+dfr-negative R isolates: Sci234 1/117 neither + 1/6 dfr-only; Oxford a
  small minority). Expected FN mode, surfaced in the artifact `fidelity_caveat`.

## Scope (honest)
- **EXPERIMENTAL_SCORED / scorer_local** — NOT in the frozen `DRUG_RULE` / `mic_tiers.DRUG_BREAKPOINTS` /
  `shipped_decoder_surface`. Renders in the SEPARATE `external_validation_report_card`, marked distinct
  from deployed-decoder cells. Promotion to the frozen deployed surface is a separate, deliberate decision.
- **Frozen files byte-unchanged** (amr_rules / mic_tiers / cohort_manifest / build_validation_report_card
  / compute_lineage_metrics) — guarded by a no-leak test.
- **Independence caveat:** both cohorts use curated acquired-gene callers (Oxford AMRFinder / Sci234
  ResFinder-style) — independent of THIS project's pipeline, not of the caller class.
- **Breakpoint [to pin]:** co-trimoxazole TMP-component S≤2 / R≥4 is well-established but should be pinned
  against EUCAST v14.0 / CLSI M100 2024 before any deployed claim.

## Artifacts
- Cells: `wiki/external_validation_sci234tmpsmx_tmpsmx20260616_2026-06-16.json` + `..._oxfordtmpsmx_...json`.
- Roll-up: `wiki/external_validation_report_card.{md,json}` (run-scoped `--run-id tmpsmx20260616`).
- Rule: `dna_decode/data/experimental_drug_rules.py`; scorer: `scripts/tmp_smx_external_validate.py`.
- Tests: `tests/test_experimental_drug_rules.py` (36) + `tests/test_tmp_smx_external_validate.py` (5, incl. frozen-leak guard) + `tests/test_build_external_validation_report.py` (branding passthrough).
