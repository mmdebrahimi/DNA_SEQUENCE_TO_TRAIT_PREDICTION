# dna-amr multi-drug validation — cef + tet + gent — 2026-06-06

> Option-C ("double down on what works"): extend the cipro-validated deterministic AMR caller to
> ceftriaxone + tetracycline + gentamicin using cached AMRFinder runs, with per-drug calibrated rules
> baked into `dna_decode/eval/amr_rules.py::DRUG_RULE`. Mechanism features, NOT embeddings (per
> `plans/AMR_embedding_niche_decision_2026-06-05.md`). **All 4 drugs now validated.**

## Result (cached BV-BRC cohorts, no Docker — `evaluate_cohort`)

| Drug | Cohort (N) | Rule | acc | sens | spec |
|---|---|---|---|---|---|
| ciprofloxacin | stage2_n150 (147) | threshold=2, broad class (QRDR point-mut) | 0.939 | 0.931 | 0.947 |
| **ceftriaxone** | gate_b (60: 26R/34S) | **threshold=1 + extended-spectrum subclass** | **0.933** | **0.962** | **0.912** |
| ceftriaxone | gate_b_mini_cef (12) | same | 1.000 | 1.000 | 1.000 |
| **tetracycline** | gate_b_mini_tet (12: 6R/6S) | threshold=1, class (acquired genes) | **0.833** | **1.000** | **0.667** |
| **gentamicin** | pooled (128: 28R/100S) | **threshold=1 + GENTAMICIN-subclass** | **0.945** | **0.893** | **0.960** |

The cef + gent results share ONE mechanism: a broad AMR class (BETA-LACTAM / AMINOGLYCOSIDE) over-calls
because it counts genes conferring resistance to OTHER members of the class. AMRFinder's **Subclass**
field is the drug-specific discriminator — refine to the drug's own subclass token (CEPHALOSPORIN/
CARBAPENEM for ceftriaxone; GENTAMICIN for gentamicin). Same one-line fix, +0.4 spec each.

## The ceftriaxone fix (the load-bearing finding)

The naive "any drug-class determinant" rule failed cef badly:

| cef rule | acc | sens | spec |
|---|---|---|---|
| broad BETA-LACTAM, threshold=1 | 0.65 | 0.962 | **0.412** (over-calls) |
| broad BETA-LACTAM, threshold=2 | 0.817 | 0.654 (under-calls) | 0.941 |
| **extended-spectrum subclass, threshold=1** | **0.933** | 0.962 | 0.912 |

**Root cause:** E. coli broadly carries intrinsic/narrow-spectrum β-lactamases (blaTEM-1, blaEC) that
confer ampicillin resistance but **not** 3rd-gen-cephalosporin resistance. AMRFinder encodes this in the
**Subclass** field: extended-spectrum determinants (ESBL blaCTX-M, AmpC blaCMY, carbapenemases blaNDM/
blaKPC) carry Subclass `CEPHALOSPORIN`/`CARBAPENEM`; narrow ones carry plain `BETA-LACTAM`. The
discriminator was clean on N=60:

| Subclass | S strains | R strains |
|---|---|---|
| plain BETA-LACTAM | 15 | 23 | ← noise (excluded) |
| CEPHALOSPORIN | 3 | 37 | ← kept |
| CARBAPENEM | 0 | 9 | ← kept |

Fix: for ceftriaxone, count only determinants whose Subclass ∋ {CEPHALOSPORIN, CARBAPENEM}
(`DRUG_RULE["ceftriaxone"]["subclass_any"]`).

## Per-drug rule config (baked in, `amr_rules.py::DRUG_RULE`)

- **ciprofloxacin** — threshold 2, no refinement (QRDR point-mutations need ≥2 hits). Unchanged.
- **ceftriaxone** — threshold 1 + extended-spectrum subclass refinement.
- **tetracycline** — threshold 1, no refinement (acquired tet genes; N=12, small — provisional).
- **gentamicin** — threshold 1 + GENTAMICIN-subclass refinement (excludes aph/aadA streptomycin-kanamycin
  genes that don't confer gentamicin-R) — N=128 acc 0.945.

`call_resistance(main_tsv, drug)` now auto-selects the per-drug threshold + refinement; pass an explicit
`resistance_threshold` to override.

## Honest caveats

- **cef sens 0.962 / spec 0.912** — 1 FN (an R strain with no detected extended-spectrum determinant —
  likely porin-loss/efflux-only resistance, which AMRFinder curated determinants don't capture) + 3 FP
  (S strains carrying a CEPHALOSPORIN-subclass gene but phenotypically susceptible — possible
  low-expression or MIC near the breakpoint).
- **tet N=12 is small** — sens 1.0 but spec 0.667 (2 FP). Provisional; a larger tet cohort should
  re-validate before any non-internal use.
- Same-source caveat as the cipro external validation: all cohorts are BV-BRC broth-microdilution;
  not a cross-lab check.
- gentamicin is unvalidated — the config entry is a mechanism-analogy default, flagged as such in output.

## Reproduce

```
uv run python -c "from dna_decode.data.cohort import load_cohort; from dna_decode.eval.amr_rules import evaluate_cohort; \
c=load_cohort('data/processed/gate_b_cohort.parquet'); \
pairs=[(s.assembly_accession,int(s.ast_labels[[k for k in s.ast_labels if k.lower()=='ceftriaxone'][0]])) for s in c.strains if any(k.lower()=='ceftriaxone' for k in s.ast_labels)]; \
print(evaluate_cohort('data/amrfinder_runs',pairs,'ceftriaxone'))"
```

Regression pinned in `tests/test_amr_rules.py` (`test_cef_cohort_opchars_regression` + the cef/tet
synthetic-determinant tests). 13 tests total.
