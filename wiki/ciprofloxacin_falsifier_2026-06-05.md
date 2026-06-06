# ciprofloxacin embedding-vs-classical falsifier (2026-06-05)

**De-confound gate:** `[DE_CONFOUNDED·promotable] primary lineages: 6 shared, matched minority 10 (R-in-shared 35, S-in-shared 10); country: matched=22/16g; year: matched=59/13g; 6 shared lineages, matched minority 10; no secondary axis aliases the label. (Screen passed — necessary, not sufficient.)`
**Cohort:** `data\processed\stage2_n150_cipro_cohort.parquet` (N=140; 67R/73S) · pooling mean · CV leave_one_accession_out
**Best NT:** NT-XGBoost 0.914 · **k-mer-XGB:** 0.943 · **gap -2.9pp**
**95% bootstrap CI on gap:** [-9.0, +2.9]pp (eff 1000/1000; paired on 140/140 strains)
**VERDICT:** FAIL (gap -2.9pp < 3)

| Variant | AUROC | AUPRC |
|---|---:|---:|
| NT-XGBoost | 0.914 | 0.888 |
| NT-logreg | 0.863 | 0.828 |
| POINT-XGB | 0.943 | 0.962 |

## Notes
- De-confound gate is a PRECONDITION (CONFOUNDED cohort → blocked, no verdict).
- CI-aware verdict: point gap >= 3pp AND bootstrap CI lower bound > 0 for a PASS (else NOISY).
- 'best classical' comparator = **POINT-XGB** (0.943); POINT-XGB present = QRDR/plasmid KNOWLEDGE baseline included (the real bar).
- Single drug + single cohort ⇒ NOT an architecture-class promotion regardless of verdict.
- per-strain scores persisted to the .scores.json sidecar (crash-recovery + within-lineage diagnostics).
- calibrate=False; verify_complete cache integrity = follow-up.