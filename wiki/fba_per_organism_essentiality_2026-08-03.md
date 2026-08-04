# Per-organism FBA essentiality validation — status (2026-08-03)

The FBA cell's engine generalizes across organisms (v0.11.0), but a *validation number* needs a
per-organism experimental essential-gene gold standard whose keys join to that GEM's gene ids. This is the
honest status across the shipped organisms.

| organism | model | status | accuracy | MCC | discrimination | gold standard |
|---|---|---|---|---|---|---|
| E. coli | iML1515 | **SCORED** | **0.954** | **0.652** | strong | Keio RB-TnSeq fitness (Bernstein 2023) |
| S. cerevisiae | iMM904 | **SCORED** | 0.824 | 0.252 | **weak** | SGD `phenotype_data.tab` inviable-null (Giaever/SGD) |
| S. aureus | iYS1720 | **LABEL_WALLED** | — | — | — | model ids are STM#### (Salmonella-style) → needs a gene-name crosswalk |
| P. aeruginosa | iJN1463 | **LABEL_WALLED** | — | — | — | needs a fetchable PAO1 Tn-seq essential set keyed by PA-number |

## The finding (honest, not flattering)

**FBA gene-essentiality accuracy does NOT transfer strongly from E. coli — it is organism/model/medium-
dependent.** E. coli iML1515 discriminates essentiality strongly (MCC 0.652); the same deterministic method
on yeast iMM904 is **weak (MCC 0.252)** despite a superficially-similar accuracy (0.824, flattered by the
15%-prevalence majority class — this is exactly why MCC, not accuracy, is the reported signal). So the
v0.11.0 "cross-organism engine generalizes" claim is now **quantified, not assumed**: the *engine* runs
everywhere; the *predictive quality* is E. coli-strong, yeast-weak, and label-walled for the two bacteria.

## Why yeast is weak + the demotion/improvement path

The yeast metric is on iMM904's **default medium** (WT growth 0.288 /h — low). Essentiality is
medium-dependent, and a mismatched default medium depresses the metric. **Named improvement path (not done —
would be a re-score, not a new claim):** set the standard yeast SD glucose-minimal medium on iMM904 and
re-run; a fair medium likely lifts the number. The default-medium result is the honest v0.

## What would close the two walls (external, not code-closable here)

- **S. aureus:** a `STM#### → S. aureus symbol/locus` crosswalk (the model carries real gene *names* like
  ArgD, so a name-keyed gold standard + a name join would work) — or a different S. aureus GEM whose ids are
  real SAUSA300 locus tags.
- **P. aeruginosa:** a fetchable PAO1/PA14 Tn-seq essential-gene set keyed by PA-number (Jacobs 2003 /
  Lee 2015), which the Chinese/Japanese essentiality servers (DEG/OGEE) host but time out from this host.

## Reproduce

```
uv run python scripts/fba_essentiality_validate.py --organism yeast        # SCORED
uv run python scripts/fba_essentiality_validate.py --organism saureus      # LABEL_WALLED (honest)
```

Harness: `scripts/fba_essentiality_validate.py` + `dna_decode/fba/essentiality_labels.py` (per-organism
sources + pure parsers). Tests: `tests/test_fba_essentiality.py`. Frozen AMR/forward surfaces byte-unchanged.
