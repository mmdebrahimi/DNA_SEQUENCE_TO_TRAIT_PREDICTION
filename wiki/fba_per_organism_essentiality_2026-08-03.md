# Per-organism FBA essentiality validation — status (2026-08-03)

The FBA cell's engine generalizes across organisms (v0.11.0), but a *validation number* needs a
per-organism experimental essential-gene gold standard whose keys join to that GEM's gene ids. This is the
honest status across the shipped organisms.

| organism | model | status | accuracy | MCC | discrimination | gold standard |
|---|---|---|---|---|---|---|
| E. coli | iML1515 | **SCORED** | **0.954** | **0.652** | strong | Keio RB-TnSeq fitness (Bernstein 2023) |
| S. cerevisiae | iMM904 | **SCORED** | 0.824 | 0.252 | **weak** | SGD `phenotype_data.tab` inviable-null (Giaever/SGD) |
| ~~S. aureus~~ | ~~iYS1720~~ | ~~LABEL_WALLED~~ | — | — | — | ⚠️ **WRONG — iYS1720 is a *Salmonella* pan-reactome.** See correction below |
| ~~P. aeruginosa~~ | ~~iJN1463~~ | ~~LABEL_WALLED~~ | — | — | — | ⚠️ **WRONG — iJN1463 is *P. putida* KT2440.** See correction below |

> ## ⚠️ CORRECTED 2026-08-07
>
> The two walled rows above were **misdiagnosed**. Verified against the BiGG Models API: `iYS1720` is a
> ***Salmonella* pan-reactome** (1262/1707 gene ids carry the *S.* Typhimurium `STM` prefix — the very
> detail this artifact read as an "ID convention" quirk) and `iJN1463` is ***Pseudomonas putida* KT2440**.
> Neither cell was label-walled: **no gold standard could ever have joined, because the model was the
> wrong organism.** Corrected status (v0.12.1):
>
> | organism | model | status | why |
> |---|---|---|---|
> | S. aureus | **iYS854** (USA300_TCH1516) | LABEL_WALLED | ids are `USA300HOU_####`; NTML is USA300 JE2 (`SAUSA300_####`) → still a crosswalk |
> | Salmonella | iYS1720 | LABEL_WALLED | ids are `STM####`; no keyed gold standard wired |
> | P. putida | iJN1463 | LABEL_WALLED | no keyed essentiality gold standard wired |
> | P. aeruginosa | **none** | **MODEL_WALLED** | BiGG has no P. aeruginosa reconstruction; the alias now refuses |
>
> The **E. coli and yeast rows are UNAFFECTED** — both models were correctly assigned, so both scored
> numbers and the "does not transfer strongly" finding stand. Full evidence:
> `wiki/fba_wrong_organism_model_bug_2026-08-07.md`.

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

> **⚠️ SUPERSEDED 2026-08-07 by the correction above.** The "different S. aureus GEM whose ids are real
> locus tags" was the right instinct — that GEM is **iYS854**, and v0.12.1 now uses it. Updated walls:

- **S. aureus (iYS854):** a `USA300HOU_#### → SAUSA300_####` ortholog crosswalk, to join the NTML /
  Nebraska transposon library (USA300 JE2). Near-1:1, but a crosswalk nonetheless.
- **P. aeruginosa:** needs a **genome-scale model**, not a label. The gold standard already exists and
  IS fetchable from this host — PLOS Comput Biol 2026 `pcbi.1013945.s011`, sheets `GOLD_84` / `GOLD_115`
  (downloaded + parsed 2026-08-07; the earlier bot-block was transient) — but it is **PA14**-keyed and
  BiGG has no P. aeruginosa reconstruction to join it to.

## Reproduce

```
uv run python scripts/fba_essentiality_validate.py --organism yeast        # SCORED
uv run python scripts/fba_essentiality_validate.py --organism saureus      # LABEL_WALLED (honest)
```

Harness: `scripts/fba_essentiality_validate.py` + `dna_decode/fba/essentiality_labels.py` (per-organism
sources + pure parsers). Tests: `tests/test_fba_essentiality.py`. Frozen AMR/forward surfaces byte-unchanged.
