# P. aeruginosa: the model wall IS breakable — but the validation is WEAK (2026-08-07)

Follow-on to `wiki/fba_wrong_organism_model_bug_2026-08-07.md`, which reclassified P. aeruginosa from
LABEL_WALLED to **MODEL_WALLED** (BiGG has no P. aeruginosa reconstruction). This memo resolves whether
that wall is breakable, and measures what breaking it would actually buy. **Measured, not asserted —
and deliberately NOT wired into the shipped registry (see the recommendation).**

## Both halves of the wall are now sourced and verified fetchable

| piece | source | verified |
|---|---|---|
| **model** | **iSD1509** (Dahal 2022/2023), BioModels `MODEL2205090001`, `iSD1509.omex` (754 KB COMBINE archive) | downloaded + parsed from this host; **1509 genes keyed `PA14_#####`**, 2021 reactions, WT growth 1.5909 /h |
| **label** | **GOLD_84 / GOLD_115** (Van Maele 2026, PLOS Comput Biol `pcbi.1013945.s011`) | downloaded + parsed; Tn-seq curated essential genes, `Old Locus tag` column is `PA14_#####` |

**The join is DIRECT — no crosswalk.** This is the join the 2026-08-03 label-wall memo wrongly claimed
existed for `iJN1463` (which is *P. putida*, and PAO1-vs-PA14-keyed besides). Overlap:

- GOLD_84: **38 / 84** gold genes present in the model (45%)
- GOLD_115: **53 / 115** gold genes present in the model (46%)

The non-overlap is expected and honest: the gold sets are whole-cell essential genes (ribosome,
replication, secretion), while the model only contains 1509 *metabolic* genes.

## The measured result — WEAK, and inverted by functional class

Single-gene deletion over all 1510 model genes, essential = growth < 1% of wild-type.
**134 / 1510 (8.9%) of model genes are FBA-essential.**

| gold set | recall (in-model subset) |
|---|---|
| GOLD_84 | **24 / 38 = 0.632** |
| GOLD_115 | **30 / 53 = 0.566** |

Split by PseudoCAP functional class (GOLD_115), the pattern is the **opposite** of the obvious hypothesis:

| functional class | recall |
|---|---|
| Cell structure and division | 14/17 = **0.824** |
| Uncharacterized molecule | 3/3 = 1.000 |
| Macromolecular synthesis | 7/15 = 0.467 |
| **Metabolic pathway** | **6/17 = 0.353** |
| Adaptation, Protection | 0/1 = 0.000 |

FBA does **worst on the metabolic class — the one class it should own.** The expectation going in was
the reverse (FBA works on metabolism, fails on cell division); the data falsified that.

**The medium confound was checked and ELIMINATED.** GOLD_115 is defined on PA14 grown in **LB (rich)**,
so a minimal-medium model would have been an unfair test. `iSD1509`'s default medium is *already rich* —
36 open exchanges including **18 amino acids** plus glucose and nucleosides. Medium and label condition
agree, so the weak metabolic recall is a genuine model-vs-reality gap, not a setup artifact. The likely
mechanism: on a rich medium the model can bypass biosynthetic genes via supplied amino acids, while the
wet-lab Tn-seq still finds those genes essential in LB.

## What CANNOT be computed (and why no accuracy/MCC appears above)

**Specificity, precision, accuracy and MCC are NOT computable on this label.** GOLD_84/GOLD_115 are
curated *reference subsets* the paper built to benchmark Tn-seq analysis methods — they are **not
exhaustive essential-gene lists**. So "not in the gold set" does **not** mean non-essential, and every
FBA-essential-but-not-in-gold call would be scored as a false positive on an unsound premise. Only
**RECALL** is sound here. This is the same shape as the carbon-growth validation (recall 1.0,
specificity externally walled) — see `wiki/fba_carbon_growth_validation_2026-08-03.md`.

## Recommendation: do NOT wire iSD1509 into the shipped registry yet

Deliberate, and the reason matters: this session's whole subject was an **over-claimed cross-organism
capability**. Shipping `--organism paeruginosa` on the back of a 0.57 recall with no computable
specificity would recreate that error in a new place. Preconditions to ship it:

1. A **non-BiGG model source** in `model.resolve_model_path` (BioModels COMBINE/`.omex` fetch + unpack) —
   currently the loader only knows the BiGG URL scheme. Real work, not a config line.
2. An **exhaustive** P. aeruginosa essentiality label (a full Tn-seq call set, not a curated subset) so
   specificity becomes computable. Then the cell can be honestly SCORED rather than recall-only.
3. A decision on whether recall ~0.6 with no specificity clears the bar to expose as a decoder cell.
   That is an acceptance-bar call, not a technical one.

Until then P. aeruginosa stays **MODEL_WALLED** in `essentiality_labels.MODEL_WALLED`, which is accurate:
the shipped tool has no model for it.

## The finding that generalizes

This is the **third** organism where FBA essentiality is weak (E. coli strong MCC 0.652; yeast weak
MCC 0.252; P. aeruginosa recall 0.57 with the metabolic class at 0.35). The "FBA essentiality does not
transfer strongly from E. coli" conclusion in `wiki/fba_per_organism_essentiality_2026-08-03.md` is
**reinforced by an independent organism, on a matched medium** — and it now has a sharper edge: the
failure concentrates in exactly the class FBA is supposed to be best at, on rich media.

## Reproduce

```bash
# model  (BioModels COMBINE archive -> model.xml inside)
curl -L "https://www.ebi.ac.uk/biomodels/model/download/MODEL2205090001?filename=iSD1509.omex" -o iSD1509.omex
# label  (PLOS supplementary; sheets GOLD_84 + GOLD_115, column "Old Locus tag")
curl -L "https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1013945.s011&type=supplementary" -o gold.xlsx
```
Neither is committed (regenerable, and the model is a third-party artifact under its own terms).
