# Supervised learned layer — the regime boundary (R1–R3 synthesis, 2026-07-12)

A user pushback ("I don't think our AI model is useless; the right techniques can optimize it") reopened the
learned-model track. Three recommendations were executed to completion under `--until-mvp`. The honest
result: **supervised learning on the free label IS a live, deployable path — but only where resistance
evolves convergently (viruses), not where it is clonally confounded (bacteria).** The blanket "learned
models 0-for-5" is superseded for the SUPERVISED case; zero-shot likelihood + frozen embeddings stay dead.

## The arc

| step | question | verdict | headline |
|---|---|---|---|
| (prior) | does a supervised full-sequence model beat the HIV catalog on its blind spot? | YES | 0.889 in-dist, **0.81 leave-study-out** (deployable) where ESM zero-shot got 0.449 |
| **R1** | is it general across the HIV RT drug panel, or an EFV quirk? | **GENERAL_RESCUE** | 8/11 pass the deployment split — ALL 5 NNRTIs + 3 higher-cutoff NRTIs; fails only on low-magnitude NRTIs (D4T/DDI/TDF) where the position-based catalog saturates |
| **R2** | integrate — fold into the catalog, or ship a complement? | **COMPLEMENT** | the catalog fold-in was tested + REJECTED (hard accessory rules trade sens for spec, net −0.006 bal-acc); the value is the continuous WEIGHTING → shipped as `hiv_supervised_complement` (offline scorer, 21 KB) |
| **R3** | does the technique transfer to bacteria (TB RIF)? | **IN_DISTRIBUTION_ONLY** | plain 5-fold blind-spot 0.66 but **leave-one-lineage-out 0.51 (chance)** — it learned lineage, not mechanism; the user's caveat CONFIRMED |

## The regime boundary (the load-bearing finding)

**Supervised sequence models rescue the curated catalog's blind spot where resistance is CONVERGENT, and
fail where it is CLONALLY CONFOUNDED.**

- **HIV (works):** a rapidly-evolving virus where the same resistance mutations arise repeatedly across the
  phylogeny, weakly linked to genetic background. A supervised model learns the mutation→resistance mapping
  and it generalizes out-of-study (0.81) and partially out-of-subtype (0.70).
- **TB rifampicin (fails):** a clonal, slowly-evolving bacterium where catalog-negative resistance is
  lineage-correlated (or via non-rpoB efflux/promoter/compensatory mechanisms invisible to point-variant
  features). The model "predicts" the blind spot by recognizing the lineage — which works in-distribution and
  vanishes to chance the moment lineages are held out.

The single methodological lesson that carries everything: **the de-confounded split is the whole ballgame.**
Plain cross-validation said the TB rescue passed (0.66); leave-one-lineage-out revealed it was 0.51. Had we
only run plain CV we'd have shipped a lineage-memorizer as a "bacterial rescue". The HIV win is trustworthy
precisely because it survived the equivalent de-confound (leave-study-out).

## What shipped

- **A deployable learned complement for HIV NNRTI** — `dna_decode/data/hiv_supervised_complement.py` +
  `data/hiv_ref/hiv_nnrti_supervised_complement.json` (offline weighted scorer; `blind_spot_risk`), validated
  general across the drug class (R1) and out-of-distribution (leave-study-out 0.81).
- **A confirmed negative for bacteria** — `wiki/tb_supervised_vs_catalog_2026-07-12.json`: the technique does
  NOT transfer to TB RIF; do not ship a supervised bacterial blind-spot layer on the strength of the HIV win.

## Honest scope

HIV NNRTI-class (trained on EFV); supervised (needs the free Stanford label to train); in-distribution to the
Stanford + CRyPTIC knowledge bases. TB tested on rpoB point-variants only (efflux/promoter/CN resistance is
out of scope and could differ). The frozen decoder surface + `hiv_amr.py` catalog are byte-unchanged
throughout — the complement is a separate additive layer.

## Artifacts / commits
`hiv_supervised_panel` (R1, 4b8e366) · `hiv_catalog_accessory_extension` + `build_hiv_complement_model` +
`hiv_supervised_complement` (R2, a54e3c0) · `tb_supervised_vs_catalog` (R3, ea4c669). Frozen surface
verify_lock OK across all.
