# DNA Decoder v0.1 cef cached-strain handoff - 2026-05-25

## Goal

Take the next smallest product step after v0.1 genome-input ciprofloxacin:

- expand to a second drug
- keep the simpler cached-strain surface first
- avoid changing both drug and input contract at the same time

## Overnight sequence executed

1. Verified existing cef substrate in repo.
   - working cohort: `data/processed/gate_b_cohort.parquet`
   - cef pool inside cohort: `50` strains
   - label balance in cef pool: `26R / 24S`

2. Checked current NT cache coverage against the cef substrate.
   - existing cipro cache covered only `35 / 67` gate-B cohort strains
   - existing cipro cache covered only `5 / 12` mini-cef strains

3. Downloaded / verified RefSeq packages for the full gate-B cohort.
   - target root:
     `C:\Users\b0652085\OneDrive - Bombardier\Apps\Stress-DNA Project\dna_decode_cache\refseq`
   - result: `67 / 67` strains downloaded or already present

4. Built a dedicated NT embedding cache for the gate-B cohort.
   - cache path:
     `C:\Users\b0652085\OneDrive - Bombardier\Apps\Stress-DNA Project\dna_decode_cache\embeddings\nt_gate_b_cohort_67.h5`
   - local model override used:
     `C:\Users\b0652085\OneDrive - Bombardier\Apps\Stress-DNA Project\local_models\nucleotide_transformer`
   - populate result:
     `312,585` new embeddings across `67` strains

5. Probed cache completeness.
   - `64` strains complete
   - `3` strains absent with `expected=0`, not partial corruption:
     - `1045010.61`
     - `562.59220`
     - `562.59226`
   - interpretation:
     these three strains expose zero extractable CDS rows under the current parser path, so the trainer naturally skips them

6. Trained a real cef cached-strain NT model.
   - command path:
     `python -m scripts.pipeline train --drug ceftriaxone --model nucleotide_transformer ...`
   - model output:
     `data/processed/models/ceftriaxone_nucleotide_transformer.pkl`
   - effective train/eval set:
     `49` usable strains
   - class balance:
     `25R / 24S`
   - CV strategy:
     `loso`
   - AUROC:
     `0.895`
   - AUPRC:
     `0.838`

7. Generated two real cef cached-strain prediction artifacts.
   - resistant example:
     - sample: `562.12960`
     - output:
       - `reports/dna_decoder_v0_1_cef_cached_example_R_2026-05-25.json`
       - `reports/dna_decoder_v0_1_cef_cached_example_R_2026-05-25.md`
     - result:
       - prediction `R`
       - probability `0.753`
       - confidence `MEDIUM`
   - susceptible example:
     - sample: `562.7572`
     - output:
       - `reports/dna_decoder_v0_1_cef_cached_example_S_2026-05-25.json`
       - `reports/dna_decoder_v0_1_cef_cached_example_S_2026-05-25.md`
     - result:
       - prediction `S`
       - probability `0.204`
       - confidence `MEDIUM`

## What this means

- cef cached-strain expansion is now empirically viable on this machine
- this is no longer just a design recommendation
- we now have:
  - cef substrate
  - cef cache
  - cef trained NT model
  - real cef prediction artifacts

## Important caveats

1. Cef examples are currently debug-mode, not canonical audit-mode.
   - reason: no cef audit merge sidecar equivalent to the cipro merge packet is wired into the current decoder flow
   - current cef examples therefore use:
     `--allow-missing-audit --no-attribution`

2. The current cef CV uses `strain_id` grouping.
   - no duplicate-accession mitigation issue surfaced during this run
   - but the same duplicate-audit discipline used for cipro should still be applied if cef becomes a promoted slice

3. Interpretability for cef is still exploratory.
   - tonight's work established predictive viability for cached-strain cef
   - it did not build a cef audit / interpretability closeout packet

## Recommended next move

Promote this into a small formal v0.1 cef cached-strain packet before widening scope again.

Recommended order:

1. run duplicate-accession audit on the cef pool
2. write a short cef cached-strain release candidate note
3. decide whether the next axis is:
   - cef genome-input
   - cef audit-aware packet
   - or a second-drug validation panel

## Suggested prompt for the other machine

```text
We now have a real cef cached-strain substrate on Precision 7780: dedicated NT cache built, ceftriaxone_nucleotide_transformer.pkl trained at AUROC 0.895 on 49 usable strains, and two real cef prediction artifacts. Please analyze the repo and recommend the smallest credible next slice to promote cef from internal viability to a formal v0.1 product packet, including whether duplicate-accession audit, audit-sidecar design, or genome-input cef should come next.
```

## Key files

- `data/processed/gate_b_cohort.parquet`
- `data/processed/models/ceftriaxone_nucleotide_transformer.pkl`
- `reports/dna_decoder_v0_1_cef_cached_example_R_2026-05-25.md`
- `reports/dna_decoder_v0_1_cef_cached_example_S_2026-05-25.md`
- `reports/dna_decoder_v0_1_genome_input_release_candidate_2026-05-25.md`
- `reports/dna_decoder_v0_1_parallel_handoff_2026-05-25.md`
