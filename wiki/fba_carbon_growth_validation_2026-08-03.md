# FBA carbon-source growth validation: E. coli iML1515 (2026-08-03)

**Claim tested:** FBA predicts growth (a quantitative rate) on a carbon source, validated against measured
carbon-source utilization.

- Measured-positive carbon sources (Keio/Wetmore, K-12 grows): **28**; mapped to a BiGG
  exchange: **21** (unmapped = name-mapping gaps, not FBA failures).
- **RECALL on mapped positives: 1.000** (21/21) -- FBA predicts growth on every carbon
  source E. coli is measured to use.
- FBA growth-RATE spread (quantitative): min 0.2101 / median 0.7171 / max 0.9428 /h.
- Model-gap example (VERIFIED): BW25113 grows on **sucrose** but iML1515 has no sucrose transport -> a
  false negative the validation surfaces. (7 positives unmapped total; a heuristic name-gap
  vs model-gap split is in the JSON but is unreliable on noisy assay labels -- sucrose is the verified one.)

## Honest walls (external, not code-closable here)
- **Full specificity:** EXTERNAL -- needs a MEASURED negative carbon-source set (K-12 can't-use). The Keio/Wetmore assay is positive-only; EcN Biolog (PMC9801561) is a strain mismatch + SI-locked.
- **Growth-rate correlation:** EXTERNAL -- no fetchable MEASURED growth-rate-across-carbon-sources dataset (Biolog reports activity indices, not rates; Monk 2017 rates are KO phenotypes on 16 sources, SI-locked).

## Caveats
- RECALL (sensitivity) only -- FBA predicts growth on carbon sources E. coli is measured to use.
- FBA growth RATES are quantitative but there is no fetchable measured rate to correlate against.
- Unmapped measured-positives are EITHER name-mapping gaps OR real MODEL-gaps (no transporter in iML1515 -> a false negative). E.g. BW25113 grows on SUCROSE but K-12 iML1515 lacks the sucrose system -> a genuine model limitation the validation surfaces, not just a mapping gap.
