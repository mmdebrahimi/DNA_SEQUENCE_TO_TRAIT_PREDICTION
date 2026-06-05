# Brainstorm — cef embedding-vs-classical falsifier (2026-06-04)

Subject: will `scripts/cef_falsifier.py`'s verdict be trustworthy? **Bottom line: NO — the dominant
problem is the SUBSTRATE, not the code.**

## Critical issues
1. **The cef cohort is geography/lineage confounded** (the study==class trap that killed pathotype,
   re-entering via cohort construction). `[grounded — verified on gate_b_cohort]`
   - cef R=26 strains / 26 MLSTs; S=34 / 34 MLSTs; **only 1 shared lineage** (1/26 R, 1/34 S).
   - Country: **R = USA 13 + blank 12**; **S = India 10 / Kenya 9 / Mali 6.** R≈USA, S≈Africa/S-Asia.
   - ⇒ an NT classifier here can predict continent/sequencing-batch, not β-lactamase mechanism. The
     prior AUROC 0.895 is suspect; not citable as "NT beats classical." MIC is a sampling-independent
     label TYPE, but THIS cohort drew R and S from disjoint sources.
   - Contrast: cipro N=147 has 6 shared R/S lineages → genuinely de-confounded.
2. **Gate not statistically meaningful at N≈59:** verdict on point estimate (`gap≥3pp`), bootstrap CI
   computed but unused; max-of-2 NT heads selected before bootstrap (optimism bias); 3pp inside the LOSO
   noise floor.
3. **Comparator too weak:** k-mer only; AMRFinder gene-presence (the strong baseline for plasmid
   β-lactamases) deferred. Beating k-mer ≠ beating best classical.

## Medium
- `EmbeddingCache.verify_complete` not called in cef_falsifier.py or stage1; no accession-uniqueness assert.

## What's solid
- Falsifier mechanics correct (stage1 reuse, leave_one_accession_out, bootstrap, clean parking).
- GPU-only transfer + cache validation sound. **Cipro N=147 remains a clean substrate.**

## Verdict (convergent)
The cef substrate is **not salvageable** for an honest verdict as constructed. Mash-clade-out /
leave-one-country-out are diagnostics that would only confirm collapse — not a rescue (1 shared MLST →
group-out degenerates). No honest contrast extractable from this exact cohort. **Block any promotable
cef verdict; rebuild a geography/lineage-balanced cef cohort** (R+S co-occurring within multiple
clades/countries) before a cef embedding verdict. Cipro is the substrate meanwhile.

## Recommended next steps
1. Block the cef verdict (run stopped; stamp non-promotable). Don't cite 0.895.
2. Make the lineage/geography de-risk a PRECONDITION gate for every cohort (the check run for cipro must
   run before building a falsifier — it was an afterthought for cef).
3. Rebuild a de-confounded cef cohort (needs cef-S USA / cef-R Africa-India strains; BV-BRC/NARMS/EuSCAPE).
4. Fold in the 3 fixes (CI-aware gate, AMRFinder gene-presence baseline, verify_complete) on the next clean run.
