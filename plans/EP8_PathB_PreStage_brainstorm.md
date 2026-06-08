# Brainstorm review — Path B / Gate G2 manifest (2026-06-08)

Subject: `plans/EP8_PathB_PreStage_Manifest.md` (the Arabidopsis flowering-time embedding test, G2).
Outcome: manifest REVISED (see its post-brainstorm header). This is the review record.

## Critical issues (fixed in the revision)
1. **[grounded] Primary curated-locus approach contradicted the niche criterion — decisive.** EP8 defines the
   niche as "No curated mechanism catalog"; the manifest's primary test hand-picked FLC/FRI/FT, injecting the
   exact domain knowledge the niche claim avoids → it tested "can the FM represent these 8 loci", not the
   thesis. **Fix:** primary → phenotype-AGNOSTIC subsample (window selection frozen, no FT labels, no
   flowering-gene enrichment); curated loci demoted to a narrowed secondary diagnostic.
2. **[inferred] Full genome-wide-primary is NOT the fix** — operationally hostile on 12 GB (~260k–520k
   windows/accession × ~1122) AND a naive global mean-pool dilutes sparse causal signal. The agnostic
   *subsample* (bounded window budget) is the least-compromised path; also keep a variance/block-pool so a
   null reflects representation not averaging.
3. **[grounded] No continuous within-lineage diagnostic exists** — `within_lineage_diagnostic.py` is binary
   within-MLST R-vs-S concordance. **Fix:** predeclared continuous within-group test (center within group,
   within-group R², min-N gate).
4. **[inferred] PASS was point-estimate-only** (+0.05 R² across ~9 unequal groups = brittle). **Fix:** PASS
   requires a paired-bootstrap/grouped-permutation CI on (embedding − best baseline) R² that EXCLUDES 0.
5. **[grounded] Not execution-ready** — no accession list / pseudogenome manifest / annotation / coord table /
   join validation / N-fraction QC / matched indel baseline. **Fix:** mandatory CPU-only §0.5 G2 dry-manifest
   gate before any GPU.

## Medium (fixed)
- **[inferred]** Baseline must be information-matched — SNPs + short indels + missingness over the SAME windows
  as the embedding (else SNP-only is unfairly weak vs a pseudogenome embedding seeing indels + N).
- **[inferred]** Do NOT residualize phenotype vs geography in the gate (latitude is partly causal for flowering
  via vernalization — removes real biology). Kinship/PCs = primary de-confounder; geography = sensitivity only.
  *(This corrected an interim proposal made during the review.)*
- **[inferred]** Freeze FT10 as sole primary endpoint, FT16 replication only.

## Deferred
- **[speculative]** FRI/FLC structural variants may be under-represented by pseudogenomes — check in the
  dry-manifest, non-blocking.
- Stronger pooling (attention/multi-scale) only if pursuing the paid-cloud full-genome route.

## What's solid (unchanged)
- Substrate choice (sampling-independent quantitative label + organism depth 1122, verified).
- Leave-one-genetic-group-out CV as the structure stress-test.
- Iron-law data pinning (phenotype committed; genotype VCF live 19.2 GB; FM fits 12 GB).
- Keep-frontier-open-on-clean-FAIL pre-commit.

## Open tradeoffs → now manifest §8 (user decides at the dry-manifest)
Primary estimand (portability vs causal-signal vs prediction); agnostic-window selection rule; GPU-hours/window
budget; whether a true-genome-wide cloud run is ever authorized (money gate).
