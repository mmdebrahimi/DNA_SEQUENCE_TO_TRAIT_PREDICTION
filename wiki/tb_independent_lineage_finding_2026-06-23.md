# TB independent number — lineage-collapse is the WRONG headline (homoplasy finding, 2026-06-23)

Built the lineage-collapse post-processor for the independent TB number (`scripts/tb_independent_lineage_collapse.py`, reusing the FROZEN `score_tb_cryptic.run_v1b` → `clonality.cluster_weighted_confusion`
+ the Napier barcode, NO Docker). Ran it on the partial cohort (562 isolates with VCFs so far). **Verify-in-batch
surfaced a genuine methodological finding: lineage-MAJORITY-collapse is inappropriate for TB AMR, and the RAW
independent number is the honest headline (with a clonality DISCLOSURE, not a demotion).**

## What ran (partial N=562)
| Drug | RAW sens / spec (TP/FP/TN/FN) | lineage-MAJORITY-collapse | lineage assign |
|---|---|---|---|
| Rifampicin | **0.885 / 0.961** (131/16/397/17) | sens **None** (R-clusters=0, S-clusters=33, discordant=9) | 9/561 unassigned |
| Isoniazid | **0.803 / 0.989** (159/4/360/39) | sens **0.0** (R-clusters=5, S-clusters=25, discordant=12) | 9/562 unassigned |

## The finding: TB resistance is HOMOPLASIC → lineage-majority-collapse erases it
Lineage assignment WORKS (≈99% assigned). But `cluster_weighted_confusion` collapses each lineage cluster to
ONE majority-label vote — correct for a CLONAL trait, **wrong for AMR.** In TB, resistance is acquired
INDEPENDENTLY many times within the same sublineage (homoplasy): a cluster is mostly-S with a few
independently-R isolates → it votes S, and the R minority is out-voted or flagged "discordant" → 0 R-clusters
→ sens collapses to None/0. This is the SAME phenomenon the in-distribution CRyPTIC baseline already showed
(lineage RIF 0.41 / INH 0.349 ≪ raw 0.916), here pushed to the degenerate extreme by (a) partial/clustered N
and (b) AMR's homoplasy.

**This is a methodological correction, not a result:** a lineage-MAJORITY collapse measures "does the rule
work on the lineage's TYPICAL phenotype," which for a homoplasic, within-lineage-acquired trait is the wrong
question. The rule operates per-ISOLATE on the isolate's OWN determinants; that is what the RAW number
measures, and it is the honest independent number for TB AMR.

## The honest reporting decision (mirrors the project's bacterial lineage-DISCLOSURE)
The project's bacterial report card already chose **DISCLOSE, do NOT hard-dedup-and-demote** (the
lineage-disclosure layer augments, never demotes, the per-isolate number). Apply the SAME principle to TB:
- **Headline = the RAW independent sens/spec** (per-isolate, the genuine genotype→independent-phenotype test).
- **Disclosure (qualitative, not a competing number):** TB R classes are clonally structured; resistance is
  homoplasic; a lineage-majority collapse is inappropriate and degenerates to ~0 by construction. Report the
  effective-lineage spread + discordant count as DISCLOSURE, not as a demoted sens.
- Do NOT publish the lineage-collapsed sens as "the number" — that would be a methodological error.

## Status of the deliverables
- **Post-processor BUILT + wired + runs** (`scripts/tb_independent_lineage_collapse.py`) — faithful reuse of
  the frozen path; useful for the DISCLOSURE counts (effective-lineage-N, discordant), NOT as the headline.
- **Headline = the RAW number from `scripts/run_tb_independent_amr_portal.py`** (the full-cohort run is in
  flight; partial N=562 raw: RIF 0.885/0.961, INH 0.803/0.989 — same range as the in-distribution baseline).
- Frozen surface byte-unchanged; no Docker used by this post-processor.

## Provenance
`scripts/tb_independent_lineage_collapse.py` (reuses `score_tb_cryptic.run_v1b` +
`clonality.cluster_weighted_confusion` + `tb_lineage` + the pinned barcode). Partial result
`wiki/tb_independent_amr_portal_lineage_collapsed.json`. Raw headline memo
`wiki/tb_independent_number_2026-06-23.md`. Cross-ref the in-distribution baseline's identical lineage
behavior: `wiki/tb_cryptic_parquet_baseline_2026-06-22.md`.
