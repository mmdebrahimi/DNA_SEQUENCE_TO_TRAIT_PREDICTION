# TB independent gold-set acquisition runbook (deliverable b, ratified E — 2026-06-17)

Goal: a small, CLEAN **independent** TB cohort to validate the frozen TB cell OUTSIDE the WHO-catalogue
build — the honest counterpart to the in-distribution `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`.

## Independence requirement (the real gate)
Independence must be from the **WHO v2 (2023) catalogue BUILD**, not merely "not CRyPTIC". WHO v2 swept
most public **pre-2023** TB WGS + pDST, so legacy cohorts (TB Portals / Walker / ReSeqTB) are largely
IN-distribution for the rule. Real independence = **post-2023 isolates** (temporal hold-out).

## Ratified route (E): hand-curate ~30 FIRST; DAR in parallel
1. **Hand-curated ~30-isolate post-2023 gold set (FIRST — the first independent number).**
   - Inclusion: collected/sequenced **after the WHO v2 cut**, public WGS (ENA/SRA), paired **reference**
     DST (broth microdilution preferred) for RIF and/or INH, confirmed NOT in the WHO v2 build.
   - Sources to mine: post-2023 national-surveillance BioProjects; recent CRyPTIC-independent papers with
     deposited reads + DST; WHO-independent reference-lab panels.
   - Per isolate: fetch reads → masked + regeno VCF vs H37Rv NC_000962.3 (mirror the CRyPTIC shapes the
     cell consumes), OR an assembly → VCF. Record the measured DST label.
   - Emit the manifest `tb_goldset` consumes (JSON list of
     `{strain_id, masked_vcf, regeno_vcf, label}`), then run
     `scripts/score_tb_independent_goldset.py` → `INDEPENDENT_VALIDATION` (small N → read the Wilson CI).
   - **Feasibility is [unverified] until attempted.** If 30 clean post-2023 not-in-WHO isolates can't be
     reached, fall through to the DAR.
2. **TB Portals / NIAID DAR (PARALLEL — the larger confirmatory follow-up).**
   - Submit the Data Access Request now (it is slow). Caveats: access-gated; **mixed-method DST** (filter
     to reference BMD); partial WHO-build overlap (filter to post-2023). Larger N once approved.
   - **Never block deliverable (b)'s first number on the DAR.**

## Honest status until the gold set lands
`scripts/score_tb_independent_goldset.py` emits `INDEPENDENT_VALIDATION_BLOCKED_NO_GOLDSET` — an honest
"no independent number yet", never a fabricated metric. The CRyPTIC baseline is reported separately and
is never relabelled as independent.
