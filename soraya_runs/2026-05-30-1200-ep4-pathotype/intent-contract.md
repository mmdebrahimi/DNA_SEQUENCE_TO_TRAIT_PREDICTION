<!-- intent-contract.md for 2026-05-30-1200-ep4-pathotype -->

# Soraya Run — Intent Contract

- **run_id:** 2026-05-30-1200-ep4-pathotype
- **mode:** --advance (self-imposed cap: --max-actions 3; money-only gate)
- **family:** EP-4 pathotype (highest-VOI accepted family; AMR family v0-shipped)
- **goal served:** advance the bounded ExPEC+EPEC(+ETEC) slice toward a working honestly-scoped predictor
- **environment constraints:**
  - cwd = rca_engine/articles → `/project-state` REFUSEs (path-gated); ledger recording falls back to direct Edit (surfaced)
  - EP-4 caller execution (VirulenceFinder/blastn/GPU) = workhorse-only; only analysis/data-prep doable here
  - sync = commit tracked dirs on main (data/ gitignored)
- **planned batch (VOI-ranked):**
  1. [auto] Resolve von Mentzer 2021 ETEC 7-lineage reference-genome INSDC accessions (recovers ETEC by-accession; Horesh ETEC = all Sanger lane IDs). Slice 2-class → 3-class.
  2. [auto] Write ETEC reference accessions to research_outputs/ (tracked → syncs).
  3. [irreversible/un-gated] commit + push.
- **emit (user-only, parked, non-blocking):** Gate B outreach send.
- **stop:** cap hit / all-blocked / money-declined / no-beneficial-action.
