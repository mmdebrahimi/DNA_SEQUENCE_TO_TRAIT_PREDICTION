# HANDOFF — EP-4 Pathotype Learned-Model Track — 2026-05-31

> Fresh-context entry point. Read this first. Covers the pivot to the learned-model track, what's built/committed, the in-flight k-mer run, environment gotchas, and the exact next move.

## ⚠️ Session/environment facts (read before doing anything)

- **cwd is WRONG.** Sessions launch in `C:\Users\Farshad\PythonProjects\rca_engine\articles`, NOT the project. `/project-state` path-gates and REFUSEs; ledger edits fall back to direct `Edit`. **Fix for a clean session:** launch from `C:\Users\Farshad\PythonProjects\dna_decode` (do NOT "resume" the old wrong-cwd session — resume restores the original cwd). All repo ops here use `git -C <dna_decode>` / absolute paths.
- **D: external drive is FLAKY / currently DISCONNECTED** (`WinError 21: device not ready`). It's the Seagate Portable. It drops out and orphans/loses work. **All this session's caches were re-routed to C:** (`data/ena_wgs/`, gitignored). C: has ~28 GB free; D: when mounted has 4.4 TB. Do NOT assume D: is available.
- **Sync channel = git on `main`** (tracked dirs). `data/` and `reports/` are gitignored → never sync via git; re-derive from public accessions. Remote: `github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`.
- **venv:** `.venv/Scripts/python.exe` (uv). torch 2.7.1+cu118, CUDA works on GTX 860M. `uv pip install` for deps (no plain pip). openpyxl was added this session.
- **Docker Desktop daemon** must be running for Bakta (start `C:\Program Files\Docker\Docker\Docker Desktop.exe`; ~30s to ready). Bakta image `oschwengers/bakta:v1.11.4` pulled. Bakta DB light at `C:/Users/Farshad/dna_decode_stage2/bakta_db/db-light`.

## The pivot (what changed this session)

User directed: **no Gate B emails, "assume ≥2 yes", build the model and runs.** So EP-4 pivoted from the deterministic VirulenceFinder **rule-resolver (v0)** to the **LEARNED-MODEL track**: genome → NT v2 100M embeddings (per-CDS, mean-pooled per strain) → XGBoost → leave-one-strain-out AUROC. Fully laptop-local: no VirulenceFinder, no workhorse, no emails.

## What's BUILT + COMMITTED (on `main`, pushed through 10ee1f7; later commits local)

- `scripts/pathotype_laptop_pipeline.py` — download(GCA)→Bakta→read `.ffn`→NT-embed→mean-pool→LOSO. Commits `3f74aa4`, `d509ccc`.
- LT/ST pipeline-validation smoke ran end-to-end: **pipeline PASS**, AUROC 0.0 = **small-N degenerate artifact** (8 strains, two-valued rank-inverted probas), NOT a model verdict; LT/ST label is circular anyway. Commit `87c1916`. Result: `research_outputs/pathotype_smoke_lt_vs_st_result_2026-05-31.{json,md}`.
- Bake-off plan: `research_outputs/pathotype_model_state_and_bakeoff_plan_2026-05-31.md`.
- `scripts/pathotype_kmer_bakeoff.py` — **UNCOMMITTED** (committing with this handoff).

## Substrate (all resolved, by-accession, no re-assembly)
- ExPEC 135 + EPEC 125 = Horesh WGS-masters (`research_outputs/horesh_bounded_slice_wgs_accession_candidates_2026-05-30.csv`). Fetch via **ENA bare URL** `https://www.ebi.ac.uk/ena/browser/api/fasta/<PREFIX01>` (NOT `?download=true&gzip=true` — that returns empty). Some Salipante rows carry SRA run IDs → filter to WGS-master shape `^[A-Z]{4}\d{8}$`.
- ETEC: 8 von Mentzer refs (LR accessions) + 558 GCA collection (`research_outputs/etec_vonmentzer_collection_gca_2026-05-30.csv`).
- **Label provenance:** ExPEC = isolation-site (independent, STRONG external validity); EPEC = DECA-curated (medium); ETEC = toxin-typed (CIRCULAR — near-conformance only).

## 🔴 IN-FLIGHT right now: k-mer bake-off arm
- `scripts/pathotype_kmer_bakeoff.py` running (background task; 2 python procs live at handoff). 24 strains (12 Salipante-ExPEC + 12 Hazen-EPEC) fetched + cached to `data/ena_wgs/` (all 24 present). k-mer k=8 LOSO crunching; **no result packet written yet.**
- Expected output: `research_outputs/pathotype_bakeoff_kmer_expec_epec_2026-05-31.json`.
- **FIRST ACTION in fresh window:** check if that JSON exists. If yes, read `auroc` + `degenerate`. If the python died with D:/other failure, just re-run `.venv/Scripts/python.exe scripts/pathotype_kmer_bakeoff.py` (fetches are cached on C:, so it's just the LOSO compute, ~few min).

## Pre-committed interpretation of the k-mer result (the confound diagnostic)
ExPEC(Salipante) vs EPEC(Hazen) is **study-confounded** (no single study has both classes). So:
- **k-mer AUROC ≈ 1.0** → contrast is trivially/batch-separable; NT lift is UNMEASURABLE on this substrate → **PIVOT** to within-study or lineage-matched labels; do not scale NT.
- **k-mer ~0.6–0.85** → real headroom → run the **windowed-NT arm** next, test if NT beats k-mer by **≥3pp** (the project's mandated bar, CLAUDE.md:211 "don't soften it").
- **≤2 distinct scores / degenerate** → small-N artifact (like the smoke) → treat AUROC as uninformative; bump N or fix splits.

## Next steps after the k-mer number lands
1. Branch per above.
2. If NT arm warranted: add a **whole-genome windowed NT** arm (tile genome, mean-pool — NO Bakta, avoids the ~40hr/250-genome annotation bottleneck) before any per-CDS scaling. Cache embeddings (HDF5) so pooling variants don't re-pay GPU.
3. Honest discipline throughout: report resolver-conformance vs external-validity separately; tet precedent (NT mean-pool failed 0.40 on distributed mobile-element signal) is a live risk for distributed virulence (EAEC/ETEC) — the k-mer control is the honest check.

## Open risks / notes
- **Confounding:** even with a good AUROC, check it isn't study/geography/assembler batch signal. Prefer lineage/ST-aware splits.
- **tet precedent:** NT-frozen mean-pooling may dilute sparse distributed virulence signal — same failure mode as tet AMR. Baseline-beat-by-≥3pp is the gate.
- Two parked workhorse/user items (now de-prioritized per "no emails, build model"): A2 (workhorse commits `pathotype_horesh_*` scripts) and Gate B outreach. Not on the model-track critical path.

## Uncommitted at handoff (being committed with this doc)
`scripts/pathotype_kmer_bakeoff.py`, `dogfood-observations/`, `soraya_runs/2026-05-31-1333-ep4-a1-vf-runtime/`, helper scripts. Stray empty dir `2026-05-30-1200-ep4-pathotype/` at repo root can be removed (leftover from a Soraya run-dir mislaunch). The `soraya_runs/.soraya-run.lock` shows deleted — fine (no active run).

## Project north star (unchanged)
"AI DNA decoder TOOL, not papers; failure-tolerant iteration; honest scope limits." Ledger: `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md`. Roadmap: `plans/Trait_Decoding_Roadmap.md`.
