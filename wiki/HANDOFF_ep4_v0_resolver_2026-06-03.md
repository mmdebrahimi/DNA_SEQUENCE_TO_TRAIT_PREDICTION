# HANDOFF — EP-4 v0 pathotype resolver — 2026-06-03

> Fresh-context entry point. Supersedes the 2026-06-02 NT-arm handoff (that thread closed: pooled NT failed). The project pivoted to — and SHIPPED — the ledger-locked v0 deterministic virulence-gene resolver. All work committed + pushed to `main` (`github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`), synced through `ab1fffa`.

## The arc (2026-06-02 → 06-03, dogfooding Soraya on the real project)
1. **Learned bake-off closed as confounded.** k-mer top-10k 0.514 → full-spectrum presence/absence 0.729 → pooled NT 0.38 (FAIL). Then found the substrate is **study==class** (ExPEC=Salipante / EPEC=Hazen), and the 0.729's top features are rare CTAG-motif k-mers = **assembler batch, not biology**. Every learned AUROC here is uninterpretable as biology.
2. **Pivoted to the confound-immune known-marker resolver** (ledger v0 design). It keys on SPECIFIC virulence genes → valid despite study==class. Validated: eae 12/12 EPEC vs 0/12 ExPEC (AUROC 1.0), textbook biology.
3. **Built the v0 CLI** (`dna_decode/pathotype/`): 23 marker clusters → 11-class decision table + abstention → FASTA-in provenance JSON.
4. **Cohort-evaluated + calibrated to H4 ship-gate** (overnight).

## CURRENT STATE — v0 resolver, shipped + H4-passing
- Package `dna_decode/pathotype/`: `markers.py` (23 clusters), `resolve.py` (11-class table + abstention + RULE_EXPEC_001 calibration), `detect.py` (k=15 dual-strand seed presence, located provenance, QC), `cli.py` (`python -m dna_decode.pathotype <fasta>`).
- Tests: `tests/test_pathotype_resolver.py` — **20/20** green (run: `.venv/Scripts/python.exe tests/test_pathotype_resolver.py`).
- **H4 metrics (24-genome cohort):** EPEC recall 1.0, ExPEC recall 0.75, confident supported-call precision **1.0**, ExPEC abstention **0.08**. Both ship-gate targets met (precision ≥0.95, abstention ≤15%).
- Marker DB cached: `data/virulencefinder_db/virulence_ecoli.fsa` (gitignored). Per-genome coverage cached: `data/pathotype_cov_cache/` → resolver re-tuning is instant (no re-detect).
- Demos: `research_outputs/pathotype_cli_demo_{AIEY01,JSIS01}.json` (EPEC→tEPEC_COMPATIBLE; ExPEC→ExPEC_COMPATIBLE).

## ONE open ledger v0 requirement (gated — needs the user)
**CGE VirulenceFinder side-by-side diff.** Blocked unattended: no BLAST/KMA/conda, uv venv has no pip-aligner. Path documented in `research_outputs/pathotype_vf_sidebyside_feasibility_2026-06-03.md`. Install BLAST+/KMA → then a `scripts/` wrapper can run real VF + diff vs our `marker_hits`.

## Recommended next (ranked)
1. VF side-by-side diff (after BLAST+/KMA install) — last v0 ledger item.
2. Per-gene ExPEC scoring (v0.1) — split SIDEROPHORES/CAPSULE clusters into per-gene presence; push ExPEC recall past 0.75 (cheap with the coverage cache).
3. Add ETEC genomes (von Mentzer, `research_outputs/etec_vonmentzer_collection_gca_2026-05-30.csv`) → exercise the ETEC arm + 3-class supported-surface metrics.
4. Break study==class (within-study / ST-matched genomes) before any learned-vs-classical claim.

## EMIT — pending user-only `/project-state` calls (Soraya can't self-invoke)
In each soraya_runs/.../recommendation.md. Net: append-decisions for (a) confound finding, (b) resolver validated, (c) v0 CLI shipped, (d) H4 cohort eval + calibration. These advance ledger Goals 3 + 5.

## Soraya dogfood (4 runs this session, in `dogfood-observations/`)
- 🔴 `rm -rf` runs un-gated under money-only default → add a `destructive-local` always-pause subclass.
- 🟠 auto-metric/auto-MVP verdicts certify the literal predicate only — keep the model-level substance check explicit (the CTAG override; "MVP-met but ExPEC under-sensitive").
- 🟡 audit-trail appended at run-end (model-discipline) not per-step — the OT1 seam.
- ✅ gate classifier + lock/run-set + de-risk-first + MVP predicates + coverage-cache-for-instant-reruns all behaved well; unattended dep-install correctly declined.

## Environment (unchanged)
cwd launches in `rca_engine/articles` — use absolute paths / `git -C <dna_decode>`. venv `.venv/Scripts/python.exe` (uv). D: online (`HF_HOME=D:\hf_cache`). C: 44 GB free. GTX 860M = display GPU (4 GB, TDR 2 s, batch≤2 — see `feedback_gtx860m_fm_embedding_oom_tdr` memory) but the resolver is CPU-only/pure-Python.
