# HANDOFF — dna_decode session end 2026-06-05 (fresh-session entry + workhorse transfer)

> Read this FIRST. Single source of truth for current state after a long session that CLOSED the AMR
> embedding question and SHIPPED the deterministic AMR decoder. Works for both a fresh laptop session and
> the workhorse. Everything is on origin/main (HEAD **a6b7a9e**); `git pull` to sync.

---

## 0. ⚠ Do this first (cwd + env)
```
cd C:\Users\Farshad\PythonProjects\dna_decode
!pwd        # must be /c/Users/Farshad/PythonProjects/dna_decode
git pull origin main
```
Environment facts that bit us this session:
- **Sync = git on `main`.** `data/` is gitignored → caches/cohorts/AMRFinder outputs do NOT sync; re-derive from public accessions. Remote: github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION.
- **Disconnects were location changes (mom's house ↔ home), NOT a flaky machine** — overnight/stationary is stable.
- **uv cache intermittently throws "device not ready"** → run python via `.venv/Scripts/python.exe` directly (bypasses uv cache) when it flares.
- **Docker Desktop does NOT share non-C: drives** — a DB mounted from D: shows EMPTY in the container. Keep Docker-read DBs on C:. The AMRFinder DB is now at real-C: `data/amrfinder_db/` (240 MB, gitignored). `dna_decode_stage2` on this host is a C:→D: junction; the 4 GB bakta_db stays on D:.
- AMRFinder ≈ 2–4 min/strain on the 860M; full N=147 backfill ≈ 45 min wall.

---

## 1. What happened this session (one paragraph)
The AMR Phase-2 embedding-vs-classical question was driven to a **decisive, evidence-complete verdict** and
then **CLOSED**: on the only de-confounded AMR substrate (cipro N=147), the frozen-NT-mean-pool embedding
(0.914) **LOSES to the QRDR-POINT knowledge baseline (0.943)**, and within-lineage it's at chance
(0.605, p=0.365) vs POINT's 1.000 — i.e. NT learned lineage, not the resistance mechanism. With cef
confounded, tet anti-predictive, and gent infeasible, **NT embeddings have no E. coli AMR niche.** Decision:
the AMR decoder is a **mechanism-feature tool**, and that tool was **shipped** (`dna-amr` via
`scripts/amr_predict.py`). User then set the go-forward **thesis = "both in parallel"** (ship deterministic
tools AND hunt the embedding substrate).

## 2. Current state (what's DONE)
- **AMR Phase-2: CLOSED.** Decision memo `plans/AMR_embedding_niche_decision_2026-06-05.md`; verdict +
  within-lineage diagnostic in `wiki/ciprofloxacin_falsifier_2026-06-05.md` + `..._within_lineage_diagnostic_2026-06-05.md`.
- **Deterministic AMR decoder SHIPPED** (lane 1): `dna_decode/eval/amr_rules.py` + `scripts/amr_predict.py`
  CLI + `tests/test_amr_rules.py` (6 tests). cipro N=147 op-chars **acc 0.85 / sens 0.96 / spec 0.75**;
  transparent calls (e.g. `R [HIGH | gyrA_S83L, gyrA_D87N, parC_S80I, parE_S458A]`).
- **Pathotype v0 decoder** already shipped earlier (`dna-pathotype`, `dna_decode/pathotype/`).
- **Reusable Phase-2 infra** (the durable value): `dna_decode/eval/cohort_deconfound.py` (de-confound gate,
  a precondition for ANY falsifier), `scripts/amr_falsifier.py` (drug-agnostic, CI-aware, `--amrfinder-runs`/
  `--skip-kmer`), `dna_decode/eval/point_baseline.py` (QRDR-POINT comparator), `scripts/within_lineage_diagnostic.py`,
  `dna_decode/eval/loso_kmer.py` (per-genome k-mer count cache; killed the N² stall).

## 3. The thesis + the two lanes (go-forward, user-ratified)
**Thesis = BOTH IN PARALLEL:** ship deterministic tools AND fund the embedding substrate-hunt.

- **Lane 1 — ship deterministic tools.** MVP done (AMR decoder). Possible follow-ups: package `dna-amr`
  as a console entry in pyproject (mirror `dna-pathotype`); README section; a tiered clinical rule to lift
  the 0.75 specificity (cipro clinical-R usually needs ≥2 QRDR hits — the data supports a "≥2 → R, 1 →
  reduced-susceptibility" refinement).
- **Lane 2 — embedding substrate-hunt (the REAL open question).** Find a phenotype where embeddings could
  beat domain knowledge: **sampling-INDEPENDENT labels AND no curated mechanism/knowledge baseline.** This
  is gated on the de-confounded-labels constraint that has blocked the project **3×** (pathotype circular
  labels → cef geography confound → cipro within-lineage). It's an `/idea-anchor` + `/research` cycle, NOT
  modeling. **Start it FRESH** (it was deferred from the end of this marathon session).

## 4. The next concrete move (recommended)
`/idea-anchor` on: *"a phenotype where a frozen-genome embedding could beat domain knowledge — i.e. has
sampling-independent (lab-measurement) labels AND no curated mechanism catalog to lose to."* That single
question decides whether the embedding architecture has ANY honest niche. If no such substrate is found,
the project is fully a deterministic genome→trait tool (which already ships).

## 5. Explicitly NOT next
- Re-tuning NT pooling / new cef-S or gent cohorts to retry AMR embeddings — the cipro result already shows
  knowledge baselines win on concentrated mechanisms; a new AMR cohort would re-confirm, not overturn.
- Any further AMR embedding work. AMR = mechanism-feature tool, settled.

## 6. Key meta-lessons (saved to cross-project memory)
- **Beating k-mer ≠ working.** Beat the DOMAIN-KNOWLEDGE baseline AND check within-lineage concordance —
  overall AUROC on a de-confounded cohort still conflates lineage with mechanism.
  (`memory/feedback_embedding_vs_knowledge_baseline_and_within_lineage.md`)
- **Sampling-defined labels → intrinsic confound** (de-confound gate is necessary, not sufficient).
  (`memory/feedback_sampling_defined_phenotype_intrinsic_confound.md`)

## 7. Workhorse notes
Nothing GPU-required is pending. The workhorse is the GPU box; lane-2 (substrate hunt) is research, not
GPU-bound. If lane 2 ever yields a de-confounded cohort needing NT embeddings, the workhorse populates the
NT cache (the one GPU step) and returns the `.h5`; the de-confound gate + falsifier + POINT baseline all
run on the laptop. Prior AMR-result handoff: `wiki/HANDOFF_amr_phase2_RESULT_2026-06-05.md` (superseded by
this doc for overall state; still the detailed AMR-result reference).

## 8. Provenance
HEAD a6b7a9e. Driven via repeated `/soraya` + 3× `/brainstorm` (each round caught a load-bearing flaw:
cef confound → de-confound gate; k-mer-only comparator → POINT baseline; overall-AUROC conflation →
within-lineage diagnostic). AMR ledger `project_state/dna-decode-2026-05-11.md` rows 65–75 + Decisions Made.
