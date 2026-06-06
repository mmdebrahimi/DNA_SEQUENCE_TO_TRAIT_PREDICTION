# HANDOFF — AMR Phase-2 falsifier RESULT (workhorse) — 2026-06-05

> Read FIRST in a fresh workhorse session. The decisive Phase-2 question is ANSWERED. This supersedes
> `wiki/HANDOFF_amr_phase2_embedding_frontier_2026-06-04.md` (the "go run the falsifier" handoff — done).
> Everything below is committed to origin/main (HEAD 61e502b); `git pull` to sync.

---

## 0. The one-line result

**On the cleanest possible substrate (cipro N=147, the only de-confounded AMR cohort), the NT-frozen-
mean-pool embedding FAILS vs the QRDR-POINT domain-knowledge baseline.** NT beats bag-of-k-mers but
loses to mechanism features, and within-lineage it's at chance → NT's headline AUROC was largely
lineage, not the resistance mechanism. Honest north-star FAIL. The Phase-2 embedding-vs-classical
question is settled for cipro.

## 1. The numbers (all leakage-safe `leave_one_accession_out`, N=140 eff, 67R/73S)

| Variant | AUROC | within-lineage concordance (p) |
|---|---:|---|
| **POINT-XGB** (QRDR/plasmid knowledge baseline) | **0.943** | **1.000 (p<0.001)** = pure mechanism |
| NT-XGBoost (the embedding) | 0.914 | 0.605 (p=0.365) = chance |
| NT-logreg | 0.863 | 0.628 (p=0.326) |
| k-mer-XGB (sequence baseline) | 0.824 | — |

- NT vs best-classical(POINT): **gap −2.9 pp, 95% bootstrap CI [−9.0, +2.9] → FAIL** (B=1000, paired 140/140).
- NT vs k-mer: +8.9 pp (the earlier "win" — but k-mer is the wrong bar).
- **Interpretation:** cipro resistance IS gyrA/parC/parE point mutations; the POINT feature set encodes
  the mechanism directly and discriminates R/S even *within* a lineage. NT-frozen-mean-pool dilutes that
  localized signal across the whole genome → its overall edge is lineage/genome-content, not biology.

Artifacts (on origin/main): `wiki/ciprofloxacin_falsifier_2026-06-05.{md,scores.json}` +
`wiki/ciprofloxacin_within_lineage_diagnostic_2026-06-05.md` + `wiki/ciprofloxacin_mechanism_audit_2026-06-05.{md,json}`.

## 2. Reusable infra built this session (committed, tested)

- `dna_decode/eval/cohort_deconfound.py` — **cohort de-confound gate** (a PRECONDITION for any falsifier):
  within-lineage R/S co-occurrence + matched-support floor + per-axis (MLST/country/year) contrast +
  3-state promotability (ADMIT/DIAGNOSTIC/BLOCK). cef gate_b → CONFOUNDED (blocked); cipro → DE_CONFOUNDED.
- `scripts/amr_falsifier.py` — drug-agnostic falsifier: de-confound gate → NT-XGBoost/NT-logreg/k-mer/POINT
  under leakage-safe CV → CI-aware verdict. Flags: `--amrfinder-runs <dir>` (adds POINT), `--skip-kmer`.
- `dna_decode/eval/point_baseline.py` — QRDR-POINT knowledge baseline from the AMRFinder cache.
- `scripts/within_lineage_diagnostic.py` — mechanism-vs-lineage discriminator (the clincher above).
- `dna_decode/eval/loso_kmer.py` — per-genome k-mer count cache (killed the N² within-fold-rebuild stall).
- `scripts/drug_mechanism_audit.py` — machine-agnostic AMRFinder DB default + `latest`-symlink→real-dir
  Docker mount fix.

## 3. Environment gotchas hit (so the workhorse doesn't repeat them)

- **Docker Desktop does NOT share non-C: drives** — a DB mounted from D: appears EMPTY in the container
  ("No valid AMRFinder database found"). Keep Docker-read DBs on a Docker-shared drive. The 240MB
  AMRFinder DB is now at real-C: `data/amrfinder_db/` (gitignored). Diagnose: `docker run -v X:/x img ls /x | wc -l` (0 = not shared).
- AMRFinder `latest` is a relative symlink not followed in-container → `_run_amrfinder` mounts `latest.resolve()`.
- `data/` is gitignored → NT cache, refseq, amrfinder_runs/db do NOT sync; re-derive from public accessions.
- AMRFinder ≈ 2–4 min/strain on the 860M laptop (full backfill ~45 min wall for 99 strains). Faster on the workhorse.

## 4. The strategic fork — pick the next epoch (NOT decided tonight)

The cipro FAIL narrows where frozen-mean-pool embeddings could still earn their keep:

| Option | Rationale | Cost |
|---|---|---|
| **A. Reconsider embeddings for AMR at all** | cipro (concentrated) FAILS vs POINT; tet (distributed) failed earlier (0.400). If neither concentrated NOR distributed AMR works, the architecture has no AMR niche → pivot AMR to a mechanism-feature tool (POINT/AMRFinder), drop embeddings here. | research/decision |
| **B. Test the one untested AMR niche** | distributed mechanism WITH a de-confoundable cohort + NO clean knowledge baseline (where embeddings *should* have an edge). Needs a new de-confounded cohort (cef was confounded; tet substrate infeasible per 2026-05-18 census). | high (cohort build) |
| **C. Embeddings → non-AMR phenotypes lacking curated catalogs** | the honest niche: phenotypes with sampling-INDEPENDENT labels but NO domain-knowledge feature set to lose to. Per the roadmap, this needs a de-confounded labeled substrate (the recurring binding constraint). | research/substrate |

**Recommendation:** A is a decision (cheap — make it), then if embeddings are kept, C is the real frontier
(but gated on de-confounded labels, the constraint that's bitten 3×). Do NOT chase B without a substrate.

## 5. Workhorse next move (if any)
Nothing is REQUIRED on the workhorse — the result is complete + committed. The workhorse is the GPU box;
the remaining questions (A = decision, C = substrate research) are not GPU-bound. If a new de-confounded
cohort gets built (B/C), the workhorse populates its NT cache (the one genuinely GPU-bound step) and
returns the .h5; the falsifier + de-confound gate + POINT baseline all run on the laptop.

## Provenance
Driven via repeated `/soraya` + 3× `/brainstorm` (each round caught a load-bearing flaw: cef confound →
de-confound gate; k-mer-only comparator → POINT baseline; overall-AUROC conflation → within-lineage
diagnostic; CI alignment bug). Cross-project lesson saved at
`~/.claude/.../memory/feedback_embedding_vs_knowledge_baseline_and_within_lineage.md`.
