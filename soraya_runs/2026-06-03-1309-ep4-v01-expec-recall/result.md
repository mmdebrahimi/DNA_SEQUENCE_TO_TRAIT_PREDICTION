# Result — minister run `2026-06-03-1309-ep4-v01-expec-recall`

**VERDICT: `mission-mvp-reached`** ✅ (maiden real `/soraya minister` run driven to terminal).

## Mission endpoints (both live-MET)
1. `test-exit-0 python -m pytest tests/test_pathotype_expec_recall.py` → **PASS** (5/5; ExPEC recall 0.917 ≥ 0.85
   AND confident-supported precision 1.0).
2. `project-state-row …ecoli-pathotype-prediction-cli-2026-05-26.md:ExPEC recall hardened` → **present**
   (Action Log row 18).

## What shipped
- `dna_decode/pathotype/markers.py` — `EXPEC_SUPPORT_GENE_PREFIXES` + `EXPEC_SUPPORT_GENE_K=1`; RULES_VERSION → v0.2.0.
- `dna_decode/pathotype/resolve.py` — `support_gene_count` kwarg + `RULE_EXPEC_002` (≥K support genes → ExPEC LOW).
- `dna_decode/pathotype/expec_score.py` — per-gene support-burden scorer (new).
- `scripts/build_pergene_support_cache.py` + committed `data/pathotype_pergene_cache/` (24 genomes, offline).
- `tests/test_pathotype_expec_recall.py` — 5 gate tests (recall, precision-invariant, LOW-tier-cap, LEE-gate, JSLG-stays-AMBIGUOUS).

## Metrics (24-genome H4 cohort)
| | before | after |
|---|---|---|
| ExPEC recall | 0.75 (9/12) | **0.917 (11/12)** |
| confident-supported precision | 1.0 | **1.0** (n=15) |
| EPEC recall | 1.0 | **1.0** (LEE-gate dominates, incl AAJW=4 support genes) |
| resolver test suite | 20/20 | **25/25** (0 reversed) |

JSLG (1-strong/0-support) left AMBIGUOUS by design (Option A — no existing-test reversal). K=1 in-sample on N=24
(documented scope-limit; LOW tier + DEC-module gate bound the out-of-cohort risk).

## Finiteness
Φ debited 2→1 on the round-1 GENERATE (durable in `lifecycle-journal.ndjson`). Resume hit `mission_met()` at
the round's first check (work done in-session) → terminal before a second generate. Promotion/`/project-init`
seeding of the family ledger was preempted by the terminal-check fast-path (the family ledger is the MEANS,
not a mission endpoint). No money/destructive action fired (pure-Python CPU-only); lease armed + disarmed.
