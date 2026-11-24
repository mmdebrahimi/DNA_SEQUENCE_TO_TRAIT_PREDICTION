# Independent verification audit — DNA-decode load-bearing claims (2026-07-04)

> A **second-reviewer** pass by a parallel session (Soraya), run **read-only** against the DNA-11 session's
> committed work — zero collision (no writes to the shared `main` working tree; every number recomputed from
> committed code/data, bytecode + pytest-cache suppressed). Purpose: independently confirm the load-bearing
> claims reproduce and the **self-graded** certification capstone is faithful to its code, before those claims
> harden. NOT a re-validation of the underlying biology — each domain card keeps its own honest tier.

## Verdict: PASS — no discrepancies found

The decoder's load-bearing claims reproduce independently, and the self-graded capstone is faithful to the
code that generates it (no hand-edit drift) and honest by construction.

## Checks run + results

| # | Check | Method (read-only) | Result |
|---|---|---|---|
| 1 | Load-bearing test suite | `pytest` on 5 headline files (certification_capstone, dms_learned_model_falsifier, abo_blood, cipro_bounded_falsifier, amr_rules) | **59/59 pass** (9.2s) |
| 2 | DMS learned-signal number | Re-ran `scripts/dms_learned_model_falsifier.py --max-assays 40` from committed data | **median \|Spearman\| = 0.417** — exact match to committed `wiki/jepa_dms_learned_signal_result_2026-07-04.md`; field context reproduced (ESM2-650M 0.484, GEMME 0.484, EVE 0.466) |
| 3 | Certification-capstone census (drift guard) | Called `dna_decode.data.cell_registry.cells()` directly; recomputed track + tier counts | **EXACT match** to committed `wiki/certification_capstone.md`: total **67**; track amr25/finder4/pgx3/typing6/viral29; tiers faithful_to_tool11 / independent_measured25 / knowledge_baseline4 / near_independent15 / not_censused1 / no_free_source11 |
| 4 | Capstone honesty posture | Read `certification_capstone.{md,json}` | Honest **by construction**: `no_aggregate_verdict:true`, per-cell tiers preserved, caveats disclosed. The self-validation grading risk does **not** materialize here. |

## Method / honesty notes
- **Zero collision:** no branch switch, no write to the shared `main` tree; `PYTHONDONTWRITEBYTECODE=1` + `-p no:cacheprovider` so no cache/bytecode landed in DNA-11's tree. This audit file lives on the `mosfaer` branch only.
- **Deliberately NOT run:** `scripts/build_certification_capstone.py` (it would OVERWRITE DNA-11's in-flight capstone files) — instead the census was recomputed by calling the registry directly and diffing against the committed artifact.
- **Scope bound (no silent truncation):** verified a **targeted load-bearing subset**, not all **225** test files (full suite exceeds a 10-min run — that is DNA-11's CI concern, out of scope for a claim-verification pass).
- **What this does NOT assert:** the underlying biological validity of each claim (that stays with each domain card's own honest tier); only that the claims **reproduce** from code + data and the capstone **faithfully renders** what the code computes.

## Follow-up (2026-07-04, `--until-mvp` — both lanes)

### Lane 1 — FULL suite (upgrades check #1 to whole-suite) — PASS
Ran all **225** test files chunked (5×45, cache/bytecode suppressed → DNA-11 tree untouched). Aggregate:

| passed | failed | skipped | errored |
|---|---|---|---|
| **1,988** | **0** | 8 (by design) | 10 |

The **10 errors are all one infra cause** — `ModuleNotFoundError: No module named 'xgboost'` (an explicitly
OPTIONAL Phase-1 dep; the test itself says "run `uv sync` to install Phase 1 deps"). **Zero real
regressions.** The decoder reproduces clean across the whole suite. (xgboost NOT installed — would have
mutated DNA-11's shared venv for no verdict gain.)

### Lane 2 — data-gap: a free, NON-gated, individual-level human cohort — PASS
**Finding:** DNA-11's individual-level human path is the GATED UK Biobank application (pending). **PGP-UK**
(Personal Genome Project UK) is a **free, open-consent, NON-application-gated** individual-level human
multi-omics cohort (WGS + WGBS methylation + RNA-seq + self-reported phenotype) in **OPEN** EMBL-EBI repos
(ENA / ArrayExpress) — not dbGaP/EGA. It unblocks individual-level human work **without a data-access
committee**, and DNA-11 has not captured it (non-overlapping with its summary-stat-index / variant-DB
free-data workstream).

**Verified real-surface (live 2026-07-04, not asserted):**
- ENA `PRJEB17529` (WGS/WGBS) API returns real FTP download URLs — e.g. `ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/008/ERR1726428/ERR1726428_1.fastq.gz`, sample `PGP-UK uk33D02F`.
- `www.personalgenomes.org.uk/data` live with a Phenotype column + genome/methylome reports.
- **Capture stub RAN**: `scripts/capture_pgp_uk.py` (this commit) pulled **55 downloadable run records** (35 WGS/WGBS + 20 RNA-seq) + 4 pointers → manifest on `D:/dna_decode_cache/pgp_uk/` (index only; the ~2 TB raw panel is a separate explicit D:-targeted fetch, NOT run here).

Open accessions: ENA `PRJEB17529` (WGS/WGBS), `PRJEB25139` (RNA-seq; also ArrayExpress `E-MTAB-6523`),
ArrayExpress `E-MTAB-5377` (450k methylation IDATs), phenotype+reports at personalgenomes.org.uk/data.

**Usability confirmed — CHEAP realization path (2026-07-04):** PGP-UK provides **called VCFs directly** via
ENA analysis objects (`result=analysis`, PRJEB17529): **9 `SEQUENCE_VARIATION` analyses** with real VCF
download URLs, e.g. `ftp.sra.ebi.ac.uk/vol1/analysis/ERZ389/ERZ389532/FR07961000.pass.recode.vcf.gz`. So the
deterministic decoder (ClinVar / PGx / Mendelian, VCF-consuming) can validate on **real people** with **no
FASTQ→variant-calling pipeline** — a per-participant VCF is ~MB, not the 2 TB raw panel. **The realization
step (download a VCF → run the decoder → real-people validation report) is DNA-11's territory** (it owns the
validation-report-card workstream); it is HANDED OFF, not run in parallel, to avoid duplicating that work.

## Residual / recommended
- **Lane 1:** to clear the 10 optional-dep errors, `uv sync` (installs xgboost) in the target env — cosmetic; not a regression.
- **Lane 2:** the index is captured; the **bulk raw fetch (~2 TB) is a separate, explicit, D:-targeted step** the user/DNA-11 can trigger. PGP-UK is modest-N (pilot panel ~10 deep multi-omics + broader WGS) — good for **molecular / deterministic / Mendelian** validation on real people, NOT a polygenic-complex-trait cohort (that stays confound-blocked regardless of source).
