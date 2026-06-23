# Session handoff — independent-validation breakthrough + TB holy grail (2026-06-23)

Continue in a fresh session from here. This session BROKE the project's binding constraint (a free,
independent, measured-phenotype label) and delivered the first genuinely-independent numbers for the
deterministic decoder across **5 bacterial cells + M. tuberculosis**, plus a new SARS-CoV-2 cell. The repo
is now **PUBLIC**. One background job is still running (the full TB number).

---

## 0. THE ONE LIVE THING — check this first
A **background run** is computing the full-cohort independent TB number. When you resume:

```bash
cd C:/Users/Farshad/PythonProjects/dna_decode
# 1. is it still running / did it finish?
ls D:/dna_decode_cache/tb_indep/vcf/*.vcf | wc -l      # target ~2,847
wc -l < D:/dna_decode_cache/tb_indep/results.jsonl     # scored isolates (60 at handoff; ~2,800 when done)
tail -20 D:/dna_decode_cache/tb_indep/full_run.log     # final RIF/INH sens/spec if finished
cat wiki/tb_independent_amr_portal_scores.json         # the aggregate (overwritten by the full run)
```
- **If finished:** the JSON + log hold the full independent RIF/INH sens/spec. Bank it (memo update + ledger).
- **If interrupted/died:** just re-run — it is CHECKPOINTED + skips fetched/aligned/scored isolates:
  ```bash
  TB_INDEP_WORK="D:/dna_decode_cache/tb_indep" uv run python -m scripts.run_tb_independent_amr_portal \
      --max 0 --work "D:/dna_decode_cache/tb_indep"
  ```
  (Needs Docker running. Transient NCBI `IncompleteRead` fetch fails are tolerated — re-run retries them.)
- **Status at handoff:** all 2,847 assemblies fetched; alignment ~332/2,847 done; scoring runs after align
  completes (the checkpoint stays at 60 until then). N=60 smoke already verified: **RIF acc 0.917
  (sens 0.900/spec 0.950), INH acc 0.883 (sens 0.854/spec 1.000)** — `wiki/tb_independent_number_2026-06-23.md`.

---

## 1. The big picture — what changed this session
**The binding constraint was LABELS, not models.** For the whole TB gold-set saga, 5 sources were all
author-request / DUA / circular. This session found the **EBI AMR Portal (MRC CABBAGE)** —
`https://ftp.ebi.ac.uk/pub/databases/amr_portal/` — a FREE, public, per-isolate, **measured-AST**,
accession-linked, genotype-paired dataset (1.71M rows). That broke the wall:

- **Independent validation of the FROZEN bacterial decoder** (E. coli / Salmonella / Klebsiella / Shigella),
  scored by reconstructing AMRFinder `main.tsv` from the AMR Portal's own AMRFinder genotype table → the
  frozen `call_resistance` rule UNCHANGED → vs measured phenotype. 21 SCORED_INDEPENDENT cells.
- **The independent TB number** (this session's holy grail) — a different path (raw assembly → variant call →
  WHO rule), because TB's rule needs raw VCFs not AMRFinder calls. N=60 done; full 2,847 in flight.
- **Calibrated registry promoted** (Salmonella + Klebsiella cipro → `INDEPENDENTLY_VALIDATED`,
  metadata-only freeze-amendment).
- Earlier in the session: a **SARS-CoV-2 Mpro cell** (the "next HIV" — free CoV-RDB label) was built+validated.

---

## 2. Key independent results (all committed, all on origin/main, public repo)
**Bacterial (frozen rule, AMR-Portal provenance-disjoint, measured AST)** — `wiki/amr_portal_independent_report_card.md`:
| Organism | best cells (acc) |
|---|---|
| Salmonella | cef 0.995, gent 0.987, tet 0.983, cipro 0.959 (calibrated, N=24,972) |
| E. coli | gent 0.987, tet 0.980, cipro 0.950, cef 0.919, mero 0.986 |
| Klebsiella | cef 0.948, gent 0.949, cipro spec 0.994 (oqxAB-exclusion confirmed; sens 0.755 = non-QRDR blind spot) |
| Shigella | cef 0.993, tet 0.958 |

**TB (independent, this session, N=60 smoke; full run in flight):** RIF 0.917, INH 0.883 — lands in the
in-distribution CRyPTIC-baseline range (RIF 0.916/0.974, INH 0.889/0.989), as a faithful test should.

**Honest rails on ALL of the above:** independent at the ACCESSION level (BioSample/GCA disjoint vs CRyPTIC +
our cohorts — an upper bound; BioSample cross-archive resolution would tighten it); measured phenotype =
non-circular; rule applied UNCHANGED; frozen surface byte-stable.

---

## 3. What's NOT done (ranked follow-ons — none blocked, all code-closable)
1. **Bank the full TB number** when the background run finishes (Section 0). Likely just memo + ledger.
2. **Lineage-collapse the TB number** — the publication-grade headline. TB R classes are clonally dominated;
   raw sens/spec is clonality-inflated. Run Mash clustering + the FROZEN `clonality.cluster_weighted_confusion`
   (the in-distribution baseline `scripts/score_tb_cryptic.py::score_cohort` already does this; the runner
   `run_tb_independent_amr_portal.py` currently emits RAW only — wire in the lineage-collapse path). Needs
   Mash (Docker, `dna_decode/eval/phylogeny.py::compute_mash_distances(..., use_docker=True)`).
3. **Tighten independence to BioSample-resolved** (Gate-0 `dna_decode/eval/biosample_resolver.py`) for a
   headline-publishable "not even the same isolate under a different accession" claim.
4. **Fold the bacterial + TB independent numbers into the standing report card surfaces** (the AMR-Portal card
   exists; consider a top-level cross-kingdom independent summary).
5. **SARS-CoV-2 v0.1 independent number** (held-out CoV-RDB studies — the free path; `wiki/sarscov2_mpro_validation_result_2026-06-23.md` is in-distribution).
6. **HIV PI/INSTI deconfounded v0.1** (datasets already on disk `data/raw/hiv/`; zero-risk quick win).

---

## 4. Files / scripts that matter (this session's deliverables)
- **`scripts/run_tb_independent_amr_portal.py`** — the TB pipeline (fetch→align→call→score, checkpointed,
  `--work` / `$TB_INDEP_WORK` for D:). Tests `tests/test_run_tb_independent.py`.
- **`scripts/amr_portal_feasibility.py`** — provenance-disjoint powering census (74 disjoint powered cells).
- **`scripts/amr_portal_score_independent.py`** — scores the FROZEN bacterial rule on disjoint isolates
  (per-organism routing; reconstructs main.tsv from the AMR Portal AMRFinder genotype table).
- **`scripts/build_amr_portal_report_card.py`** — standing independent report card.
- **`scripts/amr_portal_tb_cohort.py`** — the disjoint-TB cohort builder.
- **`dna_decode/data/sarscov2_amr.py`** + `scripts/sarscov2_caller.py` + `scripts/sarscov2_mpro_validate.py`
  — the SARS-CoV-2 Mpro cell.
- **Memos:** `wiki/ebi_amr_portal_finding_2026-06-23.md`, `amr_portal_feasibility_result_2026-06-23.md`,
  `amr_portal_independent_validation_2026-06-23.md`, `tb_independent_number_2026-06-23.md`,
  `amr_portal_tb_independent_runbook_2026-06-23.md`, `next_independent_label_cell_feasibility_2026-06-23.md`.
- **Data (gitignored, on D:):** `D:/dna_decode_cache/data files donwload/amr_portal/{phenotype,genotype}.parquet`
  (1.71M / 1.37M rows); `D:/dna_decode_cache/tb_indep/` (TB pipeline work dir: ref + asm + vcf + checkpoint).

## 5. Frozen-surface discipline (do NOT break)
`dna_decode/eval/amr_rules.py` + `dna_decode/data/calibrated_amr_rules.json` are sha-pinned in 4 places
(leak-guard test `tests/test_tb_leak_guard.py`, `dna_decode/eval/prospective_lock.py`,
`wiki/prospective_lock_manifest_2026-06-22.json`, `wiki/reproducibility_freeze_2026-06-13.md`). Any change to
those files = a deliberate, ratify-first freeze-amendment that **re-pins all 4 shas** (this session did one
metadata-only amendment, sha `d442b768`→`ece6744b`; predictions stayed byte-stable). The TB cell + AMR-Portal
scripts are NON-frozen and touch none of this.

## 6. State at handoff
- Last commit on origin/main: **`681aa92`** (ledger row 174). Working tree clean (only pre-existing untracked
  noise). Full test suite green (1639 passed, excl. the host-torch-paging `test_models_foundation.py`).
- Repo is **PUBLIC**: `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION` (folder `dna_decode` ≠
  repo name — that's why it wasn't visible before; it was also private until this session).
- Project ledger: `project_state/dna-decode-2026-05-11.md` (through row 174).
- Two machines: do work solo on the active machine, commit to main (the sync channel); user syncs ~weekly.
