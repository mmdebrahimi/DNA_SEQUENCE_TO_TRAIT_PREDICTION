# Stage 2 N=150 Prep Plan

> Resolve the three deferred Stage 2 decisions (annotation source, AMRFinderPlus integration, Databricks vs local) and ship the infrastructure needed for a Stage 2 N=150 cipro decision-gate run, so the gate runs cleanly once Stage 1 PASSes.

---

## Problem Statement

Stage 1 N=40 (running in background as of 2026-05-14 PM) is the engineering screen that gates Stage 2 spend. Per `plans/Phase2_Decision_Gate_Plan.md`, Stage 2 is the real ship gate: N=150 cohort under leave-one-Mash-clade-out CV with the Option-C threshold (≥5 pp AUROC + top-K attribution includes ≥1 of {gyrA, parC, parE}). Three decisions were deferred from Stage 1:

1. **Annotation source** — RefSeq GFF3 carried `gene=` for only ~11% of CDSs, making gene-presence INDETERMINATE_IDENTIFIER_OOV at Stage 1 (see `plans/Gene_Presence_AUROC_Bug_Fix_Plan.md` follow-up). Stage 2 needs a working gene-presence comparator. Options: Bakta re-annotation, Roary pan-genome clustering, accept-degenerate.
2. **AMRFinderPlus POINT* SNP-table baseline** — deferred from Stage 1 per Phase2_Decision_Gate D6. For ciprofloxacin, gyrA/parC/parE point mutations are the textbook resistance signal; absence of this baseline at Stage 1 means "best classical" is currently bounded by k-mer alone. At Stage 2 this becomes load-bearing — if NT can't beat the textbook SNP table, the foundation-model premise weakens.
3. **Compute: Databricks burst vs local** — local GTX 860M observed at ~17-19 min/strain in current populate; N=150 populate would take ~45-50 hours locally. User has Databricks access. Phase2_Decision_Gate_Plan referenced a "Databricks burst" mode.

This plan resolves all three and lays out the infrastructure steps so that when Stage 1 lands with `stage2_action == BURST_STAGE_2`, Stage 2 can launch without further deliberation.

## Design Decisions

### D1: Annotation source = Bakta re-annotation (defer Roary; accept-degenerate as fallback)

**Decision:** Re-annotate the N=150 cohort with Bakta. Bakta produces cross-strain stable gene symbols + locus tags, gene-presence vocabulary overlaps materially across strains (typical Bakta gene-symbol coverage on E. coli is 60-80% vs RefSeq's ~11%), and matches the original Phase 1 design (foundation models reference Bakta annotations).

**Rationale:** RefSeq's 11% gene-symbol coverage is structurally too low — the Stage 1 INDETERMINATE_IDENTIFIER_OOV result is reproducible and the gene-presence variant can't compete. Bakta's per-strain re-annotation gives both gene_symbol (cross-strain) AND locus_tag (strain-stable hierarchy). Roary adds pan-genome orthogroup clustering on top, which would give finer-grained groups, but at N=150 the Bakta gene_symbol vocabulary is already wide enough to support meaningful gene-presence features without the additional Roary clustering pass.

**Trade-off:** Bakta runtime ≈ 5-15 min per strain CPU-bound on this hardware. N=150 = ~12-37 hours wallclock. Acceptable for a one-time prep step. Considered Roary (deferred — slower at N=150, marginal gain) and accept-degenerate (rejected — Stage 2's `PIVOT_TO_BAKTA` action assumes Bakta is available; pre-investing in the annotation is the right move regardless of Stage 1 outcome).

**Trigger:** schedule Bakta re-annotation NOW (in parallel with the populate that's running) so it's ready when Stage 2 starts. CPU-bound; doesn't compete with NT populate's GPU.

### D2: AMRFinderPlus POINT* SNP-table baseline = in scope for Stage 2

**Decision:** Install AMRFinderPlus locally, run it per-strain with `--organism Escherichia --mutation_all`, parse POINT-method rows (Method = POINT / POINTP / POINTX / POINTN), build a (gene, codon-position, alt-residue) binary feature matrix, train XGBoost on it as a 5th gate-bearing Stage 2 variant.

**Rationale:** For ciprofloxacin, gyrA + parC + parE point mutations are the textbook resistance signal. If NT can't beat this baseline at N=150, the foundation-model premise fails — that's exactly the comparator Stage 2's Option-C threshold is designed to gate against. Codex's prior /research synthesis flagged this as the most important missing baseline. AMRFinderPlus is a single CLI tool with a permissive license; integration is straightforward (~1-5 hours wallclock at N=150 once installed).

**Trade-off:** AMRFinderPlus install on Windows can be finicky (the canonical install is Linux/WSL via conda); if Windows-native install fails, run via WSL2 or skip and use local annotation extraction. Considered defer-to-Stage-3 (rejected — Stage 2 is the real ship gate; this baseline must be present).

### D3: Compute = Databricks burst for N=150 NT populate; local for everything else

**Decision:** Use Databricks burst (A100 or equivalent) for the N=150 NT-v2-100M cache populate. Once the HDF5 is built, transfer it back to local; run Stage 2 LOSO + clade-out CV + attribution on local CPU. Bakta re-annotation also runs locally (CPU-only; can parallelize with the Databricks populate).

**Rationale:** Local GTX 860M at observed ~17-19 min/strain gives ~45-50 hours for N=150 populate alone — untenable for any iteration. A100 should do it in ~3-5 hours. The classifier heads (XGBoost, logreg) and bootstrap CI computation are CPU-bound and run fast on local hardware once embeddings are cached. Hybrid (cloud populate + local analysis) minimizes cloud spend while removing the wallclock bottleneck.

**Trade-off:** Adds Databricks setup overhead (cluster spin-up, auth, dependency install, HDF5 transfer). Estimate ~2-3 hours one-time setup. Considered local-only (rejected — 45-hour iteration cycle kills the project) and full-Databricks (rejected — analysis compute is small + Stage 2 will iterate on classifier head choices).

### D4: Stage 2 cohort = N=150 expanded from gate_b_cohort.parquet (67 strains) via the audit-cohort pipeline

**Decision:** Use `scripts/audit_cohort.py` (existing) to expand the cohort from the current 67-strain audited base to N=150 by relaxing the strict assembly-quality thresholds (allowing N50 ≥ 100 K, contig_count ≤ 50) and selecting from the BV-BRC AST TSV. Target balance: 75R / 75S after broth-microdilution filter.

**Rationale:** The current `gate_b_cohort.parquet` has 67 strains (61 unique MLSTs, high diversity). Audited cohort generator already exists and supports the necessary filters. Expanding from 67 to 150 should be feasible from the BV-BRC AST TSV; rough estimate from prior pilot-gate observations is N>200 cipro strains pass the broader filters, so 150 is achievable.

**Trade-off:** Loosening quality thresholds risks degrading Mash-clade-out CV quality. Acceptable because Stage 2 explicitly uses leave-one-Mash-clade-out (not LOSO) — coarse-grained clade boundaries are more robust to per-strain assembly noise than strain-level LOSO would be. Considered hand-curating 150 strains (rejected — too time-consuming; the pipeline exists for a reason).

## Implementation Plan

**Revision 2026-05-14 PM (post-/brainstorm Round 1):** 3 critical issues incorporated:
- Phase A.0 added (cohort builder script — `audit_cohort.py` only audits, doesn't build; `pipeline ingest` CLI doesn't expose N50/contig thresholds)
- Phase A.4 added (Bakta CDS-row `gene_symbol` validation before overnight scale-out — Bakta typically puts `gene=` on parent `gene` rows; CDS rows may lack it; would re-enter `INDETERMINATE_IDENTIFIER_OOV` despite re-annotation)
- Phase A.5 added (Mash preflight + threshold reconciliation — existing default is 0.02 / ~98% identity, plan originally said 99%; conformed to default)
- Threshold-relaxation direction corrected (repo defaults are N50≥50K + contig_count≤500; the plan's values were stricter, not looser)

**Revision 2026-05-14 PM (post-A.0 execution):** N=150 target revised to N=147 actual:
- Phase A.0 executed: cohort builder shipped + ran. After all filters (broth-microdilution + MLST + N50/contig + assembly_accession): R-ceiling is 72 strains. Took ALL 72 R + 75 S (MLST-balanced within each class) = **N=147 actual, 72R/75S, balance ±3 of target 75/75**.
- Diagnostic (`scripts/diagnose_bvbrc_mlst_gaps.py`) revealed the bottleneck: 35,790 of 85,114 BV-BRC genome rows lack `assembly_accession` (un-downloadable from NCBI Datasets). Among cipro AST: 491 R lack accession, 77 retain it. Of those 77, 72 pass assembly-quality + MLST filters.
- LABEL-STRATIFIED MLST-balanced selection was load-bearing: the default `build_cohort` algorithm prioritized MLST diversity without R/S stratification, leaving available R strains on the table (49R/101S ceiling with default algo vs 72R/75S with label-stratified). The builder script reuses `_mlst_balanced_selection` per-class.
- N=147 vs target N=150 is a 2% shortfall — well within engineering-screen tolerance. Stage 2 ≥5 pp Option-C threshold remains unchanged.

**Crash event 2026-05-14 (~16:03):** First N=40 NT populate (Stage 1 prerequisite, running since 12:30) crashed at ~8 hr wallclock with HDF5 file corruption (errno=13 Permission denied; superblock EOA truncated to 6472 bytes despite 132 MB of on-disk data). Root cause likely Seagate Portable USB hiccup (known failure mode per `GATE_B_REPORT.md` 2026-05-13). H5 file unrecoverable; populate restarted from scratch as of 2026-05-14 PM. Stage 1 verdict deferred by ~10 hr.

### Phase A (parallel with populate; CPU-only, no D: drive write contention)

0. **Write `scripts/build_stage2_n150_cohort.py` (NEW; cohort builder).**
   - Calls `dna_decode.data.cohort.candidates_from_bvbrc_ast(...)` + `build_cohort(...)` + `save_cohort(...)` directly.
   - CLI flags: `--ast-tsv`, `--assembly-metadata-csv`, `--drug ciprofloxacin`, `--target-total 150`, `--per-class 75`, `--n50-min 50000`, `--contig-count-max 500` (repo defaults; loosen only if N=150 isn't achievable).
   - Hard-fail if balance imbalanced beyond ±10 (e.g., 70R/80S OK, 60R/90S NOT).
   - Hard-fail if any strain lacks MLST (matches Stage 1 invariant).
   - Output: `data/processed/stage2_n150_cipro_cohort.parquet`.

1. **Install Bakta locally.**
   - `pip install bakta` OR `conda install -c bioconda bakta` if conda is preferred.
   - Download Bakta database: `bakta_db download --output D:/bakta_db --type full` (or `--type light` if disk is tight). Database is ~30-50 GB; D: drive has room.
   - Smoke-test: `bakta --db D:/bakta_db --output /tmp/bakta_test data/cache/refseq/GCF_000005845.2/genome.fna` on K-12 MG1655. Expect a `.gff3` output with `gene=` annotations populated.

2. **Install AMRFinderPlus locally.**
   - Try `pip install ncbi-amrfinderplus` first. If Windows-native fails:
     - Option A: install via WSL2 (Ubuntu); call from main project via `wsl amrfinder ...`.
     - Option B: download standalone binary from NCBI releases.
   - Download AMR database: `amrfinder --update`.
   - Smoke-test: `amrfinder -n data/cache/refseq/GCF_000005845.2/genome.fna --organism Escherichia --mutation_all`. Expect TSV output with POINT-method rows.

3. **Run the cohort builder from A.0.**
   - `uv run python scripts/build_stage2_n150_cohort.py --drug ciprofloxacin --target-total 150 --per-class 75 ...`
   - Verifies all strains have MLST (loud-fail per Stage 1 invariant).
   - Outputs `data/processed/stage2_n150_cipro_cohort.parquet`.
   - If balance can't hit 75R/75S ±10 from broth-microdilution-filtered BV-BRC pool: report the achievable balance + escalate to user (do not silently accept worse imbalance).

4. **Bakta CDS-row `gene_symbol` validation (BEFORE overnight scale-out).**
   - Run Bakta on 2-3 strains: K-12 MG1655 (`GCF_000005845.2`), one ST131 representative, one ST10 representative.
   - Parse outputs with `dna_decode.data.annotations.parse_gff3`.
   - Assert: fraction of CDS rows with non-empty `gene_symbol` ≥ 50% per strain. (Bakta with `--prodigal-tf` etc. should hit 60-80%.)
   - **If <50% CDS-row gene_symbol coverage:** Bakta puts `gene=` on parent `gene` rows, not CDS rows. Ship `parse_gff3` extension first — for CDS rows with empty `gene_symbol` but populated `Parent=` attribute, look up the parent gene row and inherit its `gene_symbol`. Then re-run validation. Only proceed to overnight Bakta scale-out after validation passes.

5. **Mash preflight (cohort clade structure).**
   - Run `dna_decode.eval.phylogeny.compute_mash_distances` + `cluster_at_threshold` against the N=150 cohort genomes (FASTAs only; no annotations needed).
   - Use the repo default threshold = 0.02 (~98% identity). This matches existing code; the plan's original "ANI ≥ 99%" was inconsistent with the default.
   - Output: `data/processed/stage2_n150_clades.parquet` with strain_id → clade_id mapping.
   - Report: total fold count, largest-fold strain count, distribution of folds with pure-R / pure-S / mixed labels. Loud-fail any strain that ends up `__unassigned__` (current `leave_one_clade_out_cv` silently buckets these — need to either patch the CV helper to raise OR pre-validate every strain has a clade).
   - **Gate:** if fold count < 8 OR any fold is pure-R/pure-S (degenerate AUROC), escalate to user before Phase B kickoff — Stage 2 CV may need threshold loosening or cohort tweaking.

### Phase B (after Stage 1 PASSes; only if `stage2_action == BURST_STAGE_2`)

4. **Bakta re-annotation of N=150 cohort.**
   - Driver script: `scripts/bakta_annotate_cohort.py` (NEW). Iterates cohort, calls Bakta CLI per strain, writes outputs to `D:/dna_decode_cache/bakta_annotations/<accession>/<...>.gff3`.
   - Expected runtime: 12-37 hours CPU-only on this hardware. Can run unattended overnight.
   - Validation: parse each output with `dna_decode.data.annotations.parse_gff3`; assert ≥50% of CDSs have non-empty `gene_symbol`.

5. **AMRFinderPlus POINT* extraction for N=150 cohort.**
   - Driver script: `scripts/amrfinder_extract_cohort.py` (NEW). Iterates cohort, calls amrfinder CLI per strain, parses POINT-method rows from TSV output, builds a `(strain_id, gene, codon_pos, alt_residue) → True` table.
   - Materialize as `data/processed/stage2_n150_amrfinder_point_mutations.parquet`.
   - Expected runtime: 1-5 hours.

6. **Databricks burst: N=150 NT cache populate.**
   - Spin up Databricks cluster (DBR-X with A100 if available, else V100).
   - Install dna_decode + dependencies via `pip install -e .` from the synced repo.
   - Transfer cohort parquet + refseq cache to Databricks via DBFS or workspace upload.
   - Run `scripts/populate_cache.py --cohort stage2_n150 --model nucleotide_transformer --cache /tmp/nt_n150_cipro.h5 --device cuda`. Expected ~3-5 hours on A100.
   - Transfer HDF5 back to local (`D:/dna_decode_cache/embeddings/nt_n150_cipro.h5`) for the actual Stage 2 analysis.
   - Verify cache integrity: open with `EmbeddingCache`, list strains, confirm count == 150.

### Phase C (Stage 2 decision gate, fully local on populated cache)

7. **Write Stage 2 runner (`scripts/stage2_n150_cipro.py`, NEW).**
   - Modeled after `scripts/stage1_n40_cipro.py` but with 5 gate-bearing variants:
     - NT-XGBoost (primary)
     - NT-logreg (sanity)
     - k-mer-XGB (classical)
     - Bakta-gene-presence-XGB (classical; NEW for Stage 2)
     - AMRFinderPlus-POINT-table-XGB (knowledge baseline; NEW for Stage 2)
   - Diagnostic-only: NT+k-mer-fusion-logreg (carried over from Stage 1).
   - **CV strategy:** `leave_one_clade_out_cv` from `dna_decode/eval/cv.py` (Mash-clade-out, NOT LOSO). At N=150 with diverse MLSTs, expect 10-20 clades.
   - **Gate formula:** `max(NT-XGBoost, NT-logreg) AUROC ≥ best_classical_AUROC + 5 pp` AND top-K attribution includes ≥1 of {gyrA, parC, parE}. This is Option-C — tighter than Stage 1's ≥3 pp.
   - **Attribution:** run `scripts/pipeline.py attribute` per strain (Tier 1-5 ISM); aggregate top-K across the cohort; check gyrA/parC/parE presence.
   - **Result packet:** `wiki/stage2_n150_cipro_<date>.md` with all 5 variant AUROCs, gate analysis, per-clade table, paired bootstrap CI, attribution table, Option-C verdict.

8. **(Mash preflight already done in Phase A.5.)** Reuse `data/processed/stage2_n150_clades.parquet` produced earlier. Pass `clade_assignments` dict directly to `leave_one_clade_out_cv`.

9. **Execute Stage 2 runner.**
   - `HF_HOME=D:/hf_cache uv run python scripts/stage2_n150_cipro.py`.
   - Expected wallclock: ~2-4 hours (4 variants × clade-out folds × classifier training, plus attribution pass).
   - Verdict feeds into the Phase 2 ship decision: PASS → Phase 1 ships per the original technical plan's Verification section; FAIL → demote NT track.

### Phase D (epilogue, regardless of verdict)

10. **Update project ledger.**
    - Stage 2 outcome → evidence rows in `project_state/dna-decode-2026-05-11.md`.
    - Resolve H12 (NT v2 100M frozen + XGBoost at N=150 LOSO Mash-clade-out gate) + H14 (SNP-table baseline) + H8 partial (phylogeny-permutation control).
    - Decisions Made entry capturing Phase 2 final verdict.

11. **Documentation refresh.**
    - Update `CLAUDE.md` Common Commands with Stage 2 runner.
    - Append Phase 2 closeout notes to `LESSONS_LEARNED.md`.
    - Update `TODOS.md` — close Stage 2 entries, open Phase 3 prep entries (multi-organism, multi-drug expansion, attribution refinement, MIC regression Phase 2 → Phase 3).

## Open Questions / Pending Decisions

These are not load-bearing for Phase A (Bakta install + AMRFinderPlus install + cohort expansion can proceed unambiguously), but need resolution before Phase B:

- **Databricks cluster sizing:** A100 vs V100 vs CPU-only. Recommendation: A100 if available, fall back to V100. CPU-only would defeat the burst-purpose.
- **Bakta DB type:** `full` (~30-50 GB) vs `light` (~5 GB). Light is fine for our use case (gene-symbol assignment); full gives more product annotations we don't currently use. Recommendation: light.
- **Stage 2 cohort N target:** strict 150 or accept whatever balance the audit pipeline produces if it falls short. Recommendation: 150 ±10. If under 140 after filters, raise the assembly-quality threshold lower; if under 130, escalate.

## Verification

- Phase A complete when: `bakta --version` works locally, `amrfinder --version` works locally (or via WSL), `data/processed/stage2_n150_cipro_cohort.parquet` exists with ≥140 strains and balanced R/S split.
- Phase B complete when: Bakta GFF3 outputs exist for all 150 strains with ≥50% gene_symbol coverage; AMRFinderPlus POINT-table parquet exists; Databricks HDF5 transferred to local with 150 strain entries.
- Phase C complete when: `wiki/stage2_n150_cipro_<date>.md` packet exists with all 5 variants, Mash-clade-out CV results, attribution table, Option-C verdict (PASS / FAIL / borderline).
- Phase D complete when: project ledger has Stage 2 evidence row + final decision row; docs are updated; TODOs are reorganized.

## Risk Flags

- **Bakta DB download is large** (~5-50 GB depending on type). Confirm D: drive has room before kickoff (currently ~4.5 TB free on Seagate Portable — fine).
- **AMRFinderPlus on Windows is finicky.** WSL2 fallback exists; if neither works, Phase B step 5 stalls. Mitigation: validate amrfinder install in Phase A step 2 — if it fails, escalate to user before Phase B kickoff.
- **Databricks setup is the wallclock bottleneck for Phase B.** Cluster spin-up + auth + dependency install can eat 2-3 hours. Validate this with a smoke run on a tiny cohort (e.g., N=5) before committing to the N=150 burst.
- **Stage 2 cohort expansion may not hit 150.** BV-BRC's broth-microdilution-filtered cipro cohort empirically has N>200 candidates; loosening assembly-quality further should yield 150 cleanly. But this is unverified at the relaxed thresholds.
- **Mash-clade-out CV at N=150 may produce too few clades.** Target ≥10. If the cohort collapses to <10 clades, the CV is degenerate (similar to the N=12 clade-only-baseline trap). Pre-flight check at Phase C step 8 mitigates this.
- **Phase B step 6 (Databricks burst) is gated on user authorization for Databricks spend.** Cost estimate at 4-hour A100 burst: ~$15-30 depending on region. Confirm with user before kickoff.

## Time + cost budget

| Phase | Wallclock | Cost | Compute |
|---|---|---|---|
| A: Bakta/AMRFinder install + cohort expansion | 2-4 hr | $0 | Local CPU |
| B.4: Bakta re-annotation (N=150) | 12-37 hr | $0 | Local CPU (overnight) |
| B.5: AMRFinderPlus extraction (N=150) | 1-5 hr | $0 | Local CPU |
| B.6: Databricks NT populate (N=150) | 3-5 hr | $15-30 | Databricks A100 |
| C: Stage 2 runner end-to-end | 2-4 hr | $0 | Local CPU |
| D: Docs + ledger | 1-2 hr | $0 | Local |
| **Total** | **21-57 hr** (mostly Bakta CPU overnight) | **$15-30** | Hybrid |

Most of the wallclock is Bakta (CPU-bound, parallelizable, overnight). Active human time is ~6-10 hr spread across setup + monitoring.
