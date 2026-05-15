# EP2 Cef + Tet Smoke Design Plan

> Evidence Packet 2 per Bellman target-state framing (2026-05-15 reset). Parallel-fire data-shape probe at 12-strain mini-cohort for ceftriaxone + tetracycline. Tests H17 (cipro-derived NT-XGBoost architecture transfer assumption). Closes Pending Decisions row 9.

---

## Problem Statement

EP1 (cipro Stage 2 N=147/N=150) cannot be the only architecture-validation gate before Stage 2 burst spend. The cipro smoke gate PASSED 2026-05-14 with NT-XGBoost AUROC 0.750 vs k-mer 0.694 (+5.6 pp) — strong evidence for QRDR-like point-mutation signal but NOT for the architecture's transfer to distributed-mechanism resistance.

**Drug mechanism heterogeneity (CLAUDE.md gotchas + literature):**
- Ciprofloxacin: QRDR point mutations (gyrA-S83L/D87N + parC-S80I/E84V) — concentrated CDS-local signal. K-mer + NT both work because the mutations sit in conserved single-gene contexts.
- Ceftriaxone: distributed across (1) plasmid β-lactamases (blaCTX-M-15 family, blaSHV, blaTEM-extended), (2) chromosomal AmpC over-expression, (3) porin loss (ompC/ompF), (4) efflux pumps (acrAB-TolC). Signal is plasmid-context + chromosomal + regulatory — k-mer may catch β-lactamase ORFs but miss regulatory + porin signal.
- Tetracycline: tet-family efflux pumps (tetA, tetB, tetC, tetD, tetK, tetL, tetM) + ribosomal-protection proteins (tetO, tetW). Most signal is mobile-element (plasmid + transposon); chromosome integration is partial.

**Success criteria for the design plan (this document):**
- Per-drug cohort selection mechanism is unambiguous (which AST labels; class balance constraints; MLST diversity).
- Smoke runner is drug-generic at the output-path layer (currently cipro-hardcoded).
- Baseline class per drug is mechanism-appropriate (k-mer alone may be inadequate for cef/tet).
- Pass/fail threshold + falsification criteria are written before runs fire.

---

## Design Decisions

### D1: Per-drug cohort = 12 strains, 6R/6S, MLST-diverse — same shape as cipro mini

**Decision:** Use `scripts/build_mini_cohort.py --drug ceftriaxone --per-class 6` and `--drug tetracycline --per-class 6` to produce two new parquets: `data/processed/gate_b_mini_cef_cohort.parquet` + `data/processed/gate_b_mini_tet_cohort.parquet`.

**Rationale:** EP2 is a falsification-tier probe, not a powered comparison. 12 strains × 6R/6S is the same noise-floor regime as cipro mini (±8.3% per-strain; ±0.19 95% CI on AUROC). Falsifying H17 doesn't require statistical power; it requires data-shape signal that's interpretable.

**Trade-off:** Some drug × MLST combinations may not have 6R/6S available in BV-BRC AST. Mitigate: relax to 5R/5S minimum; raise to 7-8 per class if AST coverage permits. Class imbalance >1.5× → defer EP2 for that drug + log as data-bottleneck finding.

### D2: Smoke runner generalized at output-path layer

**Decision:** Rename `scripts/smoke_gate_12strain_cipro.py` → `scripts/smoke_gate_12strain.py`. Make output filename + report heading + cohort description templated on `--drug`:
- `wiki/smoke_gate_12strain_<drug>_<date>.md` (was: `wiki/smoke_gate_12strain_cipro_<date>.md` hardcoded at line 341)
- `# Smoke Gate — 12-strain <drug> cohort` (was: cipro-hardcoded at line 264)
- `**Cohort:** ... (12 strains, 6R/6S <drug>)` (was: cipro-hardcoded at line 272)
- `**Drug:** <drug>` (was: hardcoded at line 273-274)

Internal logic already accepts `--drug` (lines 75, 128, 165 use `drug_lower` from the arg). Only output strings need parameterization.

**Rationale:** Codex G3 critique — current script reads "cipro-locked" from the OUTSIDE but is drug-generic on the inside. Generalization is cosmetic + naming, not logic.

**Trade-off:** The rename breaks `wiki/smoke_gate_12strain_cipro_2026-05-14.md` references in docs/CLAUDE.md/Mid-term-table. Mitigate: keep the existing artifact at its existing path; rename script + future-output convention; backfill no historical files.

### D3: Mechanism-appropriate baseline class per drug

**Decision:**
- **Cipro:** k-mer (current) + gene-presence (now non-degenerate via Bakta annotations once EP2 cef/tet are wired). K-mer is good for QRDR point mutations; gene-presence less so but informative.
- **Cef:** k-mer + gene-presence (Bakta-annotated; should detect blaCTX-M / blaSHV / blaTEM family genes when present) + AMRFinderPlus POINT* (covers β-lactamase POINT-mutation calls). Three-baseline matrix; pick the strongest classical per fold for the gap calc.
- **Tet:** k-mer + gene-presence (Bakta should detect tetA/tetB/tetM) + AMRFinderPlus tet-family calls. Same three-baseline matrix.

**Rationale:** "Best classical" framing per EP framework requires mechanism-appropriate comparator. K-mer alone is the wrong default for distributed-signal mechanisms; the smoke gate should report each baseline's individual AUROC + the max for the gap calc.

**Trade-off:** Requires AMRFinderPlus + Bakta wired into the smoke runner. AMRFinderPlus DB + Bakta DB are installed locally as of 2026-05-15 (per `wiki/stage2_install_artifact_2026-05-15.md`). Bakta annotation per strain is CPU-heavy ~5-30 min; for 12-strain × 2-drug = 24 annotations = potentially 12 hours of CPU. Cache annotations on first run; reuse on subsequent runs.

### D4: Pass/fail threshold = same 15 pp engineering heuristic, mechanism-aware caveat

**Decision:** Per-drug smoke gate verdict:
- PASS: NT-XGBoost AUROC ≥ max(k-mer, gene-presence, AMRFinder-POINT) AUROC − 15 pp.
- FAIL: NT-XGBoost AUROC < best-classical − 15 pp.
- INDETERMINATE: any baseline returns INDETERMINATE_IDENTIFIER_OOV (gene-presence on low-symbol-coverage strain set) — re-run after Bakta re-annotation.

**Rationale:** 15 pp threshold for cipro was already a loose engineering smoke bar, NOT statistically powered. Same threshold for cef/tet maintains comparability. Mechanism-aware caveat: NT is allowed to be weaker than best-classical at 15 pp on distributed-signal drugs because the architecture was developed against concentrated-signal cipro; failure to meet the 15 pp bar is informational (data-shape divergence), not a hard demote signal.

**Trade-off:** Codex Open Q: should cef/tet have mechanism-specific thresholds? Defensible alternative: 20 pp for cef (more headroom for distributed signal), 10 pp for tet (tet-family detection should be cleaner). Defer to empirical EP2 run; revisit threshold if either smoke produces a >15 pp gap.

### D5: Falsification criteria for H17 (NOT same as pass/fail)

**Decision:** H17 (cipro-derived architecture transfers to BOTH cef AND tet) is falsified if EITHER:
- (a) Cef smoke produces NT-XGBoost AUROC ≤ 0.55 (= no detectable signal above chance) AND best-classical AUROC ≥ 0.65 (= classical CAN detect signal), OR
- (b) Tet smoke produces same shape.
- (c) Inverted: NT-XGBoost AUROC < (1 − 0.65) on either drug — anti-predictive signal indicating plumbing bug, not biological signal (per anti-predictive-AUROC LESSON_LEARNED 2026-05-14).

**Rationale:** H17 is about ARCHITECTURE transfer (does NT-XGBoost pattern work?), NOT about NT being best. Falsification = NT cannot pick up signal in a regime where classical CAN. Confirmation = NT picks up signal at the mechanism's expected magnitude (≥0.65 AUROC). The 15 pp bar (D4) is a smoke heuristic; H17 falsification is a sharper biological-signal test.

**Trade-off:** Setting the falsification bar at AUROC ≤ 0.55 vs ≤ 0.60 is a judgment call. Conservative (≤0.55) avoids false H17 falsification at small N; aggressive (≤0.60) catches weak-signal cases earlier. Default: 0.55 + revisit if smoke produces 0.55-0.60 ambiguous results.

### D6: Cohort source = BV-BRC AST + assembly availability (same parser as cipro)

**Decision:** Use existing BV-BRC AST TSV + assembly metadata pipeline (`scripts/build_mini_cohort.py`) for cef + tet. Same quality filters (contig_count ≤500, N50 ≥50K, MLST diversity).

**Pre-flight verification (before EP2 fires):** Run `build_mini_cohort.py --drug ceftriaxone --per-class 6 --dry-run` (if --dry-run exists; else write to /tmp + inspect). Confirm:
- 12 strains available with 6R/6S balance.
- ≥10 distinct MLST values (avoid singleton-clade degeneracy that broke H11 at cipro mini).
- All 12 have downloadable `assembly_accession` per BV-BRC's known bottleneck (LESSON_LEARNED 2026-05-14).

Same pre-flight for `--drug tetracycline`.

**Rationale:** The cipro mini cohort selection process is empirically validated (12-strain smoke PASSED 2026-05-14). Reuse, don't redesign.

**Trade-off:** If pre-flight reveals cef or tet AST coverage is <12 R+S strains with downloadable accessions, EP2 for that drug is blocked. Mitigate: relax to 10 strains per drug + log the data-bottleneck. If <10 available, raise as a Pending Decision (cef/tet cohort recovery vs PIVOT_TO_DIFFERENT_DRUG).

---

## Implementation Plan

### Step 1: Pre-flight cohort availability (no code changes)
- Run `build_mini_cohort.py --drug ceftriaxone --per-class 6` against BV-BRC AST.
- Run `build_mini_cohort.py --drug tetracycline --per-class 6`.
- Inspect output parquet for: row count = 12, label balance 6R/6S, MLST diversity ≥10, accession-availability 12/12.
- Branch on result: PASS → continue to Step 2; FAIL → write blocker report.

### Step 2: Generalize smoke runner output strings
- Rename `scripts/smoke_gate_12strain_cipro.py` → `scripts/smoke_gate_12strain.py`.
- Replace cipro-hardcoded strings at lines 264, 272, 273-274, 333-341 with `--drug`-templated equivalents.
- Default output path: `wiki/smoke_gate_12strain_<drug>_<date>.md`.
- +2 unit tests pinning the output-path templating.

### Step 3: Wire AMRFinderPlus + Bakta into smoke runner (mechanism-appropriate baselines per D3)
- Add `--include-amrfinder` flag → calls AMRFinderPlus via `tools/docker_runner.run` per smoke strain; parses POINT-method rows per drug-specific target manifest.
- Add `--include-bakta-gene-presence` flag → calls Bakta annotation (cached) → uses gene-symbol-based gene-presence matrix (per Gene_Presence_AUROC_Bug_Fix_Plan.md fix).
- Per-drug target manifest (configurable, defaults baked in):
  - cipro: gyrA, parC, parE, qnr-family
  - cef: blaCTX-M, blaSHV, blaTEM, ompC, ompF, ampC
  - tet: tetA, tetB, tetC, tetD, tetK, tetL, tetM, tetO, tetW

### Step 4: Run EP2 — cef smoke
- `HF_HOME=D:/hf_cache uv run python scripts/smoke_gate_12strain.py --drug ceftriaxone --include-amrfinder --include-bakta-gene-presence`
- Output: `wiki/smoke_gate_12strain_cef_<date>.md`.
- Verdict: PASS / FAIL / INDETERMINATE per D4.

### Step 5: Run EP2 — tet smoke (parallel-eligible with Step 4)
- Same invocation with `--drug tetracycline`.
- Output: `wiki/smoke_gate_12strain_tet_<date>.md`.

### Step 6: H17 verdict
- Apply D5 falsification criteria across both packets.
- Append `--update-hypothesis H17 --status <falsified | confirmed | under-investigation>` to ledger.
- If falsified: file follow-up Pending Decision on whether to demote cef/tet from Phase 1 ship gate (was 3-drug; falls to cipro-only if H17 falsified).

### Step 7: EP2 result packet write-up
- Single ledger-adjacent doc: `wiki/EP2_cef_tet_verdict_<date>.md`.
- Tables: per-drug AUROC × variant; H17 status; data-shape divergence findings.
- Surface follow-up decisions: EP5 readiness; PIVOT_TO_BAKTA / ALTERNATIVE_POOLING per-drug; Stage 2 cef/tet readiness.

---

## Verification

- `build_mini_cohort.py --drug` produces valid cef + tet parquets with 6R/6S × 12-strain × ≥10-MLST × 12/12-accession.
- `scripts/smoke_gate_12strain.py --drug X` runs end-to-end on each of cipro/cef/tet (3-drug coverage at smoke fidelity).
- Output filenames + report headings vary by `--drug` arg; no cipro-hardcoded strings remain.
- AMRFinderPlus + Bakta integration produces per-strain POINT + gene-presence matrices; smoke runner parses + aggregates.
- H17 status updated per D5 criteria after both smokes complete.
- All existing tests pass: `uv run pytest tests/ -m "not slow" -q` returns prior baseline + new smoke-runner tests.

---

## Open Questions

1. **Should EP2 PASS require BOTH cef + tet, or is one-of-two acceptable?** Conservative read of H17 ("BOTH ceftriaxone AND tetracycline") → both must pass. Pragmatic read: one-of-two is informative + lets the other drug pivot to BAKTA annotation route. Lean: conservative; H17 wording binds.

2. **Bakta annotation cache strategy.** Per-strain Bakta runs are ~5-30 min CPU; 24 strains × 2 drugs = potentially 12 hr cumulative. Cache strategy: (a) cache by `(strain_id, bakta_version, db_version)`; (b) one-shot annotation pass before any smoke; (c) lazy + cached per smoke run. Lean: (b) — explicit pre-flight annotation phase, then smoke runs are fast.

3. **AMRFinderPlus POINT* parser scope.** Existing parsing in EP1 cipro path covers QRDR point mutations. For cef/tet, AMRFinderPlus also calls plasmid β-lactamases + tet-family genes via gene-presence detection, not POINT — different rows in the same `mutations.tsv` + `main.tsv`. Parser may need to read both files per drug.

4. **Mechanism-specific threshold revisit.** D4 defaulted to 15 pp uniform. If cef smoke produces 5-15 pp gap, does that constitute weak-PASS or FAIL? Defer to empirical run; revisit threshold only if smoke results sit in the ambiguous zone.

---

## Status

- **Written:** 2026-05-15 (during dwell time for Stage 1 verdict).
- **Pre-requisite:** Stage 1 verdict (informs whether to fire EP2 immediately on PASS, or pivot per FAIL-branch decisions).
- **Estimated lift:** Steps 1-2 = ~2 hours; Step 3 = ~4-6 hours (Bakta + AMRFinder integration); Steps 4-5 = 12 hours CPU + waiting; Step 6-7 = ~1 hour.
- **Blocker B1 status:** Bakta DB installed + verified 2026-05-15 (4.0 GB at `C:/Users/Farshad/dna_decode_stage2/bakta_db/db-light/`); workaround = `--entrypoint /bin/bash -c "..."` for conda-init.
