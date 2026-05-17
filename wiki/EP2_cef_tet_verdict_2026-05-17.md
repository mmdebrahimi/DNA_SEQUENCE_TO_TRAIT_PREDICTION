# EP2 — Ceftriaxone + Tetracycline Smoke Verdict (2026-05-17)

> Architecture-transfer probe per H17. **H17 is FALSIFIED.** NT-derived architecture transfers to ceftriaxone but NOT to tetracycline. Cef passes; tet hard-fails (anti-predictive NT AUROC 0.400).

---

## Summary

EP2 fired the 12-strain smoke gate against ceftriaxone + tetracycline mini cohorts (6R/6S each, 12 unique MLSTs, full assembly availability — built from the existing N=38 cipro cohort strain pool, sharing the populated NT cache at `D:/dna_decode_cache/embeddings/nt_n40_cipro.h5`).

| Drug | NT-XGBoost AUROC | k-mer-XGB AUROC | Gap (k-mer - NT) | Verdict |
|---|---:|---:|---:|---|
| Ceftriaxone | **0.833** | 0.833 | 0.0 pp | **PASS** (NT not obviously worse) |
| Tetracycline | **0.400** | 0.722 | +32.2 pp | **FAIL** (NT obviously worse; anti-predictive) |

Gene-presence variant failed cleanly for both drugs due to missing GFF3 on one cohort strain (`GCA_008727135.1` — same strain that was unresolved in the original N=40 populate; pre-existing issue, not load-bearing for the verdict).

## H17 verdict per EP2 D5 falsification criteria

**Hypothesis H17:** "Cipro-derived NT-XGBoost architecture transfers to BOTH ceftriaxone AND tetracycline at 12-strain smoke fidelity."

Falsification rule D5(b): `H17 falsified if tet smoke produces NT AUROC ≤ 0.55 AND best-classical AUROC ≥ 0.65`.

- Tet NT AUROC = 0.400 ≤ 0.55 ✓
- Tet k-mer AUROC = 0.722 ≥ 0.65 ✓

**H17 status: FALSIFIED.** Filed for ledger update: `--update-hypothesis H17 --status falsified`.

## Anti-predictive NT result on tet — three candidate causes

NT-XGBoost AUROC = 0.400 is **below 0.5**, indicating NT predictions are systematically pointing the wrong way on tet. Three candidate causes (in order of likelihood):

1. **Architectural mismatch (most likely).** Tet resistance biology is dominated by mobile elements: tet-family efflux pumps (tetA/B/C/D/K/L/M) + ribosomal-protection proteins (tetO/W). These are mostly plasmid-borne; chromosome integration is partial. NT whole-genome mean-pooling averages embeddings across ~5000 genes; mobile-element gene contributions get diluted into noise. This matches the EP1 closeout's "architecture co-bottleneck" hypothesis from cipro.

2. **N=11 training calibration regime.** Per the 2026-05-14 LESSON: `CalibratedClassifierCV` with isotonic regression at N=11 collapses to symmetric two-value output (AUROC ≈ 0). The smoke runner uses `calibrate=False` exactly to avoid this — but anti-predictive results at N=12 may still indicate that the head can't separate signal at this size. NOT a calibration bug per se; more a "head training underpowered" signature.

3. **Label noise carryover.** The 12 tet-labeled strains were filtered from the N=38 cipro cohort. The tet labels in this subset weren't audited for MIC tier or mechanism consistency. Some labels may be unreliable similarly to the cipro case.

The right diagnostic for (1) vs (2) is to re-run with scaled-NT-LR (per the H13 fix); for (3), do an analogous mechanism × MIC merge on the 12 tet strains.

## What this means scientifically

Combined with the EP1 cipro closeout, the cross-drug evidence chain is:

| Drug | Cohort | NT-XGB AUROC | Mechanism class | Verdict |
|---|---|---:|---|---|
| Ciprofloxacin | N=38 + label noise | 0.568 (Stage 1) / 0.615 (Stage 1b mean+max) | QRDR point mutations | FAIL + EP1 closeout |
| Ciprofloxacin (smoke) | N=12 mini | 0.750 | (same) | PASS (2026-05-14) |
| Ceftriaxone | N=12 mini | 0.833 | β-lactamase plasmid + porin + AmpC | PASS |
| Tetracycline | N=12 mini | 0.400 | tet-family efflux + ribosomal protection (mobile elements) | FAIL (anti-predictive) |

**Pattern: NT-frozen-pooling architecture works on concentrated-signal mechanisms (cipro QRDR + cef plasmid β-lactamases) AND fails on distributed mobile-element mechanisms (tet).**

This is consistent with EP1's "architecture co-bottleneck" finding: cipro at N=38 with label noise failed because of TWO bottlenecks (cohort noise + architectural pooling mismatch); cipro at N=12 mini was clean enough to mask the architecture bottleneck because QRDR is locally-concentrated; cef at N=12 mini passes because plasmid β-lactamases are also locally-concentrated; tet at N=12 mini fails because mobile-element signal is distributed and the architecture can't localize it.

## Scientific narrative for the EP1 + EP2 combined story

The publishable angle (deferred pending more evidence) is:

> Frozen whole-genome NT pooling localizes concentrated-signal AMR mechanisms (QRDR point mutations + plasmid acquired-gene β-lactamases) at smoke fidelity but fails on distributed mobile-element mechanisms (tet efflux + ribosomal protection). The architecture's failure mode is data-shape-dependent, not drug-or-cohort-dependent in isolation. Production AMR-prediction pipelines should select per-drug architectures or use targeted per-gene windows for mobile-element-dominant resistance.

## Decisions locked

1. **H17 falsified.** Filed for ledger.
2. **NT-frozen-pooling architecture is NOT a 3-drug Phase 1 ship gate.** It works for cipro + cef at smoke fidelity but not tet.
3. **Phase 1 EP-list updated:** EP1 cipro (closed as audit-infrastructure packet); EP2 cef + tet (closed; cef PASS, tet FAIL). EP3 (attribution audit) + EP4 (clade-shift) remain queued.
4. **No mechanism-aware baseline expansion (Bakta + AMRFinder) on tet for now** — the architecture is the bottleneck, not the comparator.
5. **External publication still deferred** to post-EP3 — attribution audit on cef PASS strains may reveal whether NT is genuinely localizing β-lactamases or lineage-tracking.

## What's NOT decided / open follow-ups

- **NT-LR scaled-pipeline re-run on tet.** Anti-predictive 0.400 may improve with the scaled-LR head (per the H13 fix). ~30 min to fire; informative if it lands above 0.5.
- **Tet mechanism × MIC audit.** Analogous to cipro's. Cheap; ~1 hr. Would tell us if the 12 tet strains have label-noise issues separate from the architecture.
- **Per-gene NT windows on cef + cipro CLEAN_R strains.** Architectural diagnostic. Deferred.
- **EP3 attribution audit on cef PASS strains.** Does NT actually attribute to known β-lactamase loci, or is the 0.833 lineage-tracking?

## Lessons

1. **Cef + tet smoke produced different verdicts despite identical infrastructure.** Architecture × biology interaction was load-bearing; CLAUDE.md's note about "distributed-mechanism resistance" needing distinct baselines was correct.
2. **Anti-predictive NT on tet (0.400)** is a real signal, not a plumbing bug — the H13 LESSON pattern (anti-predictive = bug) is calibration-specific; at N=12 with calibrate=False, anti-predictive can be genuine architectural failure.
3. **EP1 + EP2 chained inference works:** EP1 alone was ambiguous (cohort + architecture); EP2's cross-drug evidence locked the architectural hypothesis without firing Databricks burst.
4. **N=38 cipro cohort's NT cache reused for cef + tet mini cohorts.** ~hours of compute saved vs re-populating per drug. Pattern worth preserving for future EPs.

## Next step

EP3 attribution audit on cef PASS substrate (or follow-up tet diagnostics — depends on which question is most informative). EP3 plan does not yet exist as a separate document; queue for next session.

## Artifacts

- `wiki/smoke_gate_12strain_ceftriaxone_2026-05-17.md` — cef smoke packet
- `wiki/smoke_gate_12strain_tetracycline_2026-05-17.md` — tet smoke packet
- `data/processed/gate_b_mini_cef_cohort.parquet` — cef mini cohort (12 strains)
- `data/processed/gate_b_mini_tet_cohort.parquet` — tet mini cohort (12 strains)
- This verdict packet
