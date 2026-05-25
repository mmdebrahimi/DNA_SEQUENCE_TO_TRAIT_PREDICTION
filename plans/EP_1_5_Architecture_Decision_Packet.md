# EP-1.5 Architecture Decision Packet — Distributed-Mechanism Resistance

> Pre-staged decision packet for the architectural fork that gates EP-2 multi-drug expansion. Candidate architectures + scoring rubric + proof-of-concept protocol. Final pick lives in the same file once Codex runs the small-N POCs on Precision 7780.

**Status:** DRAFT 2026-05-25. POC runs are gated on EP-0 close (need GPU + populated NT cache); decision lands after POC. The decision packet itself is what fires `plans/Post_V0_EP_Ladder_Plan.md` EP-1.5.

**Anchor docs:**
- `wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md` (the cross-drug pattern that triggers EP-1.5)
- `plans/Post_V0_EP_Ladder_Plan.md` D3 (multi-drug forks by mechanism class; EP-1.5 unblocks distributed-mechanism drugs)

---

## Problem Statement

The 2026-05-17 cross-drug architectural finding established that NT v2 100M frozen + mean-pooled embeddings + XGBoost **fails on distributed mobile-element mechanisms**:

| Drug | Mechanism class | NT-XGBoost smoke AUROC | k-mer baseline | Verdict |
|---|---|---:|---:|---|
| Cipro | QRDR point mutations (concentrated) | 0.750 | 0.694 | PASS |
| Cef | Plasmid β-lactamases (concentrated) | 0.833 | 0.833 | PASS |
| Tet | tet-family efflux + ribosomal protection (distributed) | **0.400 anti-predictive** | 0.722 | FAIL |

The proposed candidate dilution mechanism: whole-genome mean-pool contributes ~1/N_genes (N≈5000) per gene; localized signal survives this dilution because the predictive feature is confined to a small index set; distributed signal does not.

EP-2 ships multi-drug AMR — **but cannot ship for distributed-mechanism drugs (tet, plus probably gent + others) without architectural change**. EP-1.5 picks which change.

**Decision required:** which architectural fix is selected for distributed-mechanism resistance?

**Decision constraint:** must be implementable in frozen-inference mode on Precision 7780 (RTX 3500 Ada, 12 GB VRAM). NO fine-tuning of NT v2 100M (the GTX 860M has CC=5.0 < 7.0 required for bitsandbytes; the Precision 7780's RTX 3500 Ada is Ada Lovelace CC=8.9 which DOES support bitsandbytes — but fine-tuning a 100M parameter LM with extension to longer windows is multi-day GPU work that overshoots the EP scope).

---

## Candidate architectures

Three serious candidates. All operate on the same N=12 tet smoke cohort (already populated; no new data work needed). All run on Precision 7780 GPU (Codex).

### Candidate A — Per-gene NT windows with locus-specific pooling

**Approach:**
- Instead of mean-pooling all N≈5000 gene embeddings into one 512-vector, KEEP per-gene embeddings.
- For each gene in a known cipro / cef / tet target catalog (per `dna_decode/data/mic_tiers.py` `DRUG_LOCI_BY_MECHANISM[drug]`), pull its embedding directly.
- Strain feature = concatenation OR attention-pooled vector over the per-locus embeddings only (size = K_loci × 512, where K_loci ≈ 10-30 per drug).
- Train XGBoost (or logistic regression with regularization) on these per-locus features.

**Pros:**
- Preserves the signal AMRFinderPlus already knows about.
- Uses the existing audit framework's mechanism catalogs directly — no new biological knowledge required.
- Implementation reuses `dna_decode/models/cache.py` per-gene access (already exists).

**Cons:**
- Locus-catalog completeness ≠ data completeness. If a strain has a novel resistance gene NOT in the catalog, the architecture misses it entirely. (Mean-pool catches it weakly; this catches it not at all.)
- Loses the "global context" signal that may matter for regulatory/epistatic mechanisms.
- For tet specifically: tet-family genes (tetA/B/M/O) are mobile-element-borne, so they appear in many contexts (chromosomal AND plasmid). Per-gene attention may still work, but signal-to-noise depends on whether AMRFinderPlus locus naming captures all tet-family variants.

**Estimated tet AUROC under this architecture:** > 0.65 (extrapolating from the cipro 0.75 + cef 0.83 results being driven by localized signal — tet-family genes are localized to their CDS contexts when present).

**Estimated implementation cost:** 2-3 days on Precision 7780 (refactor `pipeline.py` train + predict to support per-locus feature mode).

---

### Candidate B — k-mer + AMRFinder-feature fusion baseline

**Approach:**
- Concatenate (a) NT mean-pooled 512-vector (current v0) + (b) k-mer count features (already shipped: `dna_decode/eval/loso_kmer.py`) + (c) AMRFinder per-mechanism counts (per-drug from `scripts/drug_mechanism_audit.py` outputs).
- Train XGBoost on the fusion feature vector.

**Pros:**
- Re-uses existing infrastructure end-to-end (no new code modules; just orchestration).
- k-mer baselines already work on tet (smoke AUROC 0.722) — fusion can ONLY match or beat that.
- AMRFinder features encode the "is this mechanism present" signal directly, which is what the mean-pool fails to capture.

**Cons:**
- It's a baseline-stacked classical model. The "AI DNA decoder" narrative weakens: the predictive power is coming from AMRFinder + k-mer, not from the foundation model. Codex's prior /probe critique: "your wrapper has to win on phenotype utility, provenance quality, or workflow ergonomics — not on 'it uses AI.'"
- If we win on tet via this fusion, we have a fusion-stack model — defensible as engineering, weaker as differentiation vs AMRFinderPlus + a wrapper.
- Doesn't address the underlying signal-dilution problem; it routes around it.

**Estimated tet AUROC under this architecture:** ≥ 0.72 (matches k-mer floor; AMRFinder features may add 5-10 pp).

**Estimated implementation cost:** 1 day on Precision 7780 (mostly feature concatenation + retrain).

---

### Candidate C — Lightweight transformer head on top of pooled NT (per-gene attention)

**Approach:**
- Keep NT v2 100M frozen + per-gene embeddings (no fine-tuning).
- Train a SMALL transformer head (e.g., 2-4 transformer layers, 256 hidden) that attends across the per-gene embedding sequence.
- Output is a strain-level prediction; the attention mechanism learns which genes matter per drug.

**Pros:**
- The "right" architectural answer for distributed-signal data — attention over per-gene embeddings naturally learns to ignore noise + focus on signal.
- Can generalize beyond the locus catalog (unlike Candidate A) because attention learns from data which genes matter.
- Defensible architecturally: "we use foundation-model embeddings + light learned aggregation."

**Cons:**
- New code: a Pytorch transformer head + training loop + checkpoint management. Not currently in repo.
- Training a transformer head on N=12 mini-cohort is BARELY enough data; high overfit risk; would need cross-validation discipline.
- More compute: each training run is minutes, not seconds. Hyperparameter sweep multiplies.
- Risk of "small NN doesn't beat XGBoost at this scale" — a real failure mode (per the McElfresh NeurIPS 2024 finding cited in `LESSONS_LEARNED.md` 2026-05-13: GBDT > NN at small samples-to-features ratios).

**Estimated tet AUROC under this architecture:** ?? — high variance. Could be > 0.80 if the architecture catches the per-gene signal cleanly. Could be < 0.50 if the small transformer overfits to N=12 noise.

**Estimated implementation cost:** 4-7 days on Precision 7780 (new training code; hyperparam tuning; cross-validation).

---

## Selection rubric

Each candidate runs on the existing N=12 tet smoke cohort (same data Stage 1 saw). Score each on 5 axes:

| Axis | Weight | Candidate A score | Candidate B score | Candidate C score |
|---|---|---|---|---|
| 1. Tet AUROC ≥ 0.50 (clears anti-predictive floor) | HARD GATE | TBD | TBD | TBD |
| 2. Tet AUROC ≥ 0.65 (gives credible tet path) | 1 | TBD | TBD | TBD |
| 3. Does NOT regress cipro AUROC (≥ 0.75 floor) | HARD GATE | TBD | TBD | TBD |
| 4. Does NOT regress cef AUROC (≥ 0.80 floor) | HARD GATE | TBD | TBD | TBD |
| 5. Implementation cost (days) | -1 | A=2-3 days | B=1 day | C=4-7 days |
| 6. Generalizability (does it survive a novel-mechanism strain not in catalog?) | 1 | LOW | MEDIUM | HIGH |
| 7. North-star alignment (foundation-model-driven vs feature-engineering-driven) | 1 | HIGH | LOW | HIGH |
| 8. Maintenance cost (new code to keep current) | -1 | LOW (reuses existing modules) | LOW (no new code) | HIGH (new training loop) |

**Decision rule:** the chosen candidate must pass ALL hard gates (1, 3, 4). Among passing candidates, pick the one with the highest weighted soft-axis score AND that doesn't cost > 5 days to implement.

**Tie-breaker:** if two candidates are within 1 point on weighted score, prefer the LOWER implementation cost (favor shipping).

**Fallback:** if NO candidate passes all 3 hard gates:
- (a) Document the failure in this packet.
- (b) Distributed-mechanism drugs (tet, gent) are documented as OUT-OF-SCOPE for v1; v1 ships cipro + cef only.
- (c) EP-1.5 closes; EP-2 ships cef-only sub-track; tet/gent deferred to v1.5+ with a re-design EP.

---

## Proof-of-concept protocol (Codex executes)

**Cohort:** existing N=12 tet smoke cohort at `data/processed/gate_b_mini_tet_cohort.parquet`. NT cache at `D:/dna_decode_cache/embeddings/nt_n40_cipro.h5` (shared with cipro smoke per the 2026-05-17 EP2 reuse).

**Protocol for each candidate:**

1. Reuse `scripts/smoke_gate_12strain_cipro.py --drug tetracycline` plumbing for fold structure + AUROC computation.
2. Replace the NT-XGBoost feature pipeline with the candidate's feature path.
3. Run LOSO on the 12-strain cohort (calibrate=False).
4. Record AUROC + a paired bootstrap CI vs the k-mer baseline (already 0.722 on tet smoke).
5. Run cipro + cef regression smokes with the same architecture to check hard gates 3, 4.
6. Record runtime (wallclock + GPU memory).
7. Emit a per-candidate result row in `wiki/ep1_5_poc_results_<DATE>.{md,json}`.

**Total POC compute budget:** ~1-2 hours on Precision 7780 per candidate (12-strain cohort is small). Total ~3-6 hours for all 3. Far cheaper than the EP-2 multi-drug cohort populates (~4-5 hr per drug).

**POC stopping rule:** if Candidate A or B passes all 3 hard gates with tet AUROC ≥ 0.65, STOP the POC + lock the choice. Don't run Candidate C (high-variance + expensive) unless A AND B both fail.

---

## Decision template (Codex fills post-POC)

After POC runs, this section gets filled with actual numbers:

```
## Decision (2026-XX-XX, post-POC)

Selected candidate: [A | B | C | NONE — distributed-mechanism drugs deferred to v1.5+]

POC results table:
| Axis | Candidate A | Candidate B | Candidate C |
|------|-------------|-------------|-------------|
| Tet AUROC | actual | actual | actual |
| Cipro AUROC (regression) | actual | actual | actual |
| Cef AUROC (regression) | actual | actual | actual |
| Implementation days | confirmed | confirmed | confirmed |
| Hard gates passed | Y/N | Y/N | Y/N |

Rationale:
[2-3 sentences on why the chosen candidate beat the others]

Trade-off acknowledged:
[What we give up by picking this candidate]

Next: EP-2 multi-drug starts; tet sub-track uses the chosen architecture.
```

---

## Risk flags

- **R1 (MEDIUM):** N=12 cohort is statistically thin. POC results have ±0.20 AUROC noise bands. A candidate winning by < 0.10 AUROC is within noise; pick on cost in that case.
- **R2 (MEDIUM):** Candidate C (transformer head) hyperparameter tuning could balloon scope. Hard-cap at 3 hyperparam configs per candidate; if none work, document + move on.
- **R3 (LOW):** Cipro + cef regression test could surface that a winning candidate REGRESSES on the concentrated-mechanism drugs. If that happens, the architectural fork (per D3 of `Post_V0_EP_Ladder_Plan.md`) becomes load-bearing: keep current architecture for cipro + cef, use the new one ONLY for tet/gent.
- **R4 (LOW):** AMRFinderPlus feature counts (Candidate B) require running AMRFinder on the 12 tet strains. AMRFinderPlus on Precision 7780 is already validated; runtime ~95s/strain × 12 = ~20 min. Not a blocker.
- **R5 (BEHAVIORAL):** Candidate B "wins" via classical-stack feature engineering. Resist the temptation to celebrate if it wins on tet but loses the north-star alignment (Codex's /probe critique on the "moat" — winning via AMRFinder features makes us a wrapper, not a decoder). Codex should report this honestly in the decision rationale.

---

## What this packet deliberately does NOT cover

- Cef + cipro architectural rework. The cross-drug finding showed both PASS on mean-pool — no architectural change needed.
- 4th-mechanism-class falsifier (colistin / aminoglycoside). The cross-drug pattern would be RE-OPENED if a future 4th-drug test contradicts §3. Out of scope for EP-1.5.
- Eukaryotic / multimodal architecture. Different scope entirely (EP-4+ research-program territory).
- The actual EP-2 multi-drug implementation. EP-1.5 only picks the architecture; EP-2 builds + trains the models.

---

## Open questions (for Codex / user post-POC)

1. If Candidate A wins but tet AUROC is only ~0.55 (above the 0.50 floor but below the 0.65 target), do we ship tet with INDETERMINATE attribution AND a documented "this is a weak signal" caveat? Or defer tet entirely?
2. If Candidate B wins, do we re-name the project framing internally (e.g., from "AI DNA decoder" to "AI-augmented AMR predictor")? Or keep the framing + accept that some drugs are AI-decoded and some are feature-stack-decoded?
3. Should we ALSO POC a hybrid (mean-pool 512-vec + per-locus concatenation) as Candidate D? Effort is incremental on top of Candidate A; could be worth it if A wins but is close to the gate.
4. Should the POC also test a "leave-one-mechanism-class-out" cross-validation (train on cipro+cef, test on tet) to surface whether the architecture transfers? Out of scope strictly, but cheap.

---

## Bottom line

EP-1.5 is decision-only — picks the architecture, doesn't build EP-2. Recommended starting order: run Candidate B first (cheapest, 1 day, validates a floor); then Candidate A (medium cost, 2-3 days, validates the AI-decoder hypothesis); skip Candidate C unless A+B both fail.

Most likely outcome (extrapolating from cross-drug pattern): Candidate A wins on north-star alignment + ships tet at AUROC ~0.65-0.70. Candidate B wins on cost + ships at ~0.72 but degrades the project narrative. Codex makes the final call with POC data in hand.
