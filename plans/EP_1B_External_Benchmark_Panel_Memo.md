# EP-1B External Benchmark Panel Selection — Memo

> Selection criteria + diversity rubric for the 10-genome external AST-labeled E. coli panel that fires the external benchmark smoke test (`plans/Post_V0_EP_Ladder_Plan.md` EP-1B). Side-by-side comparison: dna_decode v0 vs AMRFinderPlus vs RGI vs known cipro AST label.

**Status:** DRAFT 2026-05-25. Codex picks the final 10 genomes on Precision 7780 after EP-0 close + EP-1A.

---

## Why this memo exists

EP-1B is the first external validation of v0 — moves the moat claim ("honest output discipline + audit verdict") from internal hypothesis to evidence. Without a pre-staged panel selection rubric, panel composition could quietly skew results (e.g., 10 wild-type cipro-S = trivial 10/10 agreement; 10 strong-QRDR cipro-R = trivial 10/10 agreement; either gives a misleading "we agree with AMRFinderPlus" narrative).

Pre-stages the rubric so Codex picks a panel that actually exercises the prediction surface.

## Hard panel criteria

1. **N = 10 E. coli genomes.** Not 1 (statistical noise); not 100 (over-scope for "smoke" benchmark). 10 is the floor for a credible side-by-side; can scale up in a follow-up EP.
2. **Public + freely redistributable** (same as EP-1A criterion 1).
3. **None overlap the N=147 training cohort.** Reuse EP-1A criterion 2.
4. **All have a known cipro AST phenotype + a documented MIC or category** (R / I / S). Required for the comparison column.
5. **All have downloadable genome FASTA + GFF3.** Required for v0 to run.
6. **Mechanism class diversity** (see below) — the load-bearing requirement that makes the panel non-trivial.
7. **Assembly quality threshold** matches EP-1A criterion 6.

## Mechanism class diversity (the load-bearing piece)

Pick 10 genomes spanning these classes:

| Class | Target N | Why |
|---|---|---|
| **Wild-type / no AMR mechanism, cipro-S** | 3 | Tests S prediction path + "no false alarms" behavior |
| **QRDR point mutation (gyrA-S83L, parC-S80I, etc.), cipro-R** | 3 | v0's strongest signal class; tests R agreement |
| **Plasmid quinolone resistance (qnr family OR aac(6')-Ib-cr), cipro-R or I** | 2 | v0 + AMRFinderPlus disagreement is likely here (per the 2026-05-17 cross-drug finding — mean-pool dilutes plasmid signal); high information |
| **Regulatory mutation (marR / acrR / soxR frameshift), cipro-R or I** | 1 | Edge case; tests whether v0's audit verdict catches regulatory mechanism |
| **MULTI-CLASS** (≥ 2 mechanisms co-occurring), cipro-R | 1 | Real-world clinical case; tests audit verdict propagation under mechanism stacking |

**Why this distribution beats "10 random cipro-R + cipro-S":** the 3+3+2+1+1 split forces v0 to show its strengths AND weaknesses on shared inputs, which is the actual moat-validation question. A uniform-random panel would over-sample the easy cases.

## Soft preferences (apply to each pick within its class)

- Recent deposition (2023+).
- Diverse MLST types (don't pick 10 ST131 isolates; spread across ST10, ST73, ST95, ST131, ST167, etc.).
- Geographic + source diversity (clinical / animal / environmental — not all from one outbreak).
- AMRFinderPlus + RGI calls already published or computable from the genome (saves us re-running AMRFinder).

## Candidate sources

Same set as EP-1A, with added emphasis:

| Source | Particularly useful for |
|---|---|
| NCBI Pathogen Detection | Browse-and-filter UI; per-isolate AMRFinderPlus call available; AST sometimes attached |
| NARMS | Public AST data; categorical R/I/S labels with MIC where measured |
| PATRIC / BV-BRC | Avoid (high overlap risk with N=147) |
| ATCC | Wild-type reference (ATCC 25922) for the cipro-S class |
| ResFinderFG curated panel | Pre-vetted mechanism examples |

## The comparison table (output of EP-1B execution)

Each row in `reports/external_benchmark_cipro_<DATE>.{md,json}`:

| strain | accession | MLST | known AST | known mechanism | v0 prediction | v0 calibrated_prob | v0 confidence_tier | v0 attribution_scope_confidence | v0 top-3 attribution | AMRFinderPlus calls | RGI calls | agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|

Plus:
- Per-row "agreement_with_AST" (Y/N) for v0 + AMRFinderPlus + RGI.
- Aggregate: overall accuracy of each tool vs AST.
- "v0-additional-info" analysis: cases where v0's calibrated probability or audit verdict carries information AMRFinderPlus's gene call doesn't.
- "AMRFinderPlus-additional-info" analysis: cases where AMRFinderPlus gives information v0 doesn't (e.g., specific mutation calls v0 can't make).

## EP-1B success criteria (recap from `Post_V0_EP_Ladder_Plan.md`)

- 10/10 genomes complete without errors.
- v0 prediction agrees with AST label on at least 7/10 (≈ 0.70 accuracy floor matching v0 spec).
- Side-by-side renders a "what does v0 add vs AMRFinderPlus" answer one way or the other.

## What's measured + what's NOT

**Measured:**
- v0 agreement with AST.
- AMRFinderPlus + RGI agreement with AST.
- Cases where v0 + AMRFinderPlus disagree.
- v0's confidence tier + scope confidence on agreement vs disagreement cases.

**NOT measured (out of scope for EP-1B):**
- Clinical decision-support performance (would require ethics review + clinical-grade panel).
- Multi-drug accuracy (cef / tet — EP-2 territory).
- Speed/cost benchmark (EP-2.5+ territory).
- Multi-organism (EP-3 territory).

## Selection rubric (Codex executes)

```
For each candidate genome:
  1. Does it satisfy hard criteria 1-7? If no, skip.
  2. Which mechanism class does it fit? Track per-class count.
  3. Once a class quota is hit, skip further candidates in that class.
  4. Within each class, prefer candidates with more soft preferences satisfied.
  5. Stop once N=10 + class quotas satisfied.

Failure mode: if a class quota can't be hit (e.g., can't find a regulatory-mutation
+ cipro-R genome with public AST), document the gap + substitute with an extra
genome from another class. Note the substitution in the EP-1B report.
```

## Open questions for Codex / user

1. Should EP-1B run AFTER EP-1A's same-strain parity test passes, or in parallel? If parallel: ε in the parity test could still be tightened post-1B; if sequential: 1A becomes a hard gate for 1B.
2. Is the 7/10 agreement floor for v0 the right bar, or should it be tighter given v0's 0.87 leakage-safe CV AUROC? 7/10 = 70%, but 0.87 AUROC implies ~85% accuracy at threshold 0.5. Tighten to 8/10 OR keep at 7/10 (matches v0 spec language)?
3. Does CARD/RGI add value over AMRFinderPlus alone for EP-1B? Both are gene-call tools with curated databases; they overlap heavily. ResFinder + AMRFinderPlus might be more complementary. Codex / user decides.
4. Should an LLM-based AMR predictor (if any open-source ones exist) join the comparison? Out of scope for v0 but a meaningful "AI vs AI" data point if cheap. Defer to v1+ unless trivial.

## What this memo deliberately does NOT cover

- Actual accession numbers (Codex picks).
- AMRFinderPlus / RGI invocation mechanics (existing scripts already wrap them).
- The cef / tet / gent external benchmark (EP-2 territory).
- Clinical interpretation guidelines.
