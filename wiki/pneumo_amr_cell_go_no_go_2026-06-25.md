# Pneumococcus AMR cell — GO/NO-GO (2026-06-25)

Assessed the in-hand pneumococcal measured-AST (GPS Poland cohort, already downloaded) + probed the
determinant-catalog landscape, to decide whether to build a pneumococcus measured-MIC AMR cell. **No Docker
used** (the ktype finisher is occupying it); this is the assessment + ceiling, with any caller build deferred.

## Verdict: CONDITIONAL-GO, DEFERRED (a real cell, NOT a cheap exploit)
The free measured label is real and clean — but the cell is a scoped multi-session BUILD, not the afternoon
exploit the "highest-VOI widen move" framing implied. Verify-in-batch corrected that over-optimism.

## Rec 2 (probe) — does a free pneumo-AMR determinant catalog exist? **YES**
GPS's own pipeline (MOESM4) is a free, open-source determinant caller, and it IS the catalog:
- **β-lactams:** `pbp1a`/`pbp2b`/`pbp2x` PBP-type → MIC (the CDC PBP-typing method; a regression table, NOT
  gene-presence).
- **Other drugs:** per-drug `*_Determinant` columns (erm/mef for macrolides, tet(M/O) for tetracycline,
  gyrA/parC for FQ, folA/folP for co-trimoxazole) — gene-presence/point, AMRFinder-style.
Free + curated. So gate-1 (catalog) is CLEARED.

## Rec 1 (assess the in-hand AST) — the data + the ceiling
- **Measured AST (MOESM3):** ~12 drugs, n=49–259, ALL wet-lab (clears the circularity gate G1/G3). BUT
  **method-heterogeneous** — Penicillin/Cefotaxime/Meropenem/Vanc/Chloramphenicol are ETEST(210)+broth(49)
  MIC; Erythromycin/Levofloxacin/Tetracycline/Clindamycin are agar/disc **zone diameters** (not MIC);
  Ciprofloxacin n=0.
- **No-Docker ceiling computed (GPS determinant pipeline vs measured MIC):**
  - **Penicillin @ meningitis breakpoint (S≤0.06): acc 0.969** (n=257; TP96 FP5 TN153 FN3) — clean, usable.
  - Penicillin @ non-meningitis (S≤2): acc 0.881 (FP27) — the SAME data, worse, purely from the breakpoint.
  - Meropenem @ non-meningitis: acc 0.782 (FP50, FN0) — over-call, again breakpoint-driven.
  - **Finding: β-lactam R/S is breakpoint-AMBIGUOUS** (meningitis/non-meningitis/oral give different R/S from
    one MIC). Any cell MUST fix the breakpoint context explicitly (the project's `mic_tiers` is E. coli — pneumo
    needs its own meningitis/non-meningitis/oral breakpoint set).

## Why it's CONDITIONAL/DEFERRED, not a clean GO
1. **β-lactams (the main pneumo drugs) need a PBP-typing engine the project does NOT have.** The current AMR
   engine is gene-presence AMRFinder; pneumo β-lactam resistance is PBP-mutation-MIC (the CDC regression
   table). Building that is real, multi-session work.
2. **Wrapping GPS's calls = faithful-to-tool.** The determinant catalog IS GPS/CDC's pipeline; if our cell just
   re-runs it, we validate THEIR tool, not add value ([[feedback_validate_wrapper_vs_underlying_tool]]). The
   0.969 above is GPS's-tool-vs-measured, i.e. the ceiling — NOT our decoder.
3. **Docker-gated.** Building OUR caller (AMRFinder for erm/mef/tet + a PBP-typer) needs Docker, currently busy
   with ktype; do NOT contend (the contention lesson).
4. **Single Poland cohort, modest N, mixed methods** → G4 (surveillance/cohort domination) + G8 (clonality) +
   breakpoint heterogeneity all apply. A real cell needs multi-cohort + breakpoint discipline.

## What would make it a clean GO (the scope, if pursued)
- A pneumo breakpoint set (meningitis/non-meningitis/oral) added to `mic_tiers`.
- An independent caller: AMRFinder gene-presence (macrolide/tet/FQ — tractable) + a PBP-type→MIC engine
  (β-lactams — the hard part). Then validate vs measured AST at the correct breakpoint, headlining the delta
  over naive GPS-pipeline use.
- Multi-cohort (not just Poland n=263) to clear G4/G8.
- Cleanest tractable SUB-cell to start: **macrolide (erm/mef gene-presence) vs measured** — simple determinant,
  but the measured label is zone-diameter (needs disc breakpoints) + Docker for AMRFinder.

## Rec 3 — bank + the acquisition framing (the honest frontier)
This assessment REINFORCES the negative-results map's conclusion: free public data isn't the cheap unlock it
looks like — even a "free measured label" cell (pneumo AMR) carries a real engine-build + breakpoint +
single-cohort cost. The genuine data-availability frontier remains **non-public acquisition** (a
collaborator/biobank/clinical wet-lab AST source clears the gates by construction) — a USER relationship
decision, not a code or research task.

## Bottom line
- **GO/NO-GO: CONDITIONAL-GO, DEFERRED.** Real cell (free clean measured label, determinant catalog exists,
  0.969 ceiling proven), but a scoped multi-session PBP-engine build, Docker-gated, not a tonight-exploit.
- Banked as a scoped future cell with the obstacles named above. The cheap-exploit framing was corrected by
  verify-in-batch.
- Highest real "widen" frontier = non-public acquisition (rec 3) — your call.
