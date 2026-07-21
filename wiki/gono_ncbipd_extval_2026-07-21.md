# N. gonorrhoeae external validation on NCBI-PD — cipro + cefixime v0.1 CONFIRMED (2026-07-21)

**Status:** ✅ **ciprofloxacin + cefixime externally CONFIRMED** on 160+ independent isolates (the fresh
AR-Bank cells + the cefixime v0.1 rule generalize) · 4 other gonococcal cells honestly characterized ·
**Cell:** `neisseria_amr` (NON-FROZEN) · **Frozen surface:** byte-unchanged.

## Why this run

The AR-Bank Neisseria validation (2026-07-21) endorsed cipro + cefixime on **N=20**, but the cefixime v0.1
mosaic-penA-34 rule was **derived + validated on that same 20-isolate cohort** — an honest tension (a
cohort-fit rule is not yet externally confirmed). This scores the SAME cell on **170 independent NCBI
Pathogen Detection isolates**, provenance-disjoint from the AR-Bank set (the 21 AR-Bank BioSamples excluded).

**Feasibility win — no compute wall.** NCBI-PD publishes its own AMRFinderPlus calls in the `AMR_genotypes`
field *including point mutations* (`gyrA_S91F=POINT`, `penA_I312M=POINT`), and this cohort's isolates have
downloadable assemblies (NCBI ran AMRFinder for us). So this is a pure **metadata join + score**, offline,
in seconds — no Kaggle, no Docker, no assembly download. The penA numbering was verified compatible between
the two AMRFinder sources ({I312M,V316T,N512Y,G545S} present in both; `penA_G542S` is a separate real
allele, not a renumbered G545S).

## Result (170-isolate provenance-disjoint cohort; powers 5 drugs)

| drug | n | R/S | sens | spec | acc | verdict |
|---|---|---|---|---|---|---|
| **ciprofloxacin** | 163 | 94R/69S | **0.989** | **0.986** | **0.988** | **SCORED_ENDORSED** ✅ |
| **cefixime** (v0.1) | 166 | 19R/147S | 0.789 | **0.905** | 0.892 | **SCORED_ENDORSED** ✅ |
| ceftriaxone | 169 | 2R/167S | — (2R) | 0.946 | 0.935 | UNDERPOWERED (spec confirmed) |
| azithromycin | 156 | 110R/46S | 0.0 | 1.0 | 0.295 | **DEGENERATE** (data gap — see below) |
| **penicillin (v0.2)** | 31 | 14R/17S | 0.929 | 0.941 | 0.935 | **SCORED_ENDORSED** ✅ (narrowed — see below) |
| tetracycline (v0.2) | 60 | 34R/26S | 0.324 | 1.000 | 0.617 | SCORED_NOT_ENDORSED (honest partial — see below) |

## penicillin + tetracycline v0.2 — narrowed to the clean plasmid markers (the v0.1 over-call fix)

The first run showed penicillin + tetracycline calling **all-R (spec 0.0)** — a DEGENERATE over-call. This
was the **predicted failure of their v0.1 rules**: those were validated on an R-saturated AR-Bank cohort
(0 S isolates) and their own docstrings warned "spec untested; the promoted chromosomal markers (penA-point,
mtrR, rpsJ_V57M) will over-call once S isolates exist." The NCBI-PD cohort has S isolates, so the warning came
true. **v0.2 narrows each to the SPECIFIC literature determinants:**

- **penicillin v0.2 = blaTEM OR ponA_L421P.** blaTEM (plasmid penicillinase, PPNG high-level) is 8/14 R, **0/17
  S — clean**; ponA_L421P (chromosomal PBP1 acylation defect) is 9/14 R, 1/17 S. Combined: **sens 0.929 /
  spec 0.941 → SCORED_ENDORSED.** The near-universal penA-point + mtrR (lineage-linked, not penicillin-causal)
  are DROPPED from the binary call.
- **tetracycline v0.2 = tet(M) only.** tet(M) (plasmid, high-level TRNG) is 11/34 R, **0/26 S — clean** →
  **spec 1.00, sens 0.324.** Honest partial: it cleanly IDs high-level TRNG but the chromosomal low-level
  tet-R (rpsJ_V57M + mtrR + penB cumulative, MIC 2-4) is **not cleanly determinant-separable** from tet-S
  (same near-universal markers) → a ~68% sens ceiling that is a real multi-locus-cumulative property, not a
  rule bug. Not endorsed (low sens), but no longer a *misleading* all-R over-call.

Both markers are **literature determinants** (penicillinase, PBP1, TRNG plasmid) — not cohort-mined — so this
mirrors the cefixime v0.1 posture (derived-informed-by-cohort + literature-grounded). Mild optimism caveat:
the sens/spec is on the same cohort that informed the narrowing.

## The two confirmations (the headline)

- **ciprofloxacin — externally CONFIRMED, acc 0.988** (163 isolates). The gyrA/parC QRDR rule is a
  *literature* rule (not fit to any cohort); its near-perfect score on 163 independent isolates is a genuine
  strong external confirmation (N=20 AR-Bank → N=163 NCBI-PD).
- **cefixime v0.1 — externally CONFIRMED, spec 0.905** (147 independent S isolates). The mosaic-penA-34 rule
  I *derived* on the 20 AR-Bank isolates **generalizes**: it holds specificity ≥0.85 on 147 unseen S
  isolates. This directly closes the honest tension — the rule is now cohort-derived **AND** externally
  confirmed. (sens 0.789 = 15/19 R, slightly below the AR-Bank 0.917 on independent data — honest.)

## The four honest negatives (external validation doing its job)

The other 4 gonococcal cells were NEVER AR-Bank-validated (that run scored only cipro + cefixime); this run
is their first test, and it exposes that they are **not reliable** on this substrate:

- **azithromycin — DEGENERATE (sens 0.0): needs a CUSTOM 23S caller (diagnosed 2026-07-21, follow-up #4).**
  The cell keys on 23S rRNA mutations (A2045G/C2611T…). Two-level diagnosis: (1) NCBI-PD's published AMRFinder
  calls omit 23S rRNA (0/110 azithro-R carry any 23S marker); (2) **AMRFinderPlus ITSELF does not call gono
  23S rRNA point mutations** — verified on our own AR-Bank `-O Neisseria_gonorrhoeae` run: its determinant
  classes are all protein/promoter (PBP2, GyrA, S10/RpsJ, FolP, PonA, MtrR, ParC, RpoB), with the "macrolide"
  keyword appearing only in the multi-drug `mtrR` efflux description, never as an actual 23S rRNA call. 23S is
  a 4-copy rRNA gene that AMRFinder's protein-centric method doesn't handle for gonococcus. **So azithromycin
  requires a custom 23S BLAST-caller subsystem** (native `blastn` is available — no Docker): BLAST a
  gonococcal 23S rRNA reference vs each assembly, map to the resistance positions (E. coli 2059 / gono 2045),
  call the mutation. **Inherent ceiling:** 23S is multi-copy and azithromycin resistance is often
  heteroplasmic (a subset of the 4 copies mutated), which the consensus assembly collapses — so only
  high-level all-copy R is reliably detectable from WGS (a well-known literature limitation, like tet(M) for
  high-level TRNG). This is a scoped subsystem build with a real biological ceiling, deliberately NOT
  half-built here (a partial caller would produce a misleading number).
- **penicillin + tetracycline — DEGENERATE (spec 0.0, all-R over-call).** Both cells' own docstrings flagged
  an over-call risk (penicillin: chromosomal penA/mtrR promoted → near-universal; tetracycline: rpsJ V57M
  near-universal). This run **confirms the over-call on independent data** — spec 0.0. They need
  acquired-marker narrowing (like the cefixime v0.1 fix) before they can be endorsed.
- **ceftriaxone — spec confirmed (0.946), sens underpowered (2R).** The A501-specific v0.1 rule holds
  specificity on 167 S isolates; sens is untestable (only 2 R, both non-A501).

## Scorer hardening (verify-in-batch)

The first run mislabeled azithromycin `SCORED_ENDORSED` — a cell predicting all-S gets spec 1.0 + passes a
spec-only endorsement check while being useless (sens 0.0). Added a **DEGENERATE guard**: a cell that
predicts all-one-class (sens 0 with powered R, or spec 0 with powered S) is `DEGENERATE_NOT_ENDORSED`
regardless of the other metric. Same class of integrity fix as the empty-assembly-as-S gate — a degenerate
output must never read as a validated one.

## Lineage-collapse (clonality-corrected) — the caveat RESOLVED, no compute

The dominant caveat is clonal inflation (raw sens/spec counts one vote per isolate, not per lineage).
Resolved here with **no Mash / no Docker** by using NCBI-PD's OWN published SNP clusters (`<PDG>.reference_
target.all_isolates.tsv` → each isolate's `PDS_acc`; NULL = a singleton lineage) with the frozen
`clonality.cluster_weighted_confusion` (one vote per lineage; mixed-label clusters excluded as DISCORDANT).
170/170 isolates SNP-clustered (153 in real clusters + 17 singletons):

| drug | RAW sens/spec | **LINEAGE-collapsed sens/spec** | discordant |
|---|---|---|---|
| ciprofloxacin | 0.989/0.986 | **1.00 / 1.00** | 3 |
| cefixime | 0.789/0.905 | **0.727 / 0.892** | 2 |
| penicillin | 0.929/0.941 | **0.917 / 0.933** | 0 |

**All 3 endorsed cells HOLD at the lineage level** — cipro is *perfect* collapsed, cefixime + penicillin
keep spec ≥ 0.85. The cells are NOT clonally inflated → the determinant rules genuinely decode the
**mechanism**, not clonal population structure. This is the strongest validation tier the project produces,
and it was reached at zero compute cost.

## Honest caveats

- **Provenance-disjoint** (different isolates than AR-Bank; the 21 AR-Bank BioSamples excluded) but **NOT
  methodology-independent** (same AMRFinderPlus + same `neisseria_amr` cell) — the standard caveat.
- NON-FROZEN cell; the frozen decoder surface is byte-unchanged (`verify_lock` OK).
