# Eukaryotic Cycle — Dual-Machine Coordination Plan (2026-06-07)

> How the two machines split the eukaryotic leap with **zero duplication** + **phase gates** where we
> share progress and decide before continuing. Sync channel = **git origin** (`mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`),
> pull/push only — NO Gmail-attached code (DLP signal), NO routing personal code through the Bombardier
> machine. Companion: `wiki/HANDOFF_workhorse_eukaryotic_2026-06-07.md` (workhorse's exact scope).

## Machine roles (no overlap)

| Machine | Owns | Compute | Substrate |
|---|---|---|---|
| **LAPTOP (this, GTX 860M)** | **Path A — fungal AMR (C. auris azole)** end-to-end + ALL code/catalog/cohort-prep for both paths | BLAST scan only (no GPU); "minor compute OK to move our side forward" | C. auris WGS+MIC |
| **WORKHORSE (Precision 7780, RTX 3500 Ada ~12GB)** | **Path B — Arabidopsis embedding test** GPU compute ONLY | plant DNA-FM embedding (≥12–24GB GPU) | Arabidopsis flowering-time |

**No-duplication contract:** different substrates, different methods, different machines. The laptop never
runs the GPU embedding; the workhorse never touches the fungal-AMR determinant pipeline. The laptop *prepares*
the Path-B data/plan (cohort manifest, baselines spec) so the workhorse only does the GPU-bound step.

## Phases + share-and-decide gates

Each phase ends with **SHARE (git push + a result packet) → DECIDE (review results together) → continue /
iterate-same-phase / pivot.** No phase auto-advances.

### Phase 0 — Fungal-AMR infra ready (LAPTOP, no compute) — IN PROGRESS
- ✅ Determinant catalog (`dna_decode/data/fungal_amr.py`, 7 tests).
- ▶ BLAST ERG11/FKS1 caller (`scripts/fungal_erg11_caller.py`) + smoke on 1 C. auris reference genome.
- **Gate G0:** caller correctly calls a known ERG11 mutation on a reference genome. → SHARE → decide cohort scope.

### Phase 1 — Fungal-AMR validation (LAPTOP, light compute) 
- Extract C. auris WGS+MIC from supplementaries (S. Africa 188 / India 350); fetch assemblies; build a
  clade-de-confounded fluconazole cohort (reuse `cohort_deconfound.py`).
- Run the caller; report acc/sens/spec + efflux/aneuploidy discordance.
- **Gate G1:** fungal-AMR decoder validated (or documented failure mode). → SHARE result packet → DECIDE:
  ship as `dna-decode amr --organism Candida_auris` / iterate / stop. **This is the no-compute eukaryotic win.**

### Phase 2 — Arabidopsis embedding test (WORKHORSE, GPU) — GATED on Gate G1 + compute confirm
- Laptop pre-stages: AraGWAS download manifest + baseline spec (SNP-PRS, kinship-only) + CV fold design.
- Workhorse: VRAM-fit check (plant DNA-FM on 12GB) → embed accessions → baselines → clade-stratified CV →
  within-lineage diagnostic.
- **Gate G2:** embedding R² vs baselines → the embedding-niche PASS/FAIL (first YES/YES/YES test). → SHARE → decide.
- **Money gate:** if 12GB insufficient → cloud budget decision (user) before any paid compute.

## Sync discipline (both machines)
- Pull at session start, push at session end. `main` is the channel.
- Each machine appends to the project ledger (`project_state/dna-decode-2026-05-11.md`) ONLY for its own
  phase rows — no concurrent edits to the same section (laptop = Path A rows; workhorse = Path B rows).
- Result packets: `wiki/<phase>_result_<date>.md` — the share artifact at each gate.
- Cross-machine sync check: `scripts/cross_machine_sync_check.py` after each handoff.
