# Session handoff — 2026-06-10 (dna_decode)

Start here in a fresh session. Project root: `C:\Users\Farshad\PythonProjects\dna_decode`. Ledger:
`project_state/eukaryotic-trait-decoding-cycle-2026-06-07.md` (action log through row 65). Branch `main`,
~46 commits ahead of origin (user syncs ~weekly — do NOT push). North star: **AI DNA-decoder TOOL, breadth,
failure-tolerant iteration** (not papers).

## 1. ONE thing in flight — collect it first
- **E. coli flagship cipro provenance-disjoint validation** — background run `bfzbxz61z` (was 21/60 at handoff,
  ~40 min remaining; ~60 fresh genome download + AMRFinder). It validates the **DEFAULT `DRUG_RULE`
  qrdr_point@2** (the core shipped decoder, not a calibrated-registry rule) on a fresh, leakage-verified
  (0 overlap with the N=147 tuning + gate_b held-out parquet cohorts) provenance-disjoint cohort (30R/30S).
- **To collect on resume:** check the run output (or just re-run — it's cached/restartable):
  ```
  HF_HOME=D:/hf_cache .venv/Scripts/python.exe scripts/provenance_disjoint_validate.py \
    --group Escherichia_coli_Shigella --amrfinder-organism Escherichia --registry-organism Escherichia \
    --drug ciprofloxacin --per-class 30
  ```
  It writes `wiki/provenance_disjoint_validation_escherichia_coli_shigella_cipro_2026-06-10.{md,json}`.
  Then: verify 0 leakage (script now hardened to exclude parquet flagship cohorts), commit the artifact,
  append a `/project-state --append-action` ledger row, and the provenance-disjoint strand is **DONE**
  (flagship default-rule + 2 calibrated-registry rules, 3 organisms / 2 phyla, all leakage-verified).

## 2. The session's headline result (Anchor-3 arc)
A strategic conclusion was caught as OVERCLAIMED by `/brainstorm`, corrected, reopened, and turned into a
real **provenance-disjoint independent-ish validation of the shipped decoders**:
- **Klebsiella cipro: acc/sens/spec 0.967** (n=60, leakage-free) — `wiki/provenance_disjoint_validation_klebsiella_cipro_2026-06-10.md`
- **Campylobacter cipro: 1.0** (n=40, leakage-free) — different phylum, the TUNING-boundary single-gyrA-T86I rule
- **E. coli flagship: pending** (run above)
- **Conclusion:** the deterministic AMR decoders **generalize** to provenance-disjoint (different-submitter-lab)
  cohorts. **Honest tier (never inflate):** provenance-disjoint ≠ methodology-independent (most NCBI submitters
  still use CLSI broth microdilution) ≠ external clinical validation.

## 3. The reusable capability built this session (free, no money)
- `scripts/ncbi_pd_provenance_census.py` — Stage-1: streams NCBI-PD metadata, tallies AST by submitter-class,
  reports whether a provenance-disjoint subset is POWERED (≥20/class). Organism-DEPENDENT: Klebsiella/Campylobacter
  powered; **Salmonella NOT** (4R — NARMS/CDC-dominated).
- `scripts/provenance_disjoint_validate.py` — Stage-2: select fresh disjoint cohort (excludes ecosystem submitters
  AND all prior-cohort accessions incl. parquet flagship cohorts — the leakage fix) → ensure_run (reuse cache)
  → score the DEPLOYED `call_resistance(organism=, drug=)` → artifact with tier + leakage-control fields.
- Census + sweep + pay/no-pay decision: `wiki/independent_phenotype_{label_census,source_sweep}_2026-06-10.md`,
  `wiki/ncbi_pd_provenance_census_2026-06-10.md`. **DO NOT PAY verdict stands** (independence ⊥ genome-linkability
  across all sources; the free provenance-disjoint subset is the move).

## 4. Also shipped earlier this session
- **Antiviral influenza-NA decoder = the 4th kingdom** (commit e3b5711): `dna_decode/data/antiviral_amr.py` +
  `scripts/flu_na_caller.py` + `dna-amr` CLI route. H275Y/oseltamivir, real field-isolate validated. The
  deterministic target-site method now spans bacterial/fungal/protozoan/**viral**.
- **expression_context (Acinetobacter carbapenem EXPRESSION floor)** — built, validated, **HOLD: UNDERPOWERED**
  (corrected from an initial overclaim — the independent cohort was 14/15 acquired-carbapenemase, target subset N=1).
  Override ships OFF/experimental. `wiki/provenance...`? no → `executed_plans/Expression_Context_Acinetobacter_Meropenem_Plan/`.

## 5. Open decisions / ready next moves (all user-gated)
1. **Goal-supersession (your scope decision):** the ledger Refined goal still reads the narrow eukaryotic
   Path A/B; actual work is decoder-breadth + independent-validation. Fix = `/project-state
   eukaryotic-trait-decoding-cycle-2026-06-07 --supersede-goal --new-goal "..." --supersedes-row "..."`.
   I can draft the new goal line for ratification.
2. **NEW idea-anchor drafted (Anchor 4):** "decoder-suite provenance-disjoint validation REPORT CARD" —
   generalize the proven method into a standing suite-wide capability. Full prompt in the 2026-06-10 conversation
   (offer pending to save it into `plans/idea_anchor_prompts_2026-06-10.md`).
3. **Anchor 2 — pfmdr1 directional catalog:** already drafted at `plans/idea_anchor_prompts_2026-06-10.md`
   (Anchor 2). A new build (tests opposite-direction selection in the catalog architecture).
4. **H2 (Arabidopsis embedding, G2):** the one OPEN hypothesis — user-owned on the Precision-7780 GPU / Databricks.
   Re-run `/soraya --until-mvp eukaryotic-trait-decoding-cycle-2026-06-07` when `wiki/phase2_arabidopsis_result_*.md` lands.

## 6. Key gotchas / lessons banked this session (memory)
- **Run an adversarial `/brainstorm` on a strategic/closure conclusion BEFORE committing it** — 3 overclaims this
  session were caught post-commit (`feedback_brainstorm_strategic_conclusions_before_commit`).
- **Validation leakage:** exclude tuning/prior cohorts by accession — including parquet cohorts the data/raw glob
  misses (hardened in the validator). Gate rescue-overrides on TARGET-subset positive support, not whole-class
  (`feedback_abstain_override_needs_positive_support_gate`).
- **Env:** `python3` absent on Windows → use `.venv/Scripts/python.exe`; git-bash `/tmp` ≠ Windows-python `/tmp`;
  NCBI WebSearch policy-filter trips on AMR phrasing → use WebFetch on specific URLs; Docker churn can corrupt the
  WSL2 D: mount → restartable + `run_in_background`; commit to `main` (= the sync channel), never push.
- **SAFETY (standing):** do NOT route personal code through the Bombardier/DLP machine; money spend → hard pause + ask.

## 7. Fastest resume
```
cd C:/Users/Farshad/PythonProjects/dna_decode
git log --oneline -8
/project-state eukaryotic-trait-decoding-cycle-2026-06-07      # read current state
# then collect the E. coli flagship result (§1), or pick a §5 move
```
