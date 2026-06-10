## Lens status
- /probe: applied (2026-06-10) — surfaced the genome-not-threaded plumbing gap, the cached-run provenance gap, the loose BLAST contract (`-max_target_seqs` truncation vs 22 ISAba1 copies), and the wrong-shaped validation endpoint. All four are addressed below.
- /idea-anchor: applied (2026-06-10) — ratified acceptance bar (0 S-upgrades) + v0 scope (ISAba1 to OXA-51 only).
- /brainstorm: applied (2026-06-10, Deep, 2 rounds) — 3 grounded issues + 1 inert-detector guard, all folded into the steps below (see `## Save-time amendments`).

## Problem Statement
The deterministic AMR decoder ABSTAINs on `Acinetobacter|meropenem` (registry verdict EXPRESSION_FLOOR) because resistance is driven by ISAba1 inserted upstream of the intrinsic blaOXA-51 gene (a hybrid promoter overexpresses it) — invisible to gene-PRESENCE, since every isolate carries OXA-51. Add a deterministic `expression_context` signal that reads the SAME assembly the decoder already has, detects an upstream ISAba1 to OXA-51 junction, and upgrades ABSTAIN to R — but only after the signal is validated on an INDEPENDENT Acinetobacter cohort with ZERO susceptible false-upgrades AND at least one true R rescue. Success is a perfectly-specific minority-rescue plus an honest quantification of how much of the floor stays sequence-invisible; it is NOT "carbapenem expression solved."

## Codebase Context
- `dna_decode/eval/amr_rules.py:180` `call_resistance(main_tsv, drug, resistance_threshold=None, organism=None, registry=None)` — takes only the AMRFinder main.tsv path; no genome. `:258` `_call_resistance_calibrated(p, drug, organism, cal)` has the EXPRESSION_FLOOR to `prediction "ABSTAIN"` branch (`:270`) — the exact override insertion point.
- `dna_decode/data/calibrated_amr_rules.json` — `Acinetobacter|meropenem` + `Pseudomonas_aeruginosa|meropenem` are EXPRESSION_FLOOR (n=30 each); the others are CALIBRATED.
- `dna_decode/typing/blast_caller.py:21` `call_alleles(..., with_positions=True)` — blastn helper with qlen(coverage) plus subject positions, offline-safe, reuses `dna_decode/pathotype/vf_runner.find_blastn`/`_find_makeblastdb`. BUT collapses to the single best hit per allele (`-max_target_seqs 5`) — cannot enumerate multi-copy ISAba1 (the rescued strain has 22 copies). Reuse its resolvers plus invocation idiom; the new module needs an all-hits variant.
- `scripts/independent_cohort_validate.py` — builds disjoint 15R/15S cohorts excluding cohort-1 accessions (`_exclusion_set`/`select_independent_cohort`), AMRFinder cached/restartable, labels from NCBI Pathogen Detection AST, no money. Its `TARGETS` is hardcoded to 3 CALIBRATED cipro organisms; scores acc/sens/spec and excludes ABSTAIN from `n`.
- `scripts/organism_drug_validate.py` — `ensure_run`/`_run_dir`/`latest_metadata_url`/`select_cohort` (the per-accession fetch plus AMRFinder run primitives, reused by the validator). FASTA lands in `data/raw/<slug>/refseq/<acc>/genome.fna`.
- `data/raw/acinetobacter_meropenem/` — cohort-1: 30 genomes plus `selected.tsv` (R/S). `data/isaba1_ref/ISAba1_ref.fna` (570 bp partial), `OXA51fam_ref.fna` (825 bp). `soraya_runs/2026-06-10-7qjq-expression-frontier-isaba1/isaba1_falsifier.py` — the frozen falsifier method (`6 sseqid sstart send pident length bitscore`, 400 bp strand-aware UPSTREAM-of-OXA window, pident>=85/len>=120, NO ISAba1-orientation check). Rescued accession: GCA_000692095.1 (22 ISAba1 copies).
- External tool: blastn/makeblastdb (NCBI BLAST+ 2.17.0 at `C:/Users/Farshad/ncbi-blast/bin`). Verified surface: `-outfmt "6 qseqid sseqid sstart send pident length"` (sseqid = contig, needed for the same-contig test), `-max_target_seqs`, `-evalue`, `-perc_identity` (the flags blast_caller already uses).

### Reusable-Code Survey
- `dna_decode/typing/blast_caller.py` — reuse `find_blastn`/`_find_makeblastdb` resolvers plus the makeblastdb+blastn+offline-degrade idiom (NOT `call_alleles` itself — it collapses multi-copy hits).
- `scripts/organism_drug_validate.py` — reuse `ensure_run`/`select_cohort`/`latest_metadata_url`/`_run_dir` for cohort-2 fetch (same primitives the existing validator composes).
- `scripts/independent_cohort_validate.py` — reuse the `_exclusion_set` pattern for disjointness; extend with an Acinetobacter EXPRESSION_FLOOR target plus an ABSTAIN-rescue endpoint.
- `soraya_runs/2026-06-10-7qjq-expression-frontier-isaba1/isaba1_falsifier.py` — the EXACT frozen rule (window, identity, proximity logic, NO orientation) ports verbatim into the module as the PRIMARY detector; do NOT re-tune and do NOT add constraints to the primary.
- None additional — searched: `dna_decode/typing/`, `dna_decode/eval/`, `scripts/`, `dna_decode/pathotype/vf_runner.py`.

## Pre-Change Baseline
- Current behavior: `Acinetobacter|meropenem` (and `Pseudomonas|meropenem`) resolve to `ABSTAIN` always (EXPRESSION_FLOOR), 0 rescue of intrinsic-only-R.
- Falsifier on cohort-1 (N=30, 15R/15S): junction-positive R 1/15, S 0/15 (zero false positives); recovers 1/5 intrinsic-only-R (GCA_000692095.1). `wiki/expression_frontier_isaba1_falsifier_2026-06-10.md`.
- Test suite: full suite green pre-change (967 passed at last run); the override default must keep that at 0 regressions.

## Verification Signal
- Promotion gate (ratified + brainstorm-hardened): `s_upgrades == 0` AND `r_rescues >= 1` AND `n_S >= 15`. The `r_rescues >= 1` clause guards against an inert detector (never-firing) trivially passing the zero-S-upgrade test. The artifact reports the Wilson 95% upper bound on the false-upgrade rate next to `0/n_S` (so "specific" is bounded, not asserted).
- Module unit tests pass: junction-positive to signal R; no-junction to no signal; multi-copy (>5 ISAba1) to no truncation; missing blastn to `status:unavailable` (no raise). (Orientation cases belong to the optional-refinement re-falsification, not the primary.)
- Caller: genome-absent to ABSTAIN unchanged; non-Acinetobacter / CALIBRATED organisms unaffected; override fires only when an explicit genome/accession input is given AND the registry `expression_context.enabled` is true AND verdict EXPRESSION_FLOOR AND the primary signal is True.
- Full test suite: 0 regressions.
- If Step 1 cannot build >=10R/10S independent free, documented "floor not independently validatable on free data" (an honest terminal, not a failure to hide).

## Implementation Steps

### Step 1: Build the independent Acinetobacter meropenem cohort (de-risk FIRST)
Files: scripts/build_acinetobacter_indep_cohort.py, data/raw/acinetobacter_meropenem_indep/selected.tsv
Depends on: none

**What changes:**
- New script reusing `organism_drug_validate.ensure_run`/`select_cohort`/`latest_metadata_url` plus `independent_cohort_validate._exclusion_set` to fetch up to 15R/15S Acinetobacter meropenem strains from NCBI Pathogen Detection AST, EXCLUDING the 30 cohort-1 accessions in `data/raw/acinetobacter_meropenem/selected.tsv`.
- Download assemblies (refseq cache) plus run AMRFinder cached/restartable; write `selected.tsv` (acc TAB R/S).
- HARD KILL GATE: if fewer than 15S OR fewer than 10R have downloadable assemblies free, emit `COHORT_INFEASIBLE` plus stop (the anchor binding risk; documents floor-not-validatable rather than faking a cohort). Target n_S>=15 (the promotion gate's S floor).

**Test strategy:**
- Manual/integration: run the script; assert `selected.tsv` has >=15S and >=10R, all disjoint from cohort-1 (no accession overlap). No unit test (network plus Docker fetch).

### Step 2: expression_context detector module (PRIMARY = exact frozen rule; all-hits)
Files: dna_decode/eval/expression_context.py
Depends on: none

**What changes:**
- New module `detect_is_upstream_junction(genome_fasta, *, is_ref, target_ref, upstream_bp=400, min_identity=85, min_len=120) -> dict`. Reuses `vf_runner.find_blastn`/`_find_makeblastdb`; builds an asm blastn DB; queries BOTH refs with `-outfmt "6 qseqid sseqid sstart send pident length"` (sseqid = contig, REQUIRED for the same-contig test), NO `-max_target_seqs` cap (or a high cap, e.g. 10000) so multi-copy ISAba1 is never truncated.
- PRIMARY rule = the EXACT frozen falsifier rule: `signal=True` iff an IS hit lies within `upstream_bp` 5-prime of a target hit ON THE SAME CONTIG (`sseqid`), strand-aware on the OXA hit, with NO ISAba1-orientation constraint. This is the rule that produced 1/15 R, 0/15 S — porting it verbatim keeps the validated signal identical to the falsified one.
- Returns `{status, signal: bool, evidence: {n_is_hits, n_target_hits, junction:{contig, distance_bp}, raw_hits:[...]}}`. `raw_hits` retains the BLAST rows + the extracted OXA-upstream interval (reproducibility / off-target audit at 85% id).
- Optional refinement (DEFAULT OFF): an `is_orientation` parameter that adds the ISAba1-orientation constraint. It MUST NOT be enabled in the primary path until separately re-falsified on the original 30 (must still rescue GCA_000692095.1) — guarded by a docstring contract + a test that the default path ignores orientation.
- Offline-safe: missing blastn/makeblastdb/ref to `{status:"unavailable", reason, signal:False}` (never raises). Thresholds frozen from the falsifier (400 bp, 85% id, 120 len).

**Test strategy:**
- Unit (synthetic contigs, no network where possible): (a) ISAba1 placed 200 bp upstream of OXA, same contig to signal True; (b) ISAba1 1 kb away to False; (c) ISAba1 on a DIFFERENT contig to False (same-contig test); (d) no ISAba1 to False; (e) >5 ISAba1 copies one of which is upstream to True (proves no truncation); (f) blastn-bin forced absent to status unavailable, no raise.

### Step 3: Unit tests for the detector module
Files: tests/test_expression_context.py
Depends on: Step 2

**What changes:**
- Implement the 6 cases above using committed `data/isaba1_ref/` refs plus synthetic assembly fixtures written to tmp_path. skipif-no-BLAST on the real-blastn cases; the offline-degrade case runs without BLAST. Assert `raw_hits` is populated on positive cases.

**Test strategy:**
- This step IS the test; assert signal booleans plus status plus evidence/raw_hits fields per case.

### Step 4: Thread genome into the caller plus gated ABSTAIN to R override
Files: dna_decode/eval/amr_rules.py, dna_decode/amr/cli.py, dna_decode/data/calibrated_amr_rules.json
Depends on: Step 2

**What changes:**
- `call_resistance(...)` plus `_call_resistance_calibrated(...)` gain optional `genome_fasta: Path | None = None` (backward-compatible default). The genome path is an EXPLICIT input — NEVER inferred from the `main.tsv` directory layout (auditability).
- In the EXPRESSION_FLOOR branch: iff `genome_fasta` is provided AND the registry entry carries `expression_context: {enabled: true, is_ref, target_ref, upstream_bp}` AND `detect_is_upstream_junction(...).signal` is True, return `prediction:"R"` with `rule="expression_context_v1 (ISAba1-upstream-of-OXA51 junction)"` plus the evidence (incl. raw_hits) in the record; ELSE unchanged ABSTAIN.
- Registry: add an `expression_context` block to the `Acinetobacter|meropenem` entry with `enabled: false` + `experimental: true` (ships OFF — flipped to opt-in only by Step 7 after the Step 6 gate passes; never auto-default-on).
- `cli.py` genome-mode passes `args.genome_fasta` into `call_resistance`.

**Test strategy:**
- Covered by Step 5 (regression tests).

### Step 5: Caller regression tests for the gated override
Files: tests/test_amr_rules_expression_override.py
Depends on: Step 4

**What changes:**
- Tests: (a) `enabled:false` (default registry) plus genome present to ABSTAIN unchanged; (b) `enabled:true` (test fixture registry) plus junction-positive genome to R with expression_context rule; (c) `enabled:true` plus no-junction genome to ABSTAIN; (d) genome absent to ABSTAIN regardless (no inference from main.tsv dir); (e) a CALIBRATED organism (Campylobacter) path unaffected by the new param.

**Test strategy:**
- This step IS the test; use a temp registry plus synthetic genomes (reuse Step 3 fixtures). skipif-no-BLAST on the junction cases.

### Step 6: Independent-cohort eval overlay plus artifact plus promotion gate
Files: scripts/expression_context_validate.py, wiki/expression_context_acinetobacter_validation_2026-06-10.md
Depends on: Step 1, Step 2

**What changes:**
- New script: for each Step-1 cohort strain, resolve its assembly FASTA from the refseq cache BY VERSIONED ACCESSION (e.g. GCA_000692095.1; provenance-safe, NOT run-dir adjacency — closes the cached-run mis-pair gap), run the PRIMARY `detect_is_upstream_junction`, and tabulate against the R/S labels restricted to the prior-ABSTAIN set.
- Endpoint: `abstain_rescue_rate` (upgraded-R / total-R), `r_rescues` (count of true-R upgraded), `s_upgrades / n_S`, the Wilson 95% upper bound on the false-upgrade rate, plus per-strain evidence. Emit `.md` plus `.json` retaining raw BLAST rows, the extracted upstream interval, SHA256 of the FASTA + both query refs, and the accession version (reproducible junction calls).
- Verdict: `PROMOTE` iff `s_upgrades == 0` AND `r_rescues >= 1` AND `n_S >= 15`; else `HOLD` (plus reason). Eval-only — does NOT flip the registry (that is Step 7, contingent).

**Test strategy:**
- Manual/integration on the real cohort. Unit-test the pure tabulation function (`compute_rescue_endpoint(signals, labels)`) with synthetic inputs (0-S-upgrade + >=1 R-rescue + n_S>=15 to PROMOTE; 1-S-upgrade to HOLD; 0 R-rescues to HOLD inert; n_S<15 to HOLD).

### Step 7: Promotion (contingent, opt-in only) plus docs plus ledger
Files: dna_decode/data/calibrated_amr_rules.json, CHANGELOG.md, LESSONS_LEARNED.md, project_state/eukaryotic-trait-decoding-cycle-2026-06-07.md
Depends on: Step 5, Step 6

**What changes:**
- IFF Step 6 verdict is PROMOTE: promote to EXPERIMENTAL OPT-IN — leave `expression_context.enabled:false` + `experimental:true`, and expose an explicit opt-in surface (e.g. a `--expression-context` CLI flag that turns it on per-call). DEFAULT behavior stays ABSTAIN; the R override never fires without the explicit opt-in. Default-on is DEFERRED until n_S is materially larger (e.g. n_S>=50, narrowing the Wilson bound). IF HOLD: leave OFF; document the held result.
- CHANGELOG plus LESSONS entry (incl. the Wilson bound + opt-in posture); ledger `--append-action` row recording the rescue-rate, r_rescues, s_upgrades/n_S, Wilson bound, and gate verdict.

**Test strategy:**
- Full suite green (0 regressions) after the flag state is finalized; the registry-load test covers the new `expression_context` block schema (incl. `experimental`) either way.

## Execution Preview
- Wave 0: Step 1, Step 2
- Wave 1: Step 3, Step 4, Step 6
- Wave 2: Step 5
- Wave 3: Step 7
- Total waves: 4. Max parallelism: 3 (Wave 1). Critical path: Step 2 to Step 4 to Step 5 to Step 7 (4 steps). Step 1 (cohort fetch, long-running) runs in Wave 0 and gates only Step 6.

## Risk Flags
- Cohort infeasibility (binding): the independent Acinetobacter meropenem cohort may not be buildable free (the NCBI `assembly_accession` wall). Step 1 HARD KILL GATE surfaces this before any wiring effort; the honest terminal is "floor not independently validatable," not a forced result.
- Partial ISAba1 ref under-detection [inferred]: the 570 bp partial mainly costs SENSITIVITY (missed junctions leave ABSTAIN unchanged), not specificity — acceptable for an "upgrade only when detected" signal. Off-target detection vs related IS elements at 85% id is the specificity risk; the artifact retains raw rows to audit it.
- Assembler FN ceiling [speculative]: IS elements are repetitive; short-read assemblies fragment them, so a junction across a contig break is missed regardless of the detector. Report as a sensitivity caveat; do not attribute to the feature.
- Low rescue ceiling + wide Wilson bound [grounded]: falsifier showed 1/5 FN recovered; at 0/15 S the Wilson 95% upper bound is ~20%, so the R override ships OPT-IN, never default-on, until n_S grows. Frame success as bounded-specificity + quantification, never as cracking expression resistance.
- Orientation-refinement drift [grounded]: adding ISAba1-orientation changes the falsified detector; it stays OFF in the primary and must be re-falsified on the original 30 (preserve GCA_000692095.1) before any use.

## Open Questions
None blocking. Ratified: acceptance bar = 0 S-upgrades + >=1 R-rescue; v0 scope = ISAba1 to OXA-51 only (OXA-58/ampC/efflux-regulator generalization explicitly deferred); promotion = experimental opt-in, default-on deferred to larger n_S.

## Verification
- All Verification-Signal items met: module tests green; caller default-behavior unchanged (genome-absent to ABSTAIN); 0 full-suite regressions; Step 6 endpoint emits `s_upgrades`/`r_rescues`/`abstain_rescue_rate`/Wilson-bound + reproducible raw hits; registry ships `enabled:false` + `experimental:true`.
- Promotion only on Step 6 `PROMOTE` (0 S-upgrades AND >=1 R-rescue AND n_S>=15) and only to OPT-IN. HOLD to documented, override stays off.

## Save-time amendments

Captured at: 2026-06-10
Source: `/save-plan` arguments + pre-exec `/brainstorm` (folded into the Implementation Steps above; this section is provenance for human readers — `/execute-plan` reads only `## Implementation Steps`).

- exact-falsifier-rule is the PRIMARY detector (Step 2); ISAba1-orientation is a separate default-off refinement requiring re-falsification on the original 30 (preserve GCA_000692095.1) — avoids validating a different rule than was falsified.
- gate adds `r_rescues>=1` (Step 6) — guards against an inert never-firing detector trivially passing the zero-S-upgrade test.
- gate adds Wilson 95% upper bound on the false-upgrade rate + raises the S floor to `n_S>=15` (Step 6).
- promotion = experimental OPT-IN, NOT default-on (Step 7); default-on deferred until n_S materially larger (~50).
- BLAST outfmt adds `sseqid` for the same-contig test (Step 2); genome input is explicit, never inferred from main.tsv adjacency (Step 4); validation artifact retains raw BLAST rows + extracted interval + SHA256(FASTA,refs) + accession version (Step 6).
<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
