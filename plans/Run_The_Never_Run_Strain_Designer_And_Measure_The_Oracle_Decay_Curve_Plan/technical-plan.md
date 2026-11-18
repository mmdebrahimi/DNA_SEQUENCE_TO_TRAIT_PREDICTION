## Lens status
Repo-index: absent (`repo-index-status.sh` not installed) — source-only reading, no degradation.
sentrux: available — plan is NOT structural (adds scripts + artifacts, no module boundary change), so the architecture gate is advisory only.
gh: available. codex: absent — no external review path on this host.
project-rules.md: absent — `## Project Rules Check` omitted.
DESIGN.md: absent — no UI surface in this plan.

## Problem Statement

The ranked research memo (`wiki/generator_hypotheses_ranked_researched_2026-11-18.md`) produced an ordered work list whose first two items are now grounded, and whose third is not. This plan executes items 1 and 2 and resolves item 3's substrate question instead of assuming it.

**Item 1 — run the built-but-never-run growth-coupled strain designer.** `dna_decode/fba/design.py` implements "given a PRODUCT, which gene knockouts make producing it NECESSARY for growth?", is already CLI-wired (`fba/cli.py:155,159`) and already tested (`tests/test_fba_design.py`), yet has **no result artifact in `wiki/`**. It has never been run on the real model.

**Item 2 — measure the oracle's accuracy decay against distance from its training data.** CORRECTED 2026-11-18 after reading the primary sources (`wiki/generator_hypotheses_source_notes_2026-11-18.md`): the earlier claim here, that published work ESTABLISHES the phenomenon, was **wrong**. The edit-distance result cited for it measures training-data INFLUENCE (via Kronfluence), not accuracy, and that paper never plots performance against distance at all; the other support is on SIMULATED landscapes with Pearson, counting regimes from a seed rather than distance to the nearest training sequence. So the decay is **plausible but NOT established for the quantity this step measures** — which makes the measurement MORE novel, not less: two papers cited as establishing it instead confirm it is un-done. This repo has never plotted it either.

**Item 3 — the published expression-to-growth law.** CORRECTED 2026-11-18: calling it "parameter-free" was half right. The growth-versus-expression law itself has zero free and two FIXED parameters, but the promoter-activity half adds **one fitted parameter**, overpredicted in minimal medium and needed a correction term, and its best-fit slope sits below the model's own predicted lower bound. Testing it needs promoter strength PAIRED WITH growth rate. That pairing is not verified present locally, so this plan probes substrate reachability rather than scheduling a test that may be unrunnable.

Non-goals: no new wet-lab labels; no edits to the five sha256-pinned frozen-surface files; no reopening the natural-population zero-shot embedding track.

## Codebase Context

**`dna_decode/fba/design.py`** — public surface verified by reading: `classify_coupling` (returns `OBLIGATORY` / `POSSIBLE` / `INFEASIBLE` against `FLUX_TOL = 1e-6`), `Design`, `improves_on_baseline`, `rank_designs`, `biomass_reaction`, `resolve_target`, `product_range_at_growth`, `apply_knockout`, `evaluate_knockouts`, `set_anaerobic`, `competition_ranking`, `find_coupled_designs` (enumeration) and `find_coupled_designs_milp`. The coupling test is a two-sided LP at a growth floor: `min_flux > 0` means every growing flux distribution also secretes the product.

**Two solver paths, and the cheap one has a MEASURED defect.** `pyproject.toml:140-147` records that the bounded enumeration is exhaustive only at depth 1 and is blind to designs whose members are individually unremarkable, which it says is most real ones, citing two measured failed attempts to fix it by pre-ranking. It also records that SCIP is REQUIRED not optional, because GLPK lacks indicator constraints and the straindesign big-M fallback returned `unbounded` on iML1515 (measured 2026-08-07).

**Extras.** `[fba]` is `cobra>=0.29` alone, which bundles GLPK via swiglpk. `[design]` is `cobra` plus `straindesign>=1.19` plus `pyscipopt>=5`, and pulls matplotlib.

**Environment, measured now.** `cobra` is ABSENT (`ModuleNotFoundError`), so `design.py` cannot execute today. The model IS present: `~/.cache/dna_decode/fba/iML1515.xml.gz` (620,180 bytes). `fba/model.py` resolves it from a three-tier path whose third tier is that cache directory.

**CLI flags, read from `fba/cli.py`:** `--design-target PRODUCT`, `--max-knockouts` with choices 1, 2, 3 (default 3), `--milp`, `--milp-floor` (default 0.5), `--anaerobic`, `--candidates`, `--json`.

**H2 substrate.** No protein DMS table is local — a search for dms/pear/blatem/proteingym under `data/` returns exactly one hit, `data/forward_ref/blatem_3349172526.fna`, a reference FASTA with no measured values. What IS local: `data/mpra_gse108535/` (GSE108535, 10,898 sequence-activity pairs, every sequence exactly 150 bp, five element-identity columns) and ten fine-tuned encoder checkpoints at `data/promoteratlas/finetune_checkpoints/seed_*.pt` (22,305,931 bytes each).

### Reusable-Code Survey

- **`dna_decode/glm/expression.py`** — `load_pairs`, `leave_element_out`, `FEATURES`, `spearman`, `SEQ_LENGTH = 150`. Owns the GSE108535 substrate; its docstring records that a naive join silently loses 65% of the data over a hyphen-vs-underscore mismatch, so this loader must be reused, never reimplemented.
- **`scripts/glm_selfdistill_finetuned_teacher.py`** — `load_finetuned_teacher(ckpt, pa_repo)` returns a working chunked scorer for the fine-tuned encoder, built and verified yesterday at 66 seq/s.
- **`dna_decode/fba/design.py`** plus **`fba/cli.py`** — the designer and its CLI; nothing to rebuild.
- **`tests/test_fba_design.py`**, **`tests/test_fba_compose.py`** — existing test patterns for this cell.
- Searched: `dna_decode/` (all 12 top-level packages listed), `scripts/`, `tests/`, `pyproject.toml`, `data/`. `graphify-out/GRAPH_REPORT.md` absent.

## Pre-Change Baseline

- `dna_decode/fba/design.py` result artifacts in `wiki/`: **0** (`ls wiki/fba_design*.json` returns no match).
- `cobra` importable: **no** (`ModuleNotFoundError`), so the designer's measured output today is nothing.
- Spearman-vs-distance curve for any oracle in this repo: **absent** (only incidental mentions in `wiki/generator_hypotheses_2026-11-18.md` and `wiki/yeast_cnv_attribution_result_2026-07-02.md`).
- Suite: 12 failed / 8,899 passed / 210 skipped (2026-11-18, failure set equals the standing absent-data set).
- Disk: C: 25 GB free of 476 GB (95% used).

## Verification Signal

- A committed `wiki/fba_design_*.json` whose records carry `OBLIGATORY` / `POSSIBLE` / `INFEASIBLE` from the deployed `classify_coupling`, plus the model id and growth floor.
- A committed decay artifact with Spearman per distance bin AND the per-bin n, so a bin of 3 variants cannot read as a trend.
- Both solver paths reported side by side, with the enumeration path's documented depth-1 blindness stated next to its result rather than left for a reader to rediscover.
- Suite failure SET unchanged from the 12-test baseline (count may move only by tests this plan adds).

## Implementation Steps

### Step 1: Install the minimal FBA extra and confirm the designer imports
Files: uv.lock
Depends on: none

**What changes:**
- `uv sync --extra fba` (cobra only; GLPK arrives bundled via swiglpk, no separate solver).
- Confirm `import cobra` succeeds and `dna_decode.fba.model` loads iML1515 from the cache path.
- Record the resolved cobra version and whether `uv.lock` changed — it is a frozen reproducibility unit, so a lock diff is reported, not silently accepted.

**Test strategy:**
- `uv run pytest tests/test_fba_design.py -q` — previously collectible but engine-blocked; it should now exercise the real solver.

### Step 2: Run the enumeration design path and record the first artifact
Files: wiki/fba_design_enumeration_2026-11-18.json
Depends on: Step 1

**What changes:**
- Run `dna-fba --design-target <product> --json` at default max-knockouts 3 for a small set of canonical fermentation products, writing one artifact.
- Stamp into the artifact: model id, organism, growth floor, target reaction, every design's coupling class, and the depth-1 blindness caveat verbatim from `pyproject.toml`.
- Report `OBLIGATORY` count separately from `POSSIBLE` — a `POSSIBLE` design is explicitly NOT a design.

**Test strategy:**
- Assert the artifact parses and that every coupling value is one of the three constants exported by `design.py`, re-imported rather than retyped.

### Step 3: Install the MILP design extra, the expensive arm, declared up front
Files: pyproject.toml
Depends on: Step 2

**What changes:**
- `uv sync --extra design` (adds `straindesign` plus `pyscipopt` plus matplotlib).
- CARE pre-check before installing: C: is at 95% (25 GB free). Surface the figure; do not gate on it.
- Declared as an arm at plan time so a thin Step-2 result PROMOTES this rather than cancelling it.

**Test strategy:**
- Confirm `pyscipopt` imports and that straindesign reports SCIP as the active solver, since the recorded failure mode is a silent GLPK fallback returning `unbounded`.

### Step 4: Run the MILP path and compare it against the enumeration result
Files: wiki/fba_design_milp_2026-11-18.json
Depends on: Step 3

**What changes:**
- Run `dna-fba --design-target <same products> --milp --milp-floor 0.5 --json`.
- Record the solver actually used, the status, and the design count against Step 2's.
- If the MILP path finds designs the enumeration missed, that CONFIRMS the recorded depth-1 blindness on our own model; if it does not, record that as a finding about this model rather than retrying.

**Test strategy:**
- Assert the artifact names a solver and that a non-SCIP solver with an `unbounded` status is reported as the known-failure branch, not as an empty result.

### Step 5: Build the oracle decay-curve harness
Files: scripts/oracle_decay_curve.py
Depends on: none

**What changes:**
- New script reusing `glm/expression.py` `load_pairs` plus `leave_element_out` plus `spearman`, and `glm_selfdistill_finetuned_teacher.load_finetuned_teacher` — no reimplementation of either.
- Per held-out variant, compute minimum Hamming distance to any TRAIN sequence; bin; report Spearman against measured activity per bin with per-bin n and a Wilson-style interval.
- REFUSE to emit a trend when a bin's n falls below a floor fixed in the script before any run.

**Test strategy:**
- Offline unit tests on synthetic sequences: distance arithmetic, binning edges, and that a below-floor bin refuses rather than reporting.

### Step 6: Run the decay curve and record the artifact
Files: wiki/oracle_decay_curve_2026-11-18.json
Depends on: Step 5

**What changes:**
- Run on GSE108535 with one fine-tuned checkpoint; emit per-bin Spearman plus n.
- State the scope limit in the artifact: this measures decay on an MPRA promoter substrate, NOT a protein DMS substrate, and the encoder was fine-tuned on a DIFFERENT dataset (GSE144621), so distance here is to the held-in GSE108535 split rather than to the fine-tuning corpus.

**Test strategy:**
- Assert monotonicity is REPORTED rather than assumed, and that the artifact carries per-bin n.

### Step 7: Probe whether the expression-to-growth law is testable with local data
Files: wiki/expression_growth_substrate_probe_2026-11-18.json
Depends on: none

**What changes:**
- Read-only inventory of whether any local dataset pairs promoter strength WITH growth rate or fitness.
- Emit a reachability verdict with the roots searched, so the H5 item is either unblocked or recorded as substrate-blocked with the specific missing pairing named.
- Do NOT fetch anything; this step decides whether a fetch is warranted.

**Test strategy:**
- Assert the verdict is one of a closed enum and that an empty inventory yields `substrate-blocked` rather than a vacuous pass.

### Step 8: Test pass and suite reconciliation
Files: tests/test_oracle_decay_curve.py
Depends on: Step 4, Step 6

**What changes:**
- Add the offline tests for Step 5's pure helpers.
- Run the full suite and reconcile the failure SET against the 12-test baseline, not the count.

**Test strategy:**
- `uv run pytest tests/ -q`; any new failure is diagnosed before the batch is reported complete.

## Execution Preview

- Wave 0: Step 1, Step 5, Step 7
- Wave 1: Step 2, Step 6
- Wave 2: Step 3
- Wave 3: Step 4
- Wave 4: Step 8
- Total waves: 5. Max parallelism: 3. Critical path: Step 1 to Step 2 to Step 3 to Step 4, then Step 8.

## Risk Flags

- **The cheap design path is documented to fail on most real designs.** Step 2 may legitimately return zero `OBLIGATORY` designs; that is the recorded depth-1 blindness, not a broken run. Step 3 exists so that outcome promotes the MILP arm rather than closing the item.
- **`uv.lock` is a frozen reproducibility unit** and both install steps may touch it. Report any diff explicitly rather than committing it as incidental.
- **Disk at 95%** (25 GB free). `pyscipopt` plus matplotlib are the heaviest additions here. Surfaced per the Care pre-check; not a gate.
- **[unverified] external surfaces:** the `straindesign>=1.19` API, and that it reports its active solver, are taken from the `pyproject.toml` comment rather than from the package's own help output or docs — it is not installed yet, so Step 3 must verify before Step 4 relies on it.
- H2's result will be an MPRA-substrate decay curve. It does NOT transfer to the protein DMS regime — and CORRECTED 2026-11-18, that K=2 / K=4 / K=5 pattern was measured on **simulated NK landscapes with Pearson r**, not on real protein DMS, so it is a weaker comparator than this line originally implied.
- No external review path on this host (codex absent, LM Studio unreachable), so this plan has had no second reader.

## Open Questions

- RESOLVED 2026-11-18 by draft-then-ratify, recorded in the execution state file: **lactate first**, because a known-positive makes a null result diagnostic rather than ambiguous. Succinate and ethanol are the other canonical growth-coupled benchmarks. CAVEAT added after the source-reading: the 2005 lactate validation is **UNVERIFIED from primary source** — three independent routes to that review were blocked from this host — so it remains search-summary level and must not be quoted as established. This does not weaken Step 2, which rests on verified repo facts rather than on the literature.
- Is max-knockouts 3 the right depth for the first run, given enumeration is exhaustive only at depth 1? Running depth 1 first would separate "exhaustive and found nothing" from "sampled and missed".

## Verification

- `wiki/fba_design_enumeration_*.json` and `wiki/fba_design_milp_*.json` both exist and parse, each naming its solver and growth floor.
- `wiki/oracle_decay_curve_*.json` carries per-bin Spearman AND per-bin n.
- `wiki/expression_growth_substrate_probe_*.json` carries a closed-enum verdict and the roots searched.
- Full-suite failure SET identical to the 12-test baseline, modulo tests this plan adds.
- The five frozen-surface files byte-unchanged; prospective lock re-verified.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified -->
