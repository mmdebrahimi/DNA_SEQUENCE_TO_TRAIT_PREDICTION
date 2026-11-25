# Price The Re-Opened Epigenomics Goal On Measured Compute And Measured Substrate — Technical Plan

> Turn a completed T4 capacity measurement plus a verified-reachable ENCODE API into a decidable go/no-go on the goal abandoned 2026-05-12, without building or training anything.

## Lens status
Inputs: conversation
Degradations: repo-index absent (`repo-index-status.sh` not on this host) so Step 1 research was direct file reading with no Layer-2 summaries; graphify-out absent; `project-rules.md` absent so Project Rules Check is omitted; Architecture Impact omitted as non-structural (one new script, no shared-contract change) with sentrux available and unused

## Problem Statement
On 2026-05-12 this project gave up its most ambitious goal — "a tissue-specific gene-regulatory-element activity predictor using a transformer-based DNA model trained on ENCODE / Roadmap Epigenomics multi-tissue data (Enformer-shaped)" — on LOCAL compute grounds (GTX 860M, 4 GiB). Constraint `local-gpu-vram-blocks-dna-lm-finetune` recorded that cut with a `gpu-vram-gte 8` expiry test which is measurably unable to fire: `ConstraintRegister.check_all(include_resolved=False)` returns 11 rows and EXCLUDES it because the record carries `resolved=true` while `last_checked` and `last_verdict` are both empty, and forced in with `include_resolved=True` it returns `unknown` ("no GPU visible / nvidia-smi unavailable") because the compute that actually arrived is CLOUD.

Two claims were bundled in that cut and only one has now changed. "No >=8 GB GPU is available" is FALSE and measured: `wiki/encode_sizing_t4_2026-11-24.json` records 44 of 45 configurations fitting on a real Tesla T4, including Enformer's exact 196,608 bp context at 1024 channels / 11 layers / 5,313 tracks / 134.63M params in 10.86 GiB of 14.56 at 1028.3 ms per step. "Therefore the goal is infeasible" does NOT follow and remains unmeasured, because sufficiency needs a sourced training budget and a usable substrate, neither of which exists in-repo.

This plan prices both missing halves and stops. It produces the memo that discharges `encode-roadmap-regulatory-goal-unsized-on-cloud`, a measured ENCODE substrate census, and a gate-screen verdict from the project's own rejection gates, so the re-open decision becomes a reading rather than an argument. It does NOT build, train, fine-tune, or fetch training data.

## Codebase Context
`dna_decode/eval/rejection_gates.py` exposes `screen_candidate(packet)` with the constants `PASS` / `TRIP` / `NOT_APPLICABLE` / `INSUFFICIENT_DATA` and the layer tokens `L1_AMR_RS` / `L4_FORWARD_CONTINUOUS`; `scripts/screen_candidate_gates.py` drives it, encodes hand-worked candidates for `--verify`, and takes `--packet <json>` for a new one. Gates G1 (wet-lab versus tool-derived label) and G3 (assay reading versus collection context) are JUDGMENTS that require an explicit human evidence string and REFUSE the whole verdict without one, so Step 4 cannot be fully automated and the plan says so rather than discovering it at runtime.

`dna_decode/eval/regime.py` carries a 9-entry `REGIMES` tuple of `Regime` records with `population` / `endpoint` / `method` / `verdict` / `evidence` / `artifact` / `split_unit` / `organism_transfer`. The cell this candidate lands in is `natural_organism_supervised`, verdict `REQUIRES_DECONFOUNDING`, whose recorded condition is "report WITHIN-GROUP performance against each group's own null" — not a closed negative, which is why the screen is worth running at all.

The constraint register is `soraya_runs/constraint-expiry.json` (a `schema_version` key plus a `constraints` dict keyed by cid), driven by `ConstraintRegister` in the Soraya skill's `scripts/constraint_expiry.py`. Kinds are `gpu-vram-gte` / `file-exists` / `test-exit-0` / `manual`; `check_all` runs only the read-only executable kinds and skips `resolved` records by default. The target constraint uses `file-exists wiki/encode_roadmap_cloud_sizing.md` and currently reads `still-binding`.

The Kaggle surface is verified this session: `scripts/kaggle_push_poll.py {push,status,pull}` with `--gpu` (which sets `machine_shape: NvidiaTeslaT4`, closing the P100 trap) and `--dataset`. Tests live under `tests/` with `testpaths = ["tests"]` and `addopts = "-v"`; relevant neighbours are `tests/test_cell_regime.py` and the eight `tests/test_constraint*.py` files.

ENCODE's surface was verified live rather than inherited: a search request for `type=Experiment` with `assay_title=DNase-seq`, `limit=1` and `format=json` returns HTTP 200 `application/json` with `@graph` records carrying `accession`, `assay_title`, `biosample_ontology` and `files`, and `total` = 3,582.

### Reusable-Code Survey
- `scripts/oxford_w0_probe.py` — the W0 pattern Step 1 mirrors: pin an unverified external schema and its key cardinalities BEFORE building against it, with pure summarise helpers split from the fetch so they unit-test offline.
- `scripts/essentiality_ladder_w0_probe.py` — second instance of that pattern, and the precedent for a probe whose verdict gates every downstream parser.
- `dna_decode/eval/rejection_gates.py` — the ten-gate screen Step 4 consumes directly; already carries the judgment-gate refusal behaviour.
- `scripts/screen_candidate_gates.py` — the `--packet` entry point plus the `--verify` self-check that re-derives committed hand verdicts.
- `dna_decode/eval/regime.py` — the regime map Step 4 places the candidate on; read-only here.

Also searched, no candidates: graphify-out/GRAPH_REPORT.md (absent), src/, lib/, utils/, common/, helpers/ (all absent).

## Pre-Change Baseline
- `ConstraintRegister(".").check_all(include_resolved=False)` returns **11 rows**; `encode-roadmap-regulatory-goal-unsized-on-cloud` reads **still-binding**, detail "wiki/encode_roadmap_cloud_sizing.md still absent".
- `wiki/encode_roadmap_cloud_sizing.md` is **absent** — it is the expiry target.
- Measured compute envelope, committed at `wiki/encode_sizing_t4_2026-11-24.json`: 45 configs tried, **44 fitting**; largest fitting is 196,608 bp / 1024 channels / 11 layers / 5,313 tracks / 134.63M params / **10.86 GiB** / **1028.3 ms per step** / about **31,509 steps per 9 h**; the sole failure is that same config in fp32 (CUDA OOM).
- **Zero** in-repo artifacts carry a sourced training budget for this model class. Record `grep -ril "enformer" wiki/ | wc -l` before Step 2 as the baseline count.
- ENCODE reachability, measured: HTTP 200, `total` = **3,582** DNase-seq experiments.
- `dna_decode/eval/regime.py` carries **9** regime records; none has `organism_transfer='held_out_organism'` for a natural-population long-context epigenome target.
- Full-suite baseline failure SET captured by `uv run pytest tests/ -q` BEFORE Step 1, because the standing practice here is to compare failure SETS rather than counts.

## Verification Signal
- Re-running `check_all(include_resolved=False)` after Step 3 returns the same 11 rows with the target constraint flipped from **still-binding** to **expired**.
- `wiki/encode_substrate_census_2026-11-24.json` reports a count of aligned multi-tissue signal FILES under a named assay filter, distinct from and smaller than the 3,582 experiment count, plus the distinct `biosample_ontology` term count — a measured number, not an estimate.
- `scripts/screen_candidate_gates.py --packet` returns a verdict in which no gate reads `INSUFFICIENT_DATA` that the census could have answered, and the two judgment gates either carry an explicit evidence string or are reported as refused.
- `wiki/encode_roadmap_cloud_sizing.md` states measured capacity beside a sourced budget with its citation, and names what remains unmeasured.
- `uv run pytest tests/ -q` failure SET is identical to the baseline, with new passes equal to the number of tests added.

## Implementation Steps

### Step 1: Measure the ENCODE substrate with a W0-pattern census probe
Files: scripts/encode_substrate_census.py, tests/test_encode_substrate_census.py, wiki/encode_substrate_census_2026-11-24.json
Depends on: none

**What changes:**
- New probe mirroring `scripts/oxford_w0_probe.py`, with pure summarise helpers split from the HTTP fetch so they unit-test offline against committed fixtures and never touch the network.
- Query the verified ENCODE search API for the assay families a regulatory-activity target needs, each named explicitly in the artifact rather than implied: DNase-seq, ATAC-seq, histone-mark ChIP-seq, and total RNA-seq.
- Per assay record the experiment count, the distinct `biosample_ontology` term count (the TISSUE axis, which is what "multi-tissue" means), the genome assembly values present, and the count of released aligned-signal bigWig FILES — counting files rather than experiments, because the two differ by roughly an order of magnitude and conflating them would inflate the substrate.
- Record the total download volume implied by the file list so data reachability gets a number in bytes rather than an adjective.
- Refuse and exit non-zero if any assay returns zero experiments or if `total` is absent from a response, because a census that silently reports zero is indistinguishable from a broken query.

**Test strategy:**
- Unit tests over the pure summarisers using a committed fixture of two trimmed `@graph` records: assay tally, tissue-term dedup, the file-versus-experiment distinction, and the byte-sum arithmetic.
- A non-vacuity test feeding a fixture whose `total` key is missing, asserting the probe refuses rather than reporting zero.
- The live run is manual and its artifact is the deliverable; tests never hit the network.

### Step 2: Source the published training budget for this model class
Files: wiki/encode_roadmap_training_budget_sources.json
Depends on: none

**What changes:**
- Record, with a resolvable citation per figure, the published training cost of the reference architectures (Enformer and any successor whose paper states one): hardware type and count, wall-clock, parameter count, input context, and number of output tracks.
- Every figure carries `source_url` plus a `verbatim_quote`, and a figure that cannot be sourced is recorded as `unsourced: true` carrying no number — mirroring the refusal discipline already shipped in `dna_decode/data/ecoff_catalog.py`, where a value without a locator and a quote is rejected rather than approximated.
- Compute and record the one derived ratio the decision needs, reference training compute divided by one measured T4 session, labelled `derived` and naming both inputs.

**Test strategy:**
- A schema test asserting every numeric figure carries both `source_url` and `verbatim_quote`, and that an `unsourced` entry carries no number.
- A test that the derived ratio recomputes from its two named inputs, so it cannot drift from them.

### Step 3: Write the sizing memo that discharges the compute constraint
Files: wiki/encode_roadmap_cloud_sizing.md
Depends on: Step 2

**What changes:**
- State the measured envelope from `wiki/encode_sizing_t4_2026-11-24.json` verbatim, including that the single non-fitting config is fp32 at the largest shape, so half precision is the enabling lever rather than an optimisation.
- Record the two measured levers separately: AMP as a throughput lever on this card (2683.7 ms per step in fp32 against 191.7 under AMP at 384 channels, roughly a fourteen-fold difference from the sm_75 tensor cores) and gradient checkpointing as near-neutral here (10.86 down to 9.53 GiB, slightly slower) because memory stopped binding once AMP was on.
- Place that capacity beside the Step-2 sourced budget and state the ratio, so "sufficient" becomes a reading rather than a claim.
- Carry four named non-conclusions: random-tensor inputs so no accuracy claim exists; nothing trained to convergence; 134.63M params against the reference's larger count, so this reaches the CONTEXT and not the scale; and one free T4 in one session with no paid or multi-GPU scaling.
- State explicitly that this memo's existence is what expires the constraint under `file-exists` semantics, so a pessimistic memo expires it too — the test records that the question was answered, not that the answer was favourable.

**Test strategy:**
- Manual read-back against the JSON: every number in the memo must appear in `wiki/encode_sizing_t4_2026-11-24.json` or `wiki/encode_roadmap_training_budget_sources.json`.
- `uv run python scripts/contract_number_audit.py` to confirm no new uncited figure enters the trust surface.

### Step 4: Screen the candidate through the ten rejection gates and place it on the regime map
Files: wiki/encode_roadmap_gate_screen_2026-11-24.json
Depends on: Step 1, Step 3

**What changes:**
- Build an evidence packet from the Step-1 census and the Step-3 memo, run `scripts/screen_candidate_gates.py --packet`, and record the per-gate verdict.
- Supply the two judgment gates explicitly: G1 because ENCODE signal tracks are instrument readings from a wet-lab assay rather than tool-derived calls, and G3 because the readout is an assay measurement rather than a collection-context proxy. Each gets an evidence string and an explicit assertion, or the screen's refusal is recorded as the result.
- Record which `regime.py` cell the candidate occupies (`natural_organism_supervised`, `REQUIRES_DECONFOUNDING`) and its recorded condition, WITHOUT editing `regime.py` — adding a row is a separate authority call and is deliberately out of scope.
- State the verdict's scope: a gate PASS bounds the LABEL question only, never artifact reachability and never whether the work is worth doing.

**Test strategy:**
- `uv run python scripts/screen_candidate_gates.py --verify` must still reproduce both committed hand verdicts, proving the new packet did not perturb the screen.
- A test asserting the packet omits no gate the census could answer, so an `INSUFFICIENT_DATA` cannot pass silently as a shrug.

### Step 5: Re-check the constraint register and record the outcome
Files: soraya_runs/constraint-expiry.json, project_state/glm-umbrella-2026-10-07.md
Depends on: Step 3, Step 4

**What changes:**
- Re-run `check_all(include_resolved=False)` and record the verdict flip for `encode-roadmap-regulatory-goal-unsized-on-cloud`, with the row count before and after.
- If the census or the budget ratio reveals a NEW binding constraint — a download volume beyond local disk, for instance — register it with its own firing expiry test rather than leaving the finding in prose.
- Append one AC9 Action-Log row recording what was measured, what the gate screen returned, and which questions remain open, with the class drawn from the seven-value enum.
- Leave the old `local-gpu-vram-blocks-dna-lm-finetune` record's `resolved` flag untouched, because that field asserts a user action and flipping it is not this plan's to do.

**Test strategy:**
- Read-back of the register showing the flipped verdict.
- `uv run pytest tests/ -q`, comparing the failure SET against the baseline.

## Execution Preview

Wave 0 (2 parallel):  Step 1 — ENCODE substrate census probe, Step 2 — source the published training budget
Wave 1 (1 parallel):  Step 3 — write the sizing memo
Wave 2 (1 parallel):  Step 4 — rejection-gate screen and regime placement
Wave 3 (1 parallel):  Step 5 — re-check the register and record the outcome

Critical path: Step 2 → Step 3 → Step 4 → Step 5 (4 waves)
Max parallelism: 2 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags
- Severity: high — The constraint's `file-exists` expiry fires on the memo's EXISTENCE, not on its conclusion, so a pessimistic memo expires it exactly as a favourable one does. That is intended, since the test records that the question was answered, but the HIGH-salience 2026-11-21 decisions-log lesson applies directly: a bar that flags a confound in its rationale and then enshrines that confound in its primary measurement produces a clean, unanimous and useless verdict. Step 3 therefore states the semantics inside the memo rather than leaving a future reader to infer that expiry meant approval.
- Severity: high — The 31,509-steps-per-session figure is partly a function of MY chosen model rather than of the reference architecture. The HIGH-salience 2026-11-16 lesson is that a derived statistic which is a function of the DESIGN rather than the data is worth nothing as a finding, and this repo has four recorded vacuous-control instances. Mitigation: Step 3 reports the raw per-config measurements and the sourced-budget ratio separately, never a bare sufficiency verdict, and names the 134.63M-versus-reference parameter gap.
- Severity: medium — Experiment counts and file counts differ by roughly an order of magnitude. 3,582 is an EXPERIMENT total, while a regulatory-activity target consumes aligned signal files. Step 1 counts files and tissue terms separately and labels each, so "3,582 tracks" cannot enter any downstream claim.
- Severity: medium — Two of the ten gates are judgments the screen refuses to decide. G1 and G3 need an explicit human evidence string, so Step 4 cannot be fully automated; the plan states this up front instead of discovering it as a runtime refusal.
- Severity: medium — External tool integrations. ENCODE search API verified live this session at `www.encodeproject.org/search/` with `type`, `assay_title`, `limit` and `format=json` parameters, returning HTTP 200 `application/json` whose body carries `total` and `@graph[]` with `accession` / `assay_title` / `biosample_ontology` / `files`. Kaggle verified via `scripts/kaggle_push_poll.py {push,status,pull}`, where `--gpu` sets `machine_shape: NvidiaTeslaT4`, confirmed by pulling the pushed kernel's own metadata. `ConstraintRegister` verified by `inspect.signature`: `register(cid, *, blocked, rationale, expiry_kind, expiry_arg)` and `check_all(include_resolved=False)`, kinds `gpu-vram-gte` / `file-exists` / `test-exit-0` / `manual`. The Roadmap Epigenomics download surface is **[unverified]** — only ENCODE was probed, so Step 1 must verify it or record its absence rather than assume parity.
- Severity: low — File overlaps between steps: none. Each step writes a disjoint file set, Step 5 is the only step touching two files and both are append-only surfaces, and no intra-wave overlap exists, so no dependency was added for that reason.
- Severity: low — Steps depending on external services: Steps 1 and 2 need network. Both are read-only GETs against public endpoints with no credentials, so neither reaches a money or publish gate.
- Severity: low — Steps sharing transitive imports not captured by the Files lists: Step 4 imports `dna_decode/eval/rejection_gates.py` and reads `dna_decode/eval/regime.py` without modifying either, so both are read-only dependencies rather than overlaps.
- Severity: low — No restructuring applied. Step order follows the data dependency: the budget must be sourced before capacity can be placed beside it, and the memo must exist before the register can flip.
- Severity: low — `uv run` must execute with cwd at the project root, per the standing rule against running it inside another project.

## Open Questions
- Should the sizing memo conclude with a recommendation, or present the measured envelope and the sourced ratio and stop? A memo that recommends is more useful and also harder to read neutrally later; the plan currently assumes present-and-stop.
- Does a reduced-context or fine-tune-only variant count as the original goal re-opened, or as a different goal? This is acceptance-bar shaped, since the 2026-05-12 goal names an Enformer-shaped model trained on multi-tissue data and the measured envelope reaches its context at roughly half its parameter count.
- The old `local-gpu-vram-blocks-dna-lm-finetune` record carries `resolved: true` with `last_checked` and `last_verdict` empty, marked acted-upon with no verdict recorded. Step 5 leaves it untouched because that field asserts a user action, so it needs your call.
- If the Step-1 census reports a download volume beyond this host's roughly 26 GB of free C: space, is a Kaggle-side-only data path acceptable, or does the substrate need to land locally first?

## Verification
Run these in order after execution. First, call the Soraya register and confirm `check_all(include_resolved=False)` still returns 11 rows with the target constraint reading `expired` rather than `still-binding`. Second, read `wiki/encode_substrate_census_2026-11-24.json` and confirm it carries a file count and a distinct-tissue-term count per named assay, each labelled, with the experiment total recorded separately. Third, run `uv run python scripts/screen_candidate_gates.py --verify` to confirm both committed hand verdicts still reproduce, then read the new packet's per-gate verdicts. Fourth, confirm every number in `wiki/encode_roadmap_cloud_sizing.md` appears in either the sizing JSON or the budget-sources JSON, and run `uv run python scripts/contract_number_audit.py`. Fifth, run `uv run pytest tests/ -q` and compare the failure SET, not the count, against the baseline captured before Step 1.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high,test-strategy-leak -->
