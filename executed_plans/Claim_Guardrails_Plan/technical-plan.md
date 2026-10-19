# Claim Guardrails — Technical Plan

> Stop rung-scoped and door-closing claims from escaping as global facts, by making the absent field impossible, the judge's payload factual, and the trigger harness-bound.

## Lens status
Inputs: conversation (the ranked three-item guardrail proposal from this session, with item 2 redesigned on evidence gathered during Step 1 research)
Degradations: repo-index absent (`update-index` CLI not installed) — Step 1 used direct file reading, no Layer 2 consulted; sentrux unavailable — no architecture gate; graphify absent — reuse survey used filesystem + direct reads only; `project-rules.md` absent — Project Rules Check omitted

## Problem Statement

In one session, thirteen factual claims had to be corrected. The split is diagnostic and was measured, not estimated: **every number that passed through a code-owned gate survived scrutiny** (`0.3925`, `-0.0040` with its interval, the reconcile anchors `0.3000`/`0.3968`, the `-1.000000` consistency check, the `44,106` coverage count), and **every number that passed only through prose was wrong** ("twice", "the entire file", "~0.46-0.47", "2,095 min", "5.0 hours", "can only damage one row", "one sweep covers ten seeds").

Research into the worst of those — claiming the project had no training rows left when 321,123 fragments sat in the same cache directory — found a **structural cause, not a diligence failure**: in `wiki/mpra_substrate_census_2026-10-12.json`, the field `size_lever` is present on **exactly one rung of four** (`peak_tiles`, value `"EXHAUSTED -- 46,713 is the whole file"`). The three rungs that *do* have headroom carry no such field at all. `best_measured_arm` has the same 1-of-4 shape. An absent field is not a neutral absence — it reads as "nothing to say here" — so the single rung that spoke was the only size statement available, and it said exhausted.

Two designs were considered and **one was killed by evidence during research**: per-file orphan detection in `scripts/data_inventory.py` would *not* have caught this, because the fragment files **are** referenced in code (`dna_decode/glm/genomewide.py:48-49`). The correct signal is not "unreferenced" but "under-used relative to its size" — 30,000 of 321,123 rows, 9.3%.

Scope is the three items the conversation settled on: (A) make the census symmetric and code-validated, (B) convert the external review payload from narrative to a verifiable claims list, (C) bind a trigger to a harness event rather than to intention. Non-goals: a general "review everything" pass (rejected in conversation — it returns "looks fine" on plausible-but-wrong prose and trains both parties to skip it); any change to the four named limits; any new model integration (the review tool is already wired via `/brainstorm`).

## Codebase Context

- `scripts/data_inventory.py` (421 lines) already implements the right primitives: a `Dataset` dataclass with `consumed_by` / `reported_in` / an `orphan` property, and an `is_self_reference` guard enforcing that *describing* a dataset never makes it look *consumed*. Its limitation is **granularity**: `SCAN_ROOTS` keys on directories, so `mpra` reports as one row (104.7 MB, 8 code refs, 5 artifact refs, not-orphan) and the files inside are invisible.
- `wiki/mpra_substrate_census_2026-10-12.json` carries 4 rungs. Field presence measured: `ceiling_for_a_sequence_model` and `single_replicate_agreement` on 4/4; `size_lever` and `best_measured_arm` on 1/4. It also already carries a `ceiling_rule` string warning against quoting `single_replicate_agreement` as a ceiling — a rule stated in prose that nothing enforces, and that I violated in this session.
- **No `mpra_substrate_census` builder script exists** (`ls scripts/ | grep census` returns seven other censuses, none for this). The JSON is session-authored, so a validator over the existing file is smaller and safer than a generator that would become a second source of truth.
- `wiki/glm_genomewide_oracle_2026-10-07.json` holds the fragment rung's utilization inputs: `split.subsample = 30000`, `load_stats.LB.n_parsed = 321123`, and `honest_limits[0]` already states "The model arms run on a SUBSAMPLE ... the ceiling and the condition effect use the FULL data."
- `scripts/project_status.py` (6,789 bytes) is the carrier surface: `CLAUDE.md` line 6 instructs running it **before making ANY claim about this project's scope**. Per the measured continuity lesson (a mechanism fires when a surface someone already reads carries it), the verdict belongs here rather than in a new standalone CLI.
- `dna_decode/eval/` is the established home for evaluation-discipline modules (`rejection_gates.py`, `regime.py`, `doubt.py`, `prospective_lock.py`, `wildtype_tiering.py`) — all pure, deterministic, verdict-returning, refusing rather than defaulting.
- `~/.claude/skills/soraya/scripts/l3_filehook_probe.py` establishes the safe hook-probe pattern: build an **isolated temp project** with its own `.claude/settings.json` and point the tool at it, never touching the live settings the running agent depends on. `precompact_brief.py` establishes the `--print-settings` inert-shipping pattern.
- `~/.claude/skills/brainstorm/SKILL.md` builds its payload at Step 1 (line 56) with a narrative `## Problem Statement` (line 61), and specifies the same verbatim Problem Statement for the local probes at lines 287 and 311.
- Test infrastructure: `tests/` with pytest configured in `pyproject.toml`; the house pattern for a new discipline module is a sibling `tests/test_<module>.py` with an explicit non-vacuity test.
- **`python3` does not exist on this host** (measured this session — git-bash has `python` only, at `AppData/Local/Programs/Python/Python311`). Any toolkit or hook invocation must use `python`.

### Reusable-Code Survey
- `scripts/data_inventory.py` — reuse `Dataset`, `is_self_reference`, and `_human`; its orphan concept is sound but proven insufficient for this claim class
- `dna_decode/eval/rejection_gates.py` — reuse the ordered-verdict-ladder + "judgment fields refuse rather than default" shape
- `dna_decode/eval/regime.py` — reuse the "refuse to certify when the cited artifact is missing" pattern for the census-field check
- `~/.claude/skills/soraya/scripts/l3_filehook_probe.py` — reuse the isolated-temp-project hook probe verbatim in Step 1
- `~/.claude/skills/soraya/scripts/precompact_brief.py` — reuse the `--print-settings` inert-shipping pattern for Step 7

## Pre-Change Baseline

Captured before any change; each figure was measured this session, not recalled:

- Census field symmetry: `size_lever` present on 1 of 4 rungs; `best_measured_arm` 1 of 4; `ceiling_for_a_sequence_model` 4 of 4; `single_replicate_agreement` 4 of 4.
- Fragment-rung utilization: 30,000 of 321,123 rows = 9.3%; recorded nowhere as a field, only as prose in `honest_limits`.
- `wiki/data_inventory.md` granularity: directory-level; `mpra` = 104.7 MB / 8 code refs / 5 artifact refs / not-orphan. Zero mentions of `frag-rLP5` or `GSE144621` in the rendered output.
- Orphan detection's reach on this class: zero — the fragment files are referenced at `dna_decode/glm/genomewide.py:48-49`.
- Claim-error count this session: 13, of which the external review caught 2 (the union-of-draws arithmetic, and two artifacts I had not found), 5 are execution facts no reviewer can reach, and 1 was missed with the file open.
- Stop-hook payload schema on this build: **never captured**. Live hooks registered in `~/.claude/settings.json`: PreToolUse 1, SessionStart 2, Stop 2, SessionEnd 1, PreCompact 1.
- `/brainstorm` payload today: narrative Problem Statement, no enumerated claims list, no instruction to verify stated background.

## Verification Signal

- **Non-vacuity, the load-bearing one:** the new validator, run against the **unmodified** `wiki/mpra_substrate_census_2026-10-12.json`, must return a non-`ok` verdict naming `size_lever` and the three rungs missing it. A validator that passes today's file has not been shown to detect anything.
- After Step 4, the same validator returns `ok` on the amended census, and all 4 rungs carry a utilization figure.
- A test asserts that a field present on some rungs and absent on others is an **error**, not a pass — and that test fails if the asymmetry rule is removed.
- A test asserts the ceiling guard refuses `single_replicate_agreement` (0.4583) presented as a ceiling where `ceiling_for_a_sequence_model` (0.7928) exists for the same rung.
- `uv run python scripts/project_status.py` prints the symmetry/utilization line, verified by running it and reading the output (not by test alone).
- `/brainstorm`'s payload, on its next real invocation, contains a numbered `## Load-Bearing Claims` block.
- Step 1 returns either a captured Stop payload showing whether assistant text is visible, **or** an explicit `cannot-see-assistant-text` verdict. Either outcome is a result; the design of Step 7 branches on it rather than assuming.

## Implementation Steps

### Step 1: Capture the real Stop-hook payload in an isolated probe
Files: C:/Users/Farshad/.claude/skills/soraya/scripts/stop_payload_probe.py
Depends on: none

**What changes:**
- New probe modelled directly on `l3_filehook_probe.py`: build an isolated temp project directory with its own `.claude/settings.json` registering a `Stop` command-hook that dumps its stdin payload to a file, then run a short throwaway session against that cwd.
- Record the captured payload's top-level keys to a JSON artifact, and answer exactly one question: **does the Stop payload expose the assistant's message text?**
- Emit a three-value verdict: `sees-assistant-text` / `cannot-see-assistant-text` / `indeterminate:<reason>`. Never default to the permissive value.
- **Touches no live settings file.** The isolated-temp-project pattern is the whole safety argument, and it is prior art, not invention.

**Test strategy:**
- Manual: run the probe, inspect the captured payload artifact directly.
- Unit: the verdict function is pure over a payload dict — pin all three verdicts, including that an unrecognised payload shape yields `indeterminate`, never `cannot-see`.

### Step 2: Build the substrate-symmetry validator
Files: dna_decode/eval/substrate_symmetry.py, tests/test_substrate_symmetry.py
Depends on: none

**What changes:**
- New pure module, no I/O, operating on an already-parsed census dict.
- `check_symmetry(census) -> SymmetryReport`: for each field in a declared SYMMETRY_REQUIRED set (`size_lever`, `best_measured_arm`, `ceiling_for_a_sequence_model`), report any field present on ≥1 rung and absent on ≥1 rung. Ordered verdict ladder, first match wins: `indeterminate:no-rungs` → `indeterminate:single-rung` (symmetry is undefined over one rung — refuse rather than pass) → `asymmetric:<field>` → `ok`.
- `utilization(rung) -> UtilizationVerdict`: from `rows_available` and `largest_n_measured`, return `exhausted` (≥95% used), `under-used:<pct>` (<50%), `partial:<pct>`, or `unknown` — and **`unknown` never counts as exhausted** (the fail-safe direction, mirroring the cost-quota rail where an unparseable cell never counts as expensive).
- `ceiling_guard(rung) -> list[str]`: refuse any value presented as a ceiling that equals the rung's `single_replicate_agreement` while a distinct `ceiling_for_a_sequence_model` exists. This enforces, in code, the `ceiling_rule` the census already states in prose.
- The module **reads no files and writes none** — it is the verdict layer; callers supply the dict.

**Test strategy:**
- Unit: every rung of the verdict ladder reachable, each pinned by its own test.
- Unit: non-vacuity — a synthetic census with a field on 1 of 3 rungs returns `asymmetric` naming the two missing rungs; the same census with the field on 3 of 3 returns `ok`.
- Unit: `unknown` utilization does not satisfy `exhausted`; proven by asserting the negative.
- Unit: `ceiling_guard` fires on 0.4583-as-ceiling when 0.7928 exists, and stays silent when only one value is present.
- Unit: a test that fails if `SYMMETRY_REQUIRED` is emptied, so the rule cannot be silently disabled.

### Step 3: Run the validator on the live census and record the refusal
Files: wiki/substrate_symmetry_2026-10-13.json
Depends on: Step 2

**What changes:**
- Run `check_symmetry` + `utilization` + `ceiling_guard` against the **unmodified** live census and write the result as a dated artifact.
- Expected and required: a non-`ok` verdict naming `size_lever` as present on `peak_tiles` only, and `best_measured_arm` likewise. **If it returns `ok`, Step 2 is wrong and Step 4 must not proceed** — a validator that passes the file we know is asymmetric has not been shown to work.
- Record the measured baseline figures in the artifact so the before/after is auditable rather than asserted.

**Test strategy:**
- This step *is* the verification: run it and read the output (verify-in-batch discipline — a passing test suite is necessary, not sufficient).
- Pin the real-file refusal as a test that skips when the census is absent, so it cannot silently vacuum-pass in a fresh clone.

### Step 4: Make the census symmetric by populating the missing per-rung fields
Files: wiki/mpra_substrate_census_2026-10-12.json
Depends on: Step 3

**What changes:**
- Add `size_lever` and `largest_n_measured` to all four rungs. Fragment rungs: `largest_n_measured: 30000`, sourced from `wiki/glm_genomewide_oracle_2026-10-07.json` `split.subsample`, giving `under-used:9.3%`. Peak tiles keep `EXHAUSTED` but **re-worded to name its own rung** so the string cannot read globally.
- Add `best_measured_arm` to the fragment rungs from the oracle artifact (kmer4, 0.1313, 16.6% of that rung's ceiling) and mark `sigma70_designed` explicitly as `not_computable` with the reason already recorded in the census (`not_established`: no ceiling derivable, the two shipped tables being one monotone transform apart).
- Amend **in place** with a dated `corrections_2026_10_13` entry rather than writing a new dated census — the existing file is already cited by `wiki/promcalc_vs_learned_2026-10-12.md`, and a second census would create exactly the two-producers-one-key trap this repo has already paid for.
- Re-run Step 3's validator afterwards; it must return `ok`.

**Test strategy:**
- Re-run the validator and confirm `ok` plus four populated rungs.
- A test asserting the census parses and that every rung carries every `SYMMETRY_REQUIRED` field — this is the regression guard that stops a future rung being added silently asymmetric.

### Step 5: Surface the verdict where a pre-claim check already runs
Files: scripts/project_status.py
Depends on: Step 2

**What changes:**
- Add a substrate section to the existing status output: per rung, rows available, largest N measured, utilization verdict, and that rung's own ceiling.
- Emit the symmetry verdict as a visible warning line when non-`ok`, and the `ceiling_guard` findings when non-empty.
- **Fail soft and loud:** a missing or unparseable census prints `substrate: unavailable (<reason>)` and never blocks the rest of the status output; it also never prints a clean-looking silence.
- No change to any existing figure the script already derives.

**Test strategy:**
- Run `uv run python scripts/project_status.py` and read the output — confirm the rung table appears and that the figures match Step 3's artifact.
- Unit: the renderer is pure over a report dict; pin the `ok`, `asymmetric`, and `unavailable` renderings.
- Unit: a test that the script still exits 0 and still prints its pre-existing scope figures when the census file is deleted.

### Step 6: Require a load-bearing-claims list in the external review payload
Files: C:/Users/Farshad/.claude/skills/brainstorm/SKILL.md
Depends on: none

**What changes:**
- Add a required `## Load-Bearing Claims` section to the Step 1 payload template (and to the probe payload spec at the two places it restates the Problem Statement): a numbered list of the factual claims the conclusion rests on, each with where it came from (file path, artifact, command output, or `recalled — unverified`).
- Add to the Review Instructions: verify each numbered claim against the repo and report any that cannot be confirmed — explicitly licensing the reviewer to challenge the *framing*, not only the proposal.
- Record the rationale inline: measured this session, the review caught the claims I surfaced as uncertain and could not catch the ones I stated as settled background (the 0.46 ceiling never appeared in the payload at all).

**Test strategy:**
- Manual: the next real `/brainstorm` invocation emits the section; confirm by reading the payload that gets sent.
- No automated test — this is a prose contract in a skill file, and claiming otherwise would overstate it.

### Step 7: Ship a closing-claim check, inert, branching on Step 1's verdict
Files: C:/Users/Farshad/.claude/skills/soraya/scripts/closing_claim_check.py
Depends on: Step 1, Step 2

**What changes:**
- If Step 1 returned `sees-assistant-text`: a `Stop`-hook script that scans the turn's assistant text for door-closing language (`exhausted`, `nothing available`, `no longer`, `only`, `cannot`, `already closed`) and, on a hit, prints the substrate utilization table so it lands in the next turn's context.
- If Step 1 returned `cannot-see-assistant-text` or `indeterminate`: fall back to the one trigger that is already known to inject — `SessionStart` — and print the utilization table there. Weaker (it fires before the claim, not after) but it fires without being chosen, which is the entire point.
- **Ships INERT.** Arming requires editing the live `settings.json` that the running agent executes under, which is a recorded standing hazard. The script therefore provides `--print-settings` emitting the exact patch for the operator to apply, and arming is an explicit user decision, not part of this plan.
- Always exit 0, each section independently guarded, output capped well below the 10,000-char harness limit past which output is silently diverted to a file nobody reads.

**Test strategy:**
- Unit: the scanner is pure over a text string — pin hits, misses, and that a quoted/negated occurrence ("not exhausted") is counted as a hit anyway (over-firing is the safe direction here, and that choice is pinned so it is visible rather than accidental).
- Unit: `--print-settings` emits parseable JSON and modifies nothing on disk.
- Manual: run with a sample payload on stdin and confirm the rendered block; confirm the live `settings.json` is byte-unchanged afterwards by hash, not by `git status`.

## Execution Preview

Wave 0 (3 parallel):  Step 1 — Capture the real Stop-hook payload in an isolated probe, Step 2 — Build the substrate-symmetry validator, Step 6 — Require a load-bearing-claims list in the external review payload
Wave 1 (3 parallel):  Step 3 — Run the validator on the live census and record the refusal, Step 5 — Surface the verdict where a pre-claim check already runs, Step 7 — Ship a closing-claim check, inert, branching on Step 1's verdict
Wave 2 (1 parallel):  Step 4 — Make the census symmetric by populating the missing per-rung fields

Critical path: Step 2 → Step 3 → Step 4 (3 waves)
Max parallelism: 3 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — **Step 7's premise is unverified until Step 1 runs.** Whether a `Stop` hook can see assistant text has never been captured on this build. The dependency is explicit and the step branches on the verdict rather than assuming it; if the answer is `cannot-see`, the strongest version of this guardrail is not available and the fallback is genuinely weaker. Do not describe Step 7 as closing the hole.
- Severity: high — **This whole plan is advisory, not enforcement, and must not be described otherwise.** Nothing sits between the model and the message it emits. Steps 2-5 make a wrong substrate claim *detectable and visible*; they cannot make it *impossible*. Step 6 is a prose contract in a skill file. Step 7 injects context; it cannot compel attention.
- Severity: medium — **Step 4 edits an artifact already cited elsewhere.** `wiki/promcalc_vs_learned_2026-10-12.md` references the census. Mitigation: amend in place with a dated corrections entry (the pattern used earlier today), never a parallel census file.
- Severity: medium — **Step 6 edits the live skill file this session is running under.** Lower risk than a live hook edit (prose is re-read per invocation, not executed mid-turn), but the change should be verified by reading the file back and by the next real invocation, not assumed.
- Severity: medium — **`SYMMETRY_REQUIRED` is a hand-declared set**, so a future field can be added to one rung without being listed and the asymmetry check will not see it. This is the same declared-not-detected residual the dispatch-result rail carries, and it is named rather than papered over. A test pins that emptying the set fails.
- Severity: medium — **External tool surfaces.** `python` **verified** present (git-bash, `AppData/Local/Programs/Python/Python311`, 3.11.5) and **`python3` verified ABSENT** — the technical-plan toolkit's documented `python3` invocation must be run as `python` on this host. `codex` **verified** on PATH. `sentrux` **verified absent** → no architecture gate. Hook payload schema for `Stop`: **[unverified]** — Step 1 exists to verify it; no downstream step assumes it.
- Severity: low — **Over-firing by design in Step 7's scanner.** A negated occurrence still counts as a hit. Chosen deliberately (a false alarm costs one glance; a missed closing claim cost three hours today) and pinned by test so it reads as a decision rather than a bug.
- Severity: low — No restructuring was applied. Step ordering fell out of the dependency graph; no steps were split or merged. Intra-wave file-overlap check: Wave 0 touches three disjoint files, Wave 1 three disjoint files, Wave 2 one file. No overlaps, so no dependencies were added for that reason.

## Open Questions

- **Should Step 7 be armed at all, or shipped inert and left?** Arming edits the live `settings.json` the running agent executes under. The plan ships it inert by default and emits the patch; whether to apply it is a user decision, and a reasonable answer is "not yet — see whether Steps 2-5 alone change the error rate."
- **Should the symmetry rule generalize beyond this census?** Other censuses exist (`f3_substrate_census.py`, `ncbi_pd_provenance_census.py`, five more). The same 1-of-N asymmetry may well be present in them, but that is unmeasured and this plan deliberately does not assert it. Scoping to the one measured instance first is the honest order; a sweep is a separate, cheap follow-up.
- **Is the ≥95% / <50% utilization banding right?** Both numbers are asserted, not derived. The only measured anchors are 9.3% (fragments, clearly under-used) and ~94% of parsed rows in use (tiles, clearly exhausted), which the bands separate — but any threshold between those two points would also separate them, so the specific values carry no evidence.

## Verification

1. `uv run pytest tests/test_substrate_symmetry.py -q` — all green, including the non-vacuity and the `SYMMETRY_REQUIRED`-cannot-be-emptied tests.
2. Run the validator on the **unmodified** census; confirm it **refuses** and names `size_lever` plus the three missing rungs. This is the gate on Step 4.
3. Apply Step 4; re-run; confirm `ok` and four populated rungs.
4. `uv run python scripts/project_status.py` — read the output; confirm the rung table, the utilization figures matching Step 3's artifact, and that pre-existing scope figures are unchanged.
5. Delete the census temporarily; confirm `project_status.py` still exits 0, prints `substrate: unavailable`, and does not print a clean-looking silence. Restore.
6. `uv run pytest tests/ -q` — confirm the full suite is still green and the count has only grown by the new tests.
7. Read `~/.claude/skills/brainstorm/SKILL.md` back; confirm the claims-list section is present in all three payload locations.
8. Run Step 1's probe; record the verdict; confirm `~/.claude/settings.json` is byte-identical afterwards by hash.
9. Confirm Step 7 is inert: `--print-settings` emits a patch and nothing on disk changed.

## Save-time amendments

Captured at: 2026-10-13
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads ONLY
`## Implementation Steps` for executable work; amendments here are NOT executable instructions. If an
amendment changes step contracts (file lists, dependencies, structure), re-run `/technical-plan` before
`/execute-plan`.

- item 2 redesigned from per-file orphan detection to census field-symmetry after research showed the fragment files are code-referenced
- Stop-hook feasibility demoted to a probe step

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high,test-strategy-leak -->
