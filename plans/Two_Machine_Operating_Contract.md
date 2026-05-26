# Two-Machine Operating Contract — dna_decode

> Durable operating contract between Codex on Precision 7780 and Claude on the
> GTX 860M laptop. Supersedes the ad-hoc "work-where-it-fits" pattern that
> caused 4 cross-machine drift incidents in 7 days (2026-05-20 → 2026-05-26).

**Status:** active (2026-05-26).
**Source:** ratified from Codex's 2026-05-26 proposal at
`Downloads/dna_decoder_two_machine_division_of_labor_2026-05-26.md` plus 5
amendments derived from the same-day /brainstorm round.
**Falsification triggers:** §6. First production test = cef audit-aware
closeout integration.

---

## 1. Two lanes (locked)

| Lane | Machine | Function | Examples |
|---|---|---|---|
| **Execution** | Codex / Precision 7780 (RTX 3500 Ada, full Docker, AMRFinder DB) | Heavy compute, real validations, runner mechanics | Model training, AMRFinder runs, full-panel validations, cache populates, genome-input real smokes, release-candidate generation |
| **Discovery + planning + contract-lock conversion** | Claude / GTX 860M (4 GiB Maxwell, no Docker) | Scoping, ideation, feasibility, problem-framing, planning artifact authorship, **converting locked-design parameters from prose to executable regression locks** | Pathotype scoping, EnteroBase feasibility, EP scoping memos, technical plans, `/brainstorm` review, bundle integration with contract-lock conversion |

**Hard rules:**

- Discovery lane does NOT start implementation. Discovery may produce
  *feasibility-grade* scripts (e.g., `scripts/preflight_runnable.py`,
  `scripts/slim_ast_csv_for_bundle.py`) but no production model / runner code.
- Execution lane does NOT silently break locked contracts (enforced by §3
  gate item 4).
- Discovery-lane review of execution bundles = **contract-lock conversion**
  (every locked parameter → test or `KNOWN_DIVERGENCE_TARGETS` marker), NOT
  manual visual inspection.

## 2. Handoff gate (workflow state, not a permanent third lane)

A **handoff gate** is the workflow STATE triggered by any cross-machine
artifact transfer:

- Codex pushes to `origin/main`
- Codex drops a Gmail-able bundle in `Downloads/`
- Claude commits a planning revision that locks design parameters

Either machine may trigger the gate; Claude usually performs the gate checks
because Claude owns the contract-lock conversion responsibility.

### 2.1 The 5 gate checks (all must pass)

1. **Sender push status.** Sender ran `git push origin main` with no
   local-ahead commits. Verified by sender's `git status --short --branch`
   showing zero `ahead/behind`.
2. **Receiver pull.** Receiver ran `git pull --ff-only origin main`.
3. **Sync check.** Receiver ran `scripts/cross_machine_sync_check.py` and the
   verdict is `IN-SYNC` (all `KNOWN_DIVERGENCE_TARGETS` markers PRESENT).
4. **Locked-parameter coverage.** Any newly-locked design parameter in the
   incoming artifact has a regression test path OR `KNOWN_DIVERGENCE_TARGETS`
   marker in the SAME commit that locks it. See §4 (`Contract Locks` section
   spec).
5. **Contract tests run.** Relevant contract tests run on receiver side:
   `uv run pytest tests/test_*_contract.py tests/test_cross_machine_sync_check.py -q`.

A bundle that fails ANY of the 5 checks is **provisional**, not accepted.
Work that depends on a provisional artifact is itself provisional.

### 2.2 Stop condition on drift

If `scripts/cross_machine_sync_check.py` reports drift after a handoff: **hard
pause.** No integration / no new execution / no planning revision depending
on the artifact until classified and resolved:

| Drift class | Resolution |
|---|---|
| Local-ahead unpushed | Push the local commits |
| Local-behind | `git pull --ff-only` |
| Dirty working tree | Commit or stash before integrating |
| Missing marker (LIKELY_STALE) | Sender pushes; receiver pulls; re-verify |
| Downloads-only artifact | Either commit + push the artifact, or mark provisional with a manifest reason (§5) |
| Non-FF divergence | Explicit user review required; no automatic resolution |

## 3. Artifact ownership (source-of-truth model)

**NOT author ownership** (the division-of-labor proposal itself originated on
Precision 7780 — author can be either side). The right rule is who is
source-of-truth for each artifact class.

| Artifact class | Source-of-truth | Cross-lane action allowed |
|---|---|---|
| Planning contracts (`plans/*.md`, locked-design memos, `wiki/decoder_v0_*`, this contract doc) | **Claude** | Codex may PROPOSE updates via `Downloads/` or a `plans/proposals/` draft path, marked `provisional` until Claude integrates |
| Execution artifacts (release packets, validation JSONs, audit outputs, model pickles, runner outputs) | **Codex** | Claude does NOT modify; Claude may quote/cite |
| Repo infrastructure (`scripts/cross_machine_sync_check.py`, `tests/test_*_contract.py`, this contract's regression test) | **Claude** | Codex may extend (add new markers / tests) but should not weaken |
| Cross-lane artifacts (this contract, coordination plans) | **Provisional until integrated by source-of-truth holder.** Defaults: planning contracts → Claude; execution artifacts → Codex |

## 4. `Contract Locks` section spec (mandatory in every locked-design memo)

Every memo at `plans/*.md` or `wiki/*.md` that locks design parameters MUST
include a section titled `## Contract Locks`. The section lists each locked
parameter and its enforcement target. Three valid enforcement targets:

| Target | Format | Example |
|---|---|---|
| Regression test | `path/to/test_*.py::test_function` | `tests/test_drug_mechanism_phenotype_merge_contract.py::test_default_suspend_threshold_matches_cipro_calibration` |
| Sync marker | `KNOWN_DIVERGENCE_TARGETS: (path, marker)` | `(scripts/drug_mechanism_phenotype_merge.py, default=0.40)` |
| Not lock-bearing | `not lock-bearing — rationale` | `prose context; no enforceable behavior` |

**Empty enforcement target = unacceptable.** A parameter described as locked
in prose but with no regression test or marker is a future drift incident
waiting to fire (per the 2026-05-26 cef bundle integration finding).

**Single-commit rule:** the regression test or marker MUST land in the SAME
commit that locks the parameter. Historical splits (e.g., `cd76d6c` integration
+ `716d214` markers) are grandfathered; the future norm is single-commit.
Emergency split is allowed but the intermediate state is explicitly
provisional.

## 5. Handoff manifest (for non-repo artifacts)

`git push origin main` alone does NOT cover all artifacts. Large /
gitignored / Gmail-transferred artifacts (model pickles, AST CSVs, full-panel
validation JSONs) need a **handoff manifest**.

Manifest format — a section at the top of the handoff `Downloads/` doc:

```
## Handoff Manifest
| File | Repo-tracked? | Origin commit | Status |
|---|---|---|---|
| scripts/cef_mic_audit.py | yes | cd76d6c | committed |
| data/processed/models/ceftriaxone_nucleotide_transformer.pkl | no (gitignored — too large) | n/a | gmail-transferred 2026-05-26 |
| wiki/ceftriaxone_mic_audit_2026-05-26.json | yes | (pending Codex run) | not-yet-committed |
```

Status enum: `committed` / `gmail-transferred YYYY-MM-DD` / `not-yet-committed` / `provisional — reason`.

**24-hr alarm:** any contract-bearing artifact in `Downloads/` that is not
`committed` within 24 hr fires an alarm. Better: do not accept immediately
unless the manifest justifies why it can't be committed (size, license,
secrets).

## 6. Falsification triggers (how we know the contract failed)

The contract has FAILED in production if ANY of these occurs:

1. A machine begins work from a cross-machine artifact whose producing commit
   is not on `origin/main`.
2. A locked parameter changes or appears in prose without a same-commit
   regression test or divergence marker.
3. `scripts/cross_machine_sync_check.py` reports drift after handoff but work
   continues anyway.
4. A `Downloads/` artifact is used as source-of-truth for more than 24 hr
   without a repo commit or explicit `provisional` label.
5. A regression test has to be added LATER because a lock was discovered only
   during bundle review (i.e., the `Contract Locks` section was skipped or
   incomplete at memo time).

If any trigger fires, revise the contract before resuming cross-machine work.

## 7. First production test

The first real test of this contract is the cef audit-aware closeout
integration:

- Codex pulls origin (gets this contract + the 2 contract fixes + 5 regression tests + slim AST CSV bundle)
- Codex runs cef MIC audit + mechanism × MIC merge + 4 canonical predicts + release packet update (Artifacts 2-5)
- Codex `git push origin main`
- Claude runs `scripts/cross_machine_sync_check.py` — should report 10/10 markers PRESENT (5 v0/v0.1 + 5 cef audit-aware)
- Claude runs the 5 handoff-gate checks
- If all 5 pass cleanly: contract works in practice
- If any falsification trigger from §6 fires: revise contract before continuing

## 8. Cef closeout temporary exception

Today (2026-05-26) Codex on Precision 7780 is blocked on the BV-BRC AST CSV
(present only on the GTX 860M laptop). The slim AST CSV transfer
(`Downloads/dna_decoder_cef_ast_unblock_2026-05-26.zip`) is a one-off
**handoff exception**, not the new norm. After Codex runs Artifacts 2-5 and
pushes, default lanes resume:

- Codex returns to execution-only
- Claude returns to discovery (pathotype `/idea-anchor` candidate sentence,
  EnteroBase feasibility probe, EP-4 scoping)

## Contract Locks

| Parameter | Enforcement target |
|---|---|
| 2 lanes (Codex execution / Claude discovery+planning+contract-lock conversion) | `tests/test_two_machine_operating_contract.py::test_two_lanes_locked` |
| 5 handoff-gate checks present in §2.1 | `tests/test_two_machine_operating_contract.py::test_handoff_gate_5_checks` |
| `Contract Locks` section spec in §4 | `tests/test_two_machine_operating_contract.py::test_contract_locks_section_spec` |
| 5 falsification triggers in §6 | `tests/test_two_machine_operating_contract.py::test_five_falsification_triggers` |
| Single-commit rule for lock + enforcement target | `tests/test_two_machine_operating_contract.py::test_single_commit_rule_documented` |

Any future edit that weakens or removes these sections fails the regression
test and surfaces in CI before merge.
