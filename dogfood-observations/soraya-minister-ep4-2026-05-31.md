# Soraya Minister Dogfood — EP-4 A1 Run — 2026-05-31

## Invocation
- flags: `--advance --max-actions 6 --max-minutes 45 --gate-irreversible`
- scope: **STRICTLY EP-4 (ecoli-pathotype family)**; target **A1 only** (install VF + blastn + VF DB on D:); A2 parked (workhorse/user).
- safety override: treat EVERY dep-install / installer / pkg-manager / archive-extract / DB-download / outside-repo-write / delete / move / credential / remote-mutation / push / publish / send as IRREVERSIBLE → **pause + ask before EACH**; classify each shell segment separately; never infer approval; push only on separate explicit approval.
- run dir: `soraya_runs/2026-05-31-1333-ep4-a1-vf-runtime/`

## EP-4-only guarantee — honest statement
There is **no umbrella ledger**, so `--advance`'s mechanical family-ranking (step 4) cannot *enforce* family isolation. EP-4-only here is guaranteed by **manual scope construction**: the batch contains only A1 pathotype-runtime install steps; no AMR/NT or other-family ledger is read or mutated. Not refusing — the guarantee holds by construction, and the per-action approval gate keeps you in full control. Limitation reported per instruction.

## Initial EP-4 candidate actions (A1 only; A2+ parked)
A1 = make the laptop able to run the pathotype caller locally: VF python pkg + BLAST+ binary + VF DB on D: (C: is full, 28 GB; D: 4.4 TB free). Recon: `virulencefinder` NOT installed; `blastn` ABSENT; VF DB absent.

## Previewed batch — classifier verdict vs widened safety verdict

| # | action | classifier | **widened (binding this run)** |
|---|---|---|---|
| 1 | `uv pip install virulencefinder` (+deps) → venv | irreversible | **PAUSE+ASK** |
| 2 | download NCBI BLAST+ 2.17.0 Win zip → D: + extract | **auto** ⚠ | **PAUSE+ASK** |
| 3 | `git clone` VirulenceFinder DB (CGE bitbucket) → D: | **auto** ⚠ | **PAUSE+ASK** |
| 4 | `makeblastdb` index VF DB on D: | **auto** ⚠ | **PAUSE+ASK** |
| 5 | smoke: download 1 test genome + run VF locally | **auto** ⚠ | **PAUSE+ASK** |
| 6 | write runtime-contract memo (in-cwd) + git commit (NO push) | auto | run (reversible, in-cwd); **push → PAUSE** |

### ⚠ Classifier under-catch finding (validates the override)
`action_gate.classify_action` returned **`auto`** for steps 2–5 — real downloads + outside-cwd (D:) writes + archive extraction + DB indexing. Only step 1 (pip) was flagged `irreversible`. **The widened gate is load-bearing**: without the user override, steps 2–5 would have run un-gated under `--gate-irreversible` (which only pauses on what the classifier labels irreversible). Dogfood verdict: the classifier needs `curl/wget`, `unzip/tar`, `git clone`, `makeblastdb`, `datasets download`, and explicit `D:/`-target (outside-cwd) patterns added to its irreversible set.

## Approval requests + responses
- [PENDING] Action 1 (`uv pip install virulencefinder`) — awaiting user approval.
- (2–5 will each be requested individually before execution.)

## Actions actually executed
- (none yet — paused at action 1 for approval)

## Actions parked / blocked
- **A2** (workhorse commits `pathotype_horesh_*` scripts to git) — PARKED, user/workhorse-only, not attempted from laptop.
- **A5** (Gate B send) — out of A1 scope; user-only.

## Remaining EP-4 work at stop time / terminal reason / stop-short / ascent counterfactual
- (to be filled at run stop)
