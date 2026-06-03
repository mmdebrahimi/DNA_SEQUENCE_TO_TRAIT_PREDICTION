# Dogfood — Soraya `--advance` on EP-4 (2026-06-02)

Run: `soraya_runs/2026-06-02-1500-ep4-interpret-confound/`. First `--advance` flight on a live research task (prior dogfood was `minister` mode, 2026-05-31). Objective: advance the DNA-decoder project while observing Soraya behavior. Observations are objective; interpretation last.

## What worked
- **`action_gate.classify_action` is solid + fast.** Correctly tagged commit=auto, push=irreversible, script write/run=auto, `rm -rf`=irreversible. Compound-command splitting works. The money-only-default vs `--gate-irreversible` distinction behaved exactly per the table.
- **`run_dir` lock + `open_run`** worked first try: `acquire_lock`→True, created the 5-file run set (intent/audit/result/recommendation/approvals). `make_run_id` accepted an injected `now` (no wall-clock dependence — good for reproducibility).
- **VOI ranking led to the right move.** "De-risk the dominant confound before building" outranked "build more model" — and it paid off: the run discovered a project-invalidating structural confound (study==class) that would have silently poisoned any further learned-model work.
- **Emit rule clean.** `/project-state` (user-only here in practice) was emitted with the exact command + parked; the run continued and finished its other branches without halting.

## Gaps / friction (the useful part)
1. **🔴 `rm -rf data/nt_windows` classifies `irreversible` but runs UN-GATED under the money-only default.** A destructive recursive delete executes with no pause unless the operator remembered `--gate-irreversible`. The SKILL.md residual-risk note acknowledges money-only posture, but on a maiden flight this is the sharpest edge: data loss is not money but is irreversible. Consider promoting `rm -rf` / `git reset --hard` / `git clean` to a NARROW always-pause subclass even in money-only mode (a "destructive-local" tier between auto and money).
2. **🟠 Auto-verdict logic ignored feature IDENTITY.** My interpretation script classified the signal `SPARSE → hopeful` from concentration math alone. The correct read was the OPPOSITE (batch artifact) and required inspecting WHAT the top k-mers were (rare CTAG/restriction-site motifs, partial presence). Soraya (the model layer) caught it; the code did not. Lesson for any Soraya-driven analysis: a "shape" verdict from a metric is not a "meaning" verdict — the model must inspect the actual features before trusting an auto-label. This is an OT1-honesty-class reminder: code computes, model must still judge.
3. **🟡 No code-owned audit append.** Per SKILL.md AC9 the run should append to `audit-trail.md` after each step; in practice I wrote it once at the end (model-discipline, not enforced). On a longer run this is where steps get silently skipped — the handoff's OT1 caveat is real and observed.
4. **🟡 Cap/preview ceremony was light.** I previewed the batch in the intent contract but did not re-print caps + gate mode at execution time per step 6. Fine at 4 steps; would matter at the 8–25 range.

## Net
`--advance` did genuine, high-value work (caught a confound that invalidates a whole experimental arm) with the safety-critical pieces (gate classifier, lock) behaving correctly. The two real risks are both at the model-discipline seam, exactly as the handoff's OT1 caveat predicted: (a) destructive-local deletes run un-gated, (b) auto-verdicts can be wrong unless the model overrides on substance. Recommend the `destructive-local` gate subclass as the top follow-up.
