**`scripts/glm_tile_hybrid_arm.py`** (exists, LF, 314 added / 17 removed vs HEAD) — the single source
file this plan edits. Verified structure:

- `PREREGISTERED["stacking"]` — frozen bar; keys `question`, `primary_metric`, `comparator`, `margin`
  (0.02), `min_seeds_winning` (8), `mirrors`, `predicted_outcome`, `negative_scope`, `verdicts` (4
  keys: `STACKS`, `WEAK_STACK`, `NO_STACK`, `INDETERMINATE_STACKING_NOT_MEASURED`). Echoed verbatim
  into the artifact, so an added verdict key ships automatically.
- `stacking_verdict_for(p, pre)` at `:298` — three branches. `STACKS` requires
  `median >= margin AND wins >= min_seeds_winning`; `WEAK_STACK` requires only `median > 0`;
  `NO_STACK` is the fallback.
- `stacking_verdict(arms, pre)` at `:314` — builds `per` over `superset_combos()`, returns
  `INDETERMINATE_STACKING_NOT_MEASURED` when `per` is empty, then
  `for want in ("STACKS", "WEAK_STACK"): hits = [...]; if hits: return (want, "; ".join(... for c in sorted(hits)), per)`
  with a `NO_STACK` fallback that joins over `sorted(per)`. **The hit path joins over `hits`, not over
  `per`** — that is defect 2, located exactly.
- `main()` at `:415` — `--in-dir` defaults to `DEFAULT_OUT`, `--seeds` 10, `--out` None.
  `honest_limits` is a 6-item literal list inside the artifact dict at `:537`, so adding a limit
  requires regenerating the artifact. `orthogonality_matrix` at `:523` sorts by `-median` (stable, so
  exact ties fall back to `channel_pairs()` order — checked, deterministic, no flake).
- `render_md()` at `:344` — already emits the stacking section with a per-combo table that includes
  losers, the `predicted_outcome`/`negative_scope` quotes, and a `_skipped_combos` block worded "a
  missing measurement, not a negative one". No change needed for defect 2 on the markdown side.

**`tests/test_glm_tile_hybrid_arm.py`** (exists, LF, 63 tests). The assertions this plan moves:

- `:399` `assert arm.stacking_verdict_for(_p(s["margin"], s["min_seeds_winning"] - 1))[0] == "WEAK_STACK"`
  — median at the margin, 7 of 10 seeds. **This is the one existing assertion the change breaks**, and
  breaking it is the point: that case is a gain that cleared the effect-size bar without reproducing.
- `:398` `_p(margin/2, min_seeds_winning)` → `WEAK_STACK` — survives (8 wins clears the new condition).
- `:400` `_p(-0.01, 2)` → `NO_STACK` — survives.
- `:251` `set(s["verdicts"]) >= {"STACKS", "WEAK_STACK", "NO_STACK"}` — superset check, survives an
  added key.
- `:415` asserts every per-combo verdict is in the frozen `verdicts` map — survives once the key is
  registered there.
- `:807-829` read the SHIPPED artifact and its markdown, so they must run after the re-run, not before.
- `_p(median, wins, n=10)` helper synthesises `ci95` as `[median-0.01, median+0.01]`, so a CI-based
  condition cannot be tested through it without changing the helper. This is why the narrowed
  condition keys on `wins`, not on the CI.

**`tests/test_tile_chain_substrate_wall_amendment5.py`** (exists, **CRLF**, 14 tests) — not edited by
this plan; listed because it shares the baseline measurement and its EOL differs from every other file
here.

**`wiki/glm_tile_stack_arm_2026-10-26.md`** (exists, untracked, 8,007 B). Section headers verified:
`# ...` (line 1, the headline), `## THE PREDICTION WAS FALSIFIED, AND THAT IS THE RESULT` (9),
`## WHY THE VERDICT IS STILL WEAK_STACK...` (48), `## WHAT THE ORIGINAL DECLINE COST` (68),
`## THE 2-CHANNEL QUESTION IS UNMOVED` (80), `## The statistic that did not previously exist` (97),
`## Honest limits` (104). The reduction/superset fact appears at line 114 inside the `negative_scope`
quote only.

**`data/promcalc/tile_predictions/`** — 10 `seed-NN.json.gz` plus `manifest.json`, verified present.
Gitignored, so no worktree can see them; this is why Execution Preview recommends sequential.

**External tool surfaces** — none introduced. The only CLIs are `uv run python`, `uv run pytest` and
`git`, all already in use across this repo. `--in-dir` / `--seeds` / `--out` were read from
`add_argument` calls in the file itself, not inherited from documentation.

### Reusable-Code Survey

- **`PREREGISTERED["stacking"]["min_seeds_winning"]` (8) — REUSED as the narrowed `WEAK_STACK`
  condition.** This is the single most important reuse in the plan: inventing a second seed threshold
  after seeing the result would be a new post-hoc parameter, which is the thing being avoided. The
  frozen number already encodes "how many seeds count as reproducible" for the `STACKS` branch; the
  new branch asks the same question and takes the same answer.
- **The `INDETERMINATE_*` verdict prefix — REUSED, not invented.** The file already carries
  `INDETERMINATE_MISSING_SCORER`, `INDETERMINATE_RECONCILE_FAILED`, `INDETERMINATE_SHORT_ARM`,
  `INDETERMINATE_STACKING_NOT_MEASURED`. The new key joins that family and inherits its meaning: a
  state in which the arm declines to answer rather than answering wrongly.
- **`render_md`'s existing per-combo table — REUSED as the model for the detail string.** It already
  lists every combo with its own verdict; defect 2 is the JSON path failing to do what the markdown
  path already does. No new rendering logic.
- **`sorted(per)` in the `NO_STACK` fallback — REUSED for the hit path.** The correct iteration is
  already written four lines below the defect.
- None of the injected reuse signals applied — searched: `graphify-out/GRAPH_REPORT.md` (absent),
  `src/utils`, `src/lib`, `src/common`, `src/helpers`, `utils`, `lib`, `common`, `helpers` (all
  absent; this repo's shared code lives under `dna_decode/` and `scripts/`).
