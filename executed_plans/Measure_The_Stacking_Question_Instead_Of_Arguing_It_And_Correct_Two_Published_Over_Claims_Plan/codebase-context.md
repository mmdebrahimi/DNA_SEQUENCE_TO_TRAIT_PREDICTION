- **`scripts/glm_tile_hybrid_arm.py`** — `COMBOS` currently holds exactly two entries,
  `hybrid_conv_maxtx: ("conv","max_tx")` and `hybrid_conv_maxtx_gc: ("conv","max_tx","gc")`.
  `runs` is built **generically** as `{k: [] for k in (*COMBOS, "max_tx", "conv",
  "control_shuffled_conv")}`, and the per-combo spearman loop is `for name, names in COMBOS.items()`
  guarded by `if all(n in preds for n in names)`. So a new COMBOS entry **does** get scored
  automatically. Four things around it are **hardcoded** and do not generalise:
  - `_paired_hybrid_vs_comparator` and `_paired_hybrid3_vs_comparator` are two literal lines, the
    second guarded by `if runs["hybrid_conv_maxtx_gc"]`. A third combo gets a median and **no paired
    statistic and no CI** — a figure, not a measurement.
  - `verdict(arms, pre)` reads **only** `arms["_paired_hybrid_vs_comparator"]`. A new combo's numbers
    would be judged by nothing, and the headline would stay the 2-channel's.
  - the `missing` guard is `[s for s in ("conv", "max_tx") if ...]`. A combo needing `decomposition`
    is unguarded: the `if all(n in preds ...)` skip means a table lacking it yields a **shorter**
    `runs[name]` list, and `if vals` then accepts it and records `"n": <short>` with no flag.
  - `orthogonality` is a single `orthogonality(preds["conv"], preds["max_tx"])` call — the matrix
    above exists nowhere in the artifact.
- **`PREREGISTERED`** — `comparator: "max_tx"`, `margin: 0.02`, `min_seeds_winning: 8`,
  `n_seeds_expected: 10`, `mirrors: "scripts/promcalc_decomposition_arm.py PREREGISTERED (margin
  0.02, 8/10 seeds)"`, plus the `precondition` block and a `verdicts` map. **Every paired statistic
  compares to `max_tx`**, so "does a third channel ADD over two channels" is a question the script
  currently cannot express.
- **`tests/test_glm_tile_hybrid_arm.py`** — 36 tests. Two constrain this plan directly:
  `test_frozen_bar_mirrors_the_decomposition_arm` asserts `margin == 0.02`,
  `min_seeds_winning == 8`, `comparator == "max_tx"` (so the stacking bar must be a **new key**, not
  an edit to those three), and `test_the_precondition_is_reported_not_used_as_a_gate` asserts
  `precondition["measured_first"] is True` and `"not_a_gate" in precondition`.
  `test_combine_of_three_modalities_also_holds_the_order` already proves `combine()` takes three
  tables, so no combine work is needed for 3 or 4 channels.
- **The substrate is ready.** All 10 seed tables at `data/promcalc/tile_predictions/` carry exactly
  `['conv','decomposition','gc','max_tx']` — verified per-seed, so the short-list hazard above is
  latent here rather than live, and the guard is defence for later.
- **`render_md`** is annotated in-file as *"the human-facing surface. A caveat carried only in the
  JSON is not a disclosure."* — any new verdict must reach the markdown, not only the JSON.

### Reusable-Code Survey

- **`combine(preds, names)`** — already generic over channel count and already tested for three
  modalities. Reused unchanged for both new combos; nothing to write.
- **`paired(a, b)`** (imported from `promcalc_decomposition_arm`) — reused for every new paired
  statistic, including the new 2-channel-as-comparator one. It is what produces the CI, so reusing it
  keeps every CI in the artifact computed identically.
- **`orthogonality(a, b)`** — already implemented and already tested
  (`test_orthogonality_separates_redundant_from_independent`). The matrix is a loop over it, not new math.
- **`PREREGISTERED["mirrors"]` precedent** — the stacking bar reuses the same 0.02 / 8-of-10 shape and
  cites the same sibling, so the two "does more information help?" results stay directly comparable.
  This is the derivation for the new bar rather than a fresh assertion.
- **`reconcile_problems` + `require_real`** — reused unchanged as the refusal path; the new short-arm
  guard is added alongside them so all refusals live in one place.
- Searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/utils`, `src/lib`, `utils`, `lib` (absent),
  `scripts/glm_tile_*.py`, `dna_decode/glm/*.py`, `tests/test_glm_tile_*.py`.
