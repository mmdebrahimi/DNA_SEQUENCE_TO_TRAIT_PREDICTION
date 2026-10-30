- **G-A is the owning family, found only after a vocabulary correction.** A grep for the script's name
  `fragment_pyramid` returns **0 files** across `project_state/` and `features/`. The family talks about
  "pyramid arm" and "multi-resolution" instead: `features/glm-learned-representation/feature.md:51,58`
  cite `wiki/glm_pyramid_arm_2026-10-11.{md,json}` and line 63 carries a grounded claim about
  `multiresolution-pyramid-4-decades-UNRUN`. Searching in the vocabulary of the thing being built rather
  than the thing being touched nearly produced a false "no family owns this" finding.
- **The arm's code path.** `scripts/glm_fragment_pyramid_arm.py` imports `build_masked_pyramid` /
  `codes_to_onehot` from `dna_decode/glm/masked_pyramid.py` (line 383) and fits via `fit_masked_pyramid`
  (line 367); the comparator `kmer4_features` (line 229) goes through `_ridge_fit_predict` (line 469), a
  closed-form ridge with no sklearn dependency. Only `dna_decode/glm/masked_pyramid.py`,
  `scripts/glm_fragment_pyramid_arm.py` and `tests/test_masked_pyramid.py` touch the pyramid module.
- **G-A's ledger is live but silent on this rung.** 11 Action-Log rows, most recent dated 2026-10-29.
  Greps for `kaggle`, `INVALID_NULL`, `0.1120`, `0.0517` and `learned_contribution` all return **0 hits**.
- **G-C's ledger already carries the mechanism, not the boundary.** Row 2 records
  `partial_pooling +0.2033 vs two_model +0.2032 = gain +0.0001 against a 0.02 bar`, `verdict: NO_GAIN`,
  plus three defects found en route (the degenerate full-interaction arm, the silent pre-scaler pooling
  no-op, a mis-specified test path) and a flat pooling curve. Row 3 confirms all 6 MVP criteria on disk.
  `eukaryot` appears **0 times** in that ledger.
- **The sigma70 claim already exists, in one place.** `scripts/promcalc_decomposition_arm.py:14` reads
  that the decomposition "gives a learner the mechanism instead of asking it to rediscover sigma70 physics
  from 44k rows." `eukaryot` appears **0 times** across `scripts/promcalc_*.py`. The sibling
  `scripts/promcalc_decompose.py` docstring is long and detailed (isolated interpreter, BLAS pinning) but
  makes no organism claim — so the consequence belongs beside the existing claim, not in a second copy.
- **Line endings.** All three target files are pure CRLF (257 / 199 / 447 CRLF, 0 bare LF). The repo has a
  recorded history of CRLF being misreported, so every edit asserts endings UNCHANGED rather than absolute.
- **Test infrastructure.** `tests/test_promcalc_decomposition_arm.py` and `tests/test_masked_pyramid.py`
  exist and are the relevant suites; ledger files have no test coverage by design.

### Reusable-Code Survey

- **`scripts/promcalc_decomposition_arm.py:14`** — the existing sigma70 statement. Step 3 extends it in
  place rather than writing a second organism-scope note elsewhere.
- **G-A / G-C `## Action Log` tables with their `<!-- project-state:end:action-log -->` markers** — the
  sanctioned append point reused by all three steps; no new section or schema is introduced.
- **`wiki/glm_fragment_pyramid_arm_kaggle_{full,ext}_*.json` + `wiki/fragment_null_ext_summary_2026-10-29.json`**
  — every figure Step 1 records is read back out of these artifacts rather than retyped from prose.
- None — searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/utils`, `src/lib`, `src/common`,
  `src/helpers`, `utils`, `lib`, `common`, `helpers` (none exist in this repo).
