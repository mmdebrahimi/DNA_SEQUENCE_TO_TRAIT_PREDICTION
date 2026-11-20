`dna_decode/fba/design.py` (NOT one of the frozen five; editable). Public surface relevant here:
`FLUX_TOL = 1e-6`, `classify_coupling`, `resolve_target(model, target)`,
`product_range_at_growth(model, target_rxn, growth_frac, wt_growth)`, `apply_knockout(model, ident,
level)`, `evaluate_knockouts(model, target_rxn, knockouts, growth_frac, wt_growth, level)`,
`set_anaerobic(model)`, `competition_ranking(model, target_rxn, growth_frac, processes)`,
`find_coupled_designs(...)`, `find_coupled_designs_milp(...)`.

`find_coupled_designs` accepts the candidate pool as `gene_ids: list[str] | None` — verified
signature, so a pool CONSTRUCTOR composes with it without touching it.

`dna_decode/fba/model.py::load_model(None)` yields iML1515 at 2712 reactions / 1516 genes (verified
this session; the long-MILP driver asserts it).

`tests/test_fba_design.py` — 277 lines, 21 tests, all passing. Conventions: `pytest.importorskip("cobra")`
inside each real-model test, `@pytest.mark.slow` on model-loading tests, pools written inline as
`gene_ids=[...]` (the literature-succinate gate and the two lactate tests added this session both do
this — the pool lives in the TEST, not the module, which is the established precedent).

### Reusable-Code Survey

- **`competition_ranking` (design.py:226)** — the existing pool-RANKING helper. A pool-CONSTRUCTION
  helper is its sibling and belongs in the same module; its docstring also records two measured failed
  attempts to fix the pre-ranking, which is why this plan does not attempt a third.
- **`product_range_at_growth` (design.py:146)** — returns the (min, max) product window at a growth
  floor. This is exactly the primitive the flux-role test needs; the new helper must call it rather
  than re-deriving an FBA range.
- **`apply_knockout` (design.py:167)** — applies a reaction- or gene-level knockout non-mutatingly
  inside a `with model:` block. Reuse for the single-deletion probe.
- **`evaluate_knockouts` (design.py:183)** — the module's own verifier; reuse for validation, never a
  hand-rolled coupling check.
- **`set_anaerobic` (design.py:212)** — the condition setter. Reuse; do not set `EX_o2_e` by hand.
- Searched: `dna_decode/fba/*.py`, `graphify-out/GRAPH_REPORT.md` (absent), `src/utils src/lib
  src/common utils lib common` (absent). No existing pool-construction helper found.

Wiki replay — 2 HIGH-salience priors in this module, both load-bearing here:

- **`[plan_file: FBA_Conditional_Essentiality_Remediation_Plan/]`** — "the plan's premise was REFUTED
  by its own Step 8… the pre-committed rule fired `IS_A_SOLVER_ARTIFACT` and the label was WRONG."
  A zero / non-optimal result in this exact module has already once been mislabelled until an
  unplanned probe found the mechanism. Consequence for this plan: Step 3 pre-registers **branches**,
  including one where the fixed rule still finds nothing and that is a GENUINE negative. It must not
  pre-commit to the fix working. The same entry records an artifact that "had been invalid JSON since
  it was written" — every step here re-parses its own JSON.
- **`[plan_file: Oxford_Cohort_External_Revalidation_Plan/technical-plan.md]`** — "read the EXACT
  invocation params from the frozen run's committed artifacts and reuse them verbatim — never
  re-derive or guess them, or the new numbers stop being comparable." Consequence: Step 3 reads
  target list, condition, growth floor, depth and level from the committed
  `wiki/fba_pool_restriction_generality_2026-11-19.json`, not from this plan's prose.
