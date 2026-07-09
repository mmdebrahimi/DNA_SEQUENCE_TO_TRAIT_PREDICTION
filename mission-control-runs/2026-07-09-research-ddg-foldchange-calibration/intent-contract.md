<!-- intent-contract: 1.0 -->
# Intent Contract — /research (by-hand under Soraya, L1)

- **run-id:** 2026-07-09-research-ddg-foldchange-calibration
- **level:** L1 · **departments:** Research · **skill:** /research · **caller:** user (via /hypothesise #2 follow-up)
- **verbatim-user-input:** "How well do computed protein ΔΔG values (binding free-energy change from FoldX/Rosetta/AlphaFold-based predictors, and protein stability ΔΔG) correlate with and predict measured antimicrobial/antiviral drug-resistance phenotype (fold-change in IC50/MIC), and what calibration methods exist to map computed ΔΔG to measured fold-change? Include accuracy/RMSE of ddG predictors, whether resistance is active-site-binding vs allosteric/regulatory, and any published ΔΔG→MIC or ΔΔG→fold-change regression/calibration."
- **decomposition:** web search + synthesis → intake validation → followup queue update.
- **verification criteria:** ≥5 audit-grade rows OR honest-gap declared; supported memo with ≥1 row + audit floor; unsupported memo exists; followup queue touched.
- **out of scope:** no writes outside research_outputs/ + this run-dir; no rules/wiki/code edits; no fabrication to hit row counts; no auto-promotion.
- **escalation conditions (one FIRED):** search policy-blocked on the antimicrobial-resistance+MIC angle (2 refused queries) → handled by reframing to the thermodynamics/oncology-kinase angle + declaring the antimicrobial-MIC gap honestly (Decision #4). Not a hard halt: ≥5 audit-grade rows obtained via the open angle.
- **caps:** token 15% / wall-clock 30m / tool-calls 100 / uncertainty 5. Actuals: ~8 tool calls, ~10 min, uncertainty=1 (the policy-wall gap).
