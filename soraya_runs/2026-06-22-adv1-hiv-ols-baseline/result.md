# result — 2026-06-22-adv1-hiv-ols-baseline

**Verdict:** advanced — 1 batch executed, cap not hit. Stop = no-further-beneficial-executor-action
without a USER strategic decision (see recommendation.md).

**What shipped:** the OLS underlying-tool baseline for the HIV PI/INSTI/CAI cells (the wrapper-vs-tool
discipline the 5-class expansion shipped without). HIV report card OLS columns populated for all 14 cells.

**Findings (new):** PI OLS recovers +0.06..+0.13 balacc over the high-sens/low-spec catalog (v0.1 headroom,
like NRTI); INSTI catalog competitive (+-0.07, ties/beats OLS on EVG/BIC); CAI catalog BEATS OLS +0.112
(OLS overfits the tiny resistance-enriched n=140).

**Gates:** all `auto`. No money/destructive/irreversible. Frozen AMR surface byte-unchanged. 1568 tests pass.
