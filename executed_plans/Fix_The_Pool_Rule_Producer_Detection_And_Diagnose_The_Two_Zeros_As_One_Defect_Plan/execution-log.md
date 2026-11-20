# Execution Log — Fix_The_Pool_Rule_Producer_Detection_And_Diagnose_The_Two_Zeros_As_One_Defect_Plan
Date: 2026-11-19
Waves: sequential, 5 steps (max parallelism 1 — every step threaded the same two files)
Files changed: dna_decode/fba/design.py, scripts/fba_pool_rule_probe.py, tests/test_fba_design.py, tests/test_fba_pool_rule_probe.py, wiki/fba_pool_rule_diagnosis_2026-11-19.json, wiki/fba_pool_rule_fix_2026-11-19.json, wiki/fba_pool_restriction_generality_2026-11-19.json, wiki/decisions-log.md
Sentrux verdict: n/a — sentrux not installed on this host
Commit: 02dadb3f (branch promoteratlas-exposure-2026-11-01; no merge to main)

## What the plan established

The two zeros (acetate, formate) were ONE defect, not two puzzles. Producer detection was
done by stoichiometric SIGN, which is invalid for a reversible reaction written in the
consuming direction — ACKr / PTAr / PTA2 are absent from the sign-detected producers of
ac_c, and LDH_D from those of lac__D_c.

The first proposed replacement was ALSO wrong: a flux-role test at 90% of WILDTYPE growth
is a pure growth confound, returning one identical load-bearing set across all five
unrelated products, exactly equal to the set of deletions that cannot reach the growth
floor. The corrected method judges the product range at each mutant's OWN maximum,
mirroring evaluate_knockouts.

After the correction all five pools change (13->14, 12->13, 12->14, 10->14, 13->14) and
every design count is unchanged (63, 30, 2, 0, 0). So acetate and formate are a GENUINE
negative rather than an artifact of the rule.

## Deviation from the plan, recorded not forced

Step 4 asked for a test asserting the flux-role arm returns a NON-EMPTY load-bearing set
for acetate. Under the corrected method it returns an EMPTY set for acetate, measured.
That plan line was written under the premise Step 1 refuted, so asserting it would have
pinned the bug instead of the fix. The assertion was adapted and the reason is in the
test file.

## Test outcome

Full suite after the plan: 12 failed, 9012 passed, 191 skipped (17m09s). All 12 failures
sit in four documented absent-data families (cyp4f2/GeT-RM 7, essentiality-ladder 3,
README-quickstart 1, pneumo cps DB 1); none is a code defect and none touches a file this
plan changed. The essentiality family fell from 5 to 3. New tests all green: 26 in
test_fba_design.py, 22 in test_fba_pool_rule_probe.py, 16 in test_decay_curve_confound.py.

## Quoting correction this plan produced

Quote the guaranteed FLUX, never the knockout set. Four designs tie at 13.42495 for
D-lactate, so best_kos is settled by sort order over equal keys and varied across runs.
