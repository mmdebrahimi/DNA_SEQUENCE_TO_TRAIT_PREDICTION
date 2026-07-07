# Independent verification — HLA tag-SNP validation (DNA-11, 2026-07-06)

Second-reviewer (mosfaer) READ-ONLY drift-guard of DNA-11's HLA pharmacogenomics claims. Recomputed each
tag's sens/spec/ppv from its committed confusion matrix (`wiki/hla_{b5701,a3101,b5801}_validation_2026-07-06.json`,
`sample_concordance` block) — no re-run of the caller.

## Verdict: PASS — headlines faithful + honest demotion confirmed
| tag (drug) | recomputed | claimed | ships? | verdict |
|---|---|---|---|---|
| **HLA-B*57:01** (abacavir) | sens **0.9792**, spec 0.9924, ppv 0.8545, n=1103 (tp47/fp8/tn1047/fn1) | identical | **YES** | EXACT MATCH — deployable |
| HLA-A*31:01 (carbamazepine) | **sens 0.0** (tp0/fn74 — detects 0 of 74 true carriers) | 0.0 | no | correctly DEMOTED |
| HLA-B*58:01 (allopurinol) | sens 0.609, spec 0.824, **ppv 0.176** (82% of carrier-calls false) | identical | no | correctly DEMOTED |

## Integrity (the load-bearing check)
"Only B*57:01 ships" is **faithful to the data**, not an overclaim: B*57:01 clears a clinical bar
(sens 0.979 / ppv 0.85), while A*31:01 (useless: sens 0) and B*58:01 (unsafe: ppv 0.18 → too many false
carriers) are honestly demoted. This is the self-grading-risk the certification capstone flagged — it held.

## Method / discipline
- Read-only; DNA-11's `main` tree byte-unchanged (verified). Recompute = pure confusion-matrix arithmetic
  over the committed JSON (drift-guard), NOT a re-run of the HLA caller (which needs the 1000G CRAMs).
- What this does NOT re-verify: the tag→truth calling itself (that ran against 1000G HLA truth in-repo);
  only that each reported headline == its own committed confusion matrix, and the deploy/demote gate is
  applied honestly.
- Committed via `soraya_worktree_commit.sh` (v1.18) — the R3-tested non-colliding helper (dogfood).
