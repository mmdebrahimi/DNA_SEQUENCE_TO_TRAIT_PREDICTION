# Semantic contract-number check — right VALUE is not right QUANTITY (2026-09-28)

Closes the standing residual on `scripts/contract_number_audit.py`: *"whole-number matching is not
SEMANTIC — a cited number can match the right VALUE in the wrong QUANTITY (a sensitivity matching a
specificity of equal value)."* Script `scripts/contract_number_semantics.py`, artifact
`wiki/contract_number_semantics_2026-09-28.json`, 19 tests.

## Result

| | |
|---|---|
| numbers extracted | 180 |
| **measurable** | **16 (8.9%)** |
| `field_confirmed` | 15 |
| **`field_mismatch`** | **1** |
| `md_only_or_pooled` | 1 |
| `label_unlabeled` (not measured) | 147 |
| `no_numeric_field` (not measured) | 16 |

**No quantity-mismatch defect was found in any contract.** The one surviving flag is a
citation-coverage lead, not an error (below).

**DO NOT read that as "the contracts are semantically verified."** Only **8.9%** of cited numbers are
measurable by this check at all, so the clean result is **weakly powered**. Reporting the measurable
fraction before the verdict is the point of the artifact.

## The finding that actually matters: the checker had the bug, not the contracts

The first run flagged **12 mismatches. Ten were the checker's own fault** — I boundary-guarded the
*field* side (`path_is_consistent` is token-aware) and left the *prose* side a raw substring:

| prose word | mis-read as | flags |
|---|---|---|
| `J`**`acc`**`ard` | accuracy | 3 |
| `alphamis`**`sens`**`e` | sensitivity | 2 |
| `corr`**`oborating`** | correlation | 1 |
| `dis`**`cov`**`ery` | coverage | 1 |
| `subst`**`rate`** | rate | 1 |
| trailing-label theft (`2.44 (AMRFinder 2.46), Jaccard 0.4126` → 2.46 called a Jaccard) | — | 2 |

**This is the parent audit's own defect, reproduced one layer up, by me, while building the tool meant
to deepen it** — strict on one side, loose on the other, producing agreement (here: disagreement) it
had not earned. Every one of those 10 numbers was in fact sitting at its correct field. Fixing it took
mismatches **12 → 2**. Pinned case-by-case in `tests/test_contract_number_semantics.py`, in both
directions (the guard must reject `J-acc-ard` **and** still recognise `Jaccard` as a jaccard).

## Two false leads, adjudicated rather than averaged

- **`finder:any:forward` `corr 0.49`** — its own artifact, `wiki/esm_at_scale_2026-07-17`, ships as
  **`.md` only**. This check pools a cell's cited JSONs, so a number whose home artifact has no field
  structure gets adjudicated against unrelated siblings. Now classified `md_only_or_pooled`
  (not measurable), not a mismatch.
- **`typing:Salmonella:salmserovar` `coverage 0.705`** — the one remaining flag, and it is a
  **citation-coverage** lead: the coverage value genuinely exists at
  `wiki/salmserovar_threshold_tradeoff_2026-09-04.json:selective_classification.deployed.coverage =
  0.705`, which this cell does **not** cite; the artifacts it *does* cite carry 0.705 as
  `ours.accuracy`. Incidentally this confirms the numeric coincidence is real — that cell's accuracy
  (0.7050) and its coverage (0.705) are equal by chance, exactly the collision the check was built for.

## A guard against the failure this run actually had

The very first run printed `180 numbers | 0 measurable | MISMATCH 0` — arithmetically a clean bill of
health, actually a check running on nothing (`_load_artifact_text` concatenates brace-expanded
siblings to build a search haystack, so every `json.loads` died on `Extra data`). `main()` now
**REFUSES a verdict (exit 3)** below a 5% measurable fraction, because that is a plumbing signature,
not a finding. Pinned by an integration test asserting `n_measurable > 0`.

## Honest limits

1. **91% of numbers are not measurable** — 147 carry no quantity word near them (`52/64 = 0.8125`,
   `EUR 0.2773`), 16 are absent from any cited JSON. Label extraction is a **window heuristic**; a
   quantity named far from its number, or in a table header, is counted as NOT measured, never passed.
2. **Trailing labels are deliberately unused.** Measured across this corpus the trailing form produced
   ≥1 wrong label and 0 correct ones, so `0.833 recall` reads as unlabelled. Loses coverage, removes a
   systematic error.
3. **`field_confirmed` is weaker than it sounds** — a label-consistent field carries that value; it
   does not establish that the contract's *sentence* faithfully summarises the artifact.
4. **Artifact pooling** is unresolved in general: contracts rarely say which cited artifact each number
   came from.
5. Not a gate: exit 0 on a clean run, and it never edits a contract.

## Reproduce

```bash
uv run python scripts/contract_number_semantics.py   # read-only; exit 3 if under-powered
uv run pytest tests/test_contract_number_semantics.py -q
```
