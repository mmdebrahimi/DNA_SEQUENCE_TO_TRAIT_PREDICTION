# Two of three shipped HIV deployability numbers were unsourced — measured, and one of them is unmeasurable

**2026-09-30.** `scripts/build_hiv_complement_model.py` ships three supervised blind-spot complement models
and stamps each with a `deployability.leave_study_out_blindspot_auroc`. Its docstring called all three
**"deployment-validated"** and quoted `0.81 / 0.89 / 0.89`.

Only the first was measured. `scripts/hiv_supervised_deployability.py` hardcodes `DRUG = "EFV"` and
`NNRTI_DataSet.Full.txt`, so it covered **RT alone**; the PI and INSTI figures were **literals in a config
dict** with no artifact behind them. This measures every gene — which also closes the largest honest limit
on that morning's GLM headroom result, where all 7 passing drugs shared **one gene**.

## Result

| class | gene | n | catalog called-positive | blind spot n | R | S | **measured** | was shipped | status |
|---|---|---|---|---|---|---|---|---|---|
| NNRTI | RT | 4,222 | 0.554 | 1,883 | 201 | 1,682 | **0.8102** | 0.81 | SCORED (anchor) |
| INSTI | IN | 1,690 | 0.501 | 844 | 59 | 785 | **0.8923** | 0.89 | SCORED |
| PI | PR | 3,532 | **0.742** | 910 | **3** | 907 | **withheld** | 0.89 | **UNSCOREABLE** |

Three different outcomes, and two of them are not what "unsourced number" leads you to expect:

**1. NNRTI was already right.** The `0.81` is the real 0.8102 rounded. It is also the **anchor**: re-measured
here it reproduces the committed `wiki/hiv_supervised_deployability_2026-07-12.json` **exactly** — same
AUROC, same n=1883, same R=201 — and the deployed `call_hiv_observed` agrees with the original majors set on
**4222/4222** isolates, so the comparator swap is validated rather than assumed. A missed anchor makes the
script exit 3 and withhold the other genes.

**2. INSTI's unsourced 0.89 turns out ACCURATE — 0.8923.** This is a **confirmation, not a catch**, and it
should not be dressed up as one. It passes the same pre-registered bar (AUROC ≥ 0.65 AND > mutation-burden
0.2623 AND shuffled null 0.4958 < 0.55) across **69 studies**. It is the **second gene** for the supervised
blind-spot result: integrase, its own catalog, its own dataset.

**3. PI's 0.89 is not merely unsourced — it is not measurable on this substrate.** The deployed
position-based PI catalog calls **74.2%** of isolates resistant, which **empties its own blind spot**: only
**3 of 910** catalog-negative isolates are resistant, far below a 10-per-class floor. The field now carries
`None` plus the reason. This is a property of the **incumbent over-calling**, not a failure of the
complement — and it is exactly what this cell's own `THRESHOLD_NEVER_FIRES` verdict in
`wiki/supervised_complement_nonvacuity_2026-09-30.json` had been recording **while quoting 0.89 beside it**.

## What the number was, and was not, reaching

Stated precisely, because "shipped a wrong number to users" would overstate it:

- It **was** in three git-tracked files under `data/hiv_ref/`, returned by the public
  `hiv_supervised_complement.model_info()`, restated in a builder docstring and a `doubt.py` comment, and
  **published** into the non-vacuity artifact and two wiki memos.
- It was **not** rendered to a CLI user for PI or INSTI. The complement signal is reached only through
  `target_site_doubt`, which registers **mutant-level** catalogs; PI and INSTI are position-based and are
  deliberately excluded. Of the registered drugs only the five NNRTI ones map to a complement model — and
  that one was correct.

## The root cause is closed, not just the values

The values were literals in `CLASSES`, so **any rebuild would have restored them**. `deployability_block()`
now reads the measured artifact, **raises** if it is absent, and **raises if the RT anchor did not bind** —
stamping a number from a harness that cannot reproduce a known result is how an unsourced figure ships. A
test asserts the builder carries no literal and that a rebuild reproduces the committed blocks byte-for-byte.
The historical values survive only as `SUPERSEDED_LITERALS`, explicitly labelled history, used to report what
was superseded and by how much.

The fix **propagated on its own** through the derived non-vacuity artifact (it reads `model_info()`), which is
the check that this is a single-source correction rather than a local patch.

## Scope

- **This changed a disclosure, not a model.** Weights and intercepts are byte-identical in all three files; a
  test pins that `deployability` is the only field that moved, and every CLI call is unchanged.
- **IN-DISTRIBUTION** to the Stanford knowledge base, ~96% subtype B. Leave-one-study-out is the
  out-of-distribution split, **not** an independent cohort.
- A linear model over one-hot substitution tokens is the **weakest** supervised family member, so a pass is
  a **floor** on that family.
- Three genes is still three genes of one virus — this does not speak to the organism-level wall.
- PI's 3 resistant isolates are **cutoff-dependent** (fold ≥ 3 illustrative, no per-drug PI clinical cutoff
  is sourced in-repo), though the blind spot stays tiny under any nearby choice.

## Reproduce

```
uv run python scripts/hiv_deployability_per_gene.py --self-check   # reads the shipped literals, no data
uv run python scripts/hiv_deployability_per_gene.py               # 3 genes; exit 3 if the anchor misses
```

Artifact `wiki/hiv_deployability_per_gene_2026-09-30.json`; 12 tests
`tests/test_hiv_deployability_per_gene.py`. Frozen AMR surface byte-unchanged; prospective lock (v2,
2026-08-31) re-verified `ok=True`.
