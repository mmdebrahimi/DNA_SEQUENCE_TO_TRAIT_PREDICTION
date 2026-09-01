# Gentamicin v2 — the frozen rule was revised, and re-locked

**User-authorized 2026-08-31.** The AMR decoder's gentamicin rule now counts 16S rRNA
methyltransferases. The reproducibility freeze of 2026-06-13 and the prospective lock of 2026-06-22 are
**retired, not broken** — a new lock supersedes them and the old manifest is preserved.

---

## What changed, exactly

| | v1 (retired) | v2 (deployed) |
|---|---|---|
| rule | `subclass_any={"GENTAMICIN"}` | same, **plus** `symbol_rescue=^(rmt[A-H]\d*\|npmA\d*)$` |
| lock | `prospective_lock_manifest_2026-06-22.json`, cutoff 2026-06-13 | `prospective_lock_manifest_2026-08-31.json`, cutoff **2026-08-31** |
| `amr_rules.py` sha256 | `a983bf28…` | `8007129f…` |

**Not `armA`.** AMRFinder files `armA` under Subclass `GENTAMICIN` (24/24 rows), so the frozen rule
always counted it. The gap was `rmt*`/`npmA` — 134/134 rows filed under the generic `AMINOGLYCOSIDE`.
Including `armA` in the change would overstate what it does.

## It required an ENGINE change, not a config edit

The rule engine had exactly two refinements — `subclass_any` and `gene_prefixes` — and **both narrow,
composing as AND**. The `rmt` fix needs a **union** (`Subclass == GENTAMICIN` **OR** `symbol matches
rmt*`), which was not expressible. That is the same structural limit that forced TMP-SMX into a
scorer-local overlay rather than the frozen rule.

`cipro_determinants_from_main` therefore gained `symbol_rescue` — the **first widening refinement**.

**Its safety is structural, not conventional:** the rescue is evaluated *inside* the drug class gate, so
it can only re-admit a determinant already relevant to that drug. It cannot widen the rule beyond the
drug's own AMRFinder class set, and a test pins exactly that. Default `None` ⇒ byte-identical behaviour
for the other five drugs.

## The justification, re-derived against the deployed rule

`scripts/gentamicin_v2_validate.py` → `wiki/gentamicin_v2_validation_2026-08-31.json`. It scores the
**real `call_resistance`**, not a re-implementation, on the committed leakage-gated census:

| cell | n | v1 acc/sens/spec | v2 acc/sens/spec | Δsens |
|---|---|---|---|---|
| **E. coli × gentamicin** | 131 (65R/66S) | 0.756 / **0.523** / 0.985 | 0.939 / **0.892** / 0.985 | **+0.369** |
| Klebsiella × gentamicin | 43 (19R/24S) | 0.930 / 1.000 / 0.875 | 0.930 / 1.000 / 0.875 | 0.000 |

24 calls changed; 37 rmt carriers, **all R, zero S**. Klebsiella is unchanged because only 2 carriers
are present and both were already called R.

**Why the tuning cohorts disagree, and why that supports the change.** On the three cohorts whose
numbers are published on the report card, v2 buys **+0.013 — one isolate**. Those cohorts contain
almost no `rmt` carriers, so they *cannot* show the effect. That is the source-concentration finding
restated: v1's headline sens 0.893 came from a population lacking the determinant it is blind to. Three
independent measurements agree on what happens where `rmt` actually occurs — disjoint **0.523**,
post-lock prospective **0.429**, label hunt PPV **62/63**.

## The honest limit, unchanged

**Zero S-labelled `rmt` carriers exist in any dataset checked** — not in the 131 disjoint isolates, the
150 local labelled ones, or 63 publicly-labelled carriers. So "specificity unchanged at 0.985" is an
**absence of counter-examples, not a bound**. The over-call risk of the rescue is **untested, not
measured to be zero**. A rule cannot be shown safe on a population containing none of the thing it
newly counts.

## What this COST: the prospective evidence is spent

The 2026-08-24 accrual (63 post-lock isolates: cipro acc 0.967, gentamicin sens 0.429) scored the **v1**
decoder. Those isolates are **pre-lock for v2**, so they cannot be re-scored into evidence — the
prospective clock restarts at 2026-08-31 for **both** cells.

This is enforced in code, not just stated. `prospective_lock_verified` records that the lock held *when
the artifact was written* and is never re-checked, so after a revision every prior score keeps saying
`True` while describing a retired rule — the exact failure the lock exists to prevent, arriving through
the back door. `build_prospective_block` now re-verifies each artifact's stamped hashes against the
**live** surface and returns `superseded_by_surface_change` with the numbers **withheld**.

It applies to **ciprofloxacin too**, whose rule did not change but whose pinned surface did. Behavioural
sameness is an argument; hash-pinning exists so evidence never rests on one.

## Guards: re-pointed, never relaxed

Five tests failed on the change — all of them "the frozen surface is untouched" assertions, and nothing
else. Each was re-pointed at the new invariant:

- `test_tb_leak_guard` re-pinned to `8007129f…`, with the reason recorded — it still guarantees TB work
  never touches the AMR surface; the surface moved by a deliberate AMR revision.
- `test_gentamicin_rmt_candidate`'s "candidate is not deployed" guard is **superseded** by one asserting
  the opposite invariant: *if* the rescue ships, a newer lock manifest must exist, pin the live file,
  carry a later cutoff, and name what it retires. An undocumented surface edit now fails there.
- A new test asserts the **v1 manifest is preserved** — it is the only record of what was locked between
  2026-06-13 and 2026-08-31.
- Two lock-date assertions were **derived** rather than re-hardcoded, so lock v3 will not break them.
- `gentamicin_rmt_candidate.py`'s own control caught the switch live (149/150, refusing to report) and
  was **re-pointed at the deployed rule**, not weakened.

## Reproduce

```
uv run python scripts/gentamicin_v2_validate.py      # the justification, vs the deployed rule
uv run python -m scripts.prospective_lock_validate   # verify_lock against the v2 manifest
uv run python scripts/gentamicin_rmt_candidate.py    # before/after on the tuning cohorts
```
