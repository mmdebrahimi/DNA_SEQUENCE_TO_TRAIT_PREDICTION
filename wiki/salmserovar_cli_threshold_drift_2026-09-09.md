# The shipped `dna-salmserovar` was running a threshold that had been replaced five days earlier

A validated change reached the library and never reached the command. Measured, not asserted:
on the 200-isolate wet-lab cohort, with today's caller, the **shipped CLI default scored 0.6000 while
the validated constant scores 0.8250** — a **−0.225** accuracy gap and **67 abstentions instead of 21**.

## What happened

`dna_decode/salmserovar/runner.py` defines `SEROVAR_COVERAGE_THRESHOLD`. On **2026-09-04** it was
lowered **80.0 → 40.0** against a pre-registered bar, measured at +35 net correct calls, and shipped.

`dna_decode/salmserovar/cli.py` declared `--coverage` with `default=80.0`, written at the cell's first
commit (`e1c65e1`) and **never updated**. Because `main()` passes `args.coverage` to `call_serovar`
explicitly, the CLI's stale copy **overrode** the constant on every invocation.

Nothing caught it because **every validation bypassed the seam**. The threshold sweep, the SeqSero2
comparison and today's O-procedure port all call `call_serovar` directly, which uses the module
default. The function was measured; the command was not.

## The size of it

`scripts/salmserovar_o_fix_rerun.py --coverage 80 --skip-bar`, same 200 isolates, same labels, same
caller — only the threshold differs:

| | shipped CLI (cov 80) | validated constant (cov 40) |
|---|---|---|
| wet-lab accuracy | **0.6000** | **0.8250** |
| hit / miss / no_call | 120 / 13 / **67** | 165 / 14 / 21 |
| O unresolved | 38 | 15 |
| H1 unresolved | 17 | 0 |
| H2 unresolved | 73 | 54 |

Against the current caller, **42 isolates go `hit → no_call`** and one goes `hit → miss` at coverage 80.
The damage is not confined to the O axis it was tuned for: H1 loses 17 calls and H2 loses 19, because
the same cut governs every axis.

`--skip-bar` is used deliberately: a probe at a non-shipped threshold is not a candidate for adoption,
and scoring it against the O-port bar — whose baseline was measured at the shipped threshold — would
compare two different experiments and produce a verdict that means nothing.

## The fix, and its scope

Defaults are now **derived** from the module constants rather than restated:

```python
ap.add_argument("--coverage", type=float, default=SEROVAR_COVERAGE_THRESHOLD, ...)
```

**Nine typing cells were checked; exactly one had drifted.** Every sibling
(`serotype` / `pneumoserotype` / `ktype` / `plasmid` / `resfinder` / `disinfinder` / `pointfinder` /
`mlst`) still matched its constant — because their constants were never revised, so a restated copy
never had the chance to diverge. salmserovar is the only cell whose threshold was later *changed*.

The other six CLIs with both flags were converted to derive anyway. That edit is a **pure refactor**,
and it is proven rather than asserted: `PRE_FIX_DEFAULTS` in the guard records each cell's previous
literal and a test fails if any moved.

## The guard

`tests/test_typing_cli_threshold_drift.py` (22 tests). It **discovers** the typing CLIs from the
package tree rather than taking a hand-written list — a hand-enumerated list beside the data that
defines it is exactly what drifts — and resolves each default **through the real argparse parser**, so
it cannot pass on a source-text match while the shipped command does something else. Two clauses:

1. every CLI default **equals** its module constant;
2. no default is a **hardcoded literal** — equality today is not enough, since a restated number passes
   the first clause and drifts tomorrow. Referencing the constant makes the property true by
   construction.

**Proven non-vacuous:** reintroducing `default=80.0` fails three distinct tests, and reverting restores
all 22.

## Why this class is worth a standing guard

This is the same failure as `tests/test_advertised_commands.py`, one layer down. That guard exists
because a command a user is *told* to run is a promise, and it was resolving five commands that did not
parse. This one exists because a **default a user gets** is also a promise. Both are seams between what
was validated and what ships, and neither report card nor registry test looks at them.

The precondition for the failure is narrow and recognisable: **a constant restated at a call site, plus
a later revision of the constant, plus validation that exercises the function rather than the
command.** All three were present here for five days.

## Honest limits

- The 0.6000 figure is **today's caller** at coverage 80, not a reconstruction of what users got before
  the O port. It sizes the threshold gap on the current code; it is not a historical accuracy series.
- One cohort (200 isolates, per-serovar cap 12), so these are per-isolate accuracies on a deliberately
  diverse mix, not population-weighted rates.
- The guard covers `--identity` / `--coverage` on typing CLIs. Other cells may restate other
  constants (window sizes, k, thresholds elsewhere); that was **not** swept.
- `pointfinder` and `mlst` expose no such flag, so they are covered only by the discovery check
  finding nothing to compare.

## Reproduce

```bash
uv run pytest tests/test_typing_cli_threshold_drift.py -q
uv run python scripts/salmserovar_o_fix_rerun.py --coverage 80 --skip-bar \
  --checkpoint D:/dna_decode_cache/ss2_checkpoint/ofix_cov80.jsonl \
  --out D:/dna_decode_cache/salmserovar_cov80_probe.json
```

Frozen AMR surface byte-unchanged — typing cell.
