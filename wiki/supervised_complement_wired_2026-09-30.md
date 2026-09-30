# The supervised blind-spot complement was shipped and unreachable — now wired into the doubt layer

**2026-09-30.** `dna_decode/data/hiv_supervised_complement.py` shipped 2026-07-12 with measured
deployability and was reachable from **nowhere**. This wires it into the L2 doubt layer so a decoder call
that sits in the catalog's blind spot says so. Three of my own errors are recorded below, because two of
them produced confident wrong conclusions I had to retract mid-run.

## What was already there

| | |
|---|---|
| leave-one-**study**-out blind-spot AUROC | **NNRTI 0.81 · PI 0.89 (LPV) · INSTI 0.89 (RAL)** |
| `hiv_supervised_vs_catalog` | `SUPERVISED_RESCUES_BLINDSPOT` |
| `hiv_supervised_panel` | `GENERAL_RESCUE` — 8 of 11 drugs pass a pre-registered bar |
| `hiv_supervised_deployability` | `DEPLOYABLE_HOLDS_OOD` — train subtype B, holds on 291 non-B |
| reachable from | **nothing except its own builder.** 0 mentions in `cell_registry`, no CLI route, not in the doubt layer |

Third **"shipped but unreachable"** instance in this repo after the HCMV contracts and the `clinvar`/`hla`
routing gap. The capability inventory and the reachability surface are maintained separately, which is what
let me assert the opposite of a shipped, measured result twice in conversation.

## What the wiring does

A **third** signal in `target_site_doubt`, appended and never merged. The two existing deterministic signals
are both structurally silent on the case that matters: `position_novelty` fires only at **catalogued**
positions, and the completeness screen fires only on units whose purity was **measured** (it needs carriers).
A novel substitution at an un-catalogued position with few carriers gets nothing from either.

Three states, never collapsed — **not-measured** (class outside `SUPPORTED_CLASSES`: NRTI and CAI are not
covered) / **not-assessable** (no observed substitutions on this path) / **assessed**.

`target_site_doubt` gained an optional `call` kwarg, and the CLI passes `call.prediction`, so on an already-R
call the signal reports that the blind-spot question does not arise rather than a spurious strong doubt.

## Verified on a real isolate

A Stanford isolate carrying **`K103S` alone**. The deployed catalog calls it **SUSCEPTIBLE** — it carries
103**N**, not 103**S** — while the measured wet-lab fold is **log10 +0.833, about 6.8×: genuinely
resistant.** Before this wiring nothing escalated it.

```
$ dna-amr --drug efavirenz --observed RT:K103S
CALL: S  [screen | 0 determinant(s)]
  DOUBT [strong]: the shipped supervised blind-spot complement scores this genotype 0.717
  (at or above its 0.5 threshold for NNRTI) ... It is a RANKING complement, NOT a rule (the binary
  fold-in was tested and rejected at -0.006 balanced accuracy), and it is IN-DISTRIBUTION to the
  Stanford knowledge base rather than independent validation
```

**The call is unchanged.** L2 qualifies; it never alters.

## Is it non-vacuous? Measured, after three broken attempts

On the real cohort, partitioned by the deployed catalog (**catalog-R 1,057 / catalog-S 1,109** — a
non-vacuous filter, which is exactly what the first two attempts lacked):

**The blind spot = 1,109 catalog-susceptible isolates, 52 truly resistant (base rate 0.047).** At the
shipped 0.5 threshold the complement flags **48** of them:

| | |
|---|---|
| precision | **0.729** |
| recall | **0.673** |
| **enrichment over the 0.047 base rate** | **15.55×** |

So the shipped threshold is **good**. Blind-spot AUROC here is 0.957 but that is **in-sample** — these
isolates are the training distribution. **The honest number is the shipped leave-one-study-out 0.81.**

## The same genotype rendered two different ways depending on the caller

**Caught by an existing test, not by me.** `tests/test_doubt_not_measured_rendering.py::test_the_only_honest_silence_is_assessed_and_quiet`
failed in the full suite: it asserts that `efavirenz` + **`K103N`** renders **no doubt line**, because K103N is
a catalogued major DRM the catalog itself calls **RESISTANT** — "we checked and found nothing" is exactly
what silence truthfully means there, and its docstring says *"if this ever returns a line the renderer has
become noise."*

The first wiring treated an **unsupplied `call` as "assume susceptible" and fired.** The CLI passes
`call.prediction` and so stayed correctly quiet; a library caller passing nothing got a **STRONG blind-spot
doubt on the very same genotype.** One function, two answers, decided by what the caller happened to know.

`target_site_doubt` now **resolves the call from the same deployed catalog the call came from**
(`_deployed_call` → `call_hiv_observed`), so both paths agree. Two properties are deliberate:

- **An unresolvable call is NOT-ASSESSABLE, not a fire.** The fallback must not reinstate the bug in a new
  place: with no call and no way to derive one, the signal has no question to ask. That is the third state
  already in the vocabulary — not a doubt, and not a clean bill.
- **A supplied call outranks a derived one.** A caller asserting `S` on a genotype the catalog calls `R` *is*
  a blind-spot claim, so the signal still fires there.

The motivating case is unaffected: **`K103S` fires STRONG from either path**, and `V179F` still fires from the
*purity* screen rather than the complement (which scores it 0.180) — the "they fail on different cases"
property holds. Both new guards were verified **non-vacuous by reconstructing the pre-fix behaviour**: with
the resolver stubbed out the divergence test fails, and with the assume-susceptible branch restored the
already-R case fires again.

**All four of my tests passed the whole time** — every one of them supplies an explicit `call`, so none
exercised the path the rest of the suite exercises. That is the tell: a new signal's own tests agreeing with
it says little; the existing suite is what holds the contract.

## Three errors of mine, all mid-run

1. **I echoed the call value into the disclosure and the existing guard refused it.** `assert_no_call`
   checks **values, not just keys**, and rejected `call_supplied: "R"` — correctly, because echoing a call
   into its own disclosure makes the disclosure call-shaped. The CLI was broken until fixed. Evidence now
   records *whether* a call exists (`call_was_supplied`, `blind_spot_question_applies`), never its value.
2. **I documented the token normalisation as LOAD-BEARING. It is not.** Measured:
   `blind_spot_risk({"K103N"})` and `{"103N"}` both return **0.9585** against a 0.1416 baseline — the
   complement normalises internally. Kept as defensive redundancy with a **test pinning the two forms agree**,
   so if that ever changes it surfaces as a failure rather than a silently quiet signal.
3. **Two successive broken scans, each producing a confident wrong conclusion.** The first had
   `except Exception: continue`, which swallowed an `AttributeError` on **all 2,168 rows** (the function
   wants a dict, not a list) and returned "0 catalog-susceptible isolates". From that I concluded the 0.5
   threshold never fires and that I had misused the module against its own docstring — **both unfounded**.
   The second passed the dict but with `103N` tokens where the catalog needs `K103N`, so *every* isolate
   landed in the blind-spot set and a 0.996 AUROC was really the full cohort.
   **The fix that mattered was not the format — it was asserting the FILTER is non-vacuous before reading
   any number off it.** Both failures put everything in one bucket, a signature sitting in plain view.

## Honest limits

- **IN-DISTRIBUTION to the Stanford knowledge base, not independent validation.** The complement needed the
  free PhenoSense label to train.
- **A RANKING complement, NOT a rule.** Folding it into the catalog as a binary rule was tested and rejected
  at −0.006 balanced accuracy; the value is the continuous weighting. The doubt tier uses the module's own
  `is_flagged` / `DEFAULT_THRESHOLD` rather than a threshold I invented.
- **A quiet complement is NOT reassurance**, and the signal says so: it scores **0.180** on `179F`, the one
  confirmed gap, which the purity screen flags decisively at 15/15 carriers, p = 8.8e-06. The learned ranker
  and the purity screen fail on **different** cases, so neither supersedes the other.
- **NRTI and CAI are not covered** — `SUPPORTED_CLASSES` is NNRTI / PI / INSTI, reported as not-measured.
- **The PI and INSTI arms were unmeasured when this memo was written; they are now measured, and the two
  differ** — `wiki/supervised_complement_per_class_nonvacuity_2026-09-30.md`. INSTI is non-vacuous
  (**15.99×** enrichment, precision 0.824 / recall 0.737). **PI's threshold NEVER fires**: its blind spot
  holds 2 truly-resistant isolates out of 614 because the position-based PI catalog over-calls so
  thoroughly (catalog-R fraction 0.660) that the susceptible bucket is nearly empty. Do NOT claim the PI
  complement adds call-time value. That memo also **reproduces this one's NNRTI headline** from an
  independently written script.
- **No new registry cell.** The complement makes no independent call and has no CLI route of its own, so it
  belongs in the doubt block of the 26 existing HIV cells rather than as a 116th cell. The evidence-tier
  surface is unchanged.
- The enrichment figures above are **single-drug (EFV), in-sample**. The shipped per-class leave-study-out
  AUROCs remain the deployable claim.
