# v2: both fixes worked, and one of them introduced a new failure mode (2026-08-27)

Follow-up to `wiki/staleness_auditor_result_2026-08-27.md`. The v1 run passed at 3/5 TP / 0 FP, and I said
both misses looked fixable and that the fixes had to be **tested, not assumed**. They were. Both landed on
their predicted targets — **and the P3 fix created a false positive I did not predict.**

## The numbers

| run | TP | FN | **FP** | TN | unparseable | passes bar? |
|---|---:|---:|---:|---:|---:|---|
| mechanical screen | 0 | 5 | 5 | 0 | — | no |
| v1 (cap 1200, no facts) | 3 | 2 | **0** | 5 | 1 | yes |
| **v2 (cap 2500, + ARTIFACT FACTS)** | **5** | **0** | **1** | 4 | **0** | yes |

**Recall 3/5 → 5/5, precision 5/5 → 4/5.** Both runs pass the pre-registered bar (`>=3 TP, <=1 FP`), and
v2 sits exactly at the FP ceiling — one more and it would have failed.

## Each fix hit exactly its predicted target — which is the part worth trusting

The two fixes were kept **separable and labelled**, and pushed as a distinct kernel, so movement is
attributable rather than inferred from a blended change:

- **P4 (`unclear` → `stale`).** Predicted cause: token truncation, not judgment. Evidence was
  unambiguous — at cap 1200 exactly one item blew it (5207 chars, `</think>` never closed) while the other
  nine closed at ~1600–3000. Raised to 2500; it now answers, and **unparseable dropped 1 → 0**.
- **P3 (`supported` → `stale`).** Predicted cause: structural, not textual. `browser.py` calls *itself*
  "the deferred v1 graphical browser" and contains **zero** occurrences of "SHIPPED", so a bigger excerpt
  would have supplied more of the same misleading prose. Adding an `ARTIFACT FACTS` line (existence +
  implemented-code count) flipped it.

Two predictions, two confirmations, no reliance on "the score went up so it worked".

## The new false positive is my fix over-generalising, not model error

**N5** — claim: *"Tet + gent dropped from Phase 2 candidate list. 4th-mechanism-class falsifier (gent)
substrate also infeasible"*, artifact: `scripts/bvbrc_strict_mic_4drug_census.py`.

The model reasoned: *"the artifact's code includes Gentamicin breakpoints and strict-MIC classification
logic, directly contradicting the claim … and implying ongoing work"*. Given the rule I wrote, that is
**correct reasoning from a bad premise**. My `ALSO CRITICAL` instruction said, in effect, *implemented code
exists → the claim is stale* — without scoping it.

**The scope it needed:** that inference is valid for a **capability** claim ("X is not built", "the browser
is deferred"), where code existing refutes it. It is invalid for a **finding** claim ("the gent substrate
is infeasible"), where the script is the instrument that *produced* the finding — its existence **supports**
the claim. The census script is exactly that: it ran, and its output is the infeasibility verdict.

So the fix traded a structural blind spot for a category error. Both are real; the second is narrower and
now named.

## What I am NOT doing, and why

**Not patching the rule and re-running to chase 5/5 + 0 FP.** Three reasons, and the first is decisive:

1. **The benchmark is 10 items and I have now run it twice.** Tuning a prompt against it a third time is
   fitting to a 10-item test set. The pre-registered bar exists precisely to stop that; another round of
   tuning would make a PASS mean "I adjusted until it passed", which is the opposite of a falsifier.
2. The capability-vs-finding distinction is a **real, general** improvement — but it should be validated on
   items the prompt has never seen, i.e. on the full corpus with human adjudication, not on these ten.
3. v2 already clears the bar, and the deliverable was never a perfect score. It was: *does a semantic
   auditor beat a mechanical screen that scored 0/5 and 5/5?* Both runs answer yes.

The capability-vs-finding rule is therefore **recorded as the next change to make, and deliberately NOT
made against this benchmark.**

## Standing verdict

**F2 remains PASS**, on a bar committed before the first run. What v2 adds is that the ceiling is higher
than v1 suggested (5/5 recall is reachable) and that reaching it costs precision unless the
capability-vs-finding distinction is drawn.

**Still licenses:** a flagging pass over the full corpus (542 memos / ~525k tokens), every flag adjudicated
against its artifact by hand.
**Still does not license:** editing any doc on the model's say-so. At 1-in-5 false positives on
capability-shaped claims, the adjudication step is not optional — it *is* the product.

## Reproduce

```bash
uv run python scripts/kaggle_staleness_auditor.py --emit-kernel
uv run python scripts/kaggle_push_poll.py push scripts/kaggle/staleness_auditor_kernel.py <slug> --gpu
uv run python scripts/kaggle_push_poll.py pull <slug> <dir>
uv run python scripts/kaggle_staleness_auditor.py --score <dir>/results.json
```
v1 kernel `emanueleebrahimi/dna-staleness-auditor`, v2 `…-v2` (both private, free tier, deletable).
