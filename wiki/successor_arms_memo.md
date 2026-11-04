# The three successor arms

2 of 3 arms collected. Ceiling verdict: **blocked:untried-arm**.

## Decisions needed -- only you can settle these

**1. Whether the ceiling question is resolved.** The recomputed verdict is `blocked:untried-arm` -- 1 arm(s) declared but never run: finetuned-promoteratlas-encoder-150bp

   `met` would say the declared ARM SET is exhausted: every arm declared against this criterion has been run. It would **not** say the 0.941 replicate-derived ceiling has been reached -- the best arm sits at roughly 54% of it. Those are different claims, and accepting the first as the second is the over-claim this arm structure exists to prevent. Nothing here acts on the verdict.

**2. What to do about the arm(s) still outstanding.** The criterion cannot read `met` while one is unrun, and an arm may not be killed because a sibling measured -- an arm dies only by its own executed evidence.

   - `finetuned-promoteratlas-encoder-150bp` -- job 'finetune_arm' is NOT_STARTED, not COMPLETE -- an artifact is not readable as the arm's result until its run terminated cleanly

## Corrections -- beliefs worth updating

**The leader on the ledger's `best` field moved, and it is not a win.** The as-given arm is the new point-argmax by about +0.0010 over the incumbent, while its own paired interval is [-0.0068, +0.0004] against a margin of 0.02 -- an interval containing zero, several times the margin by which it leads. Its verdict is `ASGIVEN_TIES`. Read the two as joint leaders; `best` is a bare argmax and cannot say this.

**The heads arm is a clean interpretable negative, not a failed run.** Its shuffled-predictor null passed, which is what separates "these heads do not transfer" from "this measurement says nothing". A negative with a passing null is a result.

## What each arm measured

### `asgiven-position-resolved-promoteratlas-150bp`

- verdict: **ASGIVEN_TIES**
- the interval neither clears the margin +0.0200 nor lies below zero, so no gain of at least the margin is demonstrated and no loss is either
- metric: 0.5086 (incumbent 0.5076; 54% of the 0.941 replicate ceiling)
- paired gain interval [-0.0068, +0.0004] against the declared margin 0.02
- comparator: the frozen encoder's 0.5076 on byte-identical splits

Declared limits, verbatim and complete:

- NOT a blind pre-registration -- see NOT_A_BLIND_PRE_REGISTRATION. This arm's evidential weight is lower than the positional arm's whatever it returns.
- The orientation convention remains a 1-BIT INFERENCE from a performance difference between two rungs, not established biology. This arm tests the better-performing side of that inference; it does not establish which orientation `seq` is given in.
- No promoter anchor exists on this substrate, so even an ADDS verdict would say that position bins relative to an ARBITRARY 150 bp window help -- not that a promoter register was recovered.
- The encoder is FROZEN; this bounds nothing about a fine-tuned one.
- 689 duplicate sequences are unchanged and inherited by BOTH sides of the paired comparison.
- The rung ceiling is 0.941 and replicate-derived. No tile-rung number may be read as a fraction of it without saying so.

### `finetuned-promoteratlas-encoder-150bp`

Not readable. job 'finetune_arm' is NOT_STARTED, not COMPLETE -- an artifact is not readable as the arm's result until its run terminated cleanly

### `shipped-expression-heads-zeroshot-150bp`

- verdict: **HEADS_NO_GAIN**
- the whole interval [-0.4131, -0.3861] lies below zero, so a LOSS is demonstrated
- metric: 0.1178 (incumbent 0.5076; 13% of the 0.941 replicate ceiling)
- paired gain interval [-0.4131, -0.3861] against the declared margin 0.02
- comparator: the frozen encoder's 0.5076 on byte-identical splits
- this arm carries INELIGIBLE secondaries (-0.2795, -0.2539, -0.0817, 0.1754, 0.2087). Its own bar declares they cannot be the arm's result; they are not quoted as findings.

Declared limits, verbatim and complete:

- ZERO-SHOT TRANSFER, not head quality: all five heads were trained on other libraries, so a null bounds transfer to the GSE144621 tile rung and says nothing about the heads as predictors of the libraries they were fitted on.
- GEOMETRY MISMATCH: the heads expect 86 positions and the substrate is 150, so every number here is a windowed application and not a native one. 65 windows per tile.
- THE tslpred-kosuri2013 HEAD'S ENDPOINT DOES NOT MATCH THIS RUNG. It predicts TRANSLATION; promoter activity is a transcription quantity. It is included under the standing instruction to leave no stones unturned, and its inclusion WIDENS every head's interval via the alpha/5 correction rather than helping the arm.
- CIRCULARITY IS SCOPED, NOT ABSENT: the urtecho2018 head is circular on the sigma70 rung (GSE108535) and must never be a comparator there. It is not circular on this substrate, and that scoping is recorded so nobody carries it to the wrong rung.
- THE INCUMBENT IS A FITTED RIDGE AND THIS ARM IS NOT FITTED AT ALL. The comparison is therefore asymmetric in tuning budget, and it is asymmetric AGAINST this arm. A loss bounds the pair, not the heads alone.
- 0.941 IS THE REPLICATE CEILING, so the incumbent 0.5076 is 54%% of what is achievable. An arm that ties the incumbent has not approached the ceiling.

---

Rendered by `scripts/successor_arms_memo.py` from the committed artifacts and the live arms ledger. Every number is read from them; none is recomputed here.
