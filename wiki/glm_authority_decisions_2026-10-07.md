# GLM portfolio — the three authority decisions, delegated and recorded

**Date:** 2026-10-07 · **Delegated by the user** ("use your best judgement and give me the option that gives
us the highest output, even if it takes a lot longer... a system that works well is way more important to me
than saving some work"), with the added constraint relaxation that **free compute (Kaggle, additional free
accounts) is available and paid compute is acceptable where there is real added value.**

These three were parked as acceptance-bar / direction calls in `project_state/glm-umbrella-2026-10-07.md`.
They are now decided. Each is recorded with the reasoning, because a delegated decision with no recorded
reasoning is indistinguishable from an arbitrary one.

---

## D1 — A recorded negative IS an acceptable terminal, but ONLY behind validity predicates

**The question in plain terms:** if an experiment concludes "this doesn't work", does that count as finishing
the sub-project, or must a sub-project keep going until it produces a win?

**Decision: a negative terminates the family — but the MVP bar is UPGRADED so that a negative only counts
when the experiment was provably valid.**

The original bar was: module exists · named test exits 0 · verdict artifact exists · ledger row written. The
adversarial review found the hole: **a broken run that writes a null-result artifact satisfies all four
predicates.** So the bar as written made a *fake* negative reachable, which is worse than either extreme.

The fix is a third option, strictly better than both candidates I had framed:

| option | consequence |
|---|---|
| require a WIN | negatives get buried or dead ideas get flogged; a measured closed negative becomes unreachable, which this repo's whole method depends on |
| accept ANY artifact (original) | a broken run passes as science |
| **require VALIDITY, then accept either outcome (adopted)** | honest answers, fast, and a broken run cannot pass |

**Two predicates are therefore ADDED to every GLM family's MVP bar:**
- `null_clean == true` — the pre-registered shuffled-label / permuted-control arm must come out flat. If a
  control that *should* show nothing shows something, the pipeline fabricates signal and nothing is graded.
- `registered_protocol_used == true` — the run used the frozen seed list and frozen training protocol, not a
  post-hoc variant.

**Why this is the highest-output answer rather than the cheapest:** the project's entire track record is
built on recorded negatives (0-for-5 on zero-shot embeddings, the closed Arabidopsis arm, the ESM2 scale
regression). A bar that cannot terminate on a negative would have blocked every one of those, and each saved
months. What was missing was not permissiveness but *proof that the negative was earned*.

---

## D2 — Sequencing REVISED: triage first, then G-A and G-C in parallel, with compute spent on rigour

**The question in plain terms:** in what order, and how many at once?

**Decision: three changes to the original "2 parallel / 1 hardware-blocked / 2 downstream".**

**(a) G-E's reference-vs-reference triage runs FIRST.** The review identified a real contradiction in my own
artifacts — the decomposition says "G-E is last" (decompose:85) while the umbrella says "G-E's candidate 2 is
unblocked today" (umbrella:23). It is not scheduling taste: that measurement tests a premise **both G-A and
G-C rest on** — that designed-grid and genomic sequence are meaningfully different distributions. If they are
barely distinguishable, G-A's entire framing ("fixed features work on the grid and fail on the genome")
needs re-reading before a line of encoder code is written. It costs hours. It goes first.

**(b) Compute buys RIGOUR, not SEARCH.** This is the load-bearing distinction given the relaxed constraint.
The review's sharpest finding was that a 0.05 margin on a single fixed split is a *search surface*. More
compute spent on architecture hunting would make that worse. Spent on rigour it closes it:
- **10 split seeds** (not 1) for every leave-peak-out / position-blocked comparison, verdict on the median,
  seed dispersion reported — this is what actually kills the search-artifact risk
- **the full 297,599 shared fragments** for G-C, not a 30,000 subsample
- **the full 44,106 tiles** for G-A, plus the grid arm, plus the cross-substrate arm
- **the GC baseline re-derived over the same seed list**, so baseline and candidate share the partition set
- free GPU (Kaggle T4, already used in this project) for the encoder arms; CPU remains the fallback and the
  determinism reference

**(c) G-A and G-C then run in parallel, as planned** — they remain logically independent (representation vs
conditioning) and the WIP question resolves as: **two active BUILD families, plus sub-day triage
measurements that touch disjoint files.** The review flagged this as possibly rationalising a contradiction;
it is now an explicit rule rather than an implicit exception, which is the honest version.

---

## D3 — If generated sequence matches NEITHER reference distribution: CONSTRAIN THE GENERATOR

**The question in plain terms:** if the model's invented DNA resembles neither the lab-designed sequences nor
real genomic sequence, do we go find new wet-lab data to score it against, or force the model to produce
more realistic DNA?

**Decision: constrain the generator. And the reason is principled, not economic.**

A generator emitting sequence unlike anything natural is a generator whose proposals are **unlikely to
function in a real cell**. "Neither" is therefore primarily evidence that the *generator* is wrong — not
evidence that the oracle's substrate is wrong. Going to find a new oracle substrate so we can score
unnatural sequence is bending the measurement to fit a broken proposal.

Three supporting reasons:
1. **It is aligned with the actual goal.** The loop exists to propose edits that work in a real genome.
   Sequence off every natural distribution is, by definition, not a useful edit.
2. **It is cheaper and reversible** — conditional prompting is already wired (`generate_conditional`,
   `--context`), so the first attempt costs nothing new.
3. **The alternative walks into the project's known hard wall.** A new oracle substrate means new measured
   labels, and the binding constraint across this entire repo has been LABELS, not models (recorded
   repeatedly, and the reason the `F-E` label-acquisition family exists).

**The one condition under which this reverses:** if the generator is constrained to a natural distribution
and *then* produces only sequence that is natural-looking but functionally inert, the problem has moved and
a new substrate becomes the live question. That is a measurable future state, not a reason to hedge now.

---

## Consequences applied immediately

| from | change |
|---|---|
| D1 | `null_clean` + `registered_protocol_used` added to G-A and G-C MVP bars |
| D2(a) | G-E triage promoted ahead of both build families |
| D2(b) | multi-seed + full-data + re-derived baselines written into both technical plans; free GPU named for the heavy arms |
| D2(c) | WIP rule stated explicitly: two active build families + sub-day disjoint triage |
| D3 | recorded in the G-E ledger so the fork does not get re-litigated when it fires |

Plus the four review findings that are not authority calls and are simply fixes: G-C's model class must carry
`sequence × medium` interactions (an indicator alone can only learn a global offset, so the comparison was
rigged against it); the fragment key must be coordinate-based; the encoder must lazy-import torch; and the
Conv1d determinism claim must be scoped to CPU.
