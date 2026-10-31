# Brainstorm — the learned-teacher self-distillation probe (2026-10-31)

<!-- doc-status: LIVE as of 2026-10-31. Review of the probe design in
     scripts/glm_selfdistill_learned_teacher.py against the frozen bar
     wiki/glm_selfdistill_learned_teacher_acceptance_bar.json and the result
     wiki/glm_selfdistill_learned_teacher_2026-10-31.json.
     OUTCOME: the run's dissolve headline was WITHDRAWN as a result of this review.
     The ridge-anchor replication was NOT withdrawn. -->

> **Note: the external reviewer was unavailable for this session.** `codex` is not installed on this
> host (no binary on either PATH, no `~/.codex`), and the local probe endpoint was unreachable, so
> both auxiliary sections are absent. Analysis is single-model only, which increases the risk of
> anchoring bias — and that risk is elevated here because the design under review is my own. Consider
> re-running this review once an independent reviewer is available.

## Subject

Does the self-distillation ban's **mechanism** survive a learned teacher? The ban
(`wiki/alphagenome_lessons_actionable_2026-10-07.md` §3) reads *"forbidden until a learned
representation LANDS, then only with its real-data validation half,"* and rests on a measurement:
a ridge teacher's predictions lie in the span of its own features, so a linear student on the same
basis reproduced it at Spearman agreement **1.0000** exactly.

The probe was built to convert part of an open **user-authority** question — whether `HYBRID_ADDS`
(a conv encoder that loses alone but lifts a biophysical scalar when rank-averaged, +0.0339 paired,
10/10 seeds) satisfies *"a learned representation LANDS"* — into a measurement, without ratifying it.

## 1. Critical issues

### C1 — The conv arm does not measure what its name claims. The headline is withdrawn. `[grounded]`

Comparing how well each student reproduces the **conv** teacher:

| teacher | linear student (600 feat) | GC student (1 feat) | richer student better? |
|---|---|---|---|
| ridge | **1.0000** (9/9) | 0.7339–0.7502 | yes, decisively |
| conv | 0.3843–0.4475 | **0.4111–0.4567** | **no — 0 of 9** |

**In all 9 conv cells a single GC number reproduces the conv teacher better than 600 positional
features do.** The ridge rows prove the comparison is meaningful: there the richer student wins by
the expected margin.

A strictly more expressive student performing strictly worse at matching the teacher is the
signature of a **sample-limited or overfitting fit**, not of a teacher carrying irreducible
nonlinear structure. The ~0.40 agreement therefore cannot support "the conv teacher is not
recoverable," and `MECHANISM_DISSOLVES_WITH_LEARNED_TEACHER` rested on exactly that.

### C2 — The arms are mis-assigned relative to the teacher. `[grounded]`

In the 2026-10-07 probe the students were *same-basis* / *richer* / *poorer* **relative to their
teacher**; `student_linear` was same-basis by construction, which is what made 1.0000 meaningful.

Against a conv teacher the same student is the **poorer** one and structurally occupies the
control's role. Mirroring `agreement_floor = 0.99` across both arms does not make them
commensurable — it applies a recoverability test to a pairing that cannot satisfy it a priori.
This was flagged as a tautology risk before the run; C1 is the evidence that it bit.

**Honest scope of what ran:** the specific 1.0000-recoverability property does not hold for a conv
teacher *under an under-capacity student*. Whether distillation has anything to gain is **untested**.

## 2. Medium issues

| # | Issue | Label |
|---|---|---|
| M1 | Teacher labels and student features depend on an **unasserted coincidence** — `perturb(...)` is called twice (line 197 → labels at 207; line 217 → features at 218, paired at 226) and matches only because it is deterministic. Divergence would collapse agreement and look **exactly like the reported verdict**. | `[grounded]` |
| M2 | `student_mlp` hit `ConvergenceWarning: Maximum iterations (400) reached` repeatedly; its ridge-teacher agreement swings **0.2375 → 0.9501**. The information-bound arm is unreliable and the artifact does not say so. | `[grounded]` |
| M3 | Perturbed sample is **half** the original scale (8,000 rows vs 16,346). Plausibly feeds C1 directly. | `[grounded]` |
| M4 | The reconcile gate checks **seed 0 only**. All three seeds in fact match published per-seed conv values (0.3235/0.3526/0.3293 vs 0.324/0.353/0.329), so the gate is weaker than the available evidence. | `[grounded]` |
| M5 | `DRIFT_TOL`'s stated rationale cites ~7.15e-07 expected benign drift; the run measured **0.000e+00 in all 9 cells**. Tolerance still safe, documented basis wrong. | `[grounded]` |

## 3. Deferred

- Whether a **rank-average teacher is distillable at all**, given a rank depends on the candidate
  pool rather than on the sequence. Correctly flagged open; untested.
- Uniform random point perturbation does not resemble a generator's candidate distribution —
  already an honest limit in the bar.
- The unreproduced **1.230** drift event and the ≤1e-9 pair. Now *stranger*, since normal drift is
  exactly 0.0. `fit_encoder` was checked directly and does **not** jitter between refits, so the
  published tile-arm numbers are unaffected.

## 4. What's solid

- **The ridge anchor is a real result and is NOT withdrawn.** Agreement 1.0000 at 9/9 cells
  independently replicates the 2026-10-07 degeneracy on a **different substrate** (peak tiles; the
  original GSE108535 element rung is absent from the cache). It does not depend on the conv arm.
- **The shuffled-label control works as derived.** Agreement −0.0163 → 0.0939, rho −0.0054 → 0.0344,
  against a gate of 0.20 set *before* re-running (expected 0; null SE 1/√11197 ≈ 0.0095).
- **Excluding controls from the beats arm was load-bearing, not hygiene.** `student_gc` measured
  d_rho **+0.0309** against the ridge teacher, past the 0.02 margin; without the exclusion the run
  would have reported `DISTILLATION_HELPS_NOW` — the outcome that reopens a banned method —
  manufactured entirely by a control arm.
- Check order (vacuity → anchor → beats → headline), per-seed checkpointing with process isolation,
  and a determinism guard that caught the one genuine anomaly.
- `blocked_by` untouched; the verdict was declared non-ratifying from the outset.

## 5. Open tradeoffs

**Re-run with a conv-capacity student, or withdraw the conv arm entirely?**
Re-running costs a conv student per teacher × rate × seed — roughly 4× the conv fits at ~100 s each
on CPU with ~2.9 GB free. Withdrawing means reporting only the anchor replication: sound, but
answers less of the original question.

**Was the mid-run bar amendment legitimate?**
*For:* the v1 control was unachievable by construction, the three headline thresholds were untouched,
and the record plus a pre-commitment survive in `superseded_v1`.
*Against:* four thresholds in one session were set before asking what the arithmetic or the substrate
permitted. The amendment is defensible; the process that required it is the actual finding.

## 6. Recommended next step

**Withdraw the dissolve headline. Do not cite this run as evidence on the authority question.**
(Done — `headline_withdrawn` block added to the result artifact; the `verdict` field was left
unedited so the record of what the frozen bar mechanically returned survives.)

Then, in order:

1. Add a **conv-capacity student** as the agreement arm; demote the linear student to a control.
2. Assert the perturbed-set correspondence (M1).
3. Raise `COPIES` to 2 (M3).
4. Extend the reconcile gate to all three seeds (M4).
5. Correct M5's rationale; state M2 in the artifact.
6. **Re-freeze with the capacity analysis done first** — the one lesson this session repeated four
   times.

The anchor replication stands on its own and can be reported now.

## The reusable lesson

**A richer student that reproduces a teacher *worse* than a poorer student is a signature that the
agreement arm is sample-limited, not that the teacher is irreducible.** Check the monotonicity of
student capacity against agreement before reading any recoverability number — and set every
threshold only after asking what the arithmetic and the substrate permit. Four thresholds in this
one session were specified before that question was asked: a control that could not collapse
(GC on a rung where GC scores 0.3000), a drift tolerance below the float32 floor (1e-9), a beat
margin a control could clear (+0.0309 > 0.02), and a memory layout that got the job killed.

`Depth: deep · Rounds: 0 (external reviewer not installed) · Stop reason: Claude-only fallback;
grounded claims verified against the artifact and source.`
