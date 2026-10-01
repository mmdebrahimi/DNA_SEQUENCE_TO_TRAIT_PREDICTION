# Architecture review — specialists + arbiter, and the one piece to change

**2026-10-01.** Review of the user's proposed Phase 5 architecture: one big "Jev-style" local decision model
as arbiter, several small task-specific models as specialists, all specialists report, the arbiter accepts
or rejects; known data checked by comparison; unknown data checked by a "ChatGPT + Claude consensus" for
plausibility. Grounded in the user's own `decision-model-local` session (2026-09-30) plus this repo's
measurements.

---

## 0. Naming, because it changes the tool choice

**"JEV" is Jev from TypeSafe AI** (Sept 2026) — a *decision* model: text + typed question (choice /
boolean / rubric score) → answer with a probability, no prose, 70–500 ms. Closed, API-only, paid.
**It is not Meta's JEPA**, which is an embedding-space predictive (world) model for video/robotics. The
user's own session established this.

Consequence: the proposed arbiter is the **Jev/decision-model** shape, not the JEPA shape — and that is the
*correct* instinct for an accept/reject role. Open local equivalents already exist and were measured:
`tev1:0.8b` **90% accuracy at 0.21 s** (fits 4 GB), `tev1:4b` 91% at 7.4 s, Nimble 9B untested (timed out
on the 860M, left for the 1080 Ti).

---

## 1. Three things to KEEP — each independently supported by measurement

**(a) Several small specialists, not one monolith.** Supported twice: this repo measured ESM2 **regressing**
with scale (650M 0.484 → 3B 0.467 → 15B 0.438), and the user's session measured `tev1:4b` beating
`tev1:0.8b` by 1 point = **2 questions of 180 = noise** at 35× the latency. Bigger is not the lever.

**(b) Separating "predict" from "decide whether to trust the prediction."** This is already a validated
pattern in this codebase — `dna_decode/eval/doubt.py` (the L2 doubt layer) qualifies a call without
competing with it, and its three-state discipline (not-applicable / not-assessable / assessed) exists
because reporting "no doubt" for an unassessed case is a false clean bill of health. The proposed arbiter is
a generalisation of a layer that already works here.

**(c) Confidence-gated acceptance.** The user's session measured it: accepting only ≥90%-confidence answers
lifted accuracy on those answers **90% → 96%**. Calibrated abstention works. This is the mechanism that
should drive accept/reject — not an opinion.

---

## 2. The one piece to CHANGE — and the user's own data is the argument

> **"Escalating the doubtful ones to a bigger version of the same model fixed nothing. The two models got
> 12 of the same 18 questions wrong. Doubtful cases need a DIFFERENT KIND of model."**
> — the user's `decision-model-local` session, 2026-09-30

That is the orthogonality principle, measured independently on the user's own hardware. It is the same
finding this repo got from the modality sweep: ESM2+GEMME lifts (+0.022, 84%) because the two disagree in
different places, while ESM2+AlphaMissense does **not** lift (−0.002) because both are sequence-based and
therefore redundant.

**Applied to the proposal: a "ChatGPT + Claude consensus" is the most correlated pair available.** Two
frontier LLMs trained on overlapping web + literature corpora fail in the *same* places. By the 12-of-18
result, their agreement adds close to nothing — and agreement between correlated judges *feels* like
confirmation, which makes it worse than a single judge.

**Two further reasons it cannot be the validation layer:**

1. **It is most confident where it is least informed.** For an organism with no wet-lab characterisation,
   "is this plausible given what is known in nature" is answered from literature that by definition does
   not cover that organism.
2. **It is an opinion, not a measurement.** This repo has the precedent on record: a search tool volunteered
   "the commonly cited 0.06 mg/L ciprofloxacin ECOFF" and it was **refused**, because a wrong constant
   produces *a fully self-consistent evaluation that is entirely wrong, with nothing downstream able to
   catch it.* The rejection-gate set already encodes this as **G1 (circular / tool-derived label)**.

**The fix is a demotion, not a deletion:** the LLM panel is **discovery-tier** — generate hypotheses,
triage, propose what to test, flag surprises for a human. It must never be **audit-tier** (a validation
signal that can retire a question). That tier distinction already exists in this project.

---

## 3. What to use instead for "data we don't know" — three computable substitutes

**(a) Laws, not opinions.** Checks that are *computable* and can fail: mass balance (the FBA cell already
does this), single-nucleotide reachability of a codon change, reference-translation integrity (the HIV cell
already self-checks this at every DRM position), does the protein fold, is the edit synonymous. A physical
constraint can refute a prediction; a plausibility opinion cannot.

**(b) Held-out-organism calibration — this is the real answer.** The honest way to know how far to trust a
prediction on an organism with no labels: take organisms you *do* have labels for, **pretend you don't**,
and measure error as a function of phylogenetic distance. Then for a genuinely new organism the system can
say *"at this distance from anything I was validated on, my measured error is X."* That converts "unknown"
from a hand-wave into a **number**, and it is exactly families **F1** (transfer benchmark) and **F3**
(de-confounding) of the Phase 5 decomposition.

**(c) Specialist disagreement as a label-free uncertainty signal.** Needs no labels at all: when orthogonal
specialists disagree, that is real information about difficulty. This is the strongest part of the proposed
design and it should be the arbiter's primary input.

---

## 4. Two traps in the "for the data we know, we compare" half

**(a) In-sample flattery.** If the arbiter is trained or tuned on the known data, checking it against the
known data will look excellent and mean nothing. The user's session already hit this exact wall: the public
decision datasets are "probably in the models' training data, so the accuracy is likely flattering", and
many open models score 75–84% on public items but **24–31% on hidden** ones (Nimble fell to 18.7). The
arbiter needs its own held-out set, split **by organism / clade**, never randomly.

**(b) Label noise is not model error.** The session found ~4 of 12 shared "errors" were **wrong labels in
the dataset** (a baseball game labelled "World"). This repo has the same lesson recorded as
*high-sens/low-spec → suspect the label*. An arbiter tuned to agree with noisy labels learns the noise.

---

## 5. Hardware: the 1080 Ti is a real upgrade, with one specific limit

Pascal, compute capability 6.1, **11 GB**. Per this repo's recorded GPU tiering: Pascal "works but no Tensor
Cores / no bitsandbytes" — so **no 4-bit bitsandbytes quantization**, and fp16 gets no tensor-core speedup.

Practical consequence (arithmetic, not a benchmark): a 9B model at fp16 needs ~18 GB and **will not fit**;
the same model as a GGUF quant via Ollama/llama.cpp (which does run on Pascal) is ~5 GB at Q4 or ~9.5 GB at
Q8 and **fits**. A 12B at Q4 is ~7 GB and fits. So **Nimble 9B and the open 12B class are both reachable on
that box via GGUF** — which is the test the user's session explicitly left for the 1080 Ti. Unverified on
the real hardware; it is a sizing estimate, not a measured result.

---

## 6. Net verdict

**The skeleton is right and two details need changing.** Keep: small orthogonal specialists, a separate
trust layer, confidence-gated acceptance. Change: (1) the arbiter must be a **different family** from the
specialists — the user's own 12-of-18 result rules out a bigger sibling; (2) the LLM panel moves from
**validator** to **hypothesis generator**, with held-out-organism calibration taking over the job of
judging the unknown cases.

**The single highest-VOI next measurement** is unchanged from the decomposition: **F1**, the transfer
benchmark with held-out-organism splits. It is what tells the arbiter how much to trust anything on an
organism it has never seen — and without it, every other number in this architecture is uninterpretable.
