# dna_decode → GLM: state, plan, and what the AlphaGenome discussion actually taught us

**Date:** 2026-10-08 · **Purpose:** a self-contained brief for brainstorming (incl. with ChatGPT, which has
no context on this repo). Every number below is either **MEASURED** (from a committed artifact, path given)
or flagged **ASSERTED**. Where a figure has changed since it was last written down, the change is named.

---

## 0. What this project is — and the two framings that keep going wrong

`dna_decode` is a **published** (`v0.13.1`, PyPI) multi-kingdom **deterministic** genotype→phenotype decoder
CLI: **50 console entrypoints, 46 CLI traits, 128 registered evidence cells** across 8 tracks (typing 35,
viral 34, AMR 25, pgx 14, TB 12, finder 6, hla 1, mendelian 1). Evidence tiers: 33 `independent_measured`,
23 `near_independent`, 50 `knowledge_baseline`, 11 `no_free_source`, 10 `faithful_to_tool`, 1 `not_censused`.

Two framings to get right before any discussion, because both have cost months:

1. **It is NOT an E. coli AMR project.** AMR is 1 of 8 tracks. A separate AMR "report card" exists with 27
   cells / 10 scored — quoting that as the project's validated surface understates it ~4×.
2. **The GLM is the PRODUCT; the deterministic catalog was a BRIDGE.** The catalog was built so there would
   be *any* number to compare against. **Never benchmark a generative engine against a lookup table** — a
   lookup table cannot generate, so that comparison measures the wrong quantity. Every negative in this
   repo's corpus is about **discriminative zero-shot scoring on natural populations**, which is *not*
   evidence about **generation**.

**North star (user, verbatim):** name a trait → get the genome edits. Four steps: (1) identify the organism
— *shipped* as `dna-identify`; (2) which sections jointly make a phenotype; (3) act in bacteria; (4) climb
to animals.

### The regime map — the single most useful orientation artifact

Encoded as code (`dna_decode/eval/regime.py`), not prose, so a proposal can be screened before it is built:

| population | endpoint | method | verdict |
|---|---|---|---|
| natural | organism | **zero-shot** | **CLOSED NEGATIVE** (0-for-5, de-confounded) |
| natural | organism | supervised | `REQUIRES_DECONFOUNDING` — *not* closed |
| **constructed** | molecular | learned | **WORKS** — ρ 0.35–0.76 (DMS); genome-edit path 0.761 on TEM-1 |
| **constructed** | organism (per-condition) | learned | **WORKS** — FBA/Keio MCC 0.70–0.74 |
| constructed | organism (condition-**switch**) | learned | **OPEN** |

**The discriminating variable is POPULATION DESIGN, not organism complexity.** A yeast segregant cross
decoded 12/12 traits at r 0.46–0.80 — the project's first clean organism-level g→p positive — because a
single-cross panel randomises ancestry by construction. Natural populations cannot.

**8 of 9 regimes have NO held-out-ORGANISM evidence.** That is the project's biggest structural gap.

---

## 1. Where we are RIGHT NOW

The current work is the **GLM portfolio**: 5 families under one umbrella
(`project_state/glm-umbrella-2026-10-07.md`), namespaced `G-*`.

**Terminal condition:** a user names a target molecular phenotype and receives a **ranked set of concrete
genome edits**, each scored by a validated oracle.

| family | what it asks | status as of today |
|---|---|---|
| **G-A** Learned representation | does a learned encoder beat fixed features on *genomic* sequence? | **ANSWERED — negative** (today) |
| **G-B** Self-distillation | perturb-and-relabel to escape the measured-data restriction | **BLOCKED by measurement**, and today's G-A result *keeps* it blocked |
| **G-C** Condition conditioning | one (sequence, medium) model vs two per-medium models | **ANSWERED — negative** (today) |
| **G-D** Generator fine-tuning | LoRA on our own E. coli promoters | **THE UN-PULLED LEVER** — hardware-gated |
| **G-E** Distribution alignment | do the loop's two halves share a distribution? | downstream of G-A + G-D |

**Critical path G-A → G-E.** G-D joins from the side. G-C was the hedge (the only family whose motivating
measurement was already positive).

### The two halves of the loop, and which one is the bottleneck

The generative loop is **propose → score**. Measured 2026-10-07: the scoring half converts ρ 0.59 into
**79% of achievable selection value at k=10** (vs GC's 13%), verdict `ORACLE_READY_AS_SCORING_HEAD`. So
**the scoring half is adequate and the GENERATOR is the bottleneck.** That finding is what makes G-D the
highest-value family, not G-A.

**Measured generator failure (and it is a ZERO-SHOT failure):** GENERator-1.2B-prokaryote at default
sampling scores distinguishability **0.7706** against a 3-mer Markov chain's **0.6679** (lower = better).
A 1.2B model is *more detectable than a trivial Markov chain*. Two levers: sampling (in flight) and
**fine-tuning (never attempted)**.

---

## 2. What we did today — and both answers are negative

Executed a 9-step plan (`plans/GLM_G_A_And_G_C_V2_Build_Plan/`), commits `e2aaf54`…`de23ee6`. Full suite
**6,238 passed / 0 failed**. Both families reached **VALID** negatives — clean nulls and registered
protocols, stamped machine-readably so a *broken* run could not have satisfied either bar.

### G-A: a real conv encoder does not reliably beat ONE GC number

A stop gate ran first and predicted failure: adding POSITION to composition on real genomic promoter tiles
buys **−0.0068** against a +0.02 bar (`wiki/glm_tile_headroom_2026-10-08.json`). **The user decided to build
the encoder anyway.** That was the right call, and the reason is the gate's own honest limits: its proxy was
positional GC in 10 bins under a **linear** ridge, so it could see neither a precise-offset motif nor a
composition×position interaction. An encoder escapes both.

Result (`wiki/glm_encoder_gate_result_2026-10-08.md`), 44,106 peak tiles, leave-peak-out, 10 registered seeds:

| arm | median Spearman | AUROC (active) | paired vs `gc` |
|---|---|---|---|
| `gc` — one GC number (2 params) | **+0.3112** | 0.6672 | — |
| `encoder-multi` (6/20/75, 13,217 params) | +0.3234 | 0.6774 | **6/10 wins, p=0.754, CI95 [−0.0099, +0.0349]** |
| `encoder-w20` | +0.3163 | 0.6739 | 5/10, p=1.000 — a coin flip |
| `encoder-shuffled` (null) | −0.0240 | 0.4892 | 0/10, p=0.002, CI **excludes** zero |

Frozen verdict `WEAK` (margin +0.0122 vs a 0.05 bar). **The paired test says that WEAK is not
distinguishable from zero.** Do **not** quote +0.0122 as the encoder beating the baseline — the first two
seeds looked like a clean win (+0.024, +0.038) and ten seeds plus a paired test are what resolved it.

**Why this is a stronger result than the gate alone:** the encoder is demonstrably *not* broken — it learns
a planted fixed-offset motif to 0.84 against a computed ceiling of 1.0, its null collapses 10/10, and both
framings agree. So the negative is no longer a prediction from a weak proxy; it is a **measurement from the
model class the proxy stood in for**.

**Scope:** closes THIS encoder at THIS size (13k params) on THIS substrate. Capacity has already *hurt*
here — 6-mers (4,096 features) scored **worst** of every arm at 0.2022 — so a larger encoder is a different
experiment. Says nothing about G-B.

### G-C: conditioning on growth medium adds nothing measurable

`wiki/glm_condition_conditioning_result_2026-10-08.md`. 297,868 coordinate-keyed fragments shared between LB
and M9. Primary arm (partial pooling) **+0.2033** vs two independent per-medium models **+0.2032** — a gain
of **+0.0001** against a 0.02 bar. Null clean (−0.0006).

**This is a null about CONDITIONING, not about condition.** The substrate genuinely carries
condition-specific signal and that was measured first: noise-matched cross-condition agreement 0.4062 vs
within-condition 0.4676, disattenuated true cross-condition **0.7887** — so **~38% of the real variance is
condition-specific**. Something is there; a shared *linear* model over 3-mer interactions does not capture it.

---

## 3. The steps we are going to take

Ordered by value, with the blocker type named — because the type determines who can unblock it.

### Step 1 — G-D: LoRA fine-tune the generator (**the highest-value move, and the only un-pulled lever**)

- **Why it is first:** the scoring half is measured adequate (79% of selection value); the generator is
  measured broken *in zero-shot*. Fine-tuning is **the model author's own prescription** for exactly this
  failure mode ("the better approach is to fine-tune on your own data using the open weights rather than
  consuming the predictions as-is").
- **CPU-doable now, no GPU needed:** build the disjoint train/eval window split
  (`dna_decode/glm/finetune_corpus.py`) over ~4,000 E. coli promoter windows.
- **Then:** LoRA at fp32 on GENERator-1.2B-prokaryote; re-run the distinguishability falsifier.
- **Blocker: HARDWARE, and it may already have expired.** The May-2026 pivot away from learned models was
  made on a GTX 860M (4 GiB Maxwell). Since then: a **GTX 1070 (8 GB)** exists on the LAN box, and **free
  Kaggle T4s (16 GB)** have been in use since 2026-07-09. GENERator-1.2B is documented as fitting a T4.
  **ASSERTED, needs checking:** that LoRA-at-fp32 fits 8 GB. **MEASURED:** the T4 route works for this repo.
- **Sequencing constraint:** the sampling sweep must REPORT FIRST. Sampling and fine-tuning are two
  generator levers; run together, neither can be attributed.

### Step 2 — decide whether G-A and G-C are CLOSED or DEFERRED (**user authority**)

Both have valid recorded negatives. The remaining question is purely an acceptance-bar one:
- **G-C:** `NO_GAIN` at n=2 media rules out conditioning *as tested on two growth media with 3-mer
  interactions* — not conditioning in general, and specifically **not** a model whose *representation* is
  modulated by condition.
- **G-A:** rules out *this* encoder at *this* size on *this* substrate.

Executor default is to leave both as recorded negatives and start nothing new.

### Step 3 — G-B stays blocked, and today's result is why

G-B's unblock condition was literally "until a learned representation lands". **It did not land.** So the
self-distillation ban stands — and it is a *measured* ban, not a preference (see §4, lesson 3).

### Step 4 — G-E needs an oracle validated on the matching distribution

Cannot produce a sound number until G-A succeeds or closes, because on genome-like sequence **GC 0.3000
beats every learned feature**. One of G-E's candidates is unblocked today; the arm is otherwise downstream.

### Step 5 — the structural gap worth a separate conversation

**8 of 9 regimes have no held-out-organism evidence**, and the one regime with an organism-level positive
(the yeast segregant cross) points at **constructed population design** as the lever. Nothing in the GLM
portfolio currently attacks that. This is the most promising *unframed* direction.

---

## 4. Lessons from the AlphaGenome / ChatGPT discussion

> **UPDATE 2026-10-09 — there were TWO ChatGPT inputs, and the second one is the plan's actual origin.**
> (a) A *circulated summary* of AlphaGenome, processed 2026-10-07 into the two artifacts named below.
> (b) A **deep architecture discussion** — the transformer's Query/Key/Value mechanics plus 16 proposed
> improvements and a "hierarchical biological world model" framing. **That second one is where G-A and G-C
> came from**, and the ledgers prove it verbatim: G-A's originating goal is *"a LEARNED **multi-resolution**
> encoder plus a trivial head"* (its **priority #1**) and G-C's is *"the bacterial analogue of
> **AlphaGenome's cell-state conditioning**"* (its **priority #3**). It is now archived verbatim at
> `wiki/refs/chatgpt_alphagenome_architecture_discussion_2026-10-09.md` and screened item-by-item against
> our measurements at `wiki/chatgpt_architecture_discussion_lessons_2026-10-09.md`. **Read the screen
> alongside this section** — in particular that both of today's negatives tested the *weakest form* of the
> proposals they came from, and that the 3 proposals it got right were already shipped here.

**Provenance of (a), stated honestly:** it was processed on 2026-10-07 against **primary sources** — the
Pushmeet Kohli interview transcript (28,293 chars) and a DeepMind research-scientist technical talk (52,417
chars), both read in full — into `wiki/alphagenome_integration_analysis_2026-10-07.md` and
`wiki/alphagenome_lessons_actionable_2026-10-07.md`. What follows is from those artifacts, with statuses
**updated for today's results**.

### First: two claims in that ChatGPT summary were WRONG

Worth leading with, because it is the reusable lesson about the source itself:

1. **"published Oct 2024" / a "Nature 2024" citation.** Both wrong. AlphaGenome was announced **25 Jun
   2025** and published in *Nature* **28 Jan 2026** (vol 649, pp 1206–1218). An Oct-2024 Nature citation most
   likely belongs to **Enformer or AlphaMissense**. Verified against primary sources; do not propagate.
2. **The architecture description was the half that does not apply to us** (see the split below) — not
   false, but misleading by omission for a bacterial project.

**Reusable:** a circulated LLM summary is a *lead*, not a source. The cheap fix that worked here was fetching
the actual interview transcripts, which also corrected two of our *own* statements.

### The one-line split

> **The ARCHITECTURE does not transfer. The TRAINING METHOD does.**

The 1 Mb context, U-Net compress/expand, transformer tower and 2D contact branch all exist for **long-range
enhancer–promoter interaction across hundreds of kb** — and a bacterial promoter is **100–300 bp**. Its 11
output heads predict splicing, histone marks and chromatin contacts, **which a bacterium does not have**.
It is human/mouse, and its licence **forbids training on its outputs**. So AlphaGenome cannot be integrated
into the bacterial GLM, and the reason is structural, not effort.

*(Where it could legitimately integrate: the HUMAN clinical track. The AlphaGenome Atlas precomputes ~9
billion SNVs with a calibrated impact score — but it is a **model, not a label**, so using it as ground truth
would be the circular-label failure this repo already gates against. Legitimately a **comparator** or a
**doubt-layer signal**, never a label.)*

### The six lessons, with statuses updated to today

| # | Lesson | Status (2026-10-07) | **Status NOW (2026-10-08)** |
|---|---|---|---|
| 1 | Validate on the **downstream task shape**, not the upstream metric | DONE, measured | unchanged — and it is what identified the generator as the bottleneck |
| 2 | Put a **trivial head on a LEARNED representation** | queued | **RUN. Negative** — the learned half does not beat a single GC number on genomic sequence |
| 3 | **Self-distillation** on perturbed inputs | measured degenerate, blocked until #2 | **STILL BLOCKED, now more firmly** — #2 ran and did not land |
| 4 | **Fine-tune on your own data**; don't consume zero-shot | "the un-pulled lever" | **STILL UN-PULLED** — now the clear #1 move |
| 5 | **DNA-only input** is a deliberate philosophy, not an omission | adopted, zero cost | unchanged, and independently corroborated |
| 6 | **"Prerequisite, not sufficient"** | honesty rail | unchanged — the sharpest line in either source |

**Lesson 2 in detail (the big change).** AlphaGenome puts the capacity in the *representation* and keeps the
heads simple. Our oracle was the **trivial head without the learned representation** — ridge on fixed
features — so its numbers were a *floor*. Today we built the learned half and it did **not** clear that
floor reliably. That reframes the oracle's ρ 0.59 from "a floor awaiting a better representation" to "close
to what this substrate supports".

**Lesson 3 in detail (the most transferable idea, and we still cannot use it).** The method: take the
pre-trained model, **randomly perturb the input, train on the model's own outputs as labels**. Because labels
are *predicted* rather than measured, the held-out-region restriction lifts and you can train across the whole
genome. It solves exactly our problem — measured data exists only for sequences someone assayed, so there is
no direct signal for what an *edit* does. **But it is measured FORBIDDEN here:**

| arm | ρ vs measured | agreement with teacher |
|---|---|---|
| teacher (ridge, real labels) | +0.5963 | — |
| student_linear (same basis, teacher labels) | +0.5963 | **1.0000** |
| student_mlp (richer basis) | 0.5933–0.5992 | 0.999 |
| student_gc (**non-vacuity control**) | +0.1736 | 0.232 |

**Why exactly:** a ridge teacher's predictions lie *in the span of its own features*, so a linear student with
the same basis recovers the same function no matter where in sequence space you perturb. Agreement is
**1.0000 — exactly, not approximately.** There is nothing for the perturbation to regularise until the student
has capacity the teacher lacks. The MLP straddling the teacher (0.5933–0.5992) **measures the information
bound**: distillation cannot create information the teacher lacks. The GC control collapsing 0.596 → 0.174
proves the probe can detect a difference at all.

**And the honesty half that must travel with it:** at DeepMind the distilled model is validated on *real
measured* benchmarks. **Adopting distil-for-smoothness without validate-on-real-data is circular.**

**Lesson 6 in detail.** The author is explicit that AlphaGenome does **not** predict disease or height: it
predicts processes at the level of the cell and nucleus, which is a **prerequisite** for a genetic change
having any effect. *Prerequisite* — necessary, not sufficient. That is exactly our north star's step (2)→(3)
seam: the expression oracle predicts **expression**, which is a prerequisite for a **trait**. ρ 0.59 and 79%
selection headroom are **not trait claims**.

### Smaller facts worth keeping

- DNA methylation is **not** among the 11 outputs.
- AlphaGenome does **not** beat Enformer on everything — on CAGE at 128 bp the author says probably within
  statistical noise, a corrective to the "25 of 26" headline.
- Resolution is **1 bp AND 128 bp** (2048 bp contacts), not uniformly base-pair.
- The "multi-species training improves human performance" claim has **no ablation in the paper** —
  interview-only.
- Their own best variant-effect numbers are in a range where **this problem class is simply hard**; our
  forward cell sits at 0.35–0.76 on protein variants.

---

## 5. Methodological lessons from today — the transferable half

These are worth more than either verdict, and all four were caught *by instruments, not by insight*.

1. **A difference of MEDIANS is not a paired comparison.** G-A's frozen rule compared per-arm medians and
   returned a positive `WEAK`. Paired on identical splits, the same numbers give 6/10 wins, p=0.754, CI
   including zero. **And check the null's own CI excludes zero** (ours: p=0.002) — otherwise "includes zero"
   might just mean the test is underpowered for everything.
2. **A per-feature weight applied BEFORE a standardiser does nothing.** G-C's pooling scale was applied
   before `StandardScaler`, which divides out any constant column factor — measured r **0.9999970** against
   the unscaled model, identical to 7 decimals. All three conditioned arms were secretly **one model**. The
   broken version has **no symptom**: every arm runs and returns a plausible number.
3. **Compute a metric's ceiling under a PERFECT predictor before setting a bar on it.** A planted-motif test
   demanded 0.8 where tie structure capped a *perfect* detector at **0.8661** — asking for 92% of an
   unreachable maximum. The model's "failure" at 0.71 was 82% of what was achievable. **Third mis-specified
   bar in one plan**, same root cause each time: *a threshold frozen before the mechanism analysis that says
   what a correct result can achieve.*
4. **Run the FULL suite before claiming a step is done.** A `U+2016` character in a module docstring would
   have crashed `--help` on a legacy console; it sat green for two commits because only targeted test files
   were being run. A repo-wide guard cannot fire from a targeted run.

Plus one from an independent audit of the code we had just shipped: **five fail-open paths**, including a
guard that attested protocol compliance from a *missing* field, and a vacuous `assert … or True` that could
never fail.

---

## 6. Good questions to brainstorm

Framed so an outside reader can engage without repo context:

1. **Is G-D (fine-tuning) really the right #1?** The argument: scoring is measured adequate, generation is
   measured broken *in zero-shot*, and fine-tuning is the author's own prescription. The counter: both
   learned-representation attempts on this substrate just came back negative — is there a reason to expect
   a fine-tuned *generator* to behave differently from a learned *encoder*?
2. **Does the G-A negative generalise to "learned representations don't help on bacterial promoters", or is
   it specific to 13k parameters on 44k tiles?** Capacity has already hurt here once (6-mers worst of all
   arms). What would a principled next size be, and what would make it a different experiment rather than a
   bigger version of the same one?
3. **Is the substrate the problem rather than the model?** GC alone gets 0.3112 against a derived noise
   ceiling of 0.941 — **33% of achievable**, and the best learned arm only reaches 34%. Is bacterial
   promoter strength from 150 bp simply mostly *not* in the sequence at this resolution — and if so, what
   is the right target instead? (Note the tile set is 14,913 active / 29,193 inactive, so a Spearman over
   all of it may be substantially ranking noise among inactive tiles; that is why the gate reports AUROC
   on the active column as a second framing, and both framings agreed.)
4. **The condition-switch regime is OPEN.** FBA/Keio works per-condition (MCC 0.70–0.74) but the
   conditional-switch question is unanswered, and the measured bottleneck there is that expression
   ∩ condition coverage is 11 of 28 conditions with 60% of samples on glucose. Is there a better substrate?
5. **Constructed population design is the one organism-level lever that worked** (yeast cross, 12/12 traits,
   r 0.46–0.80). Nothing in the GLM portfolio attacks it. What would the bacterial analogue of a segregant
   cross be?
6. **Where does the "prerequisite, not sufficient" seam actually bite us?** We predict expression; the north
   star wants traits. What is the honest shape of the step from one to the other?

---

## 7. Honest limits on everything above

- Both of today's verdicts are **negatives scoped to the method, size and substrate tested** — not to
  learned models in general, and explicitly not to G-B or G-D.
- All thresholds in both gates (0.02, 0.05 margins; 0.05 null bar) are **ASSERTED**, not derived. The noise
  ceilings that contextualise them (0.941 tiles; 0.7925 LB / 0.8040 M9 fragments) **are** derived.
- **Read every tile number against 0.941, never against 1.0.**
- G-C tested **two** growth media. A null at n=2 conditions is weak evidence against conditioning in general.
- The G-A encoder's early stopping fires at 8–16 of 40 epochs on every seed, so the result is not an
  under-training artifact — but the inner validation split, not the epoch cap, is choosing the model.
- Peak-blocked splits hold out whole peaks, but a tile can still share a little sequence across a peak
  boundary.
- **Anything in this document that is a count or a scope claim should be re-derived** with
  `uv run python scripts/project_status.py` rather than trusted from this prose. A written count is not
  evidence — that rule exists here because a written count has been wrong before.
