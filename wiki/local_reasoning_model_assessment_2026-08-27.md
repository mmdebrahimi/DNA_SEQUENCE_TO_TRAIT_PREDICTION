# Would adding a reasoning model (qwen9b-class) improve the answers? (2026-08-27)

> ## ⚠️ CORRECTED same-day — the first version scoped this to LOCAL hardware and was wrong to
>
> User: *"we can run it on the online free resource(s) that we have and are available."* Correct, and it
> changes the recommendation. This project already uses **free Kaggle T4 GPU** (creds live at
> `~/.kaggle/access_token`; 6+ kernel scripts in `scripts/`) and Databricks. A T4 has **16 GB VRAM**, so
> the 4 GB local-VRAM ceiling I built the analysis on is **irrelevant for batch work** — a 9B (5.0 GB Q4)
> fits trivially, and so does a 14B.
>
> What survives: the interactive-review slot is still filled by Codex, and **F3/F4 stay closed on
> evidence, not compute**. What flips: **F2 (batch work) goes from "no current workload" to the top
> recommendation** — because there IS a 525k-token workload that targets the dominant error class, and it
> is batch-shaped, which is exactly what free GPU is good for. The corrected reason F1 stays blocked is
> **latency, not VRAM**.
>
> The original local-hardware analysis is kept below (§1) because it is still the right answer for
> *interactive* use, and because the reasoning trail matters.

**Verdict: your instinct is right. A second reasoning model already catches real errors — that slot is
filled by Codex. The new, genuinely open opportunity is a BATCH semantic-staleness audit on free GPU,
which targets the error class that actually dominates (10 of 11 this session) and has a ready-made
pass/fail benchmark.**

Every number below is measured today, not recalled.

---

## 0. The compute picture, corrected — two tiers, and they answer different questions

| tier | VRAM | 9B Q4 (5.0 GB) fits? | latency | good for |
|---|---|---|---|---|
| **Local** (GTX 860M) | 4.0 GB | **NO** | interactive | nothing here — see §1 |
| **Kaggle T4** (free, creds live) | **16 GB** | **YES**, easily | **batch**: push → queue → boot → run → poll | **bulk/offline work** |

The T4 removes the capacity objection completely. It does **not** remove the *latency* objection: a Kaggle
kernel is a batch job with no inbound endpoint — you cannot call into a running kernel, so a cold round
trip (queue + boot + ~5 GB model pull + inference) is tens of minutes. Codex returns a `/brainstorm`
critique in ~1–2 minutes. **So free GPU rescues batch work and cannot rescue in-the-loop review.**

That single distinction is what re-sorts the four families below.

## 1. Local hardware: the specific ask does not fit, and neither does the configured default

*(Still the correct answer for INTERACTIVE use, which is why it stands. Superseded only for batch — see §0.)*

| model | Q4 weights | fits 4.0 GB VRAM? | fits 15.9 GB RAM? |
|---|---:|---|---|
| Qwen3-4B (dense) | 2.2 GB | **YES** | YES |
| Qwen3-8B (dense) | 4.5 GB | NO | YES |
| **a 9B dense** | **5.0 GB** | **NO** | YES (CPU only) |
| Qwen3-14B (dense) | 7.8 GB | NO | YES (CPU only) |
| `qwen3.6-35b-a3b` (MoE) | 19.6 GB | NO | **NO** |

Host: **GTX 860M, 4.0 GB VRAM, compute capability 5.0** (Maxwell, 2014), 15.9 GB RAM, C: 8.1 GB free
(97% full), D: 4.2 TB free.

Three consequences:

- **A 9B runs on CPU only here** (5.0 GB > 4.0 GB VRAM). That matters more than usual for a *reasoning*
  model: reasoning models buy quality by emitting long chains of thought, so a slow token rate is taxed
  per-answer far harder than for a non-reasoning model of the same size.
- **~4B dense is the GPU ceiling** (2.2 GB + context inside 4.0 GB). CC 5.0 also predates flash-attention
  kernels and is below the CC ≥ 7.0 the repo already records as the bitsandbytes floor.
- **The existing integration has never been runnable here.** `~/.claude/scripts/lmstudio-auditor.sh`
  defaults to `LMSTUDIO_AUDITOR_MODEL=qwen3.6-35b-a3b`, which needs **19.6 GB — more RAM than this
  machine has**. Its MoE design activates only ~3B params/token (fast *if* resident), but all 35B of
  weights must be held. And no runtime is installed at all: `lms`, `ollama`, `llama-server`, `llama-cli`
  are all absent, with no model directory anywhere. `/d/models` holds only `smolvla_base` (865 MB).

**Disk is NOT the constraint** — I assumed it would be (C: at 97%) and was wrong: D: has 4.2 TB. Any
install must be pointed at D: explicitly, or it will default into the 8.1 GB left on C:.

## 2. The slot already exists, is filled, and worked

`/brainstorm` already routes every review through **Codex** — a frontier model — with the local model as
an *auxiliary second voice*. This session that machinery produced three substantive critiques that changed
the work: the subset p-values are not inferential (nested, sharing variants), "made it causal" overstates
what a definitional permutation shows, and Control B did not leave the gain "intact".

So the answer to "would a second reasoning model help?" is **yes, and you already have one.** The open
question is only whether a *weaker, local* one adds anything on top.

## 3. What actually caught errors this session — the load-bearing evidence

Classifying every substantive error/finding across this session's 14 commits by **how it was caught**:

| caught by | count | examples |
|---|---:|---|
| **Executing / deriving / probing** | **10** | H2's wrong-η² (ran the partial correlation); the colour-cell curation wall (derived from catalogs); `note` vs `notes` field bug (read the dataclass); two regressions in the final suite pass; my own test forcing a false gate on donkey/roe deer; DEG reachable (curl probe); C. auris not label-walled (read artifacts on disk); NT test guarded on the wrong precondition |
| **A reasoning model** | **1** | Codex's `/brainstorm` critique — the wrong-η² attribution, the biggest correction of the session |

**~91% of what went wrong was an unverified assertion, and the fix was running something.** No amount of
additional reasoning capacity addresses that class: the failure is not "thought about it badly", it is
"asserted without checking". A local 9B would sit on the 1-in-11 branch that a frontier model already
covers better.

This is also the honest counterweight to the raw count: that single reasoning-caught error was the most
consequential one. The slot has real value. It is just already occupied by something stronger.

## 4. Decomposition — four candidate families, three recommended against

| # | family | what it would do | verdict (corrected) |
|---|---|---|---|
| **F2** | **Batch semantic-staleness audit** (free GPU) | read each claim + the artifact it cites; flag claims the artifact no longer supports | **TOP CANDIDATE — build the falsifier** (was "no workload"; that was wrong) |
| F1 | Adversarial-review diversity / offline voice | a third opinion alongside Codex + me | **still blocked — but on LATENCY, not VRAM.** Kaggle has no inbound endpoint, so it cannot be in-the-loop. Local is capped at ~4B. Marginal: Codex already fills it and worked |
| F3 | **Decoder-side prediction** (an LLM inside `dna_decode`) | predict phenotype / score variants | **closed negative — do not build.** Evidence-based, so free GPU changes nothing |
| F4 | Curation assistance | draft the 40 unrecorded colour-cell causal variants | **no** — the fabrication hazard by construction; each locus needs OMIA/literature sourcing and verification regardless, so the model removes none of the cost. Compute-independent |

### Why F2 flipped — the workload is real and it targets the dominant error class

I called this "no current workload." Measured, that was wrong:

- **The corpus:** `wiki/` holds **542 memos / ~525,000 tokens**, plus `CLAUDE.md` at ~22k tokens loaded
  into *every* session. Too large to hand-audit; batch-shaped by nature.
- **The target:** §3 shows ~91% of this session's errors were **stale or unverified semantic claims**.
- **The gap is documented, not speculative.** `tests/test_claude_md_citations.py` says in its own
  docstring: *"It checks that cited FILES exist. It cannot check that a cited file still says what
  CLAUDE.md claims it says — the ProSST case would NOT have been caught here… Reading the artifact before
  repeating a claim remains model discipline."* That is precisely the job.
- **Independence matters here.** The stale claims were *mine*. A checker that is not me is worth more than
  me re-reading my own work.
- **Safe failure profile.** Every flag is adjudicated against the artifact, so a false positive costs one
  check; it cannot install a false belief.

### The honest counter-evidence — a mechanical version of this already failed

A proximity-based screen ("a deferral marker within N chars of an existing repo path") was built and run:
**5 hits, 100% false positives**, and was deliberately not shipped, because a guard with that FP rate gets
disabled. So the open question is narrow and real: **does a semantic model do better than the mechanical
rule that failed?** That is exactly what the benchmark in §5 decides — it is not assumed.

### Why F3 is closed, specifically

The repo already ran this experiment with a *better-suited* model. **ESM2-650M scored AUROC 0.454 —
below chance — against the curated catalog's 0.926** on HIV RT drug-resistance mutations, and the
mechanism was measured across three unrelated pathogens: resistance is reached via chemically
**conservative** substitutions at ordinarily-conserved sites, so every exchangeability/likelihood scorer
calls them benign. A 1992 BLOSUM62 matrix reproduces the same blindness. Scale made it worse (3B and 15B
regress below 650M).

A general reasoning LLM has *weaker* sequence priors than ESM2, not stronger — it has no evolutionary or
structural prior over protein sequences at all. The regime map is explicit: curated-catalog → deterministic
wins; organism-polygenic → neither; molecular-property → learned wins only when fitness-aligned.
**Drug resistance inverts the signal.** Adding an LLM to the decoder is the closed arm, re-entered.

## 5. The plan -- F2 on free GPU, with a pre-registered benchmark that already exists

The benchmark is the good news: **both classes of ground truth are already in the repo**, so this is a
falsifiable experiment, not a hopeful build.

**Ground-truth POSITIVES** -- known stale claims, each pinned as a regression test in
`tests/test_claude_md_citations.py`, plus one found today:

| # | the stale claim | why mechanically invisible |
|---|---|---|
| 1 | ProSST "the real forward pass is deferred to a Kaggle run" | it had RUN locally on CPU the same day; the cited file existed the whole time |
| 2 | TB "PENDING DATA RUNS (BLOCKED-gated by design)" | both blockers had been resolved -- by its own sub-bullets |
| 3 | genome-map "a visual browser deferred" | the browser had SHIPPED; the same file said so two bullets later |
| 4 | BV-BRC "strict-MIC 3-drug census ... deferred" | it ran as a FOUR-drug census; stale on count *and* status |
| 5 | C. auris `no_free_source` (found today) | a POWERED result exists -- 12 isolates, sens 1.00 |

**Ground-truth NEGATIVES** -- the 5 hits the mechanical proximity screen produced, all confirmed false
positives (marker and path in the same prose region but referring to different things; 3 of the 5 were
inside *correction* text explaining a fix).

| step | action | gate |
|---|---|---|
| 1 | Build the claim-to-artifact extraction pass locally: for each `(claim, cited artifact)` pair in `CLAUDE.md`, emit the claim sentence + the cited file's relevant section. Pure text, no model | auto |
| 2 | Write a Kaggle kernel (mirroring the 6 existing `scripts/kaggle_*.py`) that runs a **9B-14B reasoning model** over those pairs, emitting `{claim, verdict: supported/stale/unclear, evidence}` | auto |
| 3 | **Run it on the 10-item benchmark FIRST**, not the full corpus | auto |
| 4 | **Pre-registered pass condition: >=3 of 5 true positives AND <=1 of 5 false positives.** The mechanical screen scored 0/5 TP and 5/5 FP -- that is the bar to beat, and it is a low one | auto |
| 5 | **PASS** -> run the full 542-memo / 525k-token corpus; every flag adjudicated against its artifact before any doc is edited. **FAIL** -> close F2, record the negative beside the mechanical screen's, and keep artifact-reading as model discipline | auto |

**Why the pass condition sits there:** the value is catching real staleness, and a false positive costs one
adjudication. <=1/5 FP keeps signal-to-noise above the threshold at which a guard gets ignored -- the
failure mode that killed the mechanical version.

**Cost:** free (Kaggle T4, creds live). Two kernel runs. Fully reversible. No money.

**Kaggle gotchas already recorded, applied here:** pin `machine_shape: NvidiaTeslaT4` (`enable_gpu: true`
alone provisions a P100 that current Kaggle torch cannot run in fp16); max 2 concurrent kernels; set
`PYTHONUTF8=1` or logs return 0-byte.

### F1, for completeness

Not worth building now. Free GPU does not rescue it (batch latency, no inbound endpoint), local is capped
at ~4B, and Codex already fills the slot and demonstrably worked. If you want the offline voice anyway,
the original plan is in git history -- but it is a convenience, not a capability gain.

## 6. The variant worth more than the one asked about

If the goal is *catching my errors*, the higher-value move is not a weaker local model — it is a **second
frontier voice with a different lineage**. Codex already provides one. A third distinct frontier reviewer
would add more diversity per unit effort than a 4B, and needs no install, no VRAM, and no tok/s
measurement. The local model's genuine and *only* edge is **working with no network**.

## 7. Recommendation (corrected)

1. **Build the F2 falsifier on Kaggle** (§5). It is free, batch-shaped, targets the error class that
   actually dominates, has a ready-made 10-item benchmark, and a pre-registered pass condition so a null
   result closes the family instead of being rationalised. **This is the answer to your question.**
2. **Do not put an LLM in the decoder** (F3) -- a measured closed negative (ESM2 0.454, below chance, vs
   the catalog's 0.926), and a general LLM has weaker sequence priors than the protein model that failed.
   Free GPU does not change this; it was never a compute problem.
3. **Do not build F1 (interactive review) now.** Codex fills that slot and worked this session; free GPU
   cannot rescue it because Kaggle is batch with no inbound endpoint, and local is capped at ~4B.
4. **Do not use a model to curate the 40 colour loci** (F4). Each locus needs sourcing and verification
   regardless, so the model removes none of the cost and adds a plausible-looking-but-wrong failure mode.
5. **Standing principle, unchanged by any of this:** ~91% of this session's errors were unverified
   assertions caught by *running* things. Anything that makes checking cheaper beats anything that makes
   thinking longer -- which is exactly why F2 (a checker) outranks F1 (a thinker).