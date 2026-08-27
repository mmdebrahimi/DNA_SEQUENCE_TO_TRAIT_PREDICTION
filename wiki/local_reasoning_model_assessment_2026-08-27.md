# Would adding a local reasoning model (qwen9b-class) improve the answers? (2026-08-27)

**Short verdict: your instinct is right and already implemented — but not by a local 9B.**
The "second reasoning model" slot exists, is filled by Codex via `/brainstorm`, and demonstrably worked
this session. A local 9B cannot fit this host, and the error class it would target is not the one that
actually bites. One narrow variant is worth a cheap bounded test; three others are recommended against.

Every number below is measured on this host today, not recalled.

---

## 1. Hardware: the specific ask does not fit, and neither does the configured default

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

| # | family | what it would do | verdict |
|---|---|---|---|
| **F1** | **Adversarial-review diversity / offline voice** | a third opinion alongside Codex + me; works with no network | **worth a bounded test** — the only live candidate |
| F2 | Bulk triage | thousands of cheap classifications where frontier calls are wasteful | **no current workload** — nothing in the project needs volume judgment today, and any output would still need sourced verification, which is the expensive half |
| F3 | **Decoder-side prediction** (an LLM inside `dna_decode`) | predict phenotype / score variants | **closed negative — do not build** |
| F4 | Curation assistance | draft the 40 unrecorded colour-cell causal variants | **no** — this is precisely the fabrication hazard; each locus needs OMIA/literature sourcing and verification regardless, so the model removes none of the cost and adds a plausible-looking-but-wrong failure mode |

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

## 5. The plan for F1 — a bounded test with a decisive backtest

The point is not "install a model and see if it feels smart". It is a **retrospective backtest against a
case whose ground truth is known**: this session's wrong-η² error, which Codex caught and I had missed.

| step | action | gate |
|---|---|---|
| 1 | Install a runtime (ollama or LM Studio), **with the model dir forced to D:** — `OLLAMA_MODELS=D:/models/ollama`. C: has 8.1 GB left and a default install will eat it | dep-install → `irreversible`/reversible-outward, runs un-gated; Care resource check done |
| 2 | Pull **Qwen3-4B** (2.2 GB, GPU-resident) *and* an 8B (4.5 GB, CPU) — the two feasible tiers | auto |
| 3 | Measure real tok/s on each. A reasoning model emitting 1–3k thinking tokens at <5 tok/s is ~10 min/answer — that number decides usability, and I will not estimate it | auto |
| 4 | Point `LMSTUDIO_AUDITOR_MODEL` (or an ollama OpenAI-compatible endpoint) at whichever passes step 3 | auto |
| 5 | **The decisive test:** feed it the H2 memo *as it stood before the correction* and ask whether the η² attribution is sound. Codex caught this; I did not. Score: does the local model flag the score-side η² confound, the nested-subset p-values, or the definitional Control A? | auto |
| 6 | Verdict — **PASS** only if it catches ≥1 of the three real defects. Anything less means it adds noise to a slot Codex already covers | auto |

**Pre-registered so the result cannot be rationalised:** if the local model catches none of the three, F1
is closed and the auditor hook should be documented as unusable on this host rather than left as dead
config. If it catches one or more, it earns the auxiliary-voice slot in `/brainstorm`.

**Cost:** ~3 GB download, ~1 hour. Fully reversible. No money.

## 6. The variant worth more than the one asked about

If the goal is *catching my errors*, the higher-value move is not a weaker local model — it is a **second
frontier voice with a different lineage**. Codex already provides one. A third distinct frontier reviewer
would add more diversity per unit effort than a 4B, and needs no install, no VRAM, and no tok/s
measurement. The local model's genuine and *only* edge is **working with no network**.

## 7. Recommendation

1. **Do not add a 9B.** It does not fit the GPU, and the reasoning slot is already filled by something
   stronger that demonstrably worked this session.
2. **Do not put an LLM in the decoder** (F3) — that is a measured closed negative, and a general LLM is a
   worse fit than the protein model that already failed.
3. **If you want the local voice anyway**, run the §5 plan with a 4B and hold it to the §5 backtest. It is
   cheap, bounded, and pre-registered — but expect it to be a diversity/offline convenience, not a
   capability gain.
4. **The highest-leverage version of your instinct** is more verification, not more reasoning: ~91% of this
   session's errors were unverified assertions caught by running things. Anything that makes checking
   cheaper beats anything that makes thinking longer.
