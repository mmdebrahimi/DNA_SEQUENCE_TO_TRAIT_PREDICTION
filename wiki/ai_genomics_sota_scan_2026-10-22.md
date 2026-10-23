# AI-genomics state-of-the-art scan — 2026-10-22

Scan run while the promcalc stage-B regeneration was mid-flight (18,000/43,611 tiles at the time of
writing). **Read-only.** Nothing in the frozen surface, the running job, or the pre-registered
promcalc bar was touched. Scope: the state of the art in **generative** and **representation-learning**
genomics, pointed at this project's own open levers rather than surveyed for its own sake.

> **AMENDED 2026-10-22 after an adversarial review of this memo.** Five claims below were
> wrong or overstated and are corrected in place, each marked **[AMENDED]**: the zero-shot vs
> supervised scope of F1; decision D2 (withdrawn -- it named a model this repo does not pin);
> the tier on the NAR "unasked question"; the nine-component split; and the "freezing is not
> the bottleneck" reading. The corrections are kept visible rather than silently rewritten.

Two tiers throughout, and they are not interchangeable:

- **MEASURED HERE** — I reproduced the number myself in this session. One item qualifies.
- **REPORTED** — read off a primary source (arXiv HTML, journal page). Not independently checked.

---

## DECISIONS NEEDED (only the user can settle these)

### D1. Adopt the NAR 112,955-promoter library as a second promcalc substrate?

This is a **scope decision on a frozen, pre-registered plan**, which is why it is a question and not
something I did. The substrate is real, free, MIT-licensed, and I verified it (see M1). It would be a
*new arm*, not a change to the running comparison.

What makes it worth asking about: the running arm tests whether thermodynamic decomposition adds over
GC on **44,106 × 150 bp** MPRA tiles. This library is **112,955 × 29 bp** cores where the *randomized
variable is the 17 bp spacer and the two hexamers* — exactly the inputs to `dG_35`, `dG_10`,
`dG_spacer`, `dG_ext10`. So it tests a **subset** of the decomposition, at 2.56× the n, with a
16,227-fold label range.

**It is NOT a drop-in replacement, and the reason is structural.** 29 bp = 6 (−35) + 17 (spacer) +
6 (−10), measured. There is no room upstream of −35 for the **UP element** and none downstream of −10
for the **discriminator** or **ITR**. Our stage B calls `c.run(seq, TSS_range=[0, len(seq)])` on 150 bp
tiles and gets ~73 candidate TSSs each; a 29 bp input cannot host the full motif architecture.

**[AMENDED]** The nine `dG_*` components split **4 / 3 / 2**, and the first version of this memo
said "roughly four computable and three absent", which silently dropped two:

- **4 clean:** `dG_10`, `dG_35`, `dG_ext10`, `dG_spacer`
- **3 structurally absent:** `dG_UP`, `dG_disc`, `dG_ITR`
- **2 composites that become partially-defined:** `dG_total`, `dG_bind` -- they *sum over* the
  absent terms, and `dG_total` is the headline feature. Under the fixed-context construction
  below the absent terms are constant, so they are absorbed by the intercept and contribute
  nothing; the arm would test four components plus two contaminated composites.

If you want it, the clean construction is: **embed every core in one fixed flanking context**, then run
the calculator. The context contribution is then identical across all 112,955 sequences and cannot
confound the comparison. The NAR authors themselves tested 8 designs across 6 contexts and report
R = 0.93 across them (REPORTED), so context effects are modest in their hands.

### D2. ~~Re-download the pinned GENERator weights?~~ **WITHDRAWN [AMENDED]**

**There is nothing to re-download and nothing to decide.** This decision is kept visible rather
than deleted, because the record should show the claim was made and retracted.

It was wrong twice over. It named `GENERator-prokaryote-1.2b-base`, but the code pins
`GENERator-**v2**-prokaryote-1.2b-base` at five sites (`dna_decode/glm/generate.py` ×2, both
Kaggle sweep kernels, `tests/test_glm_generate.py`) — a different HuggingFace repo id, and the
non-v2 id does not resolve publicly (`Invalid username or password` from the model API). And the
pin is already current, measured rather than argued:

| | revision |
|---|---|
| local snapshot | `8b2f768b0d293953518ff91d34600f9322ef1f94` |
| upstream `GenerTeam/GENERator-v2-prokaryote-1.2b-base` HEAD | `8b2f768b0d293953518ff91d34600f9322ef1f94` |

Identical, with upstream `lastModified 2026-07-05` — already inside the July window the original
note worried about. The whole question dissolved on one grep and one API call, which is what
should have happened before it was ever put to a reader.

---

## CORRECTIONS

### C1. My own near-miss, caught before it shipped

My first pass at the NAR label range printed **25,114,130-fold**. That was an artifact: the column-B
header cell carries `t="s"` and stores a *shared-string index* (112956) inside `<v>`, so a regex that
reads every `<v>` as a number treats the header as the largest measurement. Excluding `t="s"` gives
**16,227-fold**, which reproduces the paper's 16,226 and confirms the file's identity. The tell was
that the "max" equalled the row count exactly.

Same class as the PEAR G6 near-miss (a normalizer measured against itself) and the HIV `L234L`
self-to-self entries: **a degeneracy or range figure has to be asked what the extreme rows actually
are, not just how large they are.**

### C2. Nothing here corrects a project conclusion — it corroborates two

Stated plainly because the opposite would matter more: this scan found **no external result that
overturns** the GLM-arm negatives. **[AMENDED]** It found two that point the same way *on a
different readout* — both are zero-shot likelihood probes, not supervised-representation tests,
so they are corroboration of direction and **not** replications of our gate. The first version of
this line said "independently reproduce them", which overstates it in exactly the way this repo
has already paid for once (see F1).

---

## M1 — MEASURED HERE: the NAR core-promoter library is public and is what it claims

Source: *Deep learning guided programmable design of E. coli core promoters from sequence architecture
to strength control*, Nucleic Acids Research 53(16) gkaf863 (2025).

The paper's Data Availability statement deposits only **BioProject PRJNA1189185** (raw NGS for library
L4) and a Zenodo code DOI — it does **not** say the 112,955 measured strengths are released. They are,
in the code repo. Verified end to end:

| property | measured value |
|---|---|
| file | `Transformer/data/source data.xlsx`, 3,379,617 bytes, MIT licence |
| sheet dimension | `A1:B112956` → header + **112,955 data rows**, rows contiguous |
| columns | `Promoter_x` (sequence), `strength` (measured) |
| sequences | 112,955 strings, **every one exactly 29 bp**, pure ACGT |
| architecture | 6 + 17 + 6 — hexamer, 17 bp spacer, hexamer |
| −10 box | varies: TATAAT 49,870 (44.2%), TATACT 6,255, TATATT 5,001, TATAGT 4,842, … |
| −35 box | diverse — modal hexamer TAGACT at 950 (0.8%) |
| strength | min 0.00449771 · median 1.3208 · max 72.9837 |
| fold range | **16,227×** — reproduces the paper's 16,226× |
| GC fraction | mean 0.3814, sd 0.0947, range 0.000–0.759 |

Two consequences worth keeping:

1. **The paper's own description is looser than its data.** It describes a *fixed* 6 bp −10 box
   (TATAAT); in the released library that holds for 44% of sequences. The −10 box is a varying
   position, which makes the substrate *better* for our purposes, not worse.
2. **GC genuinely varies** (sd 0.095, spanning 0.00–0.76). That matters specifically for us: the
   promcalc question is "does decomposition beat a GC baseline", and a substrate where GC were
   near-constant could not pose it. Here it can.

Reported model numbers on this substrate, for calibration (REPORTED, not checked): Transformer with
505,531 parameters reaches held-out Pearson **R = 0.87**; CNN 0.76; LSTM 0.81. A conditional diffusion
model paired with the predictor reaches design-to-measurement **R = 0.95**.

**The gap the field left open, and it is ours:** the authors used an RNAP-binding-energy matrix
(Phillips et al.) to *design* the −10/−35 boxes and report **R = −0.75** between binding energy and
measured strength — but, so far as the parts of the paper I read show, they **never benchmarked
that biophysical model as a predictor against their Transformer.** That would make it the promcalc
arm's question, on their substrate, unasked.

**[AMENDED] Tier: `[inferred]`, not a finding.** I read characters 0–200,000 of a 209,274-character
page; **the final 9,274 were never read**, and a data-availability or comparison note is exactly the
kind of thing that lives in a tail. The first version of this memo asserted it as fact and the
phrase reached a commit message. Reading the tail would settle it; this memo does not.

---

## F1 — Two independent 2026 results point the same way as our encoder-gate negative, on a
different readout **[AMENDED]**

**The scope rail first, because the original heading got it wrong.** This section was titled
"…results *reproduce* this project's encoder-gate negative". They do not, and the distinction is
one this repo has already paid three months for. **Both external results are ZERO-SHOT
LOG-LIKELIHOOD probes:** the Mechanistic Invariance Test's CSS/SCR are LL-derived (it reports its
effects in "LL units" and an `AT-LL r` column), and Nullsettes defines a detection as *"mutant
log-likelihood significantly lower than its nonmutant original."* **Our encoder gate is not a
likelihood probe at all** — `wiki/glm_encoder_gate_result_2026-10-08.md` contains the word
"likelihood" **zero times** (measured), and the comparison is a *fitted* conv encoder against one
GC feature, judged by per-seed deltas and a paired CI.

CLAUDE.md names this exact conflation and its price: *"the unscoped version cost three months… The
body below is entirely about ZERO-SHOT likelihood."* And the repo holds a measured case where the
two readouts give **opposite verdicts** — on the same data, supervised **0.8142** vs zero-shot
**0.4750**, a PASS where zero-shot FAILED. So a zero-shot result is evidence about *what the
pretrained representation carries*, never about what a supervised head can extract from it.

What survives, and it is still worth having: two independent 2026 benchmarks find that gLM
regulatory signal is compositional, which is **consistent in direction** with our supervised
negative and is **not a replication** of it.

**Our result (2026-10-08):** a conv encoder over raw sequence vs **one GC number** — median +0.0122,
paired CI includes zero, 6/10 seeds. Recorded as a valid negative; the learned representation did not
land, so the self-distillation ban stayed.

**External result 1 — the Mechanistic Invariance Test** (arXiv 2604.06549, Cheng & Zhang, 8 Apr 2026,
ICLR 2026 workshop; REPORTED from the arXiv HTML). A 650-sequence, 8-class benchmark of **E. coli σ70
promoters** with scrambled controls that hold nucleotide composition fixed while moving an element's
*position*. Five gLMs across autoregressive, masked and state-space architectures:

| model | params | CSS | SCR |
|---|---|---|---|
| HyenaDNA | 6.6M | 0.63 | 0.48 |
| Evo2-1B | 1B | 0.60 | 0.46 |
| NT-500M | 500M | 0.54 | 0.40 |
| GROVER | 117M | 0.52 | 0.52 |
| Caduceus | 256M | 0.49 | 0.42 |
| *random baseline* | — | *0.50* | *0.46* |
| **position-aware PWM** | **~100** | **1.00** | **0.98** |

Apparent positional sensitivity correlates with **AT content at r = 0.78–0.96**. A ~100-parameter PWM
scores perfectly. The authors' reading: **misaligned inductive bias, not insufficient capacity** — and
larger models showed *stronger* compositional bias.

*The honest caveats, which the paper itself supplies:* the headline **46× compositional-over-positional
ratio is not like-for-like** (21.0 LL units across an AT sweep vs 0.46 LL units across a position
sweep — different manipulations, one model). The **strand-blindness claim is weaker than it reads**:
HyenaDNA's 44% forward-vs-revcomp accuracy has the authors' own binomial p = 0.24, so it supports
"no better than chance", not "prefers the wrong strand". And the PWM's hardcoded positions invite a
"true by construction" objection — the authors answered it with a position-free variant (CSS 1.00,
SCR 0.92), which is the right non-vacuity check and it holds.

**External result 2 — Nullsettes** (arXiv 2506.10271v3; REPORTED). Twelve gLMs / 28 variants
(Evo 1/1.5, Evo2-7B, **GENERator-3B**, GENERanno-0.5B, NT variants, DNABERT-2, Caduceus, GPN, gLM2,
METAGENE-1, PDLLM) scored on synthetic expression cassettes whose control elements are translocated so
that no flanking sequence can rescue expression. **9 of 12 fall below 50% detection on at least one
dataset.** Detection depends strongly on the *unmutated* sequence's likelihood, and the required
threshold scales linearly with length. Their conclusion: *"Current gLMs recognize what natural
expression cassettes look like, but not how they work."* Scale did not reliably help within families.

**What the three share, stated at the width the evidence supports [AMENDED].** The shared finding
is about **what the pretrained representation carries**: on regulatory sequence the signal is
compositional (GC/AT), and position/mechanism is largely absent. The MIT measured that as
"apparent sensitivity is AT-content correlation"; Nullsettes as "detection tracks likelihood, not
function"; ours as "a *fitted* conv encoder does not beat one GC number". The first two are
statements about a **frozen model's likelihood**; ours is a statement about a **trained
representation**. They are consistent, and they are not the same experiment — the original
"three teams, three benchmarks, one mechanism" compressed that away. What the agreement does
explain is why a GC baseline is hard to beat here.

**It also matches the project's scale findings.** ESM2 peaks at 650M and regresses at 3B/15B here; the
MIT reports larger models with stronger compositional bias; Nullsettes reports no consistent gain with
size. Three more agreements.

---

## F2 — The one asymmetry where this project is *ahead* of the published work

**Frozen but Not Always Accessible** (arXiv 2608.05329; REPORTED) probes frozen gLM embeddings —
DNABERT-2 117M, NT 500M, HyenaDNA, **GENERATOR-v2 1.2B**, Omni-DNA — with LR / RBF-SVM / shallow MLP.
Result relevant to us: on **promoter tasks, frozen embeddings recover 95–100% of fine-tuned
performance** (splice sites much worse, 60–88%).

**[AMENDED] What that does and does not license.** The first version concluded it "retires 'we
should have fine-tuned' as an explanation for our encoder-gate result". That over-reads on two
axes. The paper's promoter tasks are **classification** — its own limitations say evaluation is
classification-only — while our gate is **expression regression**. And our conv encoder was
**trained from scratch, not a frozen pretrained LM**, so frozen-vs-finetuned is not the axis our
gate sits on at all. The honest reading is narrower: on promoter-shaped *classification*, frozen
embeddings lose little to fine-tuning.

But the paper states plainly that it **does not compare against simple baselines such as GC content or
k-mer features**, and therefore "cannot say whether frozen embeddings outperform hand-built features."

**That is precisely the comparison this project ran.** Our encoder gate is a paired,
multi-seed, CI-reported test of a learned representation against one GC number. The field's own
representation-analysis paper names that comparison as the thing it did not do. Our negative is
therefore not a lagging result — it is a measurement the published benchmark is missing.

Corroborating the same direction from outside genomics-proper: **Specialized Foundation Models Struggle
to Beat Supervised Baselines** (arXiv 2411.02796) reports a 10.5M-parameter supervised workflow (DASHA,
no pretraining) beating all FMs on aggregate across the 18 Nucleotide Transformer tasks — though its FM
numbers are *copied from the NT paper's tables* rather than re-run, and the margin is driven largely by
histone-modification tasks. And *Nature Methods* (2025) reports deep-learning gene-perturbation
prediction still not outperforming simple linear baselines.

**One counter-example, recorded so this is not a one-sided file:** an organism-level-holdout
virulence-factor benchmark found LoRA-fine-tuned DNABERT-2 at mean F1 0.18 / AUC 0.54, statistically
indistinguishable from k-mer logistic regression at AUC 0.51 — *but* protein-level representations
outperformed nucleotide models on the same task. So the failure is specific to the **nucleotide**
modality, not to learned representations as such. (ScienceDirect blocked a direct fetch; these numbers
are from the search summary and are the least-verified in this memo.)

---

## F3 — Generative DNA: where the state of the art actually is

The project's measured bottleneck is the **generator**, not the scorer (2026-10-07: the oracle converts
ρ 0.59 into 79% of achievable selection value at k=10). So this is the part that matters most.

**Nobody has published a generative foundation model that beats a task-specific designer on
wet-lab-measured regulatory activity.** The sharpest datum (REPORTED): **GENERator**'s designed
enhancers, validated by UMI-STARR-seq, beat the strongest *natural* sequence by >2× — and still came in
**45% below the top DREAM-challenge-designed sequence.** A generative genomic foundation model loses to
a supervised task-specific designer on the designer's own metric. That is the same shape as F1/F2 one
level up.

What is actually working, by class:

| approach | status | relevance here |
|---|---|---|
| **Small unconditional diffusion** — DiPro (~May 2026): 50 bp *E. coli* promoters from noise, >90% called real by an independent recognizer, wet-lab validated driving a lysis protein with phase-specific tunable lysis | works, small | **Most T4/CPU-feasible generative route found.** 50 bp and no auxiliary discriminator. |
| **Conditional diffusion + predictor** — the NAR platform above: design-to-measurement R = 0.95 | works, validated | Pairs a generator with a *validated oracle* — our own design-scored-by-oracle shape |
| **Multinomial diffusion / PromoterDiff** — E. coli promoter design, in-vivo activity for a majority of designs | works | corroborates the diffusion-for-prokaryotic-promoters class |
| **RL fine-tuning of an AR gLM** — ICLR 2025, HyenaDNA + RL for yeast promoters / human enhancers | computational only | frames offline-MBO vs active learning; **no wet-lab validation** (2 incidental mentions in the full text I extracted) |
| **Dirichlet flow matching / GFlowNets / score-conditioned generators** | method literature | the general design-under-oracle toolkit |
| **Scale leaders** — Evo 2 (40B, 1 Mb context, 9.3T bp, *Nature* Mar 2026); a 20B variant released Feb 2026 fitting a single H100 | unreachable here | **still not a T4.** CLAUDE.md's hosted-API-only note stands. |
| **AlphaGenome** (*Nature*, 28 Jan 2026) — SOTA on 22/24 track and 25/26 variant-effect tasks, API-available | scoring, human | confirms the standing read: training method transfers, architecture does not |

**The actionable shape:** every generative success above is **small, task-specific, and paired with a
validated oracle** — diffusion over a short window, scored by a predictor fit on real measurements. That
is the architecture this project already argued for from the oracle-selection result, and it is what the
field's wet-lab-validated wins look like. The un-pulled lever is *not* a bigger pretrained model.

**One concrete model note:** GENERanno-0.5B matched Evo2-7B on Nullsettes at **14× fewer parameters and
12× less pretraining data** (REPORTED), with prokaryote weights on HuggingFace
(`GenerTeam/GENERanno-prokaryote-0.5b-base`, 8 kb context, single-nucleotide resolution). If a learned
representation is ever revisited, that is a materially cheaper starting point than the 1.2B pin.

---

## What this scan does NOT support

- **It does not reopen the GLM encoder gate.** Two external results agreeing with our negative is
  reason to trust it *more*, not to re-run it.
- **It does not predict the promcalc verdict.** The pre-registered bar is unchanged and the reconcile
  gate is still the validity anchor. `gc_seed0` already passed at exactly 0.3000.
- **It is a scan, not a systematic review.** ~12 searches and ~12 fetches, biased toward prokaryotic
  promoters and generative design because that is where this project's open questions are. Eukaryotic
  regulatory genomics, single-cell, protein structure and clinical variant interpretation were touched
  only where they bore on those questions.
- **Three fetches were blocked** (ScienceDirect 403, PMC cookie wall ×2, Nature IDP redirect), so the
  virulence-factor and *Nature Communications* benchmark numbers rest on search summaries rather than
  primary text. They are flagged inline.
- **The MIT and Nullsettes numbers are single, un-replicated studies.** The MIT is a workshop paper; a
  third-party review flagged that its scrambled controls may not isolate cleanly. Its direction agrees
  with our own independent measurement, which is why it is load-bearing here — not its authority.

---

## Sources

Primary, fetched: [arXiv 2604.06549 HTML](https://arxiv.org/html/2604.06549) ·
[arXiv 2506.10271v3](https://arxiv.org/html/2506.10271v3) ·
[arXiv 2608.05329](https://arxiv.org/html/2608.05329) ·
[NAR gkaf863](https://academic.oup.com/nar/article/53/16/gkaf863/8246948) ·
[Zenodo 15737101](https://zenodo.org/records/15737101) ·
[YuDengLAB repo](https://github.com/YuDengLAB/Intelligent_design_of_Escherichia_coli_core_promoters) ·
[arXiv 2503.07981](https://arxiv.org/pdf/2503.07981) (full text extracted locally)

Reported via search: [arXiv 2411.02796](https://arxiv.org/pdf/2411.02796) ·
[GENERator](https://arxiv.org/html/2502.07272) ·
[GENERanno weights](https://huggingface.co/GenerTeam/GENERanno-prokaryote-0.5b-base) ·
[GENERator repo](https://github.com/GenerTeam/GENERator) ·
[DiPro](https://pubmed.ncbi.nlm.nih.gov/42212299/) ·
[DNA-Diffusion, Nat Genet](https://www.nature.com/articles/s41588-025-02441-6) ·
[AlphaGenome, Nature](https://www.nature.com/articles/s41594-026-01763-1) ·
[Evo 2 one year later](https://arcinstitute.org/news/evo-2-one-year-later) ·
[NABench](https://arxiv.org/pdf/2511.02888) ·
[Nat Commun DNA-FM benchmark](https://www.nature.com/articles/s41467-025-65823-8) ·
[Nat Methods perturbation baselines](https://www.nature.com/articles/s41592-025-02772-6) ·
[virulence-factor organism holdout](https://www.sciencedirect.com/science/article/abs/pii/S1476927126004548) ·
[multinomial diffusion promoters](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11550136/)
