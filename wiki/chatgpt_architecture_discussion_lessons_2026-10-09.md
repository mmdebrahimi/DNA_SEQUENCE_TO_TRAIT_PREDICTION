# Lessons from the ChatGPT architecture discussion, screened against what this repo has MEASURED

**Date:** 2026-10-09 · **Source archived verbatim:**
`wiki/refs/chatgpt_alphagenome_architecture_discussion_2026-10-09.md` · **Screened against:**
`dna_decode/eval/regime.py`, the GLM umbrella, and the committed artifacts named inline.

This is the screen, not the source. The source is a **secondary** source (a model's proposal), and its
sibling summary already contained two factually wrong citations — so every item below is marked by what
*this repo has measured*, not by how persuasive the proposal is.

---

## 0. SCOPE CORRECTION (user, 2026-10-09) — bacteria is the STARTING SUBSTRATE, not the domain

**The first version of this screen was wrong, and wrong in the direction that destroys the goal.** It
discarded ~6 of the 16 proposals as "structurally inapplicable to bacteria". The user's correction:
*"who told you we are just concentrating on bacteria? I said that is what we should start with if it is
easier."*

**The record already said this, in three places, and the screen contradicted all three:**

- **The north star is four steps and step 4 is `climb to animals`** — bacteria is **step 3**, a staging
  post. CLAUDE.md states it on the same line as the already-recorded sibling correction: ***"AMR was
  substrate #1, NOT the subject."*** The identical compression, caught once before, and I reproduced it.
- **The human track is LIVE NOW, not hypothetical** — 16 shipped cells (pgx 14, hla 1, mendelian 1). So
  diploid haplotypes, cell-state conditioning and CTCF/enhancer structure are not "someday"; they are
  applicable to code that already ships.
- **Our own measurement says organism complexity is NOT the barrier:** *"The discriminating variable is
  **population design, not organism complexity**."* And the sharpest fact — **the project's one clean
  organism-level g→p positive is in a EUKARYOTE**, the yeast segregant cross, **12/12 traits at
  r 0.46–0.80**, which is better than anything achieved in bacteria.

**So "bacteria is easier" deserves examination rather than inheritance.** The Arabidopsis eukaryotic
negative is often read as "eukaryotes are harder", but it was *zero-shot embeddings on a natural
population* and the model **learned population structure rather than the causal signal** — that is the
regime the map closes **for every organism**, not evidence about eukaryotes. Conversely the yeast positive
is a *constructed* panel. If population design is the discriminating variable, the right question is not
"which organism is simplest" but **"which organism has the best CONSTRUCTED panel available"** — and on
that criterion yeast segregants and human DMS/MPRA libraries compete with, or beat, bacteria.

**Why this error matters beyond bookkeeping.** Screening an architecture out because the *starting*
substrate lacks the mechanism is a textbook **scope cut outliving its constraint** — the failure mode this
project documents as having already cost it its most ambitious goal once (a May-2026 pivot made on a 4 GiB
GPU that was never re-opened when 16 GB became free). Had this screen stood, it would have discarded
exactly the architecture that step 4 requires, on grounds that expire the moment the substrate changes.

**Reclassified below:** the "structurally inapplicable" class is replaced by **`DEFERRED — applies at
step 4, and several apply to the human track TODAY`**, each with the staging reason named. Nothing is
screened out on organism grounds.

---

## 1. The discussion IS the plan's origin — confirmed from the ledgers, verbatim

| ChatGPT priority | Became | Originating goal (verbatim from the ledger) |
|---|---|---|
| **#1** multi-resolution + hierarchical tokens | **G-A** | *"Replace the expression oracle's fixed sequence features … with a LEARNED multi-resolution encoder plus a trivial head, and measure whether it transfers to real genomic sequence."* |
| **#3** (2nd half) cell-state conditioning | **G-C** | *"Build the bacterial analogue of **AlphaGenome's cell-state conditioning**: one oracle conditioned on growth medium, tested against two separate per-medium oracles."* |
| **#2** regulatory graph + sparse distal attention | — | **never built** |

So two of the three stated priorities were turned into project families, and **both reported NEGATIVE on
2026-10-08**. That is the headline — and the next section is why it does not mean what it looks like.

---

## 2. BOTH negatives tested the WEAKEST FORM of the proposal. Read this before re-proposing either.

This is the most important thing in this document, in both directions: it is the honest caveat on the
negatives **and** the reason not to simply re-run them bigger.

| | ChatGPT proposed | What we actually tested | Gap |
|---|---|---|---|
| **#1** | a pyramid **1 bp → 16 bp → 128 bp → 1 kb → 10 kb** with information travelling between levels, over a **1 Mb** context | a **13,217-parameter** conv at 3 kernel widths (6/20/75) over a **150 bp** tile | our substrate spans **~1 decade**; the proposal spans **four**. There is no 1 kb or 10 kb level to build on a 150 bp sequence. |
| **#3** | a **cell-state embedding modulating a learned representation** | a **linear interaction term** (`seq × medium`) over 3-mer frequencies | the artifact's own honest limit says it: *"not the same as a model whose **representation** is modulated by condition. A null here does not foreclose that."* |

### CORRECTION (peer review, Skill_Development session, 2026-10-09) — one leg of this was CIRCULAR

The first version of this section argued: *"a single GC number takes 33% of the 0.941 ceiling, **and the
best learned arm reached 34%** … so the honest consequence is change the substrate, not re-run it bigger."*

**The peer caught a real circularity, and I verified it before conceding.** The "best learned arm at 34%"
**IS** the 13,217-parameter conv — `encoder-multi`, median 0.3234, 0.3234/0.941 = **34.4%**,
`n_parameters` 13217, read from `wiki/glm_learned_representation_2026-10-08.json`. So §2 declared that arm
too weak to be a test, and the conclusion was then computed **from that same arm**. Split into the two
claims it was smuggling together:

| claim | status |
|---|---|
| A trivial feature (GC) captures **33% of the ceiling** (0.3112/0.941), so the substrate may carry little learnable structure | **SUPPORTED** — independent of any learned arm |
| The best learned arm reached only 34%, **so more capacity will not help** | **NOT SUPPORTED** — that arm was the proxy. **Struck.** |

**The 4-decade pyramid is UNRUN, not refuted.** It is now a `planned` arm in
`soraya_runs/arms-ledger.json`, so that ceiling criterion reports **`blocked:untried-arm`** and the weak
arm structurally cannot close it.

**But the conclusion "change the substrate" survives on a third leg that depends on NO arm at all — and
this is the part the peer's framing missed.** It is **dimensional, not inferential**: the proposal
specifies levels at **1 kb and 10 kb**, and a **150 bp** tile cannot host either — not because a bigger
model would fail, but because those levels *do not exist* in a 150 bp sequence. So the honest statement is
stronger than "run it or label it unrun": **on this substrate the proposed arm is not merely unrun, it is
UNRUNNABLE.** Running it *requires* a substrate with ≥3 decades of range, which makes substrate change a
**precondition of the experiment**, not a competing alternative to it.

(The remaining capacity-hurts data point — 6-mers, 4,096 features, **worst of every arm at 0.2022**,
`wiki/glm_genomewide_oracle_2026-10-07.json` — concerns a *fixed* feature set rather than a learned
representation with more capacity, so it is weak evidence here and is no longer load-bearing.)

---

## 3. The sixteen proposals, screened

Five verdict classes. The two that matter most are **ALREADY SHIPPED** (we are further along than the
proposal assumes) and **MEASURED NEGATIVE** (do not re-propose without new information).

### DEFERRED — absent in the step-3 substrate, REQUIRED at step 4, and several are live on the human track TODAY — ~6 of 16

**Not screened out.** A bacterium genuinely lacks a spliceosome, nucleosomes, TADs, CTCF, a second
chromosome copy and eukaryotic cell types — so these cannot be *tested* on an *E. coli* promoter. That is a
statement about the **current substrate**, not about the proposal's value, and the north star's step 4 is
exactly where each becomes necessary. The `apply NOW?` column is the one that matters: **four of the six
are applicable to code this repo already ships.**

| # | Proposal | Absent in bacteria because | Apply NOW? |
|---|---|---|---|
| 3 | hierarchical attention through **TAD / chromatin-domain** levels | no TADs — the L4–L5 levels do not exist | step 4. Human Hi-C/TAD data is public and abundant |
| 4 | sparse **distal** attention over 1–10 Mb | bacterial cis-regulation is overwhelmingly promoter-**proximal**; a 300 kb enhancer–promoter bridge is not a bacterial mechanism | step 4 — and it is the *defining* human regulatory mechanism, so it is unavoidable there |
| 5 | regulatory graph with **enhancer / CTCF** nodes | no enhancers in that sense, no CTCF | step 4. **Bacterial analogue exists and is untried:** the regulon / σ-factor / TF→operon network is a real graph, just *trans* rather than 3D-cis |
| 6 | 3D structure as a latent variable | no Hi-C-scale TAD analogue to bottleneck through (nucleoid supercoiling is real but not that) | step 4 |
| 10 | diploid haplotypes | bacteria are **haploid** | **YES, TODAY** — the human clinical track is 16 shipped cells (pgx 14, hla 1, mendelian 1), and phase / compound-het / allele-specific effects are exactly what pgx and ClinVar interpretation turn on |
| 16 | mixture-of-experts over splicing / neuronal / developmental experts | those domains do not exist here; and at 13k parameters MoE is premature by orders of magnitude | step 4. The *scale* objection stands independently of organism |

**Reusable, and it is the opposite of what the first version of this screen said:** when a proposal does
not fit the substrate you happen to have started on, the honest label is **DEFERRED with the staging
reason**, never *inapplicable*. "Inapplicable" silently converts a temporary substrate choice into a
permanent architectural exclusion — and if the starting substrate was chosen for *convenience*, that
exclusion has no expiry and will quietly delete the capability the real goal needs.

### ALREADY SHIPPED here — 3 of 16, and this is the strongest convergence

| # | Proposal | Status in this repo |
|---|---|---|
| **8** | use **perturbation** data, because observational tracks give correlation not causation | **This is already the project's measured winning regime.** `eval/regime.py`: natural-population zero-shot is a **closed negative** (0-for-5, de-confounded), while **constructed** variation **WORKS** — DMS ρ 0.35–0.76, FBA/Keio MCC 0.70–0.74. Our substrate *is* an MPRA (GSE144621). ChatGPT reasoned to this from theory; we measured it. |
| **9** | train on `ref + alt → Δ`, making variant effect first-class | **Shipped** as `dna-forward` with a genome-edit path (ρ **0.761** on TEM-1 ampicillin fitness, 1,715 variants). Independently, the project's own R8 rail says the same thing: *a differential readout on the same background cancels a confound present in both arms.* |
| **15** | add **uncertainty**; the model should be able to say "I have not seen enough biology like this" | **Shipped** as the **L2 doubt layer** (`dna_decode/eval/doubt.py`) — and it is unusual: the project's own survey notes most of the field, including Mykrobe, *predicts SUSCEPTIBLE when no known mutation is found*. The doubt block structurally **cannot contain a call** (`assert_no_call` runs on its own output). |

**So the three proposals a model reasons to from first principles are the three this project had already
built and validated.** That is real external corroboration of the architecture choices — worth more than
any of the untried suggestions.

### MEASURED NEGATIVE here — 3 of 16 (was 4; #12 moved out on peer review, see its row)

Do not re-propose these without new information; each has a committed artifact.

| # | Proposal | What we measured |
|---|---|---|
| **12** | self-supervised DNA pretraining, **then** functional supervision | **RE-SCOPED 2026-10-09 (peer review) — this was a CLASS-level kill from ONE checkpoint and it belongs in the NOT-REFUTED column.** What was measured: **GENERator-1.2B-prokaryote ONLY, at default sampling, ZERO-SHOT with NO fine-tuning** — distinguishability **0.7706** vs a 3-mer Markov chain's **0.6679** (lower is better; `wiki/glm_generator_falsifier_REAL_2026-10-07.json`, note the `_REAL_` suffix). That is **stage 1 of a 3-stage recipe, used alone**; stages 2–3 (functional + perturbational supervision) were never run. **Revival condition: any fine-tuned checkpoint scored on the same falsifier — which IS family G-D.** So this negative **argues FOR G-D**, it does not argue against pretraining. Structurally the same error as the June "0-for-4 frozen embeddings → *a MODEL ceiling*" filing, which this repo records as its most expensive scope error. |
| **13** | feed **evolutionary conservation** explicitly | **Measured structurally blind for our phenotype.** Resistance-mutation sites sit at ESM entropy percentile **0.494 — exactly average**, refuting the textbook "the site is poorly conserved" account; and **BLOSUM62, which has never seen an HIV sequence, ranks the same DRMs as well as ESM2** (4.0 vs 4.5 of 19). Johnston 2024 (*PNAS*, TrpB) independently: the most-fit variants carry a substitution *"nearly absent in natural sequences — a result that conservation-based predictions would not capture."* **Conservation-keyed methods are blind to the best variants.** |
| **11** | model variant **combinations** (epistasis) | **Measured to lose, on a DIFFERENT substrate — and the peer's objection to this one does not hold.** On **HIV-1 RT**, pairwise and non-linear estimators both lost to a linear one-hot model **4/4, two with CI below zero** — the linear model is the measured **ceiling, not a floor**. Corroborated by an independent 2026-07-11 epistasis negative. The peer argued this kill is confounded because *"the substrate is ceiling-bound near 34%"* — but that 34%/0.941 figure is the **MPRA peak-tile** substrate, whereas this result is **Stanford PhenoSense fold-change on HIV RT** (`data/raw/hiv/*_DataSet.txt`). Two unrelated substrates; the tile ceiling does not apply. **Scope (narrowed anyway): one virus, one protein family, one label type, clinical-isolate co-occurrence. Revival condition: an epistasis model tested on a CONSTRUCTED multi-variant panel (combinatorial mutagenesis) rather than on co-occurrence in natural isolates.** |
| **1** | multi-resolution tokens | **Tested as G-A → `WEAK`, not distinguishable from zero** (6/10 seeds, p=0.754, CI [−0.0099, +0.0349]). But see §2 — this tested a 1-decade version of a 4-decade proposal. |

### GENUINELY OPEN and worth the conversation — 3 of 16

| # | Proposal | Why it is live |
|---|---|---|
| **7** | cell-state conditioning **of the representation** | G-C's negative explicitly does **not** foreclose it, in the artifact's own words. The bacterial analogue of cell state is growth condition, and the substrate genuinely carries the signal: noise-matched cross-condition agreement **0.4063** vs within-condition **0.4677**, disattenuated true cross-condition **0.7887** — **~38% of real variance is condition-specific**. Something is there that a linear interaction term did not capture. |
| **2** | **learned biological tokens** rather than spatial windows | The bacterial vocabulary is different but real: promoter, operator, RBS, terminator, operon, regulon, σ-factor binding site. Discovering those as tokens is not an absent mechanism — it is an untried one. |
| **14** | mechanistic **intermediate** supervision | Partly present already: the repo has an FBA metabolic layer that *is* a mechanistic intermediate, and `constraints/` encodes axiomatic laws. The proposal's version — latent intermediates supervised wherever data exists — is untried. |

---

## 4. The two sharpest transferable ideas, and neither is an architecture

**(a) "Attention weights are not biological causality."** Quoted from the source, and it generalises far
past transformers: *a high attention score means the model found this region useful when constructing the
representation* — not that the enhancer controls the promoter. This repo has the same failure mode under a
different name: a feature's importance is not a mechanism, and the project's gates G1–G10 exist because
tool-derived and sampling-defined labels look like measurements. **Add it to the standing vocabulary: an
internal weight is evidence about the MODEL, never about the biology.**

**(b) "I would not build your model as a much bigger AlphaGenome."** The strongest framing in the
discussion, and it matches this repo's own north-star seam exactly. AlphaGenome stops at
`DNA → molecular phenotype`; the goal is `genotype → phenotype`, which needs protein effects → cellular
state → tissue → physiology → development → environment. The authors themselves say molecular predictions
do not capture complex traits because **development and environment lie outside scope**.

This is the same thing `eval/regime.py` encodes and the same thing lesson 6 ("prerequisite, not
sufficient") records: **the expression oracle predicts expression, which is a PREREQUISITE for a trait.**
ρ 0.59 and 79%-of-selection-value are not trait claims. The discussion arrives at the project's own honesty
rail from the outside, which is the best possible corroboration of it.

---

## 5. What this changes about the plan

1. **G-D (fine-tune the generator) remains #1, and this discussion strengthens it.** Proposal #12's
   sequencing is *pretraining → functional supervision → perturbational supervision*. We have measured the
   pretraining half failing at zero-shot and have never run the functional-supervision half. That is
   exactly G-D.
2. **Do NOT re-run G-A bigger on the same substrate.** The range the proposal needs (1 kb–10 kb levels)
   does not exist in a 150 bp tile, and GC already takes 33% of a 0.941 ceiling. If G-A is revisited, the
   substrate must change, not the parameter count.
3. **G-C's live successor is representation-level conditioning**, not a wider pooling grid — and the
   measured 38% condition-specific variance says the target is real.
4. **Nothing here unblocks G-B.** Self-distillation is blocked on a learned representation existing, and
   proposal #9's "train on ref+alt Δ" is already shipped by a different route.
5. **New candidate worth a plan: learned biological tokens** (proposal #2) in the bacterial vocabulary —
   a never-tried item that needs no new substrate.
6. **ADDED after the §0 scope correction — the two moves that only become visible once bacteria is read as
   step 3 rather than the domain:**
   - **Substrate choice is now an open strategic question, not a settled one.** Point 2 says G-A must
     change substrate rather than parameter count, and §0 says organism complexity is not the measured
     barrier — so the candidate substrates are a **constructed eukaryotic panel** (the yeast cross is the
     only clean organism-level positive here) and the **human DMS/MPRA libraries** that the forward cell
     already scores well on, not just a bigger bacterial tile set.
   - **The human track can absorb the DEFERRED proposals today.** Diploid/phase reasoning (#10) is
     immediately relevant to 16 shipped pgx/hla/mendelian cells, and cell-state conditioning (#7) is the
     *same* question G-C asked, in the organism where the proposal was actually designed for it. Testing
     conditioning on a human cell-type substrate is not a step-4 deferral — it is a different, available
     substrate for a question we have already built the harness for.

## 6. Honest limits on this screen

- The source is a **model's proposal**, not a paper. Its sibling summary had two wrong citations. Nothing
  above treats its claims about AlphaGenome's internals as verified — only the *ideas* are screened.
- The **DEFERRED** verdicts rest on standard prokaryotic biology, not on a measurement in this repo. They
  are high-confidence but **asserted**, not derived here — and per §0 they are statements about the
  *current substrate's* reach, never about the proposal's value or about eukaryotes being harder.
- **The "start with bacteria because it is easier" premise is itself unexamined**, and §0 gives the reason
  to doubt it: the only clean organism-level g→p positive here is a *eukaryote* (yeast, 12/12 traits), the
  eukaryotic negative (Arabidopsis) is explained by *population design* rather than by organism, and 16
  human cells already ship. Nothing in this screen tests the premise — it just stops inheriting it.
- The epistasis negative (#11) is scoped to **HIV RT**; it is not a general claim about epistasis.
- The conservation negative (#13) is scoped to **pathogen target-site resistance** — Friedman 2013
  supplies a counterexample class (EGFR/ALK resistance mutations tend to be *radical*), so it must not be
  stated as a general law.
- §2's "33% of achievable" is specific to **peak tiles at 150 bp**; it is not a statement about bacterial
  promoters in general.
