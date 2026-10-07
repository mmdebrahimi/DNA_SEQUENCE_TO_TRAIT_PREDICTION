# The GLM: prior-art scan and what it changed

**Date:** 2026-10-06. **Design it checks:** `wiki/glm_architecture_design_2026-10-06.md` (written first, on
purpose). **Decomposition:** `wiki/glm_decomposition_2026-10-06.md`.

**Scan coverage, stated honestly.** Five parallel research agents were launched; **two completed, three died
on a session rate limit** (generative-genome-LM survey + the MBO/oracle/MLDE arm landed; the regulatory-design
/ expression-dataset arm, the construct-design-tools arm and the terminator arm did not). So the **expression-
dataset question — the single highest-information unknown in the decomposition — is STILL OPEN**, and nothing
below should be read as having settled it. Both completed reports tag every figure as fetched vs snippet vs
unverified; the unverified ones are not repeated here as fact.

---

## 1. The model choice is settled, and it is not Evo 2

| Candidate | Verdict |
|---|---|
| **GENERator-1.2B-prokaryote** (`GenerTeam/GENERator-v2-prokaryote-*`) | **ADOPT.** MIT licence, decoder-only *generative*, ~2.4 GB in bf16 → fits the free T4 and plausibly the local 4 GB card, prokaryote weights refreshed 2026-07-06, **and it is wet-lab validated for exactly our use case**: ~12,026 oligos by UMI-STARR-seq, designed enhancers **+35% over the strongest natural** housekeeping-promoter enhancer and **>2× natural** for developmental. It is the only candidate that is simultaneously licence-clean, hardware-feasible, prokaryotic and experimentally validated at regulatory design. |
| **Evo 2** | **Ruled out locally; usable only via the hosted API.** 7B is ~14 GB bf16 against ~15 GB of T4, before KV cache; and the official path needs FlashAttention, which does not support Turing (sm_75). Perversely the **1B is worse**: it requires FP8 via Transformer Engine, i.e. a Hopper GPU. A pure-PyTorch community port (~2 GB) fits a T4 but **does not fix the FP8-numerics reason Arc flagged that checkpoint**, so its likelihoods are untrustworthy until cross-checked. Apache-2.0 and the hosted NIM endpoint remain a fine comparison arm. |
| **NTv3_generative** | **Cannot depend on it**, though it is the closest published analogue to the north star — masked discrete diffusion conditioned on species **and an activity level 0–4**. It is **gated and non-commercial**. |
| **megaDNA** | Keep as a **negative control**, not a generator: it is the one model in this space with a published independent falsification. |

## 2. The architecture is validated prior art — cite it, do not claim it

Evo 2's chromatin-accessibility result **is** this design's L1→L2→L3 loop: generate 128-bp chunks, score
against a target with an **Enformer + Borzoi oracle ensemble**, beam-search keep-best, continue — then
synthesize and measure by ATAC-seq, **experimental AUROC 0.92–0.95**. Design quality scaled log-linearly with
beam width. So propose→score→select with an oracle committee is **published, wet-lab-validated practice**, and
claiming novelty for the shape would be an over-claim.

**What the scan did NOT find is the organism-level link:** no instance of *generative sequence model proposes
an edit → a genome-scale metabolic model scores the organism-level phenotype*. The nearest neighbours are
OptKnock-style MILP strain design — **which this repo already has** (`fba/design.py`, reproducing the
literature anaerobic-succinate design). That composition is the genuinely unclaimed contribution here, and it
should be stated once, with evidence, not as a boast. Caveat: three agents died, so "not found" is weaker
than it would otherwise be.

## 3. Four design changes forced by the scan

1. **Never use generator likelihood as the oracle.** gLM2's BGC decoder redesigned a chimeric PKS to a 9.4×
   titre improvement — and reports that **neither gLM2 likelihood nor AlphaFold confidence predicted which
   designs worked**. That kills the cheapest tempting shortcut. (An opposing 2026 phage preprint reports
   likelihood *does* predict viability; treat the question as open and do not build on either.)
2. **The naturalness falsifier must be far harder than I wrote it.** The design's bar was "constraint-pass
   rate above a 3-mer Markov baseline". Two independent 2026 audits show that is too weak: generated
   sequences stay **trivially separable** from natural ones — megaDNA at **93% sens / 97.9% spec**, and a CNN
   at **AUROC 0.97 (eukaryote) / 0.82 (prokaryote)** against Evo 2 *and* megaDNA, with long-range
   organisation the named failure. **Add a held-out discriminator check.**
3. **Report a committee, never a single oracle number.** *Overconfident Oracles* measured that oracles
   differing **only by training seed** produce conflicting method rankings — on GFP, **three different
   methods are "SOTA" across five seeds** — and that GFN-AL ranks **1st under ESM-1b and 10th under
   Design-Bench/TAPE**. Worse for us: the TFBind-8 oracle's largest errors are **at the extremes of the score
   distribution**, which is exactly where an optimizer pushes.
4. **Adopt GILC's disjoint-oracle protocol as the minimum hygiene bar:** train two oracles on
   **non-overlapping** data subsets; one guides the search, the other only scores. This is the concrete form
   of the design's "held-out oracle" guard, and it is already standard practice in the strongest papers.

## 4. The mandatory baseline panel — where most published work in this space dies

Every GLM number must ship beside these, because the literature says they win more often than not:

- **A site-independent model, which cannot represent epistasis by construction, beats EVERY protein language
  model at mutational depth 2, 4 and 5+** (ProteinGym Table A9). ESM-1b falls **0.384 → 0.149** from depth 1
  to 4, −61%.
- **Mutation count — a single integer — reaches Spearman 0.45 on GFP; BLOSUM62 reaches 0.50**, against 0.60
  for a fine-tuned protein LM. Three-quarters of the headline is recoverable from an integer.
- **Plain linear regression ties ESM exactly on TAPE fluorescence, 0.68 vs 0.68.**
- **Random splits inflate results ~3×:** GB1 1-vs-rest best 0.32 → 0.92 on a random split, where an
  *untrained* ESM scores 0.79 — higher than anything achieves on the honest split. Three independent
  measurements of this artifact (FLIP, ProteinGym, PRIME R² 0.90 → 0.23).
- **One-hot encoding matched or beat the best pLM in nearly all settings** on mutation-dense landscapes,
  where PLMs sit at **ROC-AUC 0.5–0.6, i.e. chance** (Talpir & Fleishman 2026).
- **Randomly-initialised models scale in downstream performance in similar proportions to pretrained ones**,
  dissolving the apparent benefit of pretraining scale (Li & Lu 2024). Peer-reviewed corroboration: PLM
  fitness performance **declines beyond a certain size** (Hou 2026, *Nat Comput Sci*).

**Consequence: biologically-motivated splits only, never random; and no GLM result is reportable without
mutation-count, BLOSUM62, one-hot-ridge, site-independent and a shuffled null beside it.**

## 5. The finding that cuts in the project's favour — and it is the strongest argument for building this

The same literature that is brutal about *scoring* contains a measured, specific argument for *generative
design*, and it has to be reported because it is the pro-GLM evidence:

> **Johnston et al. 2024, *PNAS*** — a combinatorially complete **160,000-variant** landscape across four
> residues of an enzyme active site (TrpB). Verbatim: significant epistasis and many local optima *"prevent
> simulated directed evolution approaches from efficiently reaching the global optimum"*, and **"the most-fit
> TrpB variants contain a substitution that is nearly absent in natural TrpB sequences — a result that
> conservation-based predictions would not capture."**

That is the whole case in one sentence. **The best variants are evolutionarily unprecedented, so any method
keyed on natural conservation is structurally blind to them** — and that blindness is shared by a curated
catalog, by every evolutionary-likelihood scorer, and by BLOSUM62. It is the same shape as this repo's own
measured result that resistance is reached by chemically *conservative* substitutions that every
exchangeability scorer calls benign.

Corroboration that design-side ML genuinely pays when the landscape is hard: on GB1, **92% of the landscape is
holes**, greedy directed evolution reaches the global optimum **1.2%** of the time, MLDE with random training
**8.8%**, and informed-training-set MLDE **up to 92% — a 77-fold improvement**.

**The reconciliation, which is the point:** the negative results are about **scoring existing variants**; the
positive results are about **designing new ones**. Those are different tasks. The literature is not ambiguous
once they are kept apart — which is exactly the distinction that `orientation error 0` exists to protect.

## 6. Named open questions

1. **Is there a free, measured, large bacterial regulatory-sequence→expression dataset?** STILL OPEN — the
   agent that would have answered it died on the rate limit. This gates whether the expression oracle has a
   real supervised target or falls back to biophysics, and therefore gates the honest strength of the whole
   claim. **Highest-information next action.**
2. **There is no standard benchmark for generative DNA design.** The nearest usable protocol is
   **Design-Bench TF Bind 8/10**, whose ground truth is an exhaustively measured table — small (8-mers) but
   pre-registered and complete, which is the shape this repo's falsifiers already demand.
3. **Nobody has drawn a reward-overoptimisation curve for a DNA or protein oracle.** The RLHF
   overoptimisation literature and the biological-design literature are citation-disjoint in everything
   fetched. That is a real seam, and it is directly relevant to the oracle-hacking risk this design carries.
4. **No offline-MBO method has a verified wet-lab hit rate**, and a 2026 learnability result states there is
   *"a regime in which no offline method can avoid over-optimistic extrapolation"*. The one quantified
   in-silico→wet-lab attrition datapoint found: all 20 synthesized peptides had predicted MIC < 30 µmol/L and
   *"several showed weak or no detectable experimental activity"*.
