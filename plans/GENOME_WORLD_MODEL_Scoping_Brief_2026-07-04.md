# Genome World Model — Requirements, Architecture & Deep-Research Element List (2026-07-04)

> **Purpose.** Close the loop on *what it would take* to build a "genome world model" — before committing
> build compute — so parallel deep-research (farmed to different LLMs) can proceed **while the imputation
> path is built**. Every section ends with **farmable ELEMENTS (E#)**: self-contained deep-research prompts.
> **Status:** SCOPING (class-e, no build). Anchors on the project's embedding-track evidence
> (`plans/AMR_embedding_niche_decision_2026-06-05.md`, `wiki/arabidopsis_af001_stop_vs_escalate_decision_2026-06-21.md`,
> `wiki/celegans_ben1_result_2026-07-02.md`, `plans/EP8_Arabidopsis_Embedding_Test.md`,
> `plans/Hybrid_Learned_Deterministic_Decoder_Plan.md`).
> **Honest prior (do not bury):** frozen general-DNA-FM embeddings are **0-for-4** on the project's
> de-confounded phenotype tests; the *fair* plant-FM shot (EP8) is unrun, compute-gated. A genome WM must be
> designed to **earn its keep on the polygenic residual + generative/design axis**, not to beat mechanism
> scans head-on. This brief is how it *could* work, honestly — not a promise it will.

---

## 1. Disambiguation — "world model", CLIP, and masking are THREE different things

The request blends three distinct ML paradigms. Naming them precisely is the first design decision, because
they need different data, objectives, and eval:

| Paradigm | "Obfuscate & predict" analog | What it learns | Genomic instance |
|---|---|---|---|
| **Masked modeling (MLM / MAE)** | BERT / masked-image-modeling: mask tokens, predict them | contextual **representations** | DNABERT-2, Nucleotide Transformer, ESM-2 |
| **Autoregressive (AR)** | next-token prediction | representations **+ generation/likelihood** | Evo, Evo2, ProtGPT2 |
| **Contrastive multimodal (CLIP)** | align TWO modalities in a shared space (NOT masking) | **cross-modal** retrieval / zero-shot transfer | seq↔expression, seq↔structure, genotype↔image |
| **World model (JEPA / latent-predictive)** | predict **masked/future EMBEDDINGS**, not tokens | a **predictive latent model of the system** | I-JEPA analog on genome; genotype→molecular-state dynamics |

**What "a full-out genome world model" most usefully means here:** a *long-context, multimodal representation
+ predictive* model that (a) is pretrained self-supervised on sequence (masked and/or AR), (b) is **aligned
(CLIP-style) to molecular/organismal phenotype modalities**, (c) optionally carries a **JEPA latent-prediction**
objective enabling *counterfactual/in-silico perturbation* ("what if this allele flips?"), and (d) is **fused
with dna_decode's deterministic mechanism decoders** so the learned part only carries the *residual* the
mechanism scan can't. Each of (a)-(d) is an independent design axis, decomposed below.

> **E1 — Objective taxonomy for genome.** For each of {MLM, AR, CLIP, JEPA, diffusion}, what does the 2024-2026
> literature show it buys on genomic downstream tasks (variant-effect, expression, regulatory, phenotype)?
> Which combine well (e.g. MLM+AR like Evo2)? *Good answer:* a table objective→downstream-strength→cost, with
> ≥2 papers each. *Sources:* Evo/Evo2, Caduceus, ESM-3, JEPA papers, HyenaDNA, DNABERT-2.

---

## 2. SOTA sweep (early 2026 — VERIFY currency via farmed deep-research; flagged where my cutoff bites)

### 2a. DNA / genome foundation models
| Model | Backbone | Max context | Tokenization | Objective | Note |
|---|---|---|---|---|---|
| DNABERT-2 | Transformer | **up to ~10 kb** (BPE/GUE) | BPE | MLM | efficient; context depends on BPE compression |
| Nucleotide Transformer | Transformer | **~6 kb** (1000×6-mer) | 6-mer | MLM | multispecies (up to 2.5B); *frozen NT is the project's 0-for-4 negative* |
| HyenaDNA | Hyena (implicit conv) | up to **1 Mb** | single-nt | AR | long-context, single-nucleotide |
| Caduceus | **bi-directional Mamba (SSM)**, RC-equivariant | ~131 kb | single-nt | MLM | reverse-complement equivariance is genome-native |
| Enformer / Borzoi | Transformer / U-Net | ~200 kb / ~500 kb | bp | supervised (tracks) | regulatory/**expression from sequence** — the strongest "function from DNA" line |
| Evo-1 | **StripedHyena** (SSM+attn hybrid) | 131 kb | single-nt | AR | ~300 B tokens, prokaryotic; generative |
| **Evo2** | StripedHyena 2 | up to **1 Mb** | single-nt | AR (generative) | **~9.3 T nucleotides / ~128k genomes**, 1B/7B/**40B**; prokaryote→**eukaryote**; **generation + zero-shot variant effect** (the current frontier) |
| GENA-LM | Transformer | ~4.5 kb (BPE) | BPE | MLM | (wired in `foundation.py`) |
| AgroNT / **PlantCaduceus** | Transformer / Mamba | — | — | MLM | plant-specific — the EP8 target (≥24 GB GPU) |

**One architecture axis (not THE only one): CONTEXT LENGTH.** A gene+regulatory context is 10⁴–10⁶ bp; attention
is O(n²) → SSM (Mamba/Caduceus), Hyena/StripedHyena (Evo), **and dense-attention + context-extension (e.g.
Gene42)** are the long-context contenders. **BUT** (adversarial-review correction, §11): long context is *not* the
dominant bottleneck for phenotype signal — **causal labeling, population structure, assay noise, cell-type/
regulatory context, epistasis, and structural variation** bound it more. Don't over-index on context length.

### 2b. Protein FMs (relevant — you have `proteingym`)
ESM-2 (≤15B, MLM), **ESM-3** (multimodal *generative* — sequence+structure+function tokens), ProtT5, AlphaFold2/3
(structure), ProteinMPNN (inverse folding), ProtST (protein↔text CLIP). Protein is the *most mature* modality for
the "world model" idea (ESM-3 already is one).

### 2c. Expression / single-cell / functional
Geneformer, scGPT, scFoundation (rank-based single-cell); Enformer/Borzoi (bulk expression). Relevant for a
CLIP *phenotype* modality and for `depmap`/`gdsc` functional data.

### 2d. Multimodal / CLIP-style in biology
seq↔expression, genotype↔histopathology-image, protein↔text (ProtST), sequence↔structure. The CLIP recipe:
two encoders + a contrastive (InfoNCE) loss aligning paired samples. **The pairing choice is the whole game.**

> **E2 — SOTA currency check.** Refresh 2a-2d to *today*: newest genome/protein FMs, their context length,
> params, pretraining corpus size (tokens), and reported downstream SOTA. Flag any that beat Evo2/ESM-3/Caduceus.
> *Good answer:* updated tables + 1-line "what's new since early-2026" each.
> **E3 — Backbone for genome-length context.** Quantced comparison Mamba/Caduceus vs Hyena/StripedHyena vs
> attention+FlashAttention: max context achieved, throughput (tokens/s/GPU), memory scaling, RC-equivariance.
> Which is the best backbone for a solo/small-cluster genome WM? *Sources:* Caduceus, Evo2, HyenaDNA papers.

---

## 3. Data requirements — modality × granularity × scale

### 3a. The granularity ladder (your "any granularity data is available")
`nucleotide → k-mer/codon → gene → transcript/isoform → protein → variant/allele → haplotype → pathway →
cell-state → organismal phenotype`. Data DEPTH is very uneven across rungs:

| Rung | Data at depth? | Best public source | For which paradigm |
|---|---|---|---|
| nucleotide/sequence | **abundant** | RefSeq/GenBank, multispecies (Zoonomia 240 mammals), pangenome (HPRC) | self-supervised pretraining |
| variant/allele | abundant | gnomAD, 1000G, ClinVar (effect labels) | variant-effect eval |
| protein | abundant | UniRef/UniProt (250M+), PDB, **ProteinGym (DMS)** | protein FM + DMS eval |
| expression | good | GTEx (eQTL), ENCODE, Roadmap | CLIP phenotype modality |
| organismal phenotype (paired w/ genotype) | **SCARCE + confounded** ← the binding constraint | UKB / AoU / FinnGen / **OpenSNP** / **PGP** / model-organism GWAS | CLIP head + supervised eval |

**The bottleneck is NOT sequence data** (public genomes are plentiful for pretraining) — it is a **de-confounded,
sampling-independent, paired genotype↔phenotype label**, which has blocked the project 3× (pathotype circular
labels, cef geography confound, cipro within-lineage). Any WM's *supervised/CLIP* value is gated on this.

### 3b. Corpora by paradigm (what to actually assemble)
- **Self-supervised sequence pretraining:** reference + multispecies genomes (evolutionary signal is where NT/Evo
  get zero-shot variant effect); population variation (1000G/gnomAD/HPRC) for allele diversity.
- **AR/generative:** same, framed as next-token; Evo2 shows prokaryote→eukaryote transfer.
- **CLIP pairs (pick the pairing per goal):** genotype↔phenotype (biobanks incl. the ones you HAVE: FinnGen,
  OpenSNP, PGP), genotype↔expression (GTEx), sequence↔structure (PDB/AlphaFold DB), genome↔epigenome (ENCODE).
- **Protein WM:** UniRef + PDB + ProteinGym (you have proteingym cached).

### 3c. Scale numbers (order-of-magnitude — grounded in published FMs; farmed-research should tighten)
- Tokens: **Evo ≈ 300 B tokens** (prokaryotic); NT ≈ ~10¹² bp multispecies; **ESM-2 ≈ 65 M sequences** (~tens of B
  tokens). Chinchilla-optimal ≈ **20 tokens/param**.
- ⇒ *if* Chinchilla held: a **1 B-param** WM ≈ 20 B tokens min. **CAVEAT (§11):** 20:1 is a language-model prior,
  **not established for low-alphabet / repetitive / phylogenetically-redundant genomic corpora** — treat as a
  sanity check only; genomic entropy + dedup + evolutionary redundancy could shift it substantially. A 100 M-param
  domain model (~2 B tokens) is *feasible to assemble from one organism's population genomes*.
- Storage: human WGS ≈ 100 GB/genome raw, ~1–3 GB variant-only; a population corpus = **TBs** (⚠ your C: is at
  97% / 8 GB free — a corpus lives on D: or cloud, not C:).

> **E4 — Corpus assembly spec.** For a chosen target (e.g. plant flowering-time, or human PRS trait), specify the
> exact pretraining corpus: which genomes/populations, token count, dedup strategy, licenses, storage footprint,
> download plan. *Good answer:* a concrete manifest w/ URLs + sizes.
> **E5 — Scaling law for genome.** What is the empirical params↔tokens↔compute scaling for genomic FMs (is
> Chinchilla ~20:1 right for DNA, or does low sequence entropy change it)? *Sources:* Evo2 scaling, Caduceus.
> **E6 — The de-confounded paired-label problem.** Enumerate public genotype↔phenotype datasets that are
> (a) ≥ few-hundred samples, (b) sampling-independent measured labels, (c) de-confoundable (lineage/geography/
> ancestry). Rank by "cleanest CLIP/supervised substrate." *This is the project's #1 binding constraint — treat
> as the highest-priority element.* *Sources:* AraGWAS, DGRP, biobanks, model-organism panels.

---

## 4. Compute requirements — tiered, with $ estimates

| Tier | What | Hardware | Rough cost | Feasible for you? |
|---|---|---|---|---|
| **T0** | from-scratch *full* foundation (Evo2-40B class) | 100s-1000s A100/H100·months | **$1M–$10M+** | ❌ not solo |
| **T1** | small-mid genome WM from scratch (100 M–1 B params, focused corpus) | 8–64 GPUs × weeks | ~$10k–$100k cloud | ⚠ needs a grant/budget |
| **T2** | **continued-pretraining / fine-tune an open FM** (Evo2/Caduceus/NT) on your domain | 1–8 GPUs × days | ~$100–$5k, OR workhorse RTX 3500 Ada (~12 GB, borderline) | ✅ **near-term path** |
| **T3** | frozen-FM embeddings + light heads / **CLIP head** / adapters / hybrid-residual | 1 GPU or CPU inference | ~$0–$100 | ✅ **laptop-feasible now** (cached nt.h5/dnabert2.h5 exist) |

VRAM anchors: Caduceus/PlantCaduceus **≥24 GB for training** (RTX 3090/A100 in the papers); inference much less.
The workhorse RTX 3500 Ada (~12 GB) is the borderline T2 gate (EP8's open question — needs a VRAM check).

> **E7 — Compute/$ plan per tier.** Turn T1/T2 into concrete quotes: GPU-hours for N params × M tokens, spot vs
> on-demand $, wall-clock, and the minimum viable run. *Good answer:* a costed table + "smallest experiment that
> would prove/disprove the thesis."
> **E8 — Continued-pretraining recipe.** The exact recipe to domain-adapt an open genomic FM (LoRA/full-ft,
> masking ratio, LR schedule, context length, RC-augmentation) on a target organism's genomes. *Sources:* Caduceus/
> Evo2 fine-tuning repos.

---

## 5. Hypothesized best architecture (recommendation — draft-then-ratify)

Given (a) genome length → long-context, (b) the 0-for-4 frozen-embedding evidence, (c) solo/small compute, (d) the
hybrid-with-mechanism goal:

> **RECOMMEND: a long-context SSM-hybrid backbone (Evo2 / Caduceus lineage) + multi-objective SSL (masked *and*
> autoregressive) + a CLIP contrastive head aligning sequence-embeddings to a phenotype modality + a JEPA
> latent-prediction auxiliary for counterfactual perturbation — initialized by CONTINUED PRETRAINING of an open
> FM on the target organism, and FUSED with dna_decode's deterministic decoders as a hybrid-RESIDUAL** (the
> learned part carries only the polygenic signal the determinant scan provably misses — the ben-1 result names
> exactly this residual).

Component rationale:
- **SSM/Hyena backbone** — the only way to reach gene+regulatory context (10⁴–10⁶ bp) at solo compute.
- **Masked + AR** — MLM for representations, AR for generation/likelihood + zero-shot variant effect (Evo2 recipe).
- **CLIP head** — turns the representation into a *phenotype-aligned* space; enables zero-shot on new phenotypes.
- **JEPA auxiliary** — the actual "world model" bit: predict masked-region *embeddings* → in-silico perturbation
  ("flip this allele → predicted molecular-state shift").
- **Hybrid-residual fusion** — the honest niche; embeddings-alone lose, embedding-as-residual-over-mechanism is
  untried and targets the documented gap.

> **HONEST CAVEAT (integrity rail).** Every prior embedding test lost or was mixed. This architecture is a
> *hypothesis to be falsified*, not a plan of record. Its eval MUST be pre-registered (EP8-style): beat a
> SNP-PRS + kinship baseline AND a mechanism-feature baseline under de-confounded, clade-stratified CV, with a
> within-lineage diagnostic to rule out population structure. A clean fail closes the frontier honestly.

> **E9 — Multimodal fusion architecture.** Compare fusion designs for (sequence-FM ⊕ mechanism-features ⊕
> phenotype): late-fusion CLIP vs cross-attention vs residual-boosting. Which best captures the polygenic
> residual without re-introducing confounds? *Sources:* multimodal-bio fusion papers, DepMap multimodal.
> **E10 — "World model" semantics for genome.** Precisely define what a *world model* (vs a representation model)
> buys here: latent dynamics? counterfactual genotype→phenotype? in-silico saturation mutagenesis? Which are
> evaluable, and how? *Sources:* JEPA, Evo2 generative eval, in-silico mutagenesis (Enformer).

---

## 6. THE ELEMENT LIST — farm each to a different deep-research LLM (consolidated)

Each is self-contained. Format to hand off: **"Research question E#: <question>. Deliverable: <what a good
answer contains>. Check these sources: <list>. Flag uncertainty + cite."**

**Group A — Objective & paradigm:** E1 (objective taxonomy), E10 (world-model semantics), E3 (backbone).
**Group B — Architecture & scale:** E2 (SOTA currency), E5 (scaling law), E9 (fusion architecture).
**Group C — Data (highest priority):** E6 (de-confounded paired-label datasets ← #1), E4 (corpus manifest),
E10-data (granularity-availability matrix — which rung has depth for the target trait).
**Group D — Compute & training:** E7 (compute/$ per tier), E8 (continued-pretraining recipe).
**Group E — Evaluation & integration:** E11 (benchmark suite — BEND / Genomic Benchmarks / GUE / ProteinGym /
Enformer-tracks), E12 (hybrid-with-decoder fusion + the pre-registered falsifier), E13 (reuse dna_decode assets).

> **E11 — Benchmark suite.** The standard genomic-FM eval batteries (BEND, Genomic Benchmarks, GUE, Nucleotide
> Transformer tasks, ProteinGym for protein) + which the WM must win to claim SOTA.
> **E12 — Pre-registered falsifier.** Adapt EP8's falsifier to the WM: exact metrics (r²/spearman/within-group),
> baselines (PRS, kinship, mechanism-features), CV scheme, PASS/FAIL thresholds — *defined before any run*.
> **E13 — Asset reuse.** Map dna_decode's existing assets (foundation.py FMs, cached nt.h5/dnabert2.h5,
> substrates on D:, deterministic decoders, de-confound gate) to WM build steps — what's already done.

---

## 7. Decompose — build components (project families, if/when built)

| Family | Deliverable | Depends on | Gate |
|---|---|---|---|
| **G1 corpus** | assembled pretraining + paired-modality corpus (manifest + on-disk/cloud) | E4, E6 | storage (D:/cloud) |
| **G2 backbone** | chosen backbone + continued-pretraining checkpoint | G1, E3, E8 | **compute (T2 GPU)** |
| **G3 CLIP head** | contrastive alignment to a phenotype modality | G2, E6 | paired labels |
| **G4 JEPA/WM** | latent-prediction + in-silico perturbation | G2 | compute |
| **G5 hybrid fusion** | embedding-residual fused w/ mechanism decoders | G3, existing decoders | — (laptop) |
| **G6 eval/falsifier** | pre-registered benchmark + de-confounded verdict | G5, E11, E12 | — |
| **G7 compute** | provisioned GPU/cloud (VRAM check or budget) | — | **money/cross-machine** |

**Critical path:** G7 (compute) → G2 (backbone) → G3/G4 → G5 → G6. **G5 + G6 partially runnable NOW on frozen
cached embeddings (T3)** — the cheapest thesis test, no new compute.

---

## 8. --plan — phased build (each phase gated + falsifiable)

- **Phase 0 (NOW, laptop, $0):** this scoping brief + farm E1–E13 to parallel deep-research LLMs. Build the
  imputation path in parallel (user's stated concurrent track). *Deliverable:* this doc + returned research.
- **Phase 1 (T3, laptop/1-GPU, ~$0):** frozen-FM embeddings (cached nt.h5/dnabert2.h5) + a **CLIP head + hybrid-
  residual** on an existing de-confounded substrate (reuse decoders + de-confound gate). *Falsifier:* does the
  residual beat mechanism-features + PRS on the full metric surface? Cheapest test of the whole thesis.
- **Phase 2 (T2, workhorse/cloud, ~$100–5k):** continued-pretraining an open FM (Evo2/Caduceus) on the target
  organism; run the EP8-style pre-registered falsifier. *Gate:* workhorse RTX 3500 Ada VRAM check OR cloud budget
  (**money gate → user approval**).
- **Phase 3 (T1, ~$10–100k):** small from-scratch genome WM **only if** Phase 1–2 justify. *Gate:* budget/grant.

---

## 9. Existing assets to reuse (mix with what we already have)
- `dna_decode/models/foundation.py` — Evo / DNABERT-2 / Nucleotide Transformer / GENA-LM wired (+ Mock for tests).
- Cached embeddings: `D:/dna_decode_cache/embeddings/{nt.h5, dnabert2.h5, nt_n40_cipro.h5}` → Phase-1 is unblocked.
- Substrates on `D:/dna_decode_cache/`: **human** (opensnp, pgp, finngen), **fly** (dgrp), **protein** (proteingym
  DMS), **functional** (depmap, gdsc), **bacteria/TB** (cryptic, tb_indep, ena_wgs, refseq), **worm** (celegans_ben1),
  **fungal** (fungal_g1). A ready multi-kingdom paired-data pool for CLIP/eval.
- `plans/Hybrid_Learned_Deterministic_Decoder_Plan.md` (existing hybrid plan), the deterministic decoders, the
  **de-confound gate + CI-aware falsifier + within-lineage diagnostic** (the proven honesty infra — reuse verbatim).

---

## 10. One-line synthesis
A genome world model is buildable *in tiers*; the **data for self-supervised pretraining is not the bottleneck —
de-confounded paired phenotype labels and compute are.** The honest, cheapest first test (Phase 1, laptop, frozen
embeddings + CLIP + hybrid-residual) can be run now; everything heavier is gated on a compute decision. Farm
E1–E13 to parallel LLMs; the highest-priority element is **E6 (the de-confounded paired-label problem)** — solve
that and the WM has a winnable target; leave it unsolved and it repeats the project's 0-for-4.

---

## 11. Adversarial review integration (2026-07-04) — corrections, additions, caveats, verdict

> Produced by a genomic-FM adversarial critique pass on §1–§10. Corrections above are applied inline; the
> additions/caveats/elements below EXTEND the farm-list and the design rails.

**VERDICT on the thesis: MIXED — honest, not motivated reasoning, but leans optimistic.**
- **Residual-phenotype thesis = plausible but weakly justified.** After 0-for-4, the default prior is that generic
  embeddings learn population/phylogeny + sequence grammar, *not* de-confounded mechanism. Run Phase 1 ONLY as a
  cheap **falsification** against brutal baselines (E19) — not as a plan of record.
- **Generative/design thesis = MORE credible than phenotype-prediction** (esp. proteins / regulatory / variant-
  effect / constrained design) — but needs a DIFFERENT eval frame: *"does it propose sequences/edits that satisfy
  measured molecular constraints better than existing design tools?"* not *"predict phenotype from genome."* This
  may be the **higher-VOI first win**.
- **Do NOT bundle** MLM+AR+CLIP+JEPA+residual into one architecture before ONE component beats simple baselines on
  a clean substrate. Prove one axis first.

**Correctness fixes (applied inline):** NT ≈6 kb (1000×6-mer); DNABERT-2 up to ~10 kb (GUE); Evo2 ≈9.3 T nt /
~128k genomes / 1B-7B-40B / ≤1 Mb (distinct from Evo-1 ~300 B); Caduceus 24 GB only for constrained/fine-tune
(full pretraining multi-GPU); ESM-3 largest ~98 B multimodal (verify open vs closed); Chinchilla 20:1 = weak prior
for genomic corpora.

**Missing architectures to add to the E2 sweep:** **Gene42** (dense long-context attention + context-extension —
counters "SSM dominates"); **Basenji2 / Akita / DeepSEA** (supervised regulatory baselines — mandatory
comparators, not FMs); **graph / pangenome models** (variation & haplotype graphs — SNP-only linear windows miss
SV + reference bias); **discrete-diffusion** for sequence design; **perturbation-aware models**; **RAG /
database-conditioned** (CARD / ClinVar / MaveDB / RegulonDB — may beat a monolithic WM for mechanism-aware design).

**Missing data sources:** variant-effect (**MaveDB, CAGI, saturation genome editing**, ClinVar-w/-caveats);
regulatory (**MPRA / STARR-seq / CRE-seq, ENCODE cCREs, FANTOM, eQTL Catalogue, GTEx v10**); perturbation→state
(**Perturb-seq, CRISPRi/a, DepMap/Achilles/Score, LINCS/L1000, PRISM**); clean model-organism panels (**AraGWAS/
1001G, DGRP, CeNDR, BXD/CC/DO mice, maize NAM, rice 3K**); structural-variation / pangenome (**HPRC, gnomAD-SV,
1000G-SV, microbial pangenomes**).

**New farm-list elements (add to §6):**
- **E14 — Tokenization & coordinate representation.** single-nt vs k-mer vs BPE vs codon/gene vs byte vs graph;
  preserving SNPs / indels / strand-symmetry / coding-frame / repeats / SVs / coordinate metadata. (Its own axis,
  not a table column.)
- **E15 — Perturbation / causal substrate.** MPRA, Perturb-seq, CRISPR screens, saturation editing, DMS, segregant
  panels — ranked by causal interpretability × sample size. (The natural substrate for a *world*-model claim.)
- **E16 — Leakage / confounding audit.** Formal tests for homology, chromosome / species / population / lineage /
  batch splits, and public-DB train/test contamination. **MANDATORY given 0-for-4.**
- **E17 — Pangenome & structural-variation axis.** Representation for non-reference alleles, indels, mobile
  elements, CNV, plasmids, recombination / HGT, graph coordinates.
- **E18 — Design-tool validation & safety.** Success = likelihood vs naturalness vs edit-distance vs wet-lab vs
  constraint-satisfaction; off-target + **biosafety filters**; novelty **without pathogenic enhancement**.
- **E19 — Baseline gauntlet + negative controls.** BLAST/HMMER/motif, GWAS/PRS/kinship-LMM, Enformer/Borzoi/
  Basenji, tree/lineage, k-mer/logistic/CNN — PLUS shuffled-label + geography-only + nearest-neighbor-sequence +
  PCA-only + held-out-clade controls. Required before ANY WM value claim.

**Blind-spot design rails:**
- CLIP is not free — **InfoNCE will learn ancestry / batch / tissue / species / ascertainment** unless there are
  many INDEPENDENT de-confounded pairs.
- **"World model" ≠ masked-embedding prediction** — it must predict INTERVENTIONS / state-transitions under
  perturbation, else it's a representation model. (Sharpen E10.)
- Separate **sequence-plausibility** from **functional design** (+ safety) in any generative claim.

**Resolve-first open questions (yours — they set the whole build):**
1. First target trait: bacterial-AMR-residual / plant-flowering / worm-ben1-background / human-PRS / protein-DMS?
2. First win = phenotype-prediction, variant-effect-ranking, or generative-design?
3. Which paired dataset has enough INDEPENDENT samples after lineage/ancestry blocking? (E6/E16)
4. Reference-windows vs sample-haplotypes vs pangenome-graph? (E17)

