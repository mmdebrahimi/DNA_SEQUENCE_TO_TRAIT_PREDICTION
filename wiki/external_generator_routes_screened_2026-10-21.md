# Two external generator routes screened and CLOSED — AlphaProtein Novo, AlphaGenome Atlas

**2026-10-21.** Both were screened because our measured bottleneck is the GENERATOR, not the scorer.
Both close, for different and structural reasons. **Do not re-propose either as a generator route for
bacterial regulatory sequence.**

> **PROVENANCE, read before citing anything here.** This screen was carried out by a sibling session
> (`scientific-research`) and **relayed** — nothing below was fetched or re-derived in this repo. The
> relay marked every answer PRIMARY / SECONDARY / NOT FOUND, and that marking is preserved verbatim
> below because it is the whole value of the record. **This memo is deliberately NOT wired into
> `wiki/negative_results_map_2026-06-13.md`, `eval/regime.py` or `cell_regime.py`.** Those are
> load-bearing gates, and a relayed claim must not silently become one. If either direction is ever
> re-proposed, re-derive A1 or B1 from the primary source first — each is a single cheap check.

## DO NOT RE-PROPOSE — the two closing facts

### A1. AlphaProtein Novo is hard structure-dependent, so it cannot reach our substrate *(PRIMARY)*

The generator's required input is a **structure file**: the README's input spec is a path to a CIF file
containing the catalytic motif, and the package co-generates protein structure and sequence
*conditioned on that motif*. The preprint states the generator has nothing else — the design process
has only the input structural catalytic motif plus a ligand pose.

**Why this is not an AlphaGenome-shaped situation.** Our recorded AlphaGenome finding is that the
*architecture* did not transfer but the *training method* did. Here **both are structure-bound**: the
architecture is a 3D structure+sequence co-diffusion model and the method *is* motif scaffolding in 3D
space. There is no method/architecture seam to exploit. A 100–300 bp bacterial promoter has no
catalytic motif, no ligand pose and no fold to scaffold, so the direction is closed by a property of
the method rather than by a budget or a dataset. **This will not expire.**

### B1. The AlphaGenome Atlas is human-hg38-only *(PRIMARY)*

It scores all ~9 billion SNVs in hg38; the underlying model was trained on **mouse and human DNA
only**, and the official FAQ states prediction quality is expected to degrade with evolutionary
distance and that this has not been formally benchmarked. No bacterial or prokaryotic coverage, and no
prokaryotic equivalent surfaced. **Irrelevant to the bacterial work.** It touches only the human
clinical cells (`pgx` / `clinvar` / `hla`).

### B2. AVI cannot be a label under our own G1 — at two layers *(PRIMARY)*

- **Features:** AVI integrates AlphaGenome and AlphaMissense predictions — both model outputs — plus
  two conservation measures. 18 input features, against >150 for CADD.
- **Labels:** it was trained to discriminate variants with filtering allele frequency above vs below
  0.1%. That is a **population-genetics proxy, not an assay readout**.

A prediction built from predictions, supervised by a frequency proxy. It can be a **comparator** or an
**input feature**; it can never be ground truth. This is an instance of the existing G1 gate, not a new
gate.

## THE ONE THING WORTH CARRYING ACROSS, AND IT IS NEGATIVE *(PRIMARY)*

**A well-resourced lab shipped state-of-the-art results with its generator grading its own proposals,
and its filters tuned on the data it then reports against.** Per the relay: the AP Novo diffusion model
can itself be used for the filtering step, the filtering uses same-group AF3 predictions, and the
filters were iteratively improved over multiple rounds of design against the same experimental results.
No guard against a model grading its own designs was found anywhere.

This is **external corroboration of a rail we already enforce** — our hard rule that a generating model
must never grade its own proposals, which exists because the forward/inverse falsifier had to be graded
on measured wet-lab DMS with calibration and selection on *disjoint position splits*. Consequence for
reading their numbers: the headline hit rates are **post-hoc-filtered**, not a clean
propose-then-independently-score result.

## Calibration for our own propose-k / assay-k rails *(PRIMARY)*

Hit rates on hard chemistry, after filtering, as reported: ~4.9% (81 designs, 4 active ≥3 SD over
background) · ~5.6% (90 designs, 5 over threshold) · ~1.6% at a stricter kinetic bar (1,536 designs, 24
reaching kcat/KM > 1,000 M⁻¹s⁻¹) · ~4.8% (188 variants, 9 over threshold). Validation is expression in
*E. coli*, purification, then activity and kinetics.

So the usable planning figure is **order-10¹ proposals per hit after filtering, not order-1** — which is
the right order of magnitude for our own propose-k / assay-k / keep-best discipline.

**Their own stated caveat, preserved:** hit-rate definitions are reaction-specific and are *not*
intended to be compared across reactions. So do not average these, and do not carry one of them into a
different context as "the" hit rate.

## What was NOT established — recorded so it is not mistaken for settled

- **The widely-quoted AVI benchmark figures (top-50 in 29.5% of cases vs CADD 12.5%) are
  SECONDARY ONLY.** Not present in the preprint abstract, and not in the first 100k characters of the
  full text read. Two outlets report them consistently as the paper's figures, so they are probably
  right — **but they are not promoted here, and must not be quoted as primary.**
  What *is* primary is the benchmark DESIGN: **n = 112** likely-pathogenic variants from solved cases,
  comparator **CADD v1.7**, and the ranking is **intra-patient** (each case's causal variant ranked
  against that same patient's background variants, not a global pool). Status: medRxiv preprint, not
  peer-reviewed.
- **AP Novo train/test split construction: NOT FOUND.** The structural-novelty cutoffs were not
  retrieved; TM-align, Foldseek and the PDB are cited, so a novelty search was clearly run, but the
  thresholds live in the Methods of the full PDF and were deliberately not guessed.
- **AVI's own split / training-data description: NOT FOUND.**
- A saturation-genome-editing comparison (AVI 0.668 vs CADD 0.647 across ten unseen screens) is
  **secondary, single-source, small margin, metric unlabelled.** Not usable.

## Two license facts that are decisions, not findings

1. **AP Novo weights are separate from the code and non-commercial.** Code is Apache 2.0; the weights
   are not in the repo and carry their own terms, free for non-commercial organisations only. The clause
   to read: commercial use of the parameters *or their outputs* is prohibited, **and outputs may not be
   used to train, fine-tune or distil models for de-novo protein or enzyme design.** Our substrate is
   DNA rather than protein/enzyme design, so it plausibly does not bind us — **the relay explicitly
   declined to rely on its own reading, and so does this memo.** A USER/legal call if the direction were
   ever reopened, which A1 says it should not be.
2. **Whether AVI scores may be used as model FEATURES is UNRESOLVED.** The FAQ says Atlas outputs are
   non-commercial and should not be used to train other ML models, but points to an unenumerated
   "Permissive Use Downloadable Artifacts" exception in the Terms of Service. Secondary coverage claims
   the static AVI SNV download *is* that permissive artifact; the relay could not confirm it (the ToS
   page returned navigation chrome on two fetches). **This is the single thing to check before building
   anything on AVI**, and it only matters for the human clinical cells.

## Net

**A closes on A1 alone** — structure-conditioned generation cannot reach an unstructured 100–300 bp
promoter, and there is no seam. **B closes on B1 for bacteria and on B2 for labels.** The only live use
anywhere is AVI-as-a-feature on the human cells, gated on the unresolved licence question.

Three further items from the same source were screened out before research and are recorded in
`NEXT.md` rather than here: the embedded-ribonucleotides paper (off-substrate, AI-free), the Anthropic
ART / "CRISPR-like" system (explicitly unvalidated — enzyme not shown active, no cas genes), and the
virtual-cell coalition (a watch item; first data ~a year out with a commercial head start, and a free
*measured* phenotype dataset is the only thing in it that would touch our label wall).
