# The instrument question, re-answered after finding the probes were measuring the wrong geometry

**Status:** `DIAGNOSTIC`. Supersedes `wiki/glm_seed_role_decomposition_2026-11-21.md`, which is kept as the
record of what was believed and when. Artifacts, all read-only, CPU-only, no torch, no GPU, no network, no
D: volume:

| artifact | what it settles |
|---|---|
| `wiki/glm_grid_reconcile_2026-11-21.json` | the geometry gate, and a floor that can actually be subtracted |
| `wiki/glm_margin_sample_size_sweep_2026-11-21.json` | whether an unresolved margin is fixable by data |
| `wiki/glm_spread_vs_sampling_noise_2026-11-21.json` | what the across-seed spread actually is |
| `wiki/glm_instrument_levers_2026-11-21.json` | which feature block does the work, and how the null's order matters |
| `wiki/glm_margin_across_organisms_2026-11-21.json` | whether any of it survives a change of organism |
| `wiki/glm_seed_role_variance_v2_2026-11-21.json` | the corrected seed-role decomposition |

Reproduce with `uv run python scripts/glm_grid_reconcile.py` and its four siblings. **Use `uv run`, not
bare `python`:** `falsifier.heldout_auroc` catches `ImportError` and returns `None`, which is also its
"insufficient data" return, so on an interpreter without sklearn every measurement here silently reads as a
data problem. Every new probe asserts the dependency on entry.

## DECISIONS NEEDED

**All four of these are locked behind the same contract, and that is the single most useful thing in this
memo.** `dna_decode/glm/knob_contract.py:243` freezes nine `COMPARISON_GEOMETRY_FIELDS`, and three of them
are exactly the levers that turned out to matter — `feature_modes`, `markov_order` and
`n_generated_per_candidate`. Changing any one invalidates comparability with every committed anchor and
with the frozen acceptance bar. So none of what follows is a free fix; each is a priced fork, and the price
is a fresh Kaggle run on a fresh slug plus the loss of the anchor chain.

**One. Should the scoring surface become `positional` instead of `both`?** This replaces the old
"retire `kmer`?" question, which rested on a measurement that was wrong. Measured seed-paired at the
committed geometry: `positional` alone scores **+0.1740** above the floor, `both` scores **+0.1376**, and
the difference `both − positional` is **−0.0364 with a 95% interval of [−0.0804, 0.0077]** — not
resolvable, and pointing the wrong way for `both`. Meanwhile `positional − kmer` is **+0.1411**, resolvable,
positive at 9 of 10 seeds. So the 256 k-mer features are 73% of the feature vector and contribute nothing
measurable; the 96 positional features carry the whole signal. My recommendation is to move to
`positional`, because it is the same or better margin on a quarter of the features, which also relieves the
p-against-n problem (`both` fits 352 features on roughly 360 training rows). It is yours because
`feature_modes` is geometry-locked.

**Two. Strike plan Steps 4 and 6?** The recommendation stands and every reason I gave for it last time was
wrong. The honest reason now: **the across-seed spread is the sampling noise of a held-out AUROC, so a
seed-role split has nothing to remove.** A held-out AUROC between `n1` and `n2` items has a known null
sampling standard deviation, `sqrt((n1 + n2 + 1) / (12·n1·n2))`. At the committed geometry the
discriminator's test half holds about 300 natural against 60 generated, for which that formula gives
**0.0409**. The committed grid's own ten nulls spread **0.0410**. Measured at 240 replicates, the
all-coupled arm gives **0.0435** (kmer, CI [0.0394, 0.0474]) and **0.0448** (both, CI [0.0404, 0.0488]) —
both intervals contain the benchmark. Step 6 also needs a gated Kaggle session, so striking it still saves
a real cost.

**Three. Raise `n_generated` from 120 to 400?** This is the one genuinely new fork, and it is the cheapest
way to make the instrument more sensitive — but it is not free. At the committed 120 the `kmer` margin is
**+0.0329** and does not resolve; at 400 on the *same* corpus it is **+0.0535 with CI [0.0338, 0.0733]** and
does. `both` widens from **+0.1376** to **+0.1570**, and at the largest geometry tested (3600 natural, 600
generated) to **+0.2004**. The cost: `n_generated_per_candidate` is geometry-locked, so this breaks the
anchor chain, and the generator side needs about 3.3× the generation on a fresh slug. My read is that this
is **not urgent**, because the quantity an arm set is scored on is the generator's distance from the null,
and that distance is already large and registerable — the generator sits above its own null in **40 of 40**
committed cells, mean distance **+0.1924** in kmer and **+0.1355** in both. So I would hold this until a
candidate actually needs the extra resolution.

**Four. What the generator ceiling criterion should optimise — still open, with corrected numbers.** Against
a floor measured in the null's own geometry, the zero-shot generator sits **+0.2561** above the floor in
`both` mode, of which **+0.1376** merely reaches parity with a 3-mer Markov chain. (The figure I gave you
before, +0.2150, was its distance above the inflated ceiling; see CORRECTIONS.) Nothing has been written to
`soraya_runs/arms-ledger.json`.

## CORRECTIONS

Several of these correct corrections I made earlier today. That is the honest shape of it, and the earlier
memo is kept unedited so the sequence stays visible.

**The two probes were measuring a geometry no committed artifact uses, and nothing caught it.** Both
hardcoded `N_GENERATED = 400`; the committed sweep kernel froze **`n = 120`**
(`scripts/kaggle/glm_paired_sweep_kernel.py:45`). Every other constant matched, so both ran clean, produced
sensible numbers and passed their own tests. The variance probe's baseline arm carried a docstring saying it
*"reproduces the status quo"* and nothing in it ever compared against the status quo: at n=400 that arm
reports a spread of **0.0140** where the runner's own ten values give **0.0410**. Corrected to n=120 the
harness reproduces the grid at **20 of 20** seed-mode cells exactly. **The reusable rule: a probe that
re-implements a runner's measurement path must reproduce that runner's committed output, per cell, before
anything built on it is read.** An off-geometry probe does not fail — it answers a different question in the
same units.

**`margin = null − ceiling` was subtracting two numbers that cannot be subtracted.** The null is measured as
600 natural eval windows against 120 generated sequences; the ceiling splits those same 600 windows 300
against 300. Different class balance, different test-set size, and a held-out AUROC's sampling variation
depends on both. The replacement is a floor in the null's own geometry: the same eval half scored against
**120 real natural windows held out of both halves** — what a perfect generator would earn. It lands at
**0.4973** and **0.5007**, which is chance, exactly where real-against-real must sit. The max-of-five
self-control ceiling sat at **0.542** for the same comparison, so it was never a floor.

**"`kmer` mode has no winnable region" is withdrawn.** Against the like-for-like floor the margin is
**+0.0329**, not −0.0121 — opposite sign. Its interval is **[−0.0131, 0.0788]** with a standard deviation
1.95× the mean, so the honest verdict is **unresolved at this replicate count, not absent**; "positive at 7
of 10 seeds" reads as a finding while being consistent with a true margin of zero. The sign also flips on
the floor *rule* alone: against the mean of the same five self-control values rather than their maximum, the
committed margin is **+0.0274**.

**And the small `kmer` margin turns out to be the correct answer rather than a defect.** A Markov chain's
sufficient statistic *is* the k-mer conditional distribution, so an order-3 chain saturates a 4-mer
discriminator by construction. Measured across null orders 1 to 5, the `kmer` margin runs **+0.3259,
+0.1675, +0.0329, +0.0175, +0.0269** — a range of **0.3084** — while `positional` runs **+0.1850, +0.1547,
+0.1740, +0.1828, +0.1506**, a range of **0.0344** with *no* order resolvably different from the committed
3. The k-mer discriminator and the Markov null are the same class of object; the positional discriminator
measures where composition sits, which no order of Markov chain models. `positional_features`' own docstring
argued this from first principles and nobody had measured it.

**My stated mechanism for the seed coupling is refuted by the arm built to test it.** I said the roles
cancel because the Markov null is fitted on the fit half of the very corpus draw it is scored against. That
is true of the runner — and it was equally true in *every* arm of the v1 probe, because each arm drew fit
and eval from one corpus call, so the matching was constant across arms and could not explain a difference
between them. The v1 probe could not have measured its own stated mechanism. Arm F breaks the matching
directly and returns an SD ratio of **1.129 with CI [0.891, 1.445]** — contains 1.

**And my replacement mechanism does not replicate either.** I then reported that stream sharing suppresses
the spread, on an arm that gave ratio **1.515, CI [1.191, 1.912]**. Two things undo it: the arm it was
compared against was a low draw at 60 replicates (**0.0352** against **0.0448** at 240), and the inflation
appears at **one of five** offset pairs — the pair that run happened to use — and at **none** in kmer mode.
Verdicts: `NO_STREAM_EFFECT__V2_RESULT_WAS_ONE_OFFSET_CHOICE` (kmer) and
`STREAM_EFFECT_IS_OFFSET_DEPENDENT__NOT_A_MECHANISM` (both).

**"Pairing helps model-vs-null" was stated more strongly than measured.** Across all ten checkpoint/null
pairs in the committed grid, a resolvable paired-delta benefit appears in **2 of 4** model-vs-null pairs in
`both` mode (median ratio 0.684) and **0 of 4** in kmer, against **1 of 6** checkpoint-vs-checkpoint pairs
(median ratio 0.928). The strongest case is real — zero-shot minus null in `both`, paired SD **0.0275**
against an unpaired **0.0643**, correlation **+0.832** with CI **[0.424, 0.959]**. But four and six
overlapping pairs drawn from five shared series cannot carry a clean rule, and the negative correlations I
noticed have intervals spanning zero, so they are **not** evidence of anti-correlation.

**The grid I read these nulls from declares its own verdict `INSTRUMENT_INVALID`.** For an unrelated reason
— a per-cell null gate that could only ever pass at seed 7, corrected in
`scripts/glm_paired_sweep_verdict.py`. The per-cell scores are not in dispute and I reproduce all twenty of
them exactly; the verdict label is. Worth stating because I built on that grid without mentioning it.

## What the instrument can and cannot do

**The winnable region survives a change of organism in `both` mode, and this needed no download.**
`data/refseq_cache/` holds 459 assemblies carrying both a genome and an annotation — 400 *Campylobacter
jejuni*, 26 *C. coli*, 17 *Klebsiella pneumoniae*, 16 *E. coli* — spanning GC **0.266 to 0.505**.
Representatives were chosen by a rule fixed in advance and applied blind to every number
(lexicographically smallest accession); with 400 candidates, a results-aware rule could have produced any
margin wanted.

| substrate | GC | kmer margin | resolves | both margin | resolves |
|---|---|---|---|---|---|
| *E. coli* K-12 MG1655 (reference) | 0.453 | +0.0329 | no | +0.1376 | yes |
| *E. coli* (draft, GCA_059754755.1) | 0.457 | +0.0374 | yes | +0.1450 | yes |
| *K. pneumoniae* (GCA_061544025.1) | 0.505 | +0.0232 | no | +0.1396 | yes |
| *C. jejuni* (GCA_060772745.1) | 0.266 | +0.0291 | no | +0.1709 | yes |
| *C. coli* (GCA_060776685.1) | 0.272 | +0.0533 | yes | +0.1871 | yes |

`both` resolves at **5 of 5**. `kmer` is **positive at all five** and resolves at two, which reproduces the
sample-size finding by an independent route. Two controls are built into the substrate list rather than
argued for: the reference and the draft *E. coli* agree to **0.0074** with overlapping intervals, so
assembly quality and the clinical-isolate sampling frame are not driving the result; and the two
*Campylobacter* congeners agree, so this is a genus differing rather than one odd genome.

**The margin-against-GC correlation is not a GC finding and is reported as such.** It comes out at
**−0.941** in `both` mode, but the two low-GC substrates are both *Campylobacter* and the three high-GC ones
are all *Enterobacteriaceae*, so it is a two-group comparison with an **effective n of 2** dressed as five
points. Reading it as a GC effect would repeat the lineage confound that has already closed four arms in
this project. Separating them needs a high-GC *Campylobacter* or a low-GC enterobacterium, which this set
does not contain.

**The instrument has power in both modes, so a floor at chance is not a blind discriminator.** Against the
uniform i.i.d. control it reaches **0.9295** (kmer) and **0.9473** (both). That bounds the ceiling of
detectability, not the target a generator must beat — i.i.d. sequence is far more different from natural
than any real generator, so high power against it does not imply power against a near-miss.

## Two smaller things found on the way

**The runner's disjointness promise is true only up to duplicate windows.** It states that the null is
fitted on the fit half and scored on the eval half "so no generator is scored against sequences it was
trained on". The *E. coli* promoter corpus holds **4340 windows over 4323 distinct sequences** — 7 sequences
repeated across 24 windows, one of them 6 times, because distinct genes can share an identical upstream
region (operons, duplicated genes, repeated elements). At the committed geometry **0 to 3 of 600** eval
windows also appear in the fit half. Measured, recorded in proportion, and far too small to move a score —
but the unqualified claim is not exactly true, and it is stamped into every artifact. This is also why the
held-out floor filters spares by **sequence**, not by index: one leaked duplicate puts the same feature
vector in both classes and makes the floor read cleaner than it is.

**A defect in my own comparison helper, found by a test and fixed.** `paired_diff` guarded resolvability
with `bool(se and ...)`, so a difference with *zero* variance — the most resolvable case there is — was
reported as not resolvable. It survived its own unit test because that test built the second arm by
subtracting a constant, leaving floating-point residue that made the standard error tiny rather than
exactly zero. Real data never hits the branch, so no reported number moved.

## Honest limits

- Everything here is measured on the **Markov null**, not the HF generator. The generator's own paired
  delta is read from the committed grid at no compute cost and is reported, but its sampling variance
  cannot be decomposed without re-running it on a GPU.
- The Mann-Whitney benchmark is the null variance for **exchangeable classes**. It applies exactly to the
  perfect-generator floor, whose true AUROC is 0.5 by construction, and only approximately to the Markov
  null. It also assumes a fixed scoring function, while the classifier is refitted per split, so it is a
  **lower bound** on expected spread rather than an exact prediction.
- Arm A2 was tested at five offset pairs. Offsets are an unbounded space and five is not a proof that none
  behaves differently.
- One representative genome per organism, so within-organism variation is unmeasured except for *E. coli*.
  The non-reference genomes are **draft** assemblies from the AMR cohorts.
- A high-order Markov chain fitted on 600 windows of 150 bp approaches memorising them, so part of the
  order effect is overfitting the fit half rather than modelling the substrate. The sweep measures the
  effect and does **not** separate those causes.
- The sample-size sweep moves fit size, train size, test size and class balance together; it shows the
  direction the margin moves with data and does not identify a sample size at which a mode becomes usable.
- Ten seeds per cell almost everywhere, so each interval is itself imprecise. Arm A is the exception at 240.
- **A winnable region is necessary, not sufficient.** All three feature sets see LOCAL composition. The
  positional block adds *where* composition sits, not long-range organisation, which is what the published
  audits name as the real failure of generated sequence. Clearing this instrument would not establish
  functional realism.
- One window length, one null family, one organism per row.

## What I did not do

I have **not** edited the plan, changed any geometry field, written to the arms ledger, or run anything on
a GPU. All four decisions above are scope or authority calls. The peer session has been sent a correction
withdrawing the figures I gave it earlier.
