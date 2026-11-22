# What this session tested, what is still open, and what is worth doing next

**Status:** `LEDGER`. A record of a single day's work on the GLM generator arm, written so a fresh context
can pick it up without re-deriving any of it. Every number here is in a committed artifact and is cited;
nothing is quoted from memory. Commits `44fe325c` → `d7bf80ba` on `promoteratlas-exposure-2026-11-01`.

## DECISIONS NEEDED

All four change a field in `COMPARISON_GEOMETRY_FIELDS` (`dna_decode/glm/knob_contract.py:243`), so each
costs comparability with every committed anchor plus a fresh Kaggle run on a fresh slug. That shared cost
is the reason they are one decision wearing four hats.

**One. Move the scoring surface from `both` to `positional`?** Measured seed-paired at the committed
geometry: `positional` alone scores **+0.1740** above the floor against `both`'s **+0.1376**, and
`both − positional` is **−0.0364, CI [−0.0804, 0.0077]** — not resolvable, and pointing the wrong way for
`both`. `positional − kmer` is **+0.1411**, resolvable, 9/10 seeds. So the 256 k-mer features are 73% of
the vector and contribute nothing measurable.

**RECOMMENDATION REVERSED 2026-11-21, by running the scorer rather than reasoning about it. Keep
`both`.** Positional-ONLY mode **fails the committed falsifier's own power gate**: its uniform control
measures **0.8611** over 5 seeds against `MIN_POWER_AUROC = 0.90`, so `run_falsifier` returns
`INSTRUMENT_UNDERPOWERED` and emits **no number at all** for that mode. I discovered this when the first
real scoring run printed `kmer` and `both` rows and silently omitted `positional`.

The mechanism is worth keeping, because it is the same shape as the k-mer saturation finding inverted:
the power gate is calibrated on UNIFORM sequence, which is an easy target that a 256-dimensional k-mer
histogram crushes (kmer 0.9317, both 0.9494). Positional features are **weaker against uniform and
stronger against the Markov null** — the target that actually matters. So the gate rejects the mode that
is best at the hard task because it is worse at the easy one.

**`both` is therefore the practical route to positional signal**, since it passes the gate at 0.9494 and
carries the positional block. The case for positional-alone was always feature economy rather than score
— its advantage over `both` was never resolvable (−0.0364, interval crossing zero) — and that does not
justify a mode the committed scorer refuses to score. Whether `MIN_POWER_AUROC` or the uniform control is
the right calibration for positional features is an instrument question and remains open; it is not
settled by preferring the answer that is convenient.

**Two. Strike plan Steps 4 and 6?** Yes, and for a reason unrelated to the one originally given. The
across-seed spread is the sampling noise of a held-out AUROC, so a seed-role split has nothing to remove.
Step 6 also needs a gated Kaggle session, so striking it saves real cost.

**Three. Raise `n_generated` from 120 to 400?** It makes the kmer margin resolvable (**+0.0535, CI
[0.0338, 0.0733]** against **+0.0329, CI [−0.0131, 0.0788]**) on the *same* corpus. Not urgent: the
generator already sits **+0.1924** above its own null in kmer mode across 40 of 40 committed cells, so
resolution is not the binding constraint. Recommendation: hold.

**Four. What the generator ceiling criterion should optimise.** Unchanged and untouched. Against a
like-for-like floor the zero-shot generator sits **+0.2561** above it in `both` mode, of which **+0.1376**
merely reaches Markov parity. Nothing has been written to `soraya_runs/arms-ledger.json`.

**A fifth, newly created by the scope correction and owned by the user:** does the generator emit
**promoter windows** or **any genomic window**? The download differs by ~300× and only promoter mode is
scoreable by the committed falsifier. `wiki/glm_corpus_spec_all_domains_2026-11-21.json`.

## WHAT WAS TESTED, AND THE ANSWER

Everything in this section is measured, with the artifact that holds it.

| # | Question | Answer | Artifact |
|---|---|---|---|
| 1 | Do the seed-role probes reproduce the runner they extend? | **No** — both ran at `n_generated=400` against the kernel's frozen `n=120`. At n=120 they reproduce it 20/20 | `glm_grid_reconcile` |
| 2 | Is `margin = null − ceiling` a valid subtraction? | **No** — 600-vs-120 against 300-vs-300. Replaced by a floor in the null's own geometry, which lands at **0.4973 / 0.5007** = chance | `glm_grid_reconcile` |
| 3 | Does `kmer` mode have a winnable region? | **Yes, +0.0329**, not −0.0121. Unresolved at the committed n, not absent | `glm_grid_reconcile` |
| 4 | Is that fixable by data? | **Yes** — resolvable at 4 of 5 larger geometries; n 120→400 on the same corpus suffices | `glm_margin_sample_size_sweep` |
| 5 | Is the coupling mechanism fit-set matching? | **No** — the arm built to test it returns ratio **1.129, CI [0.891, 1.445]**, containing 1 | `glm_seed_role_variance_v2` |
| 6 | Is it RNG stream sharing? | **No** — appears at 1 of 5 offset pairs, against an arm A that was a low draw at n=60 (0.0352 vs **0.0448** at n=240) | `glm_spread_vs_sampling_noise` |
| 7 | Then what IS the spread? | **Sampling noise.** Mann-Whitney null sd predicts **0.0409**; the grid's own nulls spread **0.0410**; arm A at 240 replicates gives 0.0435/0.0448, both intervals containing it | `glm_spread_vs_sampling_noise` |
| 8 | Does `positional` alone match `both`? | **Yes.** `both − positional` = −0.0364, interval contains zero | `glm_instrument_levers` |
| 9 | Does the null's ORDER matter? | **Decisively for kmer, not at all for positional.** kmer range **0.3084** across orders 1–5; positional **0.0344** with no order resolvably different from 3 | `glm_instrument_levers` |
| 10 | Is that structural or an E. coli quirk? | **Structural** — replicates at GC 0.266 (ratio 6.53× vs E. coli's ~9×) | `glm_instrument_levers` |
| 11 | Does the winnable region survive a change of organism? | **`both` yes, 5/5** (+0.1376 to +0.1871, GC 0.266–0.505). `kmer` positive at 5/5 but resolvable at 2 | `glm_margin_across_organisms` |
| 12 | Is the margin driven by GC or lineage? | **Confounded, and the raw correlation REVERSES under lineage control** (−0.669 → +0.294). Suggestive only: 3 of 10 phyla have replicates | `glm_margin_gc_vs_lineage` |
| 13 | How much distinct sequence is cached locally? | **30.9 Mbp** from 0.92 Gbp raw — **63× redundant**, 400 of 459 one species | `glm_corpus_distinct_content` |
| 14 | Can a 50M generator be pretrained on it? | **No — 32× short.** E. coli alone is 215× short | `glm_corpus_distinct_content` |
| 15 | Does the organism-verification gate catch anything real? | **Yes** — 15/16 verified, 1 flagged (*Caulobacter vibrioides* vs header *crescentus*), adjudicated accept on matching strain NA1000 | `glm_corpus_fetch` |
| 16 | Does the leak gate catch real leaks? | **Yes, twice.** A single-seed gate would leave **2,749 of 3,349 (82%)** scoring windows in training; and **1,389 windows** were withheld from genomes that are not the reference, via byte-identical promoters in relatives | `glm_training_corpus_manifest` |
| 17 | Does the diffusion training loop actually learn? | **Yes** — recovers ≥90% of a sequence from a fully masked input; untrained control does not | `tests/test_glm_diffusion.py` |

**Reusable instrument facts that came out of this.** A Markov chain's sufficient statistic **is** the
k-mer conditional distribution, so a k-mer discriminator is saturated by an order-k−1 chain **by
construction** — the small kmer margin is the correct answer, not a defect. The positional margin
measures something no Markov order models, which is why it is order-invariant and why it is the surface
worth scoring on. `positional_features`' own docstring argued this from first principles; nobody had
measured the order dependence that tests it.

## WHAT IS STILL OPEN

**In flight.** The 20,000-step training run on the 9.59 Mbp corpus. Checkpointed every 500 steps,
resumable.

**Never measured, and the most important gap in my own reasoning.** I chose diffusion over autoregressive
for a measured reason — positional architecture is what matters here — but I did **not** measure that
diffusion beats autoregressive on this substrate. The architecture choice is REASONED, not tested. An AR
baseline on the identical corpus, scored by the identical falsifier, is the honest comparison and it has
not been run.

**No generator number exists yet** on the corrected instrument. Everything in the table above bounds the
INSTRUMENT. The first real generator score is the falsifier run against the trained model, which is minutes
of work once training ends and has not happened.

**The GC question needs a different corpus, not a bigger one.** With one genome per phylum, phylum predicts
GC perfectly by construction (collinearity 0.983 across 10 phyla). It needs 3–4 species *within* each
phylum spanning the GC range — a different selection rule from the one used, which was correctly optimised
for corpus diversity instead.

**Cross-domain is blocked on the fork decision**, and no eukaryotic assembly exists locally to measure
(`data/arabidopsis/` holds phenotype tables, not sequence). Anything beyond bacteria must land on Server1
(767 GB free) rather than this host (**19 GB free, 97% full** — measured, against the 66 GB the project
notes still record from October).

**The HF generator's own sampling variance** cannot be decomposed without a GPU. What IS free from the
committed grid, and was done: its paired-delta behaviour, which shows a resolvable pairing benefit in only
**2 of 4** model-vs-null pairs in `both` mode and **1 of 6** checkpoint pairs — suggestive, not a rule.

**Two species are unmeasurable at the committed geometry** — *Borreliella burgdorferi* and *Chlamydia
trachomatis*, reduced-genome pathogens with fewer than the 1,320 windows needed. They remain in the
training corpus, which has no such floor.

## NEXT HIGHEST VOI, RANKED

1. **Score the trained generator with the falsifier.** Minutes of work, and it is the north-star number
   this entire arm exists to produce. Nothing else here is a generator result. Do this first.
2. **Autoregressive baseline on the identical corpus.** Closes the gap named above: it converts an
   architecture choice from reasoned to measured, and it is the actual Axis-B question. Same corpus, same
   scorer, same geometry — the only variable is the generative mechanism.
3. **The authority decisions**, which cost me nothing and gate everything downstream. Decision One
   (`positional`) has the strongest measured support of the four.
4. **Within-phylum GC replicates.** ~12 more genomes closes a limit declared twice today and makes a
   reversed correlation interpretable rather than merely withdrawn. Off the generator's critical path.
5. **Cross-domain corpus**, gated on the fork. At the measured promoter-mode yield of ~0.53 Mbp/genome,
   a 7M-param model needs hundreds of species — which is the cost to keep in front of the decision.

Deliberately NOT on this list: more seed-role work (answered), more GC phyla (the wrong axis — it needs
replicates, not breadth), and anything requiring Server1's GPU (the corpus, not compute, is the binding
constraint, and that stays true until the fork is settled).

## AMENDED 2026-11-21, after a field-literature sweep and a design review

Four changes. Two correct numbers I published, one corrects a claim I made about the field, and one is
a defect in the training design that the sweep led me to find.

**1. The scored run trains 95% on organisms it is never scored against. CRITICAL.** The falsifier scores
only against *E. coli* MG1655 eval windows, and *E. coli* is **3,463 of 63,925 corpus windows = 5.4%**.
The corpus was built for lineage diversity, which is right for a PRETRAINING corpus and wrong for a
SCORED one, and I did not notice the mismatch when I pointed the scorer at it. Left alone it yields a
poor score for a reason unrelated to the architecture, which is one step from being misattributed to
diffusion. **Fix costs nothing and is half-done:** the pooled run IS stage one of a pretrain-broad then
fine-tune-per-organism recipe. Finish it, fine-tune on *E. coli*'s own windows, then score. The leak gate
is unaffected — fine-tuning on the fit half is legitimate and the eval windows stay withheld by sequence.

**2. The `20 tokens/param` heuristic is wrong for this domain by roughly an order of magnitude.** I used
it throughout while recording in every artifact that it is imported from language-model scaling and NOT
validated for DNA. A published prokaryotic promoter generator sits at **~1.74 tok/param** (7.4M params on
~12.9M tokens) and produced in-vivo-functional sequence. Recomputed at that density our 9.59 Mbp supports
**~5.5M params, not 0.479M**, and our 813,700-param model at 11.8 tok/param is conservatively sized
rather than starved. **Every "N-times short" figure in this ledger and in the all-domains spec is soft by
about 10x.** The direction holds — promoter corpora are small, cross-domain needs many species — the
magnitudes do not. Marked here rather than silently restated.

**3. "Our positional instrument is ahead of the field's norm" is WITHDRAWN.** I generalised from one
paper's 6-mer validation metric to the field. The specialist promoter-generation work does score
positionally: two bacterial-promoter diffusion papers measure the -10/-35 inter-motif spacing. Corrected
claim: we are ALIGNED with the best promoter work on the canonical spacing, and broader only in that our
instrument measures positional composition across the whole window rather than one known motif pair.

**4. The promoter-vs-whole-genome fork is RESOLVED to promoter mode, and it costs us nothing** because
the corpus was already built that way (`natural_set_arg: "promoter"`, 63,925 x 150 bp). Recorded because
a peer message described our 9.59 Mbp as whole-genome and proposed building the extraction we already
have — the two figures are the same object, and the confusion is worth pinning so it does not recur.

**What this does NOT change.** Every instrument result stands: k-mer saturation is structural and
replicated at a second GC extreme, the positional margin resolves at 14/14 species, the leak gate caught
1,389 real leaks, the overfit test proves the loop learns, and the morning's reconciliation of two
off-geometry probes is untouched.

**One new option, with a practical obstacle.** An open, MIT-licensed, CPU-capable pretrained bacterial
promoter diffusion model exists and uses rotary positional encoding plus inpainting-based targeted
design. It generates **50 bp**; our falsifier scores **150 bp**, and `window_length_bp` is one of the
nine `COMPARISON_GEOMETRY_FIELDS` — so it is not a drop-in scoreable baseline. It is a third arm with its
own geometry, not a resolution of the autoregressive-versus-diffusion question, which still needs matched
size, matched corpus and matched scorer.

**Revised next-VOI order:** (1) finish the pretrain, fine-tune on *E. coli*, score — the design fix and
the north-star number in one move; (2) the autoregressive baseline at matched size, which now has a
published 7.4M character-level template to match rather than a strawman; (3) the authority decisions;
(4) within-phylum GC replicates; (5) evaluate the pretrained external model as a separate arm.

**Provenance:** the literature findings here are a colleague's summaries of papers I have NOT read. The
architecture templates and the transfer-learning reframe are actionable regardless, but the specific
figures are second-hand and unverified in this session — the same standard applied to the GPU reading
earlier today.

## PROCESS LESSONS, EARNED NOT ASSERTED

**I asserted a mechanism, had it refuted, asserted a replacement, had that refuted too — and only then
checked the closed-form null that made both unnecessary.** The Mann-Whitney sd took two minutes and
matched the observed spread to three decimals, twice, across a geometry change. A mechanism story is cheap
and fits almost any direction; a closed-form null is a fixed target the data either matches or does not.
Check it first. Recorded as `feedback_compare_a_spread_against_its_closed_form_null_first`.

**A guard must be shown to FIRE, and the cheapest proof is that it rejects the data that motivated
building it.** The collinearity guard rejects the four-species case it was built for. The leak gate
withheld 1,389 real windows. The overfit test has an untrained control. Conversely: a verification run
that never reaches its guard reports a false pass — my `--smoke` flag silently overwrote the setting I was
using to create a config mismatch, so the guard had nothing to detect and I nearly read that as "did not
fire".

**Do not pin a test to the mutable state of a file you do not own.** A test asserting the peer's spec was
refused for having null accessions went red an hour later because the peer resolved them — correct
behaviour producing a red that is indistinguishable from a real break.

**Never train on what you grade on, and the withhold must be the union over evaluation seeds.** Measured:
a single-seed gate leaves 82% of the scoring windows in training. Applied by SEQUENCE across every genome,
not just the reference — relatives carry byte-identical promoters, and 1,389 windows were caught that way
from genomes that are not the reference at all.

**Raw byte totals are not corpus size.** 0.92 Gbp of cached assemblies is 30.9 Mbp of distinct sequence,
because the directory was built for a different purpose where 400 strains of one species was the point.

**A constraint can quietly become the goal.** My corpus planning drifted to bacteria-only because the
cached genomes are bacterial. The user caught it. The project's own north star already said bacteria is a
stage.

**And the one the user had to say four times:** peer messages carry their authority — act on them, do not
re-ask for a go-ahead a peer already gave. I stretched an action-class carve-out to cover ordinary local
compute and bounced a cleared decision back repeatedly. The carve-out is now a closed list of five.
Recorded in `feedback_always_pass_peer_messages_through_to_claude`.
