# Primary-source reading notes — and the four citation errors the reading found in my own research pass

<!-- doc-status: LIVE as of 2026-11-18. Discharges the honest-limit stated in
     wiki/generator_hypotheses_ranked_researched_2026-11-18.md: "search-result level, not full-text;
     every quoted figure should be re-checked against the primary source before it anchors a bar."
     It is now re-checked. Four citations were wrong. One rank moves; the top recommendation does not. -->

## READ THIS FIRST — the method finding matters more than any single correction

**At least three times in the earlier research pass I quoted a SEARCH-RESULT SUMMARY as if it were
paper text.** The summaries were accurate summaries; the attributions were not. The quotes I attached
to named papers were in two cases not in those papers at all, and in one case were about a different
quantity than the one I used them to support.

This is a distinct failure from the two the ranked memo already records (asserting things were unbuilt
that were built). That one was fixed by reading the repo. **This one is only visible by reading the
source** — and it is invisible to every automated check the repo has, because a citation can resolve,
parse, and contain the quoted number while still being about the wrong thing.

## Corrections, ordered by what they cost

### C1 — H2's motivational base largely EVAPORATES (the most consequential)

The ranked memo said the distance-decay phenomenon is **established**, citing: *"the influence of a
mutant protein against its wild type decreases as a function of its edit distance to the wild type."*

The quote is verbatim and correct. **It is about training-data INFLUENCE, not accuracy.** Influence
here is computed by Kronfluence — how much a training point moves a prediction — which is a different
quantity from predictive performance. Checked against the source: the paper *"never plots or tabulates
zero-shot performance against edit distance"*, reports no correlation statistic for the trend, fits no
functional form, and its own text/figure disagree on whether the axis is edit or Hamming distance.

The second support fares no better than claimed but differently than claimed. The K-pattern quote
(*"at K = 2, the GBT model can perform reasonable positional extrapolation to 3 mutational regimes"*,
and complete failure at K=5) **is verbatim and real**, but: it is on **simulated NK landscapes**, not
real protein data; the metric is **Pearson r, not Spearman**; and the paper explicitly *"does not
measure distance from the nearest training sequence"* — it counts regimes from a seed. It also notes
all models fail at K=5, not only the tree model, and never defines a threshold for "fails".

**Net effect, and it cuts both ways.** The exact measurement H2 proposes — Spearman against distance
to the nearest training sequence, for a learned model — is now confirmed un-done **twice over**, by
two papers I had cited as establishing it. So H2 is *more* novel than claimed. But it is now
**exploratory rather than confirmatory**: the decay is plausible and partially evidenced in adjacent
quantities, NOT established for the quantity we would measure.

### C2 — H7's promotion was OVER-STATED, and this is the one real rank change

Both quotes confirmed verbatim: *"inter-VEP concordance showed no relationship with agreement to
experimental MAVE data"*, and *"We therefore propose using inter-VEP disagreement as a practical
strategy to prioritise proteins"*.

**But the strategy is PROPOSED, not tested.** Checked: the paper's analyses are measuring concordance,
checking it against MAVE agreement, and examining features associated with concordance. There is no
experiment showing that assaying low-concordance proteins yields more informative results. Its
*"should maximise"* is a reasoned expectation.

Also: my *"ten predictors and more than 13,000 proteins"* attaches the 13,000 figure to the
concordance-vs-MAVE finding. The 13,000 covers the predictions analysis; **the MAVE-comparison subset
size is not stated**. And the disagreement metric itself is not defined in the abstract.

So H7 is a **published hypothesis, not a published method** — more novel, less de-risked than I said.
It moves DOWN out of the top tier.

### C3 — H5's "parameter-free" is HALF right

The headline is confirmed: the model *"predicts a linear relationship between microbial growth rate
and heterologous protein expression"*, and the growth law has *"zero free and two fixed parameters"*.
That half is genuinely parameter-free in the fitting sense.

**But the chain is not parameter-free end to end.** The same paper adds a separate promoter-activity
model carrying **one fitted parameter** (around 0.40 per hour) to predict expression across media.
Further qualifications read directly from the source: the best-fit slope (2.33) falls **below the
model's own predicted lower bound** (2.7); in minimal medium the promoter-activity model
**overpredicted** at low growth rates and needed a correction term; it assumes steady-state balanced
growth and host-average translation rate; and the authors state it sets a **lower bound** on fitness
cost, applicable mainly to lab-adapted bacteria.

**One genuinely good find, unprompted.** The library is **synthetic constitutive promoters made by
mutagenising the −35 and −10 regions** — structurally very close to GSE108535's combinatorial element
design. That is a better substrate match than I knew. Offsetting it: only 22 eGFP plus 9 amiE
variants, and the measured axis is **growth rate**, which we do not have locally. The one source that
might have said whether a promoter-strength-paired-with-growth dataset exists (the eLife random-promoter
paper) was **unreachable** — see the blocked list.

### C4 — H4's conclusion SURVIVES, but my citation was not from the paper I cited

The invariance-versus-equivariance framing and the promoter-strand quote I attributed to the
regulatory-genomics equivariance paper **are not in it**. Checked: it *"doesn't use the terms
'invariant' or 'promoter classification'"* and *"does not discuss whether a gene's promoter lies on
the + or − strand"*. What it does offer is a general directional caveat — that some layer weights
*"are not expected to be symmetric"* — and it never mentions MPRAs, so it does not address our case.

**The argument is now properly sourced from elsewhere.** Caduceus states the premise behind RC
equivariance outright in its motivation: the two strands *"are reverse complements of each other and
that carry the same information"*. **That premise is exactly what fails for a directional MPRA insert**
whose reverse complement is a different construct with a different expected activity. So C2 holds, on
a better foundation than the one I gave it.

### C5 — H1's binder-design support is UNVERIFIED

The claim that hallucination/oracle-optimising methods outperform generative models in binder design
came from a preprint I never fetched. The paper I *did* reach (the systematic generative-model
comparison, read from its PDF) **compares generative models to each other** and to reference sets on
secondary-structure and amino-acid enrichment. It does not test oracle-guided search against
generative models. H1's re-aimed form still stands on **repo-internal** grounds — `fba/design.py`
performs bounded enumeration plus ranking and works — not on this citation.

### C6 — H6/ALDE gets STRONGER, and one of my framings was wrong

Both quotes confirmed verbatim (*"frequentist uncertainty quantification works more consistently than
typical Bayesian approaches"*, *"incorporating deep learning does not always boost performance"*).

**New and better than I said: there is a real wet-lab campaign**, improving yield from 12% to 93% over
three rounds, reaching 99% total yield. I had described ALDE as essentially simulation-based.

**New and worse than I said:** the comparison against directed evolution was **never run in the lab**
— the "around five rounds for DE" figure is an estimate. The outperformance claim is a **simulation**
result on existing combinatorial datasets. The deep-learning caveat is also weaker than its own
introduction implies (results show higher-dimensional encodings generally work better with
deep-learning models). Batch selection is *"a straightforward heuristic"* (top-q by UCB) in
simulation and Thompson sampling in the lab; batch expected improvement was tested and dropped. And
the extension to larger design spaces — *"should naturally extend"* — is explicitly **not tested**.

## What did NOT change

**H3's demotion stands and is the most solid finding in the whole set** — the no-relationship result
is verbatim and is the paper's own headline.

**H8 stays Tier 0, and its backing was never the literature.** Three independent routes to the
strain-design review were blocked, so the 2005 lactate numbers and the growth-coupling caveat remain
**search-summary level and UNVERIFIED**. This does not weaken the recommendation, because Tier 0 rests
on repo facts verified directly in code: `fba/design.py` exists, is CLI-wired, is tested, has zero
result artifacts, `cobra` is absent, and the model file is present. **The recommendation to run it is
grounded in the repo, not in the papers.**

## Depth actually achieved per source — the honest ledger

**Primary-source contact (8):**

| source | depth | outcome |
|---|---|---|
| Disagreement among VEPs (biorxiv 712765) | abstract | both quotes CONFIRMED; proposal-not-tested found |
| Promoter strength / growth rate (PLOS ONE) | full text | CONFIRMED with major qualification (C3) |
| Determinants of ML fitness performance (biorxiv 319780) | full text, first 100k chars | quote CONFIRMED; simulated, Pearson, regime-counting (C1) |
| ALDE (biorxiv 605457) | full text, first 100k chars | both quotes CONFIRMED; wet-lab result found; DE comparison never run (C6) |
| PLM fitness is a matter of preference (biorxiv 616542) | full text | **citation WRONG** — influence, not accuracy (C1) |
| RC equivariance in regulatory genomics (biorxiv 368803) | full text, first 100k chars | **quotes NOT from this paper** (C4) |
| Caduceus (arXiv PDF, text-extracted locally) | 17 pp extracted | supplies the correct premise for C2 |
| Systematic generative-model comparison (biorxiv PDF, extracted) | 30 pp extracted | does NOT support H1's claim (C5) |

**Unreachable from this host (7), each tried and named:**

- Strain-optimization review — PMC cookie wall, journal 403, Europe PMC 403. **Three routes, all blocked.** H8's numbers unverified.
- Resource-allocation / RBS-and-promoter model — PMC cookie wall on both URL forms.
- Random-promoter eLife paper — render error, twice. **This is the one that could have answered H5's substrate question.**
- Reverse-complement consistency (OpenReview) — browser-verification screen.
- ALDE journal version — publisher auth redirect (covered by the preprint instead).
- Genie 3 preprint, Nature extrapolation paper, OptCouple, METL, query-by-committee PDF, patents — not attempted after the PMC/Nature/publisher block pattern was established.

**Two method notes worth keeping.** Publisher and PMC routes are largely blocked from here while
**bioRxiv full-text and arXiv PDFs work**; prefer preprint servers. And a PDF that WebFetch returns as
a byte dump is **saved locally** and can be text-extracted with an ephemeral `uv run --with pypdf` —
no project dependency added, `pyproject.toml` and the frozen `uv.lock` untouched.

## Net effect on the ranking

| | before | after | why |
|---|---|---|---|
| H8 strain designer | Tier 0 | **Tier 0** | unchanged; rests on verified repo facts, not papers |
| H2 decay curve | 1 | **1** | holds, but motivation now exploratory not established (C1) |
| H5 nearer target | 2 (was 1) | **2** | parameter-free corrected; substrate match found (C3) |
| H7 disagreement | 3 | **5** | published PROPOSAL, not method — the one real rank change (C2) |
| H1 enumeration | 4 | 4 | support unverified; stands on repo grounds (C5) |
| H6 ALDE | 5 | **3** | real wet-lab result; better understood (C6) |
| H4 symmetry | 6 | 6 | conclusion survives, better sourced (C4) |
| H3 agreement filter | 7 | 7 | demotion confirmed verbatim — most solid in the set |

## Honest limits of THIS pass

- Three sources were read only to their first 100,000 characters; the tail was not read, and the
  omission is flagged per-source above.
- The disagreement preprint was read at **abstract level only**, and it is the source behind two
  rank decisions. Its disagreement metric is still unknown to me.
- Seven sources are unverified, one of them (the strain-design review) after three independent
  attempts. Those claims remain at search-summary level and are labelled as such wherever they appear.
- WebFetch answers come from a small model reading the page. Where it labelled its own reading an
  inference rather than a finding, I have carried that distinction through rather than collapsing it.
