# The PromoterAtlas / Promoter Calculator circularity check — one half RESOLVED, one half UNRESOLVED with the flag possibly aimed at the wrong rung

<!-- doc-status: LIVE as of 2026-10-31. Discharges Step 6 of
     plans/Reconcile_The_Three_Disagreeing_Registries_And_Re_Scope_The_Stack_Arm_Plan (which
     prescribed the filename wiki/promoteratlas_circularity_2026-10-23.md; written today instead).
     The plan says an unresolvable check is a RESULT, not a failure, and prescribes stamping
     circularity: "unresolved" onto the arm with the list of what was checked. Both are done. -->

## The question, as the repo recorded it

Two committed artifacts carry the same open flag, verbatim:

- `wiki/mpra_substrate_census_2026-10-12.json` — *"whether the Promoter Calculator's 346 parameters
  or PromoterAtlas's urtecho2018 head were fit on GSE108535, which would make either a circular
  comparator on the sigma70 rung"*
- `wiki/sigma70_loader_reconciliation_2026-10-13.json` — *"circularity UNRESOLVED: the Promoter
  Calculator's 346 params and PromoterAtlas's urtecho2018 head may both have been fit on GSE108535;
  neither was checked"*

The plan's Step 6 required settling this **before any spend decision** on
`promoteratlas-base-embeddings-frozen-150bp-padded-to-200`, because *"spending the arm before this
check risks a contaminated number that has to be retracted."*

## HALF ONE — the Promoter Calculator is CLEARED, and this is the half that mattered most

**Verdict: NOT circular with respect to GSE108535.**

The Promoter Calculator (LaFleur, Hossain & Salis 2022, *Nature Communications*; PMC9440211) was fit
on **the Salis lab's own library of 14,206 designed sigma70 promoter variants**, built by oligopool
synthesis, and validated on 17,396 promoters. Its 346 parameters are ΔG terms for RNAP/σ70 binding
and bubble formation. GSE108535 (Urtecho, Tripp, Insigne, Kim & Kosuri, *Biochemistry* 2019) is a
**different library from a different lab** and does not appear as a training source.

**Why this is the load-bearing half.** The flag named the Promoter Calculator's 346 params, and those
are the `max_tx` channel that scores **0.3968** — the best single arm on the tile rung — and a
component of both the 2-channel hybrid (**0.4254**) and the 3-channel stack (**0.4404**), the
project's best measured predictor. Nobody had checked it. Had it been circular, the three strongest
numbers on this criterion would have needed retraction.

**Residual caveat, stated rather than buried.** Both libraries are *designed sigma70 variant*
panels, so the two assays are of the same general type even though the sequences and labs differ.
That is sequence-independence, not assay-type independence.

## HALF TWO — PromoterAtlas is UNRESOLVED, and the flag may be pointed at the wrong rung

**Verdict: `unresolved`.** Stamped onto the arm rather than left in a wiki open-question list.

### What could not be established

No source found uses the name **"PromoterAtlas"**. The repo never names the weights source either —
`wiki/ceiling_arm_wall_classification_2026-10-26.md` says only *"a public-weights fetch plus disk"*
and records that `D:/hf_cache/hub` holds exactly one model, GENERator-v2, not this one. So the arm
names a model whose identity and provenance are both unpinned.

Routes attempted:

| route | outcome |
|---|---|
| search for the name plus its `urtecho2018` head | no page uses the name |
| `ecolipromoterdb.com` (the atlas site named in the Urtecho reviews) | **FETCH FAILED** — TLS cert covers only `*.netlify.app`, not the bare domain |
| the repo's own artifacts | flag restated in two places, never checked; weights source never named |

### The lead, and it points somewhere worse than the flag does

The Urtecho **genome-wide** study (bioRxiv 2020.01.04.894907; eLife reviewed preprint 92558) measured
promoter activity of **321,123 sheared genomic fragments** across both strands of the *E. coli*
genome, plus 17,635 TSS-containing promoter regions and a **tiled** library, and *"trained several
machine learning models on these datasets to ... quantitatively predict promoter function from
sequence."*

**321,123 is this repo's fragment rung exactly.** The arms ledger carries
`fragment-rung-321k-rows-size-lever-LB-M9` and the G-A ledger records *"LB 321,123 / M9 318,457"* —
i.e. **GSE144621**, which is also the peak-tile substrate the ceiling criterion is defined on
(44,106 tiles, 0.941 replicate ceiling).

So **if** the arm's model is that study's, its training data is the tile/fragment rung **itself**,
and the arm would be circular **on the TILE rung** — not on the sigma70 rung the flag names. That is
a worse problem than the one recorded, because it would make the arm **inadmissible as a ceiling
comparator** rather than merely low-merit.

**This is a LEAD, not a finding.** The year labels do not line up cleanly with anything: the repo's
flag says `urtecho2018`, GSE108535 is Urtecho 2019 (*Biochemistry*), and the genome-wide work is a
2020 preprint / 2023 eLife reviewed preprint. Which dataset `urtecho2018` denotes is still open, and
three Urtecho-family datasets are now in play rather than one.

### What would settle it

Identify the weights source first — the arm cannot be checked while the model it names is unpinned.
Then read that model card's or paper's training-data section. `ecolipromoterdb.com` needs a route
that tolerates the cert mismatch, or the eLife paper's methods read directly.

## The merit judgement, carried forward unchanged and still INFERRED

The plan already recorded a separate, negative merit inference that this check does not alter:
PromoterAtlas would be a **third sequence-learned channel** beside the `conv` already in the stack,
which is the failure mode `gc` demonstrated at 0.4217 against the hybrid's 0.4254. Every
sequence-statistics arm on this rung clusters **0.2691–0.3234** while the one **thermodynamic** arm
sits at **0.3968**, so the next channel should carry a different **PRIOR**, not a different neural
net. The orthogonality measurement that would settle it needs per-tile vectors the arm has not
produced, so this stays an inference.

## Consequence for the arm

**Do NOT spend it yet.** Two independent reasons, and the first is new:

1. `circularity: unresolved`, with a live lead that the contamination would be on the **measured
   rung**, which would void the number rather than weaken it.
2. The standing inferred merit-negative above.

Per R6 the arm **stays `planned`** — not killed. The blocker now travels with the arm.

## Honest limits

- Half two is **unresolved, not negative.** Nothing here shows the arm *is* circular.
- The Promoter Calculator clearance rests on secondary sources (the paper's PMC record and the
  authors' own conference abstract), not on a line-by-line read of its supplementary training tables.
- The 321,123 match is a **numerical coincidence plus a plausible identification**, not a
  provenance chain. It is strong enough to block a spend and not strong enough to assert.
- One prior assumption of mine was wrong and is corrected here: I recommended this arm as the next
  highest-VOI move without checking that the repo had a prescribed pre-step and a recorded merit
  inference against it. Both existed.
