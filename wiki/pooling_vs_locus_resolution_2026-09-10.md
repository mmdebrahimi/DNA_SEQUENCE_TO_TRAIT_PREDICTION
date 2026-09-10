# Locus resolution does NOT recover what pooling destroys — the survivor's deployable claim is refuted (2026-09-10)

**Verdict `REFUTED_LOCUS_IDENTITY`, computed mechanically from a bar frozen before any number was seen.**

The 2026-09-10 `/innovate` sweep's top survivor, `POOLING-dilution`, argued that mean-pooling — not
population design — is the mechanism behind the 0-for-5 de-confounded embedding failures, on the measured
grounds that PC1 explains **0.807** of pooled-embedding variance. That measurement is correlational. Its
**deployable half** was a prediction, and this is that prediction executed:

> restricting the representation to candidate **causal loci** before pooling recovers within-lineage
> signal that whole-genome pooling destroys.

Reproduce: `uv run python scripts/pooling_vs_locus_resolution.py` (offline, no new embedding compute —
the per-gene vectors were populated on Databricks in May and are simply re-aggregated).

## Result

N=123 strains (60R / 63S) over 71 MLST lineages, leave-one-lineage-out CV, both arms mean-pooled to the
same 512 dims so dimensionality is held constant:

| representation | genes pooled | AUROC |
|---|---|---|
| whole-genome pooled | ~4,900 | 0.8029 |
| **QRDR-restricted** (gyrA, gyrB, parC, parE) | 4 | **0.8386** |
| **random 4 single-copy genes** (50 draws) | 4 | mean 0.7479 · p95 0.8337 · **max 0.9020** |

**The frozen rule required BOTH** `qrdr > whole-genome` **and** `qrdr > max(null)`. The first holds; the
second fails. Per the pre-registered middle branch this is a **refutation of the deployable claim**, not a
partial win — and that branch was written down in advance precisely so the outcome could not be re-read
afterwards.

## What actually happened: arbitrary genes beat the causal ones

**A random draw of four arbitrary single-copy genes reached 0.9020 — higher than the four genes that
actually cause ciprofloxacin resistance, and higher than pooling the whole genome.** QRDR sits around the
95th percentile of the null: better than a typical random draw, comfortably inside its range.

That is what **ancestry-correlated resistance** looks like. If resistance tracks lineage, then *any*
representation that encodes lineage predicts it — and four arbitrary genes encode lineage perfectly well.
Narrowing the representation does not escape the confound; it just changes which genes carry it.

So the sweep's PC1 = 0.807 observation stands as a fact about the pooled representation, but the
conclusion drawn from it — that **locus identity** is the lever — does not survive contact with its own
control. **The control was the load-bearing arm, and it is what killed the claim.**

## This does NOT contradict the 0-for-5 record — it is a weaker test

Leave-one-lineage-**out** lets a held-out ST sit phylogenetically close to a training ST, so ancestry
signal survives it. The project's 0-for-5 rests on the strictly harder **within-lineage** question. A
whole-genome AUROC of 0.8029 here and "within-lineage = chance" there are answers to *different questions*
and are not in tension.

**And the strict test is barely answerable on this cohort**, which is why it was not attempted: of 71
lineages only **5 carry both labels**, and the largest of those is ST131 at **25R / 1S** — effectively
single-class. Recorded rather than worked around.

## Two defects found in this run, both in my own code

1. **The null ran ZERO draws and printed nothing.** It drew on `gene_id`, which is **strain-unique**
   (`cds-ETY41380.1`), so the intersection across strains was the empty set — the documented 0%-overlap
   trap, biting the control rather than the model. The first run therefore produced two AUROCs and no
   null, which per the frozen bar is *uninterpretable*. It now draws on single-copy **product** strings —
   the same cross-strain identity `resolve_qrdr` uses — so each draw is genuinely the same four genes in
   every strain. **It failed loudly by printing nothing, which is the only reason it was caught.**
2. **Two annotation vocabularies, one collision.** These GFF3s come from two sources: one names gyrA
   *"DNA gyrase subunit A"*, the other *"DNA topoisomerase (ATP-hydrolyzing) subunit A"* — because gyrase
   **is** a type-II topoisomerase. A loose `topoisomerase.*subunit A` pattern therefore matches **gyrA and
   parC** and silently conflates them. The numeral differs too (*topoisomerase 4* vs *IV*). Handling only
   the first vocabulary resolved gyrA on 78/140 strains; handling both raised it to 123, and
   `resolve_qrdr` **refuses** a strain rather than guess when a gene resolves ambiguously.

## Honest limits

- `gene_symbol` is empty on every row (GenBank GCA assemblies), so gene identity runs through `product`
  strings — the reason both vocabularies had to be handled explicitly.
- One drug, one organism, one cohort. 50 null draws, so `max` is a coarse upper edge.
- The QRDR arm carries prior knowledge the whole-genome arm would have to discover unaided. That
  asymmetry **is** the hypothesis, and it makes the refutation stronger, not weaker: the favoured arm lost.
- 16 strains were refused for unresolvable QRDR sets, so the scored set is biased toward whatever
  annotation vocabulary resolves cleanly.
- This tests the **representation** choice only. Nothing here speaks to the model, and the decisive
  separation of pooling-vs-population-design still needs the constructed-cross run
  (`scripts/yeast_bloom_fm_prep.py`, blocked on coordinate liftover).

Artifacts: `wiki/pooling_vs_locus_resolution_2026-09-10.json` · bar
`wiki/pooling_vs_locus_acceptance_bar.json` · tests `tests/test_pooling_vs_locus_resolution.py`.
