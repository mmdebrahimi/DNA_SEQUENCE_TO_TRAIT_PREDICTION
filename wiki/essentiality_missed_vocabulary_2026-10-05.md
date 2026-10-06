# Why the conserved-core decoder misses 83% of human essentials — and what that does to the ladder plan (2026-10-05)

## The question, and why it was worth answering before building anything

A 5-rung cross-organism "transfer ladder" was drafted to decide whether the decoder's E. coli→human decay
(AUROC 0.695 → 0.580) is **(A)** continuous with phylogenetic distance or **(B)** a cliff at the eukaryote
boundary. An adversarial review of that plan raised a third possibility neither covers, and asked whether a
cheaper direct test should run **first** and possibly pre-empt the whole ladder:

| mechanism | meaning |
|---|---|
| (i) | bacteria-only function — correctly absent in human (peptidoglycan, BAM complex) |
| **(ii)** | **conserved function under a BACTERIA-SPECIFIC NAME — a catalogue phrasing gap, not biology** |
| (iii) | host-specific function — genuinely outside a bacterial catalogue (proteasome, spliceosome) |

(ii) matters out of proportion: if it dominates, the honest conclusion is *"the catalogue has a phrasing
bug"*, nothing about phylogeny, and a ladder would read a vocabulary artifact as a biological law.

The decay is also not a discrimination failure. It is **silence**: **83.1% of human essential genes score
exactly zero** (E. coli 62.1%), and in human **11 of 13 core patterns have P(hit | non-essential) =
0.0000** — the patterns never misfire, they just never fire.

## Answer: NO, it does not pre-empt the ladder — and my first attempt at this said the opposite

**Verdict: `INDETERMINATE_BAR_SENSITIVE_VOCABULARY_CANNOT_ADJUDICATE`.**

| bar | frac mechanism (ii) | band |
|---|---|---|
| loose — any content word shared with a caught E. coli gene | **0.5548** | dominates |
| strict — words the decoder's REGEXES key on | **0.2580** | material |
| strict minus regex-fragment debris (`cell`, `iii`) | **0.2509** | material |

**The first version of this probe used only the loose bar, reported 0.5548, and produced
`PHRASING_GAP_DOMINATES_LADDER_IS_WRONG_INSTRUMENT`.** That verdict was retracted before publication, not
after. The loose bar is inflated >2× because the words doing the work are generic modifiers — `binding`,
`small`, `beta`, `alpha`, `cell`, `repeat`. Sharing "alpha" with E. coli says nothing about whether the
catalogue could see a gene. The honest definition of the catalogue's reach is the vocabulary its **regexes**
key on, and that lands in a different band. Three defensible tokenizations spanning a threshold crossing
means the middle of this quantity is not measurable this way, so the probe now **refuses a verdict** rather
than picking the flattering bar.

## What survives every bar

Two floors do not depend on the reach definition at all — an absent word is absent under any definition,
and the (ii)-floor words are both in the regex vocabulary and unambiguously functional:

| floor | n (of 566 missed) | share |
|---|---|---|
| mechanism (ii) — `polymerase`/`helicase`/`replication`/`division`/`topoisomerase`/`primase` | **57** | **0.1007** |
| mechanism (iii) — no E. coli vocabulary at all | **141** | **0.2491** |

So **host-specific absence is ~2.5× the demonstrable phrasing gap.** The top `ABSENT_FROM_ECOLI` words are
exactly the machinery a bacterial catalogue cannot contain: `nuclear` 34, `proteasome` 24, `splicing` 19,
`ribonucleoprotein` 19, `ubiquitin` 15, `mitochondrial` 14, `nucleolar` 11, `dead-box` 9.

## Consequence for the ladder plan

1. **The ladder is NOT invalidated.** The dominant robust mechanism is (iii) — genuinely host-specific
   function — which is a real biological signal, so a distance-vs-cliff question is meaningful.
2. **But ~10% of the human miss is a fixable phrasing gap with nothing to do with phylogeny**, and it must
   be separated out *before* any rung is scored or it will be silently misattributed to distance or to the
   domain boundary. Two patterns go to exactly 0.00 hits/1000 human genes — `DNA polymerase III|…|
   replicative` and `cell division protein Fts|divisome|…` — while human plainly has replicative
   polymerases and cell division.
3. **This is a floor, not an estimate.** 57 genes is what survives an adversarial bar; the true (ii) share
   sits somewhere between 0.10 and ~0.25 and this instrument cannot narrow it.

## Method — no curated biology, by construction

The tempting implementation is a hand-written synonym map ("POLD1 means replication"). That is the
fabrication hazard this repo refuses elsewhere: `dna_decode/data/ecoff_catalog.py` requires a `source_url`
plus `verbatim_quote` and *raises* rather than return a plausible default, because a wrong curated value
yields a fully self-consistent evaluation that is entirely wrong. A synonym map would be exactly that — an
unsourced biological claim doing the adjudicating.

So every word here comes from NCBI's own `description` / `name` fields, and the adjudication is a corpus
comparison the data performs on itself: for each content word of a missed human essential, *where does that
word live in the E. coli corpus* — inside the catalogue's regex reach, present-but-uncaught, or absent
entirely. The only biological list in the probe is the robust-(ii) floor, and a test asserts every word in
it already appears in the decoder's own patterns.

## Honest limits

- **A word is not a function.** Shared vocabulary does not prove the same complex. This measures vocabulary
  overlap — which is what mechanism (ii) *is* — and nothing more.
- Human labels are BAGEL CEGv2/NEGv1, **two curated extremes, not a genome-wide screen**, so the base rate
  0.431 is a construction artifact and is **not** comparable to E. coli's genome-wide 0.0928.
- `ECOLI_BUT_UNCAUGHT` (111 genes, 0.1961) cannot be split into (i) vs (ii) without asserting biology, so
  it ships as the joint bucket.
- Stopwords are hand-written, but they are English function words plus annotation boilerplate — not
  biological judgements. Removing `protein`/`domain` changes which words surface.
- One organism pair. Nothing here speaks to the bacterial rungs.

## Reproduce

```bash
uv run python scripts/essentiality_missed_vocabulary.py --top 25          # read-only
uv run python scripts/essentiality_missed_vocabulary.py --emit            # + wiki artifact
```

Needs the D: essentiality cache (BAGEL CEGv2/NEGv1 + `Homo_sapiens.gene_info.gz` +
`ecoli_k12_feature_table.txt.gz`); exits 2 if absent. 5 guards at
`tests/test_essentiality_missed_vocabulary.py`, including one that fails if a single bar could ever produce
the verdict again, and one that rejects a vacuously-empty robust floor.
