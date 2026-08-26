# The colour-cell family has a CURATION wall in front of its substrate wall (2026-08-26)

The animal colour/plumage family is **19 CLI cells / 65 loci**, every one shipping as `KNOWLEDGE_BASELINE`.
It has been put to a measured test **exactly once** — the dog cell against the free Darwin's Ark cohort
(`wiki/dog_coat_darwins_ark_measured_2026-07-30.md`; N=3,277 genotypes × 29M biallelic SNVs, N=1,930
owner-reported colours) — and mostly failed on **substrate**, not biology:

> black 160/161 = **0.994** · blue/grey 11/31 = 0.355 · every other base colour **unscorable**

because K^B is a 3 bp deletion, ASIP A^y/a^t a SINE insertion, MLPH d3 a frameshift, and MC1R `e` fell in an
imputation gap. An imputed biallelic-SNV panel cannot represent any of those.

This screen asks how far that generalises, **derived from the committed catalogs** rather than asserted:
`scripts/colour_cell_substrate_screen.py` classifies the causal variant each locus records in its own
catalog text. Two findings, and **the second is the bigger one**.

## Finding 1 — where a causal variant IS recorded, most of it is unreachable by a SNV panel

| class | loci | SNV-panel representable? |
|---|---:|---|
| SNV | 11 | yes |
| INDEL | 11 | **no** |
| STRUCTURAL | 3 | **no** |

**14 of the 25 loci with a recorded variant (56%) are indel/structural.** The dog result generalises: this
family rests disproportionately on variant classes a SNP array or imputed panel cannot carry.

## Finding 2 — 62% of loci do not record a causal variant AT ALL

**40 of 65 loci (62%) record no causal variant**, and **7 of 19 cells record none for any locus**:
alpaca, cattle, mouse, pig, pigeon, rabbit, sheep.

Those cells encode allele *symbols* and dominance order — `rabbit C (TYR): C full > chinchilla > Himalayan
> c albino` — with no variant to genotype. That is a stronger statement than "unvalidated":

> **They are unvalidatable as written.** You cannot score a locus whose causal variant is unspecified, so
> no cohort would help. The blocker is **curation, not data.**

Only two cells are fully SNV-tractable, and both are small: **donkey** (3/3) and **roe deer** (1/1).

| verdict | cells |
|---|---|
| `UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED` | alpaca, cattle, mouse, pig, pigeon, rabbit, sheep (7) |
| `PARTIALLY_SNV_TRACTABLE` | camel, cat, chicken, dog, fox, goat, guineapig, horse, mink (9) |
| `NO_LOCUS_SNV_TRACTABLE` | buffalo (1) |
| `FULLY_SNV_TRACTABLE` | donkey, roe deer (2) |

## What this corrects about the family's own framing

The strategic read going in was "the colour family is at a **substrate** wall." That is only half right —
the substrate wall is real (Finding 1) but it sits **behind** a curation wall (Finding 2) that blocks
two-thirds of the loci before substrate is even reachable.

It also corrects the dog cell's status. `coatcolor` ships as `KNOWLEDGE_BASELINE` with the Darwin's Ark
scoring described as "the v0.1 measured tier" — **pending framing for work that already ran**. Its one real
measured result (black 0.994, N=161) is not reported on the trust surface, and under-claiming is as much a
trust-surface falsehood as over-claiming.

## The screen's own honest limits

- It reads **what the catalog records**, not what the literature knows. `UNRECORDED` is a statement about
  the catalog, never evidence about the substrate.
- The classifier is a **text heuristic** over free-form provenance strings. It is anchored by a self-check
  against the dog cell (`--self-check`), which is the one case with measured ground truth.
- **The self-check earned its keep on the first run.** It flagged dog `A` as `UNRECORDED` against an
  expectation of `STRUCTURAL` — and the *classifier was right*: the dog ASIP entry names only the locus and
  papers ("OMIA (A locus ASIP); Dreger & Schmutz 2011; Bannasch 2021"), never the SINE. The expectation had
  encoded the literature instead of the catalog. Recorded as `catalog_gaps_vs_measured_artifact`, because
  even the most-developed colour cell fails to record one of its five causal variants.
- No cell here has been measured except dog. A blocked locus means a SNV panel cannot **represent** the
  variant; it does not by itself predict that a cell would fail.

## Consequence

**Adding colour cell #20 adds a rule that cannot be validated on any substrate.** The cheap, real move for
this family is *curation* — record causal variants for the 40 unrecorded loci — which would convert
"unvalidatable" into "unvalidated but screenable". Whether that curation is worth doing at all is a scope
call, since the substrate wall (Finding 1) still waits behind it for ~56% of what gets recorded.

## Reproduce

```bash
uv run python scripts/colour_cell_substrate_screen.py --self-check   # anchor the classifier first
uv run python scripts/colour_cell_substrate_screen.py                # no network, no GPU
```
Artifact: `wiki/colour_cell_substrate_screen_2026-08-26.json`. Frozen AMR surface untouched.
