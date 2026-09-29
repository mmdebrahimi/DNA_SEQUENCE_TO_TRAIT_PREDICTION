# A learned genotype alphabet loses to a 1992 substitution matrix — on the protein's own measured data

**2026-09-29.** The cheapest decisive test of the genotype-alphabet premise, run on PEAR CTX-M-14.
Both pre-registered bars **fail**. Scripts: `scripts/pear_genotype_alphabet.py` +
`scripts/pear_genotype_alphabet_robustness.py`. Artifacts:
`wiki/pear_genotype_alphabet{,_robustness}_2026-09-29.json`.

## A correction first, because I made the error in the message that proposed this test

I told the user *"ESM2 gets ρ = 0.352 on it; the orthogonal hybrid was never run."* **The second half is
false.** The hybrid **was** run on this exact substrate on 2026-09-02 and it **lost**:

| on PEAR CTX-M-14 / cefotaxime (pooled, n=1,513) | Spearman |
|---|---|
| ESM2-650M | **0.3523** |
| ESM2 + ProSST (orthogonal hybrid) | **0.2042** |
| ProSST alone | **−0.0404** (below chance) |
| BLOSUM62 | 0.1982 |

Asserted from memory, in the same breath as proposing a measurement. The correct framing is *stronger*,
not weaker: the incumbent is hard to beat here **and** the standard modality trick already failed on it.

## The design, and why the split is the whole experiment

**Substrate:** 2,114 constructed CTX-M-14 variants, measured cefotaxime fitness, independent lab, clears
all ten rejection gates. 1,513 missense over 264 protein positions.

**Split by PROTEIN POSITION, not by variant.** Each position carries ~19 variants, so a random variant
split leaks the position mean and any model scores well by memorising "position 104 is tolerant". Holding
out whole positions forces generalisation across the alphabet — which is exactly the claim a genotype
alphabet makes.

**Mapping verified before any correlation was read:** all 2,114 REF bases match the CDS (the script
*refuses* otherwise), and the consequence ordering came out **silent +0.0009 > missense −0.0457 >
nonsense −0.1286**. No frame-shifted mapping produces that by accident.

**The key framing:** **BLOSUM62 already *is* a learned genotype alphabet** — a 20×20 substitution matrix
learned from alignments in 1992. So "learn an alphabet" is not a new idea; the testable question is whether
learning one from *this protein's own measured DMS* beats the generic 34-year-old one.

## Result

Held-out-position Spearman. **Not comparable to the pooled 0.352** — different quantity.

| variant | stride 2 | stride 3 | stride 4 | stride 5 |
|---|---|---|---|---|
| A · BLOSUM62 (zero-shot, generic) | +0.226 | +0.254 | +0.254 | +0.308 |
| B · ESM2-650M (zero-shot PLM) | +0.352 | +0.311 | +0.411 | +0.353 |
| C · learned alphabet, pair means | +0.184 | +0.216 | +0.197 | +0.267 |
| D · learned alphabet, ridge one-hot | +0.133 | +0.246 | +0.135 | +0.210 |
| E · ESM2 **+** supervised alphabet head | +0.262 | +0.351 | +0.339 | +0.355 |

**B1 — a learned alphabet beats BLOSUM62: FALSE on all four strides. STABLE.** Trained on 1,017 measured
variants *from this very protein*, the position-free learned alphabet is **worse** than a matrix that has
never seen it. Not a coverage artifact: only **4 of 496** held-out variants fell back to a default
(0.8%), far under the 50% pre-registered invalidator.

**B3 — a supervised head on top of ESM2 beats ESM2 alone: UNSTABLE.** Wins on strides 3 and 5, loses on 2
and 4. Position-clustered bootstrap (resampling **positions**, not variants — ~19 variants share a position
and are not independent): **mean +0.0387, 95% CI [−0.0218, +0.0992], P(Δ>0) = 0.891. The CI includes
zero.** So B2 also fails, since it passed only through E.

**Had the run stopped after the first stride, +0.040 would have been published as the cornerstone
result.** It is a coin-flip.

## What this does and does not close

**Closed, and diagnosably so:** a **position-free** learned substitution alphabet carries almost no
marginal signal on this substrate. The reason is interpretable rather than mysterious — "what does X→Y do
on average" is already better estimated from millions of alignment columns than from 1,017 measurements
on one protein. More DMS data on one protein will not fix that; it is the wrong axis.

**The sharpest reading, and it reframes the goal rather than blocking it: ESM2 already *is* a genotype
language model.** A transformer over amino-acid tokens with positional context, scoring 0.31–0.41 here
zero-shot — beating every learned-alphabet variant. In the molecular regime the GLM is *already
instantiated*; bolting a learned substitution alphabet onto it adds nothing reliable.

**NOT closed — and must not be over-compressed into "GLM is dead"** (over-compressing a scoped negative
is this project's most-repeated error, three recorded instances):

- **An alphabet WITH context** was not tested. Position-free was the *hypothesis under test*, not a
  limitation of the implementation — but it also means nothing here speaks to a contextual model.
- **Multi-assay supervised training** was not tested, and it is the version that matches the real GLM
  shape. One protein cannot support it; the repo already holds the harness (ProteinGym 217 assays, MaveDB
  2,383 held-out assays at ESM2 median 0.478).
- **Nucleotide/codon-level alphabets** were not tested here at all.

## Honest limits

- **One protein, one assay, one lab.** A position-held-out result on CTX-M-14, not a general law.
- **Strides are nested partitions of the same 264 positions** — they test partition-sensitivity, not
  replication. Four rows are not four experiments.
- **The asymmetry favours the supervised side**: ESM2 sees no measured fitness while C/D/E see two-thirds
  of positions. So **C/D losing to BLOSUM62 is decisive, while E beating B would only have been
  suggestive** — and it did not survive anyway. That asymmetry was pre-registered, not invoked afterwards.
- Held out by **position**, not by structural neighbourhood; adjacent positions share local context, so
  the supervised numbers are if anything optimistic.
- BLOSUM62 is the **shipped** `variant_effect.blosum62_score`, checked non-vacuous (I→V above W→P) rather
  than hand-typed.
