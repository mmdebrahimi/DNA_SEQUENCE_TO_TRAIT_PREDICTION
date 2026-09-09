# The ktype caller lands near, but below, the wzi method's own published ceiling

`typing:Klebsiella:ktype` was the **last** cell carrying the shared "never measured" default. It looked
cohort-blocked — Klebsiella capsule typing has no free wet-lab label — until reading the caller showed
the blocker was mis-stated.

## Why a comparator existed after all

Our caller is **not** a Kaptive wrapper. It types the capsule from a **single conserved gene** (`wzi`,
BIGSdb Pasteur scheme as bundled by Kleborate) and maps the allele to a K locus. **Kaptive types the
full K locus.** Two different methods over different amounts of sequence — the same relationship as
AMRFinder to ResFinder — so Kaptive is a genuine independent comparator, and **no wet-lab cohort is
needed to run it**. The wet-lab cohort is what an `INDEPENDENT_MEASURED` tier would require; this is a
tool-tier measurement.

## The question is not "do they agree 100%"

They cannot, and the caller's own docstring says so: **wzi → K-type is ~94% predictive and NOT
one-to-one**, because isolates with distinct K types can share a wzi allele (Brisse 2013 *JCM*). A
shortfall from 100% is a property of the **method**. The measurable question is:

> does our implementation **achieve** the wzi method's ceiling, or fall short of it?

Falling well short would be an implementation defect — the thing a measurement can fix. Landing at the
ceiling means the caller is faithful and the residual is the method's limit.

## Result — 307 Klebsiella genomes

| | |
|---|---|
| genomes with a cached assembly | **307** |
| Kaptive errored / not typeable | **0 / 0** |
| our caller abstained (not an error) | **43 (14.0%)** |
| comparable | **264** |
| **agreement (strict)** | **232 / 264 = 0.8788** |
| agreement crediting ambiguous calls | 0.8902 |
| published wzi ceiling | ~0.94 |

**Verdict `NEAR_THE_WZI_METHOD_CEILING`** — about 6 points below strict, 5 lenient. Within 10 points, so
this is *consistent with* the method's own limit plus cohort composition and is **not clean evidence of
an implementation defect** — but it is also not "reaches the ceiling", and I am not claiming it does.

**Two figures, and which is the headline.** An ambiguous multi-KL string is **not** a resolved call, so
the strict figure counts it a miss and is the headline. Only **3** of 32 disagreements are of that class
— all `KL22/KL37` vs Kaptive's `KL22`, all from `wzi_37` — and in all 3 Kaptive's answer is inside ours.
That is the caller reporting the method's limit honestly rather than guessing wrong, which is why the
lenient figure is printed beside it rather than hidden.

## The disagreements are concentrated, not diffuse

| ours → Kaptive | n |
|---|---|
| KL15 → KL52 | 6 |
| KL19 → KL57 | 4 |
| KL15 → KL51 | 4 |
| KL22/KL37 → KL22 | 3 (ambiguous, see above) |
| KL163 → KL183 | 2 |
| KL18 → KL10 | 2 |

**KL15 accounts for 10 of the 29 genuine misses.** Concentration is the same tell that located the
serotype defect, so it is recorded as a **lead** — but it is *not* a diagnosis. A single wzi allele
legitimately mapping to several K loci is exactly the published non-one-to-one behaviour, so this
pattern is equally consistent with the method working as designed. Distinguishing the two would need
the per-allele mapping examined against Kaptive's locus-level evidence, which this run did not do.

## Honest limits

- **The comparator is a TOOL, not a wet-lab label.** This measures agreement with an independent
  implementation, never correctness — both could be wrong together. **The cell stays
  `FAITHFUL_TO_TOOL`**; a serology-labelled Klebsiella cohort is what a measured tier needs and none is
  free.
- Genomes are Klebsiella from this project's **AMR cohorts**, so they are **enriched for resistance**
  and their K-type distribution need not match a population.
- The **~94% ceiling is a published property of the wzi method on its own cohorts, not a constant.**
  Treating it as an exact bar would over-read it, which is why the verdict has a near-ceiling band
  rather than a single threshold.
- Kaptive rows it does not call typeable would be excluded from the denominator; here there were **none**,
  so the denominator is the full comparable set.
- Our caller's **43 abstentions are counted separately** from disagreements. Abstention is not error,
  and pooling them would understate the caller.

## Reproduce

```bash
uv tool install kaptive && kaptive db install kpsc_k
uv run python scripts/ktype_kaptive_concordance.py
```

Needs blastn + cached Klebsiella assemblies. Kaptive is installed **isolated** (`uv tool`), deliberately:
it pulls numpy 2.5.3, and installing it into the project venv would have moved numpy under the entire
test suite. Frozen AMR surface byte-unchanged — this is a typing cell.
See [`ktype_kaptive_concordance_2026-09-08.json`](ktype_kaptive_concordance_2026-09-08.json).
