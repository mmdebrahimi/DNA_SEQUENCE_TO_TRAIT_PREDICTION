# The GLM question, answered on the bar the zero-shot arm failed

**2026-09-30.** Two results in this repo looked contradictory, and which one generalises decides whether a
genotype language model has anything to learn. This puts both on the **identical** isolate set and judges
them by the **same pre-registered bar**.

| | |
|---|---|
| 2026-07-09, recorded as a CLOSED NEGATIVE | ESM2-650M **zero-shot** on the catalog's blind spot: **0.4485** against a pre-registered `>=0.65 AND > burden AND null < 0.55`. **FAILED.** |
| 2026-09-30 (today) | the shipped **supervised** complement flags that same blind spot at **15–16×** enrichment on two independent drug classes |

Same blind spot, same label, same tokens available. The methods differ in exactly one way: **zero-shot
likelihood vs supervised fitting over position-resolved substitution tokens.**

## The result

**The bar is copied verbatim from the committed 2026-07-09 artifact** — the bar the zero-shot arm was held
to and failed. Judging a new method by a fresh bar is how a negative gets laundered into a positive, so
`PREREGISTERED` is frozen in code and a self-check asserts the old zero-shot number still fails it.

| dataset | n (blind spot) | R | zero-shot ESM2 | **supervised, strictest split** | burden | null | delta |
|---|---|---|---|---|---|---|---|
| validation | 1,111 | 53 | **0.4485** FAIL | **0.8446** PASS *(ungrouped — optimistic)* | 0.4655 | 0.5101 | **+0.396** |
| full | 1,883 | 201 | **0.4750** FAIL | **0.8142** PASS *(leave-one-**study**-out)* | 0.3946 | 0.4695 | **+0.339** |

**The headline is the second row: 0.8142 under leave-one-study-out** — the model predicts studies it never
trained on. Patient-grouped is 0.8223, so removing patient leakage costs almost nothing.

**PIPELINE ANCHOR — the check that makes the comparison trustworthy.** On the non-Full set our ESM2 returns
**0.4485 against the committed 0.4485**. The subset construction therefore agrees exactly with the
2026-07-09 run; had it disagreed the script exits 3 and publishes nothing, because a head-to-head on a
differently-built subset is not a head-to-head.

**The load-bearing clause is `>=0.65`, not the burden clause.** Mutation burden is *anti*-predictive here
(0.3946 — more mutations, less resistance), so beating it is easy and is not what carries the pass. The
real hurdle is the absolute floor: zero-shot 0.475 vs supervised 0.814.

## It cannot be rediscovering the catalog, and that is structural

The subset is **catalog-NEGATIVE by construction** — no isolate carries any catalogued major DRM. So every
feature the model uses is a substitution the deployed catalog does not carry. This is headroom, not
rediscovery.

**What it actually learned corroborates three independent findings in this repo:**

| token | carriers | coef | independently flagged by |
|---|---|---|---|
| **`103S`** | 13 | **2.69** (top) | the isolate hand-verified today: catalog says S, measured wet-lab fold ~**6.8×** |
| **`179D`** | 32 | 2.08 | the doubt purity screen (`V179F`, p=8.8e-06) **and** the 2026-09-01 curation measurement (`V179D`/`V179E`) |
| `101E` / `98G` / `101Q` / `190C` / `230I` / `138Q` | 5–84 | 1.4–1.9 | — |

Position **103** and position **179** are now each reached by **four** and **three** unrelated routes
respectively — a purity screen, a multivariate-OLS curation pass, a hand-verified isolate, and this
supervised fit. None of those methods shares a code path with another.

**Not asserted:** whether `98G` / `101E` / `190C` / `230I` / `138Q` are established NNRTI accessory
positions in the literature. That is a sourcing question I have **not** done, and stating it from memory is
the fabrication hazard this project guards against.

## What this changes, and what it does not

**It does not overturn the 2026-07-09 negative — it scopes it.** That result is about **zero-shot**, which
is exactly what `eval/regime.py` already records for natural populations (`natural × molecular × zero_shot
= LOSES_TO_CATALOG`). The supervised direction was never the thing tested, and it has **~0.34 AUROC of
headroom** on the same data.

**It does not contradict the June 5-drug bacterial functional-alphabet negative either.** That probe tested
**raw-sequence k-mers over whole genomes** against a **curated determinant alphabet**, within lineage, and
found the curated alphabet wins at 0.93–1.00 while k-mer sits near chance. Different alphabet, different
arm, different question. Reconciling the two gives the actionable rule:

> **The learned layer's headroom is exactly the size of the curated catalog's incompleteness — and that
> varies by arm.** Bacterial AMR for the 5 drugs tested: the catalog already captures the signal, so there
> is nothing to learn. HIV NNRTI: 1,111 of 2,168 isolates are catalog-negative and 53 of those are truly
> resistant, so there is real signal to learn — and a supervised model over genotype tokens finds it.

**It does not beat the catalog.** The catalog scores 0.926 on the full cohort; this model works *where the
catalog is silent*. Complement, not replacement — which is also why it belongs in the L2 doubt layer rather
than in the L1 rule.

## Honest limits

- **A linear model over one-hot tokens is the WEAKEST member of the supervised family.** A pass is a
  **floor** on that family, not a ceiling — and it says **nothing** about whether attention or sequence
  context adds anything. Which is the right shape for a go/no-go: the cheap floor cleared the bar, so the
  expensive question is now worth asking.
- **In-distribution** to the Stanford knowledge base. Leave-one-study-out is the strongest available
  de-confounding here, but the data is ~96% subtype B, so **subtype** generalisation is untested.
- One drug (EFV), one gene (RT), one label cutoff (fold ≥ 3). The cutoff is the repo's uniform illustrative
  one, not a per-drug clinical breakpoint.
- 201 resistant isolates in the full set — adequate, not large.
- `n_carriers` for several top features is single-digit (`230I` n=5, `138Q` n=8); those individual
  coefficients are not individually reliable even though the aggregate AUROC is.

## Reproduce

```
uv run python scripts/glm_alphabet_headroom.py --self-check   # pure bar logic, no data
uv run python scripts/glm_alphabet_headroom.py                # both datasets + the anchor check
```

Artifact `wiki/glm_alphabet_headroom_2026-09-30.json`. Exits **2** when the gitignored Stanford data or the
cached masked-marginals are absent, and **3** when the pipeline anchor disagrees with the committed run —
it refuses to publish rather than report a mismatched comparison.
