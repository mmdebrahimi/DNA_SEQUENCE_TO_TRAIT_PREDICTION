# The GLM headroom result generalizes off one drug — 7 of 8 RT drugs, leave-one-study-out

**2026-09-30.** `wiki/glm_alphabet_headroom_2026-09-30.md` showed a supervised model over genotype tokens
clearing the pre-registered bar that zero-shot ESM2 failed — on **one drug, one gene**. That was its
largest honest limit. This closes it across the RT drugs.

It was nearly free: NNRTI and NRTI drugs **share the RT protein**, so the cached masked-marginals already
cover every drug, and there is a **committed zero-shot all-RT baseline** to compare against
(`wiki/hiv_esm_vs_catalog_allrt_2026-07-09.json` — median subset AUROC 0.4558, **0 of 6 genuine passes**).

## The result

| variant | split | median zero-shot | median supervised | supervised passes | zero-shot passes |
|---|---|---|---|---|---|
| non-Full (anchored) | ungrouped — **optimistic** | 0.4558 | 0.8084 | **5 / 5** | **0 / 5** |
| `.Full` | **leave-one-STUDY-out** | 0.5121 | **0.8141** | **7 / 8** | **0 / 8** |

**The headline is the second row: 7 of 8 RT drugs pass under leave-one-study-out, against 0 of 8 for
zero-shot on the identical subsets.**

Per-drug, `.Full`, leave-one-study-out:

| drug | class | subset n | subset R | zero-shot | supervised | delta | verdict |
|---|---|---|---|---|---|---|---|
| EFV | NNRTI | 1,883 | 201 | 0.4750 | 0.8141 | +0.339 | PASS |
| NVP | NNRTI | 1,976 | 348 | 0.4833 | 0.7244 | +0.241 | PASS |
| ETR | NNRTI | 786 | 118 | 0.4933 | 0.8788 | +0.386 | PASS |
| RPV | NNRTI | 364 | 103 | 0.4828 | 0.8357 | +0.353 | PASS |
| 3TC | NRTI | 924 | 28 | 0.6331 | 0.7742 | +0.141 | PASS |
| AZT | NRTI | 944 | 17 | 0.6077 | 0.6728 | +0.065 | PASS |
| DDI | NRTI | 914 | 19 | 0.5945 | 0.8727 | +0.278 | PASS |
| **D4T** | NRTI | 897 | 11 | 0.5121 | 0.9053 | +0.393 | **FAIL_NULL_TOO_HIGH** |
| DOR | NNRTI | 97 | 54 | 0.5937 | 0.6981 | +0.104 | PASS *(small subset — excluded from genuine)* |
| ABC / TDF | NRTI | 874 / 719 | 9 / 2 | — | — | — | UNDERPOWERED |

## Three things that make this trustworthy rather than flattering

**1. Six per-drug anchors, all binding.** For each drug the zero-shot subset AUROC is recomputed and
required to match the committed value; **6/6 reproduce** (tol 0.02). A drug whose anchor missed would be
reported `ANCHOR_MISMATCH` with its supervised number **withheld** — a head-to-head on a differently-built
subset is not a head-to-head. The subset recipe (`DRIFT`, `CUTOFF`, `subs_of`) and the deployed comparator
`call_hiv_observed` are **imported** from the committed run's own modules, not re-written.

**2. The nominal/genuine split — the trap this repo already documented.** A raw pass-count launders a
small-subset artefact: the committed run reports `nominal 1 / genuine 0`, its nominal pass being **DOR on
a 37-isolate subset with mutation-burden 0.783**. The same split is applied to the supervised arm, and
every excluded drug is **named with its n**. Excluding DOR moves zero-shot from 1/6 to **0/5 — independently
reproducing the committed run's own `genuine: 0`**, a consistency check this design did not aim for.

**3. The bar bites.** **D4T FAILS** (`FAIL_NULL_TOO_HIGH`, shuffled null 0.571 ≥ 0.55) despite the highest
raw supervised AUROC in the table (0.9053) — on 11 resistant isolates of 897, grouped CV is unstable and
the null says so. A run where every drug passed would be a bar that cannot discriminate.

The bar itself is **imported** from `glm_alphabet_headroom.py`, so it has one definition across both
scripts and cannot drift between them; a test asserts it is not re-declared here.

## Honest limits

- **Every drug shares ONE gene (RT).** This generalizes across **drugs**, not across **genes**. A second
  gene (protease, integrase) is a separate test — and the PI blind spot is nearly empty
  (`wiki/supervised_complement_per_class_nonvacuity_2026-09-30.md`), so the drug is not the only thing
  that has to be chosen.
- **A linear model over one-hot tokens is the WEAKEST supervised family member.** A pass is a **floor** on
  that family and says nothing about whether attention or context helps. That remains the open question.
- **IN-DISTRIBUTION** to the Stanford knowledge base throughout; ~96% subtype B, so subtype generalisation
  is untested.
- The `.Full` and non-Full variants are **different isolate sets** (the non-Full files carry no `RefID`, so
  that arm cannot be grouped and is labelled `OPTIMISTIC` in its own split name). The anchored numbers and
  the grouped headline are not two estimates of one quantity.
- **2 of 11 RT drugs remain underpowered** (ABC 9 resistant, TDF 2) and are reported, never scored.
- Fold ≥ 3 is the repo's uniform illustrative cutoff, not a per-drug clinical breakpoint.

## Reproduce

```
uv run python scripts/glm_headroom_allrt.py --self-check   # bar + guard provenance, no data
uv run python scripts/glm_headroom_allrt.py               # 11 drugs x 2 variants + 6 anchors
```

Artifact `wiki/glm_headroom_allrt_2026-09-30.json`; 10 tests `tests/test_glm_headroom_allrt.py`. Exits 2
without the gitignored data / ESM cache, and 3 if no anchor binds.
