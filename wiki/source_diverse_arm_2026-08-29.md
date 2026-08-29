# The source-diverse arm, and the first thing it did was retract one of my own numbers

Last run I measured source concentration and showed that 3 of 10 SCORED cells rest on a single BioProject.
The counter-measurement — scoring the frozen decoder on a source-*diverse* cohort — was made in a throwaway
snippet. This turns it into tooling, with one design constraint that turned out to matter more than the
tooling itself.

## The self-applied standard

An arm whose argument is *"your cohort was too concentrated to see this"* must not ship a concentrated
cohort of its own, or it is just a second opinion from a differently-biased sample. So the arm gates every
cell on its own bar — **≥5 BioProjects and no source above 60%** — and a cell that fails emits
`status: source_concentrated` with **no metrics at all**.

**It refused 6 of 10 cells** (7 after the effective-sources gate added later the same day — see the Update at the end; `E. coli x meropenem` moved from `underpowered` to refused).

| cell | n | sources | largest share | verdict |
|---|---:|---:|---:|---|
| E. coli × ciprofloxacin | 131 | 8 | 31% | **scored** — acc 0.977 / sens 0.958 / spec 0.988 |
| E. coli × gentamicin | 131 | 8 | 32% | **scored** — acc 0.756 / sens 0.523 / spec 0.985 |
| E. coli × ceftriaxone | 117 | 6 | 36% | underpowered (S class < 20) |
| E. coli × meropenem | 71 | 5 | 48% | **refused** (see Update) — 2.42 effective sources |
| E. coli × tetracycline | 73 | 3 | 58% | **refused** — 3 BioProjects |
| Klebsiella × ceftriaxone | 66 | 5 | **68%** | **refused** — over the 60% share bar |
| Klebsiella × ciprofloxacin | 44 | 2 | 98% | **refused** |
| Klebsiella × gentamicin | 43 | 2 | 98% | **refused** |
| Klebsiella × meropenem | 44 | 1 | 100% | **refused** |
| Klebsiella × tetracycline | 41 | 1 | 100% | **refused** |

## The retraction

**`Klebsiella × ceftriaxone` is one I published last run** as one of three "independent validation
numbers" — acc 0.924, sens 0.978, spec 0.810 on n=66. Its cohort has 5 BioProjects, which looked diverse
enough that I didn't check the share. **68% comes from one of them**, over the arm's own bar.

Reporting a number from a cohort as concentrated as the ones this work criticises is exactly the error
being exposed. It is now `source_concentrated` with no metrics, and the earlier memo carries a retraction.

The lesson is the same one the concentration measurement itself produced, applied to me: **source COUNT is
not diversity; the SHARE is.** 5 sources sounds fine until one holds two thirds. I wrote that sentence last
run and then failed to apply it to my own number in the same run.

## What survives, and it is the part that mattered

The two cells that clear the bar are exactly the two that carry the findings:

- **E. coli × gentamicin** — sens **0.523** against the frozen cell's 0.893, on 8 sources (4.25
  effective). This is the number the whole `rmt` thread rests on, and it is now produced by an arm that would have refused it had the cohort been
  concentrated.
- **E. coli × ciprofloxacin** — spec **0.988** against the frozen cell's 0.700, on 8 sources (4.32
  effective) vs the frozen cell's 2.

Nothing about the `rmt` conclusion changes. What changed is that the number now comes with a checkable
statement about the cohort that produced it.

## Namespace

Results land in `wiki/source_diverse_validation_<organism>_<drug>.json` — **never** in
`provenance_disjoint_validation_*`, which is the glob the report card's `load_scored()` reads. Writing
there would silently overwrite a frozen cell with a different number: the shared-key trap, guarded by test.

These **augment** the frozen cells. No frozen metric changes, no cell state changes.

## Honest limits

- The pool is 200 PD-labelled genomes with cached AMRFinder output; **only E. coli reaches the diversity
  bar at all**. Every Klebsiella cell in the pool is 1–5 sources. So this arm currently says something about
  E. coli and nothing about Klebsiella, which is a property of what happens to be cached here.
- Diversity is measured by BioProject. Two BioProjects from the same lab would count as two.
- `underpowered` and `source_concentrated` are different refusals and are reported separately — a cell can
  be diverse and too small, or big and too concentrated.

Reproduce: `uv run python scripts/source_diverse_validate.py` (offline, given the census + provenance
sidecars).


---

## Update (same day): the bar had a loophole, found by trying to use it

**The target was `campylobacter x ciprofloxacin`** — the most concentrated SCORED cell (100% one
BioProject, n=40) and the only flagged cell with no source-diverse replication. Sizing it first, before
spending any AMRFinder time, produced two findings.

### 1. Campylobacter cannot be meaningfully source-diversified from public data

All of NCBI-PD's Campylobacter holdings, after the leakage manifest, are **6 BioProjects** — and one holds
**79%**:

| BioProject | R | S | |
|---|---:|---:|---|
| PRJNA292664 | 617 | 2449 | |
| PRJNA292668 | 145 | 559 | |
| PRJNA560409 | 13 | 110 | ← **the source of the existing SCORED cell** |
| PRJNA562719 | 2 | 3 | |
| PRJNA287430 | 0 | 2 | |
| PRJNA239251 | 1 | 1 | |

Excluding the cell's own project leaves **2 substantial projects plus scraps of 1–3 genomes**. A 20R/20S
cohort is constructible in the letter of the bar and not in its substance. **This is a data-availability
wall, not a code wall** — no amount of building fixes it, and it is why no replication was attempted.

### 2. My own bar would have accepted that cohort — so the bar changed

A cohort of `18/18/2/1/1` passes **both** shipped rules (5 projects; largest 45% ≤ 60%) while being two
real projects wearing three tokens. Added **inverse-Simpson effective sources** (`MIN_EFFECTIVE_SOURCES =
3.0`), the same effective-N idiom the lineage layer already uses for clonality:

| cohort | projects | largest share | **effective** |
|---|---:|---:|---:|
| 18/18/2/1/1 | 5 | 45% | **2.45** |
| 12/10/8/6/4 | 5 | 30% | **4.44** |

**It fired on real data immediately.** `E. coli x meropenem` — 5 nominal sources, 48% share, **2.42
effective** — moved from `underpowered` to `source_concentrated`. Refusals went 6 → 7 of 10. The two
scored cells (E. coli cipro and gentamicin, 8 sources each) are unaffected, so **no finding changes**.

I added this rule because my own gate would have accepted a cohort I would not defend. Count and share
are both necessary and neither is sufficient.
