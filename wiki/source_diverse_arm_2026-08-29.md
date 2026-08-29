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

**It refused 6 of 10 cells.**

| cell | n | sources | largest share | verdict |
|---|---:|---:|---:|---|
| E. coli × ciprofloxacin | 131 | 8 | 31% | **scored** — acc 0.977 / sens 0.958 / spec 0.988 |
| E. coli × gentamicin | 131 | 8 | 32% | **scored** — acc 0.756 / sens 0.523 / spec 0.985 |
| E. coli × ceftriaxone | 117 | 6 | 36% | underpowered (S class < 20) |
| E. coli × meropenem | 71 | 5 | 48% | underpowered |
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

- **E. coli × gentamicin** — sens **0.523** against the frozen cell's 0.893. This is the number the whole
  `rmt` thread rests on, and it is now produced by an arm that would have refused it had the cohort been
  concentrated.
- **E. coli × ciprofloxacin** — spec **0.988** against the frozen cell's 0.700, on 8 sources vs the frozen
  cell's 2.

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
