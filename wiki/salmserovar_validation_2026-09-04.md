# The serovar caller met a wet-lab label for the first time, and lost to the tool it wraps

**`typing:Salmonella:salmserovar` shipped at `FAITHFUL_TO_TOOL` — checked against the reference method,
never against reality.** Measured now against a free, independent, wet-lab label on 200 isolates:

| comparison | hit | miss | no-call | accuracy |
|---|---|---|---|---|
| **our caller vs wet-lab label** | 99 | 42 | **59** | **0.702** |
| **in-silico incumbent vs the same labels** | 184 | 15 | 1 | **0.925** |
| our caller vs the incumbent (the old number) | 98 | 43 | 59 | 0.695 |

**Delta −0.222, with a 29.5% abstention rate.** The cell is worse than naive use of the tool it mimics.

That delta is the whole point of running both. Our 0.702 would look respectable in isolation; against an
incumbent scoring 0.925 **on the same isolates** it is a clear underperformance. This is the recorded
wrapper-vs-tool lesson producing exactly the result it was written for.

---

## Why this label is trustworthy — and it is *measured*, not asserted

Salmonella serovar is unusual and valuable: its gold standard is **slide agglutination**, a wet-lab
antisera reaction, not a computation. That is the free, independent, isolate-level label the AMR track
never had. But public `serovar` strings are a mixture — some are agglutination results, some are
SeqSero2/SISTR output pasted into a metadata field. Scoring against the latter is scoring the tool
against itself.

**The incumbent's own score is the circularity probe.** Published in-silico-vs-agglutination accuracy is
~0.95. Had these labels been copied from the tool, the tool would score ≈1.000. It scores **0.925** —
consistent with genuine wet-lab labels. Two further supports:

- **1,570 PD rows** carry a serovar while the in-silico caller returned an unresolved serotype
  (`I -:-:-`). Those labels demonstrably did not come from that tool.
- The cohort is restricted to reference public-health labs (CDC / PHE / FDA / USDA-FSIS / state health
  departments) whose routine workflow includes traditional serotyping.

**Cohort:** 200 isolates, 74 distinct serovars, 29 BioProjects, **largest-source share 0.125** — clears
this project's own 0.60 source-diversity bar comfortably, applied at construction rather than discovered
afterwards.

## Diagnosed failure modes — specific, not "it underperforms"

> **CORRECTION, same day.** The claim below that phase-2 flagellin is "the single largest defect" is
> **wrong, and it was measured wrong by counting the wrong thing**. Counting formulas that *end* in `-`
> conflates "phase 2 is genuinely absent" with "nothing resolved on the H axes at all" — `4:H?:-` has an
> empty H2 only because H1 failed first. Partitioning by the **first axis that actually failed**
> ([`salmserovar_nocall_anatomy_2026-09-04.json`](salmserovar_nocall_anatomy_2026-09-04.json)):
>
> | cause | n | share |
> |---|---|---|
> | **O antigen unresolved** | **21** | **35.6%** |
> | H1 (phase-1 flagellin) unresolved | 16 | 27.1% |
> | O:H1 valid — only H2 blocks it | 13 | 22.0% |
> | O:H1 called but pair absent from table | 8 | 13.6% |
>
> **Only 22% is reachable by a phase-2 fix**, and the obvious such fix — resolving on O:H1 alone when
> that pair is unique — has **measured headroom of zero**, because the H2-blocked formulas are precisely
> the *ambiguous* ones, which is why they need H2. The real priority is **O-antigen DB coverage**, whose
> cost is data engineering on the wzx/wzy allele set, not a code change.

Of the 59 no-calls:

- ~~**33 have an empty H2 (phase-2 flagellin).** This is the single largest defect.~~ **Superseded — see
  the correction above.** The count of trailing `-` is real; the causal attribution and the priority
  ordering were not.
- **O-antigen unresolved** (`O?:r:1,5` for Infantis, `O?:f,g:-` for Rissen) — allele coverage gaps.
- **O-antigen mis-grouped** — Typhi called `1,3,19` (should be 9,12); Enteritidis called `9,46` (should
  be 9), which resolved to the wrong serovar "Hillingdon".
- **A malformed antigen name leaking from DB construction**: `22-gene2:z:1,6`.

The H-antigen calls were frequently correct where O failed, so **the O-antigen axis and the H2 phase are
where the work is** — not the formula lookup, which resolved correctly whenever given a complete formula.

## Fairness of the comparison

Equivalence is decided by `dna_decode.salmserovar.equivalence` and applied **identically to both
callers**, so no leniency can favour either. It grants only (a) documented notation normalisation —
`Typhimurium var. 5-` ≡ `Typhimurium`, `4,[5],12:i:-` ≡ `I 4,5,12:i:-`, monophasic Typhimurium ≡ its
formula — and (b) formula→name resolution through the committed W-K-L table. It refuses fuzzy matching:
`Newport` vs `Newbrunswick` is a miss, as it should be.

`no_call` is counted **separately from `miss`** throughout. A caller that abstains is not a caller that
is wrong, and merging them would hide a 29.5% abstention behind an error rate.

## Two defects found in this work itself

- **A quote-aware parsing bug, same class as the documented `AST_phenotypes` one.** PD's `computed_types`
  is comma-separated *and* Salmonella antigenic formulas contain commas, so a naive split shredded
  `I 4,[5],12:i:-` into `I 4` and manufactured disagreements that were not real. Fixed before any number
  was reported; agreement moved 0.848 → 0.863 on inspection alone.
- **`pending` was passing the label filter** and would have scored as a serovar miss.

Both were caught by *reading the disagreements*, not by trusting the rate.

## Honest limits

- Per-isolate agglutination provenance is **unprovable**; the reference-lab filter is a judgment. Residual
  circularity is **bounded, not eliminated** — but since both callers are scored against the same labels,
  contamination inflates both, so **quote the delta with more confidence than the levels**.
- A per-serovar cap (max 12) deliberately flattens natural prevalence. This is per-isolate accuracy on a
  diverse mix, **not** population-weighted accuracy — real-world Typhimurium/Enteritidis dominance would
  change the headline in either direction.
- Reference-lab restriction biases toward clinical/regulatory isolates.
- One organism, one trait. Says nothing about the other five `FAITHFUL_TO_TOOL` typing cells except that
  the same test is now cheap to run on them.

## Reproduce

```bash
uv run python scripts/build_salmserovar_cohort.py     # -> data/salmserovar_cohort.tsv + provenance JSON
uv run python scripts/salmserovar_validate.py         # needs blastn + the built antigen DB
```

The frozen AMR surface is byte-unchanged — this is a typing cell.
