# The serotype fix survives a lineage-disjoint split — the named limit is closed

Every prior number for this fix carried the same caveat, stated in three consecutive memos: the held-out
half was held out **by isolate, not by lineage**. Two accessions can be near-identical genomes, so an
isolate-disjoint split does not prove the rule generalizes past the clones it was measured on. This run
removes that.

## The split

Isolates are grouped by **multi-locus sequence type** (E. coli Achtman 7-locus: `adk fumC gyrB icd mdh
purA recA`), and **whole STs** are assigned to TRAIN or TEST by deterministic hash. **No sequence type
appears on both sides.**

| | |
|---|---|
| distinct STs | **145** |
| TEST | 71 STs / **169 isolates** |
| TRAIN | 74 STs / 229 isolates |
| excluded (no complete 7-locus ST) | **2** |
| largest single ST | 53 isolates (13.3% of typed) |

ST was chosen over Mash for two reasons: Docker was down, and ST is the canonical, interpretable lineage
unit for this organism. The second reason is the real one — a reader can check an ST.

## Result — TEST half, 169 isolates, no shared lineage

| rule | O | H |
|---|---|---|
| coverage-only (the old deployed rule) | 0.9797 | 0.7844 |
| **identity-primary (the fix)** | 0.9730 | **0.9461** |
| **gain** | **−0.0068** | **+0.1617** |

**28 H calls genuinely differ between the two rules on the test half**, so the comparison is powered —
this is not a null produced by two rules that never disagreed. The script refuses with
`UNDERPOWERED_RULES_NEVER_DIFFERED` if that count is zero, precisely so a vacuous agreement cannot be
read as a survival.

**Verdict: `SURVIVES_LINEAGE_DISJOINT`.**

## The gain did not shrink — and that is the informative part

| split | H gain |
|---|---|
| isolate-disjoint (previous run) | +0.1057 |
| **lineage-disjoint (this run)** | **+0.1617** |

The usual pattern when a limit like this is closed is that the effect *shrinks* — which is exactly what
happened earlier in this same track when the unblinded fix went from +0.155 to +0.106 on replication.
Here it did not. **Lineage overlap was not carrying the earlier result.**

**Do not over-read that as "the fix is bigger than we thought."** The two numbers are not estimates of
one quantity: different test sets (250 isolates vs 169), different composition. The supportable claim is
**"survives lineage-disjointness"**, and nothing more.

The small **O-axis cost replicates too** (−0.0068 here, −0.0086 isolate-disjoint). A stable trade-off
across two independent splits is more informative than either number alone: the fix buys a large H gain
for a small, consistent O cost, and that is now measured twice rather than assumed.

## Honest limits

- **ST is a lineage unit, not a clonality guarantee.** Two different STs can still be related. This is
  strictly stronger than isolate-disjoint and strictly weaker than a full phylogenetic split.
- **Isolates without a complete 7-locus ST are excluded, not pooled** (2 here). Pooling them into one
  bucket would have manufactured a fake lineage; excluding them is honest but means the test half is not
  the whole cohort.
- One cohort (400 labelled E. coli), one antigen DB build.
- The largest ST holds 13.3% of typed isolates, so the split is not perfectly balanced across lineages —
  no single clone dominates, but the STs are not equal-sized either.
- **Mash was not used.** A Mash-distance split would be the stronger test and remains available when
  Docker is up.

## Reproduce

```bash
uv run python scripts/serotype_lineage_disjoint.py
```

Needs blastn + cached assemblies. Frozen AMR surface byte-unchanged — this is a typing cell.
See [`serotype_lineage_disjoint_2026-09-04.json`](serotype_lineage_disjoint_2026-09-04.json).
