# The serotype fix replicates on held-out isolates — smaller, as an unblinded choice predicts

The identity-primary fix was found by inspecting the discovery cohort's failure pattern and then
applied to those same isolates: H accuracy 0.770 → 0.926. That is an **unblinded** result — the change
was chosen after seeing what it would repair — so its effect size is exactly the kind that shrinks on
replication. **The prediction was written down before this run**, so it could fail.

## Result — 250 held-out isolates, disjoint by accession from the discovery set

| | coverage-only | identity-primary | delta |
|---|---|---|---|
| **H accuracy** | 0.8171 (201 hit / 45 miss) | **0.9228** (227 / 19) | **+0.1057** |
| **O accuracy** | 0.9612 (223 / 9) | 0.9526 (221 / 11) | **−0.0086** |
| no-call rate | — | — | **0.0000** |

**Verdict: CONFIRMED** — clears the pre-registered +0.05 bar on isolates the fix was never chosen from.

Cohort: 250 called, **47 BioProjects**, largest-source share 0.165 (clears the 0.60 diversity bar), 108
distinct O-groups, 31 discovery-overlapping accessions explicitly excluded.

## Three things worth stating plainly

**1. The gain shrank, and that confirms the caveat rather than undermining it.**
Discovery +0.155 → held-out **+0.106**, a shrinkage of 0.050. That is precisely the regression an
unblinded selection predicts. **Quote +0.106, not +0.155.** The fix is real; the discovery number was
optimistic, exactly as its own memo warned.

**2. The O axis is slightly WORSE under identity-primary** — −0.0086, two isolates. Reported rather
than buried. The net is strongly positive (**26 H misses fixed, 2 O misses introduced = 24 net calls
corrected**), but the trade-off exists and someone tuning this later should know it is not free.

**3. The secondary prediction was met exactly, and it is the mechanistic one.**
The claim was that the rule decides *which* allele wins, not *whether* one does — so resolution should
be untouched. Observed no-call delta on 250 isolates: **0.0000**. That corroborates the *mechanism*,
not merely the effect. A fix that changed resolution would have meant the story was wrong even if the
accuracy moved.

## Design choices that make the comparison fair

- **Both rules scored from ONE blastn pass per isolate.** Only the sort key differs, so no run-to-run
  variation enters the comparison.
- **Disjointness enforced two ways**: a different cohort seed (77 vs 23) *and* explicit exclusion of
  every accession in the discovery checkpoint. Overlap was measured (31), not assumed.
- **The prediction is stamped into the artifact** (`preregistered`), so a reader can check it was not
  adjusted after the numbers landed.

## Honest limits

- **Held out by ISOLATE, not by lineage.** Two distinct accessions can still be near-identical genomes,
  so this is not a phylogenetically independent replication.
- **Still no in-silico incumbent** for E. coli (PD leaves `computed_types` null), so the *absolute*
  accuracy stays weakly anchored. What this run tests is the **rule comparison**, which is internally
  controlled.
- Shared blastn parameters and DB coverage limit both rules equally — this cannot detect a defect that
  affects them both.
- O-only labels score the O axis alone, so the two axes have different denominators.

## Reproduce

```bash
uv run python scripts/serotype_heldout_confirm.py --target 250
```

Needs blastn + network. Frozen AMR surface byte-unchanged.
