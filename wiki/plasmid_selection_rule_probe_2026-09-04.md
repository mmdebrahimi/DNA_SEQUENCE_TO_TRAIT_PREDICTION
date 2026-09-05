# The last coverage-only caller is structurally inert — and no cohort was needed to show it

`plasmid` was the third caller found carrying the coverage-only selection pattern, and the only one
still untested. It had been deferred twice on the grounds that it "needs a cohort". **It did not** — not
for the first question.

## The distinction that unblocked it

Labels are needed to say which ordering is **better**. They are not needed to say whether the ordering
changes the **answer at all**. If it never does, the concern is moot and no cohort is required.

## The structural prediction — tested, not assumed

The plasmid caller reports a **set of replicon families**: every family with at least one called allele.
The coverage comparison only decides *which allele represents a family it has already decided to
report*. So the reported set should be invariant to the ordering, and only secondary fields should move.

That is a first-principles claim about the code, so it was executed rather than published.

## Result — 40 assemblies, 446 replicon calls

| | |
|---|---|
| assemblies whose **replicon set** differs between orderings | **0** |
| assemblies where a family's **best allele** moved | **12** |

**The prediction survived.** The rule genuinely fires — it changes the representative allele on 12 of 40
assemblies, so the probe demonstrably exercised both orderings — but it **never changes the primary
output**.

**Verdict: `STRUCTURALLY_INERT_FOR_THE_REPORTED_SET`. No fix, no cohort, no change.**

## The caveat that is not buried

**Secondary fields do move.** `best_allele` and the identity/coverage printed for a family change on 12
of 40 assemblies. A consumer that reads those rather than the replicon set **is** affected, even though
the set is not. That is a real if narrow exposure, and it is the reason this is "inert for the reported
set" rather than "inert".

## This closes the sweep across all four callers

| caller | selection | outcome |
|---|---|---|
| `ktype` | identity-primary | already correct — no action |
| `serotype` | coverage-only | **live defect, fixed** — H accuracy 0.770 → 0.926, replicated held-out at +0.106 |
| `pneumoserotype` | coverage-primary | probed: flips 1 call in 25 and that flip is wrong under **both** orderings — not transferred |
| `plasmid` | coverage-only | **structurally inert** for the reported set — not transferred |

One real bug out of three suspects, and the other two were closed on evidence rather than on assumption.
**The pattern-match was a lead, not a diagnosis** — which is exactly why each was measured separately
instead of propagating the fix.

## Honest limits

- This answers *whether* the ordering changes the answer, **not which ordering is better**. That still
  needs wet-lab replicon labels (PCR-based replicon typing), which are rare in public metadata — but it
  is only needed if the set actually moves, and it does not.
- Assemblies are E. coli and Salmonella borrowed from other cohorts, **not a plasmid-focused set**;
  replicon content is whatever those genomes happen to carry.
- One DB build (`enterobacteriales`). A different allele set could behave differently.
- 40 assemblies. A larger set could surface a rare case where two families' calls interact, though the
  structural argument says they cannot.

## Reproduce

```bash
uv run python scripts/plasmid_selection_rule_probe.py --limit 40
```

Needs blastn + cached assemblies. Frozen AMR surface byte-unchanged.
