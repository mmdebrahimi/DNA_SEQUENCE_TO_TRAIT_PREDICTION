<!-- doc-status: active -->
<!-- measured-on: 2026-11-28 -->
<!-- recheck-after: 2026-12-05 -->

# The regulatory-metric arm is 5/6 built, correct where built, and entirely uncommitted

A READ-ONLY verification of the working-tree arm produced by
`plans/Build_The_Metric_And_Its_Frozen_Bar_Before_Any_Regulatory_Model_Plan/`. Nothing in that arm
was modified, committed, or staged by this pass. It is checked here because a downstream plan
(`plans/EncodeMultiBiosampleRegulatorySignal/`) already treats its acceptance bar as frozen, and a
frozen bar that exists only as an untracked file cannot be shown to have been frozen.

## Ownership — not this session's work

No transcript under `~/.claude/projects/` authored these files. The only transcript naming them is
this session's, and it names them solely as read-back `git status` output: **zero** `Write` calls
match `regulatory_*`, and this session has **no entries between 2026-11-26T03:00Z and 2026-11-28**
while the files were written 2026-11-26 16:54–17:02 local (21:54–22:02 UTC). `ListAgents` shows no
live peer in this working directory. So the arm is a sibling session's, left mid-flight two days ago.

## What IS done, and verified rather than assumed

- **Step 1 binds to the existing owner instead of reimplementing it.**
  `dna_decode/eval/regulatory_metric.py` imports `group_centered_association`,
  `group_centered_spearman` and `permutation_null` from `dna_decode/deconfound/deconfound.py`.
- **The subtle argument contract is RIGHT.** `regulatory_association` passes a RAW `y` with a
  RESIDUALIZED `x_resid` to `permutation_null`. That looks like a double-centering bug and is not:
  the owner's own docstring specifies exactly that pairing — it shuffles `y` within each group and
  re-centers internally, against a fixed pre-residualized feature.
- **The fixtures are two-sided, so the metric cannot pass by always returning zero.** One fixture
  builds a pure locus main effect (pooled `abs(rho) > 0.8`, within-locus `< 0.15`); a second builds
  genuine within-locus signal (within-locus `abs(rho) > 0.5`). A conditioned statistic that is
  always near zero is as uninformative as one that is always high, and these two fixtures move in
  opposite directions.
- **The floor refuses rather than scores:** one locus at 29 rows returns
  `UNSCORABLE_GROUP_FLOOR` with `within_locus: None`.
- **38 tests pass** across the arm's six untracked test files (7.26 s, this host, this session).
- **The bar is internally consistent, and notably does NOT carry the defect that is open authority
  call 1.** `frame_verdict.CHECK_ORDER` has nine outcomes against a frozen bar naming eight; this
  bar's four `outcomes` are all reachable from its four `check_order` rows, because the final row
  maps to a PAIR rather than to a single verdict.
- **The bar was frozen BEFORE the measurement, which is the correct order and is checkable by
  mtime:** bar 16:58:01, locus-capability probe 17:18:03, ENCODE cohort probe 20:01:09, all
  2026-11-26. Its `forecast` puts `p_indeterminate` at 0.1 and the probe then found 0 loci clearing
  the 30-row floor — a forecast being wrong after the fact is what a pre-registered forecast is
  FOR, and is not evidence the bar was written to fit a result.

## What is NOT done

- **Step 4 is built but UNRUN.** `scripts/kaggle/gcsfuse_random_read_probe.py` is a manual Kaggle
  probe and no result artifact exists anywhere in `wiki/`. Its `SINGLE_EXAMPLE_MIB_PER_S = 9.01` is
  the DERIVED requirement (9.27 MiB at 1028.3 ms), never a measured throughput. Every statement
  about the data path's adequacy still rests on the sequential arithmetic alone, exactly as the
  plan's own `[unverified]` risk flag says.
- **Step 5's memo was therefore written before the verdict it consumes**, and its filename is
  `wiki/regulatory_ingestion_scope_2026-11-26.md` where the plan contracts
  `wiki/regulatory_ingestion_scope_2026-11-25.md`. The plan's risk flag required a test pinning
  every numeric in that memo to a committed artifact; no such test is present, and `9.01` is a bare
  literal in both the memo and the probe module.
- **Step 6 never ran.** Both of its contracted files — `project_state/glm-umbrella-2026-10-07.md`
  and `soraya_runs/constraint-expiry.json` — are tracked and unmodified, so neither the decision
  state nor the constraint register records this arm. `.claude/execute-plan-state/` still reads all
  six steps `pending` with a null test baseline.

## The one consequential exposure

The bar declares its own provenance (`frozen: 2026-11-26`,
`frozen_before: any human regulatory-activity model score`) while living in an untracked file. An
untracked file has no history: the freeze is unverifiable to anyone who did not measure the mtimes
above, and it can be edited to fit a later result with no trace. The mtime ordering recorded here is
precisely the evidence a commit would make durable and a `git clean` would destroy — and
`git clean -fd` is already forbidden in this repo, which is currently the only thing protecting
two days of work.

**Committing it was NOT done here and is deliberately left as an authority call**, because making
another session's acceptance bar durable is an acceptance-bar act, not housekeeping.

## Scope of this verification

Read-only, one host, one pass. It checks the arm against its OWN plan contract and against the
owner module's semantics; it does not re-derive the statistic's correctness, does not run the Kaggle
probe, and says nothing about whether the metric is the right instrument. It also does not re-open
the substrate question: `wiki/human_locus_capability_probe_2026-11-26.json` already measured
`WITHIN_LOCUS_SUBSTRATE_UNAVAILABLE` (7,321 shared elements, exactly 2 observations each, **0** loci
clearing `min_n=30`), so on the only human substrate on disk this bar can currently return only
`INDETERMINATE_GROUP_FLOOR`. That is a substrate fact already recorded in `NEXT.md`, not a defect of
the bar.
