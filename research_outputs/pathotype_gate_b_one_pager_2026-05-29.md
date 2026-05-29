# Pathotype Gate B One-Pager - 2026-05-29

## Ask

Would you pilot a local audit/QA companion for **E. coli pathotype interpretation** within 60 days?

## What it is

A local-first CLI that sits next to existing pathotype workflows and returns:

- deterministic pathotype-compatible calls
- explicit `HYBRID`, `AMBIGUOUS`, `UNCLASSIFIED` states
- provenance / audit packet
- honest output surface instead of forced overconfident labels

Current target classes:

- `EHEC_COMPATIBLE`
- `STEC_NON_LEE`
- `tEPEC_COMPATIBLE`
- `aEPEC_COMPATIBLE`
- `ETEC_COMPATIBLE`
- `EAEC_COMPATIBLE`
- `UPEC_COMPATIBLE`
- `HYBRID`
- `AMBIGUOUS`
- `UNCLASSIFIED`
- `COMMENSAL_LOW_MARKER_BURDEN`

## What it is not

- not another opaque AI classifier
- not a hosted SaaS platform
- not a replacement pitch for your whole pipeline
- not claiming to beat ECTyper broadly today

## Why this might matter

The wedge is narrow:

- existing pathotype calling can still leave awkward edge cases
- hybrids and partial-marker cases are operationally annoying
- labs often need a more honest surface than “force a single clean label”
- auditability and provenance matter more than one more black-box score

## Current proof points

- Gate A manual sanity passed on the workhorse machine
- pinned local runtime works
- raw outputs are usable
- manual decision table is computable from real outputs
- 5-strain frozen panel behaved as expected:
  - EHEC
  - ETEC
  - EAEC
  - UPEC-compatible
  - commensal

## Current maturity

- early pilot / evaluation stage
- not broad production-ready distribution
- no heavy Horesh cohort build started yet
- build is intentionally gated on real user demand

## Ideal pilot shape

- local evaluation only
- small benchmark or known-case panel
- focus on ambiguity / hybrid / abstention usefulness
- 2-6 week evaluation
- explicit go / no-go criteria

## The one question

If this were available as a local CLI with clear outputs and provenance, would you pilot it in your workflow within 60 days?

Accepted answer types:

- yes
- yes, with conditions
- no

## What would help most

If “yes, with conditions,” the most useful extra detail is:

- what exact use case would you test first?
- what current tool/output is failing you?
- what would make this worth integrating?
