# Pathotype Gate B Response Rubric - 2026-05-29

## Purpose

Score real responses consistently.

This is the real Gate B interpretation layer.

## Pass rule

- `PASS = 2+ credible yes-to-pilot within 60 days`

## Credible yes definition

A reply counts as credible `yes` if it is:

1. from a relevant target archetype
2. explicitly willing to pilot or evaluate
3. within a near-term window
4. not purely polite encouragement

## Response classes

### `yes`

Counts toward pass.

Examples:

- “Yes, I would test this.”
- “Yes, send me a small panel/example.”
- “Yes, we could evaluate this in our workflow.”

### `yes_with_conditions`

Counts toward pass **only if** conditions are concrete and plausible.

Examples:

- “Yes, if it outputs TSV/JSON and runs locally.”
- “Yes, if you can show ambiguous/hybrid examples.”

Does **not** count if vague:

- “Maybe eventually”
- “Interesting, keep me posted”

### `no`

Does not count.

Useful subtypes:

- current tooling already sufficient
- too little distinction from existing tools
- integration burden too high
- wrong user / wrong workflow

### `no_response`

Does not count.

## Strongest positive signals

Rank these highest:

1. willing to pilot on real cases
2. asks for example panel / output schema
3. names a current pain point this would address
4. names a concrete integration slot

## Strongest negative signals

Rank these highest:

1. “we already solve this with current tooling”
2. “difference from ECTyper is not enough”
3. “interesting but not worth integration effort”
4. “would not use unless part of existing workflow”

## Decision table

| Outcome | Interpretation | Next move |
|---|---|---|
| 2+ credible yes | Gate B PASS | authorize Horesh substrate build |
| 1 credible yes | HOLD | more outreach or refine wedge first |
| 0 credible yes | FAIL | stop / pivot before heavy build |

## Logging template

Use this exact row shape in a tracker:

| contact_name | org | role | archetype | date_sent | date_replied | response_class | counts_toward_pass | key_reason | notes |
|---|---|---|---|---|---|---|---|---|---|

## Interpretation guardrails

Do:

- count only explicit pilot willingness
- prefer signal quality over politeness
- log exact phrasing where possible

Do not:

- count generic encouragement as yes
- upgrade a vague maybe into yes_with_conditions
- move the goalposts after replies land
