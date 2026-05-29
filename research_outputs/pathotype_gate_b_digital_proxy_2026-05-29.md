# Pathotype Gate B Digital Proxy - 2026-05-29

## Scope

Digital-only proxy for Gate B.

Purpose:
- recreate the most relevant buyer / user / blocker personas
- estimate likely reaction to the current Phase 4 pathotype wedge
- decide whether real Gate B outreach is worth running

This is **not** Gate B pass/fail evidence.
Real Gate B still requires real external contacts.

## Product being tested

Current wedge being evaluated:

- local-first E. coli pathotype CLI
- deterministic resolver on top of pinned VirulenceFinder / ETECFinder-style calls
- honest 11-class surface with abstention / ambiguous / unclassified states
- audit-aware packet rather than black-box confidence theater
- Gate A already passed on the workhorse

Current maturity:

- substrate sanity proven
- not yet a full Horesh cohort build
- not yet a trained classifier layer
- not yet a validated production distribution channel

## SME panel selected

### Tier 1 - invoke first

1. `external-early-adopter-champion`
   - best fit for “would an early genomics / public-health adopter pilot this despite rough edges?”
2. `technical-buyer`
   - integration / reproducibility / workflow / packaging gate
3. `user-buyer-enduser`
   - day-to-day bioinformatics / microbiology workflow reality check
4. `skeptic-blocker`
   - strongest counterweight against “this is novel so somebody will use it”

### Tier 2 - supporting lenses

5. `economic-buyer`
   - only relevant in a narrow sense here because this is closer to a low-cost pilot / OSS-style evaluation than a large enterprise purchase
6. `procurement-legal`
   - lower weight because local CLI reduces vendor/data-transfer surface, but still useful for regulated/public-health deployment realism

### Cross-cutting

7. `devils-advocate`
   - final audit of the simulated panel

## Simulated persona reactions

### 1. External Early-Adopter Champion

Verdict:
- `would-pilot-with-conditions`

Likely reaction:
- interesting if framed as a **QA / audit / ambiguity-resolution layer**
- not interesting if framed as “yet another pathotype caller”
- strongest beachhead is a lab or pipeline maintainer already unhappy with opaque DEC outputs or hybrid-edge-case handling

What they would need:
- one narrow pilot
- fixed benchmark panel
- clear “better than current workflow” story
- founder/maintainer responsiveness

What makes them walk:
- vague claim of replacing ECTyper outright
- no concrete benchmark artifact
- unresolved runtime friction

### 2. Technical Buyer

Verdict:
- `requires-POC`

Likely reaction:
- deterministic local CLI is a positive
- pinned runtime and audit packet are positives
- still too early for broad adoption without:
  - reproducible packaging
  - benchmark substrate
  - clear output schema
  - workflow fit for existing nf-core / public-health pipelines

Primary concern:
- this must fit as a small auditable step, not a fragile bespoke stack

### 3. User Buyer / End-User

Verdict:
- `would-use` only if it reduces manual ambiguity review
- otherwise `would-work-around`

Likely reaction:
- if this helps explain weird hybrid / ambiguous E. coli pathotype cases, it is useful
- if it adds another command without eliminating an existing pain point, it gets skipped

Must-have:
- dead simple CLI
- TSV/JSON outputs that slot into existing notebooks/pipelines
- examples using real known strains
- no hidden cloud dependency

### 4. Skeptic / Blocker

Verdict:
- `block-pending-real-demand`

Likely reaction:
- pathotype calling itself is not greenfield
- ECTyper already occupies much of the obvious space
- biggest risk is building a technically clean tool that nobody installs

What would change their mind:
- 2+ real users explicitly saying they would pilot within 60 days
- evidence that the audit/abstention layer solves a real unserved pain

### 5. Economic Buyer

Verdict:
- `pilot-first`, very small scope only

Likely reaction:
- only compelling if positioned as a cheap way to reduce analyst time on ambiguous cases or improve confidence in pipeline outputs
- not compelling as a standalone big-budget tool

Good sign:
- local CLI / no large SaaS commitment

Bad sign:
- no quantified saved-time or decision-quality improvement yet

### 6. Procurement / Legal

Verdict:
- `approve-pending-redlines` for pilot-style use

Likely reaction:
- easier than a hosted platform because:
  - local execution
  - limited data-transfer surface
- still needs:
  - clear license posture
  - runtime provenance
  - documented external DB/source dependencies

This is not the primary blocker right now.

### 7. Devil's-Advocate audit

Verdict:
- `HOLD`

Three load-bearing claims behind the current path:

1. users care about the audit / abstention layer enough to switch behavior
2. the wedge is distinct enough from ECTyper to motivate adoption
3. a benchmark/QA framing is strong enough to create real pilots

What would falsify them:

1. target users say current outputs are already good enough
2. target users see no meaningful difference from ECTyper or existing practice
3. users say “interesting, but not worth integration effort”

Bottom line from the audit:
- simulated support exists for running real Gate B
- simulated support does **not** justify skipping real Gate B

## Synthesis

### What the digital panel says

- this is **not** a broad-market greenfield product
- it **does** have a plausible wedge
- the wedge is:
  - honest pathotype audit layer
  - ambiguity / hybrid / abstention handling
  - local deterministic reproducibility

### Best framing for real Gate B

Do **not** ask:
- “Would you install a new pathotype caller?”

Ask:
- “Would you pilot a local audit/QA companion for E. coli pathotype interpretation that gives explicit ambiguous/hybrid/unclassified outputs and provenance, instead of forcing overconfident labels?”

That is the highest-probability truthful wedge.

### Digital proxy verdict

- `PROCEED_TO_REAL_GATE_B`

Interpretation:
- real outreach is worth doing
- heavy build is still **not** justified yet
- Gate B remains the decision boundary

## Recommended real-contact archetypes

Best five contact types:

1. public-health lab lead handling E. coli pathotyping
2. GenomeTrakr / FDA-adjacent pipeline user
3. nf-core or similar workflow maintainer
4. clinical microbiology researcher working on DEC / hybrid E. coli
5. bioinformatics lead who already maintains local pathogen typing workflows

## Recommended outreach angle

Lead with:
- local deterministic CLI
- explicit ambiguity / hybrid states
- audit packet / provenance
- not another black-box model

Do not lead with:
- AI
- non-AMR expansion vision
- full future roadmap
- “better than ECTyper” unless benchmark evidence is in hand

## Recommended next move

1. other machine: build the real Gate B one-pager + outreach drafts
2. user / discovery lane: contact 5 targets
3. this machine: stay idle on heavy Phase 4 build until those responses land

## Final conclusion

Digital recreation of Gate B personas says:

- the idea survives first-pass buyer scrutiny
- the likely adoption wedge is narrow but real
- the project should advance to **real Gate B outreach**
- the project should **not** skip directly to the Horesh heavy build on digital simulation alone
