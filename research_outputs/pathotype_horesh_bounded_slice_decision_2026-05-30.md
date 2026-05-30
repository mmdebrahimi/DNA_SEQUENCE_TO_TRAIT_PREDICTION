# Pathotype / Horesh — Bounded-Slice Decision + Concord Memo — 2026-05-30

> Resolves the two-machine disagreement over the Horesh build (laptop "off critical path" vs workhorse "build it now"). Adjudicated via adversarial review. **Both stances were partly wrong; this is the agreed path.** Doubles as the handoff back to the workhorse (Precision 7780).

## Verdict in one line

**Proceed with Horesh — but bounded, provenance-split, and honestly reframed. Send Gate B now.**

## What each side got right / wrong

- **Laptop was wrong** that Horesh work is off the critical path. The ingestion + caller + materializer + manifest code is mandatory plumbing regardless of class coverage. Building it on real rows is progress, not motion.
- **Workhorse is right** to build the pipeline on available data rather than wait for the perfect cohort.
- **Workhorse overclaims** if a 3-class Horesh result is read as a "validated predictor." The tool is a **deterministic VirulenceFinder-marker resolver**, and 52.2% of Horesh labels are themselves gene-rule-derived. Resolver-vs-marker-derived-labels measures **rule-implementation fidelity, not prediction skill.**
- **Nuance that rescues ExPEC specifically:** clean ExPEC labels are assigned by **isolation source** (blood/urine), independent of virulence markers. So "do these markers track extraintestinal isolation site" on the clean subset **is** a genuine genotype→phenotype test — not circular. Circularity bites only the dropped `(predicted)` rows.

## The decision (3 moves, in parallel)

### 1. Bounded vertical slice — NOT the full cohort
- Re-pick the smoke + first cohort to **direct-WGS-accession** rows (Salipante-style; `JSIS00000000`/upec-276 is a valid first row). **Skip Kallonen lane-ID re-assembly entirely** for now.
- Set explicit **exit criteria before scaling**, recorded in the run summary:
  - accession-resolution yield (% of selected rows that fetch by accession)
  - caller completion rate
  - provenance completeness per row
  - count of *genuinely independent* labels per supported class (isolation-source or curated)
  - runtime + storage cost per N
- Build only enough substrate to run this slice. Decide full-cohort AFTER the slice reports.

### 2. Provenance-split evaluation — two separate numbers, never collapsed
- **Resolver conformance:** agreement vs the gene-rule-derived (`(predicted)`) labels. Expect high; it measures implementation fidelity, NOT skill.
- **External validity:** agreement vs genuinely independent labels (ExPEC by blood/urine isolation; EPEC/ETEC by dedicated source-study). This is the real signal.
- Report both with explicit labels. Never publish a single blended "accuracy."

### 3. Reframe the v0 claim — honest scope
- v0 is an **"auditable marker-based pathotype-compatibility resolver with abstention,"** NOT a "general E. coli pathotype predictor."
- Supported v0 classes on Horesh substrate: **ExPEC/UPEC-compatible, EPEC, ETEC** (H1-passing, H2-clearing). 
- EAEC, commensal, clean-EHEC: **documented v0 scope-limit** (Horesh N=2/2/47). Keep them as explicit abstention/UNCLASSIFIED outputs or behind a flag — do NOT present them as supported predictions. Acquire later via Whittam DECA / von Mentzer / NCBI host_disease.

## The action that gates everything — SEND GATE B NOW

Gate B outreach is prepared (`research_outputs/pathotype_gate_b_send_kit_2026-05-29.md`) and **unsent, 0 replies.** It is the cheapest, highest-leverage move and the full heavy build is gated on it (≥2 credible yes / 60 days). It also directly tests the real open question: do users want a marker-resolver-with-abstention, or do they require prediction beyond known rules? **Send it today** (user-side action; 5 named contacts verified-then-fire).

## What does NOT change

- Gate A PASS stands.
- The reusable engine (ingest→markers→resolve→abstain→audit→CLI JSON) is sound and reused.
- Heavy full-Horesh build stays gated on Gate B + a passing slice.

## For the workhorse, concretely
1. Change the manifest's first-N selection to accession-bearing source studies (Salipante etc.); leave Kallonen rows out of the smoke.
2. Add a provenance column to the run summary; emit conformance vs external-validity separately.
3. Keep EAEC/commensal/EHEC out of the *supported* surface; abstain.
4. Don't authorize the full materialize until (a) the slice exit-criteria pass AND (b) Gate B ≥2 yes.

## Open questions to resolve during the slice
1. Per-class, which exact labels are independent of every marker the resolver uses?
2. How many direct-WGS-accession rows per class after dedup + QC + lineage-leak control?
3. What user decision does v0 serve — surveillance triage / clinical interp / research annotation / pipeline step? (Gate B answers this.)
