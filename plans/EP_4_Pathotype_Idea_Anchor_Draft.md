# EP-4 Pathotype `/idea-anchor` Draft

> Pre-staged `/idea-anchor` candidate sentence + scoping context for E. coli
> pathotype prediction — the recommended first non-AMR phenotype per
> `plans/EP_4_Non_AMR_Phenotype_Candidates.md` Rank 1.
>
> NOT yet fired. Discovery-lane deliverable per `plans/Two_Machine_Operating_Contract.md` §8.

**Status:** DRAFT 2026-05-26.
**Trigger conditions:** fire `/idea-anchor` with the sentence below ONLY after:
1. Cef audit-aware closeout pushes to origin (Codex side completes Artifacts 2-5)
2. EP-1.5 architecture decision settles (if it's needed — pathotype's concentrated mechanism may skip this gate per the EP1+EP2 2026-05-17 finding)

**Authorship:** Claude (planning artifact per contract §3 source-of-truth ownership).

---

## Ready-to-paste `/idea-anchor` invocation

```
/idea-anchor E. coli pathotype prediction from a user-supplied genome assembly: an open-source CLI tool that takes FASTA + optional GFF3 and emits a multiclass pathotype call (EPEC / EHEC / ETEC / UPEC / EAEC / commensal) with audit-grade provenance, including which acquired virulence-gene clusters drove the call (stx for EHEC; LEE for EPEC; afa/papC for UPEC; ETEC enterotoxins) and a side-by-side comparison against CGE VirulenceFinder gene-call output.
```

Copy that single line into a Claude Code prompt to run `/idea-anchor`. The skill will then produce its 6-section output (Verbatim Input → Formal Rephrase → Fundamental Clarifications → Current Assumptions → Blunt Opinion → Recommended Next Step).

---

## Why this sentence

The sentence is shaped to mirror the v1 "Honest AMR resistance predictor" framing already locked in `wiki/decoder_v0_ux_and_success_criterion.md`, applied to pathotype. Specifically:

| Frame element | AMR v0 instance | Pathotype EP-4 instance |
|---|---|---|
| Input substrate | E. coli FASTA + GFF3 | Same — no new input channel |
| Output shape | Per-drug binary R/S | Per-genome multiclass pathotype |
| Provenance | Top-K gene-level ISM + Tier 1-5 catalog | Top-K virulence-gene clusters + pathotype catalog |
| Audit framework | mechanism × MIC × opacity merge | mechanism × pathotype × opacity (extension) |
| Honest output | `attribution_scope_confidence` + `suspend_gate_fired` | Multiclass extension of the same fields |
| Comparison benchmark | AMRFinderPlus side-by-side | **CGE VirulenceFinder side-by-side** — direct analog |
| Mechanism shape | Concentrated (cipro QRDR / cef β-lactamases) | **Concentrated** (pathotype = acquired-gene clusters) |

Concentrated mechanism is the load-bearing claim: it predicts that the existing
NT mean-pool + XGBoost architecture (which PASSED on cipro + cef and FAILED on
tet per the 2026-05-17 cross-drug architectural finding) should generalize.

---

## Why NOT a different candidate

Per `plans/EP_4_Non_AMR_Phenotype_Candidates.md` ranked verdicts:

- **Growth rate** — distributed mechanism (would FAIL mean-pool; analogous to tet failure mode); no clinical buyer
- **Biofilm formation** — mixed mechanism (borderline; needs EP-1.5 to settle architecture); smaller cohorts
- **Plasmid stability** — too AMR-adjacent (doesn't pay off the "jump out of AMR" framing)
- **Lactose fermentation / metabolic auxotrophy** — too easy (gene-presence alone solves it; no AI value-add)

Pathotype clears all 6 selection criteria: paired data exists (EnteroBase
>200K labeled isolates), mechanism is concentrated, measurement is categorical
+ well-defined, clinical relevance is HIGH (clinical microbiology, veterinary,
food safety, outbreak investigation), cohort size potential is N ≥ 500
trivially, and architectural compatibility is minimal (`--drug` → `--phenotype`,
binary → multiclass classifier).

---

## Pre-anticipated answers to `/idea-anchor`'s clarification questions

The `/idea-anchor` skill asks up to 3 fundamental clarification questions.
Likely ones + the user's likely position:

| Likely question | Likely user position |
|---|---|
| Who is the primary actor? | Clinical microbiologist OR public-health researcher running pathotype calls on E. coli isolates from clinical / food / veterinary specimens |
| What exact outcome should improve? | Replace manual VirulenceFinder gene-call inspection with a single CLI that emits a calibrated pathotype probability + audit citation, in the same UX as the v0/v0.1 AMR decoder |
| Is this feature-specific or system-wide? | Feature-specific — pathotype is a NEW phenotype added to the existing decoder surface, NOT a system-wide refactor. v0.2-class change. |

---

## Pre-anticipated assumptions (`/idea-anchor` Step 3)

These will land in the assumptions section regardless of what the user types:

1. EnteroBase E. coli pathotype labels are usable as ground truth at N ≥ 500
2. Pathotype is a CATEGORICAL phenotype (multiclass — EPEC/EHEC/ETEC/UPEC/EAEC/commensal). Binary commensal-vs-pathogen is a degenerate case if cohort assembly is too hard.
3. CGE VirulenceFinder is the legitimate side-by-side benchmark (direct analog of AMRFinderPlus for AMR — same shape, same downstream comparison discipline)
4. v0/v0.1 schema generalizes from `--drug` to `--phenotype` with minimal disruption (per EP-4 candidates memo "Architecture / spec implications")
5. NT mean-pool + XGBoost architecture transfers because pathotype mechanism is concentrated (per 2026-05-17 cross-drug architectural finding: PASSES on concentrated, FAILS on distributed)

---

## Blunt opinion (pre-emptive — what `/idea-anchor` Step 4 should surface)

**Strong candidate but 3 risks to name:**

1. **EnteroBase pathotype labels may be inconsistent.** E. coli pathotype is biologically fluid: strains carry pathotype-defining genes without being unambiguous pathotypes; commensal-to-pathogen transitions exist; some isolates carry multiple pathotype gene clusters simultaneously. The label noise problem from AMR (which led to the audit framework) reappears here in a different shape: not "wrong MIC" but "wrong category boundary." Mitigation: the audit framework's mechanism-class-bounded reasoning generalizes — flag opaque/borderline calls upfront.

2. **The audit framework needs multiclass extension.** v0's `noise_class` buckets (CLEAN_R / SUSPECT_S / OPAQUE_R / NOISY) were designed for binary R/S. Multiclass means more bucket categories: CLEAN_EHEC / OPAQUE_EHEC / SUSPECT_EPEC_with_EHEC_genes / etc. The structural shape generalizes but the catalog needs work.

3. **External-benchmark discipline NOT YET LIFTED on v0/v0.1.** Per `LESSONS_LEARNED.md` 2026-05-22, cipro v0 + cef v0.1 shipped without an AMRFinderPlus/RGI side-by-side benchmark. EP-4 explicitly proposes the VirulenceFinder benchmark in the candidate sentence. That's HARDER than v0/v0.1 — and the project's track record suggests the benchmark may slip to a "v1.0 follow-on" status the same way EP-1B did for AMR. Decision needs to be explicit before EP-4 ships: benchmark on day 1, or accept the same scope-limit-doc pattern as v0/v0.1?

---

## Recommended next step (`/idea-anchor` Step 6 from this draft)

When trigger conditions are met:
1. Paste the `/idea-anchor` sentence above into the Claude Code prompt
2. Run the skill's 6-section output
3. Then proceed to `/project-init` per planning-pipeline.md
4. `/project-init` will produce a refined goal hierarchy + Bellman frame; `/probe` + `/brainstorm` follow per the chain

**Time-to-fire estimate:** ~3 hr Codex compute for cef closeout + 30 min Claude planning chain = same-day capability once cef closes.

---

## What this draft is NOT

- NOT a `/idea-anchor` invocation — just a pre-positioned candidate sentence
- NOT a `/project-init` — that fires after `/idea-anchor` produces its 6-section output
- NOT an implementation plan — Phase 4 doesn't fire until cef closes + EP-1.5 settles
- NOT a commitment — the user can override the sentence with a different framing when the time comes; this is the starting point, not the lock

---

## Contract Locks

| Parameter | Enforcement target |
|---|---|
| Pathotype as Rank 1 EP-4 candidate (Candidate 2 in the candidates memo) | not lock-bearing — sourced from `plans/EP_4_Non_AMR_Phenotype_Candidates.md` which is the source-of-truth; this draft cites |
| Concentrated mechanism shape claim | not lock-bearing — sourced from `wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md` |
| CGE VirulenceFinder as the side-by-side benchmark choice | not lock-bearing — design choice; user may override at `/idea-anchor` time |

Per the two-machine contract §4, this draft does not introduce new locks; it
cites existing locks in their source-of-truth memos.

---

## Bottom line

When Codex closes cef, paste the `/idea-anchor` sentence above + run the
planning chain. The rest of the context (mechanism shape, dataset, ranking
rationale, architectural risk) lives in
`plans/EP_4_Non_AMR_Phenotype_Candidates.md`. This draft is the bridge between
that scoping memo and the actual `/idea-anchor` invocation — no new planning,
just packaged for action.
