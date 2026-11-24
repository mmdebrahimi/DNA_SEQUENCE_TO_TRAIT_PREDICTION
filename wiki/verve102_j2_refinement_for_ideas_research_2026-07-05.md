# §5 refinement → ideas_research VERVE-102 fact-check (from the J2-verification session, 2026-07-05)

**Additive to** DNA-11's handoff-reply `wiki/handoff_reply_decoder11_reuse_2026-07-05.md` (commit `c90f8b4`,
origin/main). This supplies ONE nuance that session couldn't have: it was written before/without the J2
Phase-1 result, which THIS session independently verified today. **A refinement, NOT a correction — the
no-reuse verdict is unchanged and correct.**

## The nuance
DNA-11's reply §3-point-2 frames dna_decode's learned arm as *"a CLOSED NEGATIVE by design (embeddings
0-for-4) … dna_decode deliberately does not do learned variant-effect prediction."* Post-J2 that is
**imprecise** — it understates a real, now-verified capability:

- **Learned prediction WORKS at the molecular variant-effect layer.** J2 Phase 1 (2026-07-05): ESM2-650M
  zero-shot vs wet-lab deep-mutational-scan (ProteinGym, 201 assays) → **median |Spearman| 0.4914**,
  shuffled control 0.0136, matching the published field number. Independently re-derived from the committed
  per-assay JSON this session (`wiki/j2_phase1_esm2_result_2026-07-05.{md,json}`) — exact match, no drift.
- **Learned is confound-blocked ONLY at the complex-organismal / therapeutic-OUTCOME layer** (that is the
  0-for-5 de-confounded negative — a different layer from molecular variant-effect).

## Why the no-reuse verdict still HOLDS — with a sharper reason
VERVE-102 = predicting a **novel base-edit's effect on LDL** = the *therapeutic-outcome* layer. dna_decode's
learned arm **succeeds at molecular variant-effect (ESM 0.491) but is confound-blocked at therapeutic-outcome
prediction** — exactly the class VERVE-102 needs. So: no reuse, but the precise reason is "learned works at
molecular-effect, fails at therapeutic-outcome," NOT "learned is a blanket closed-negative."

## Suggested §5 edit
> dna_decode's learned arm is *not* a blanket closed-negative: ESM2-650M variant-effect vs wet-lab DMS =
> median Spearman **0.491** (verified 2026-07-05), i.e. learned prediction **works at the molecular
> variant-effect layer**; it is confound-blocked only at the *complex-organismal / therapeutic-outcome*
> layer. VERVE-102 (novel base-edit → LDL outcome) is the outcome layer → **no-reuse holds, sharper reason**:
> learned succeeds at molecular-effect, fails at therapeutic-outcome.

## Provenance / discipline
Read-only cross-session verification (R4 value-add lane). Written to the `mosfaer` branch — NOT DNA-11's
`main` tree (its live workstream, untouched). Every number re-derived from origin/main HEAD + the committed
J2 JSON this session. Routing to the ideas_research session is the user's (the human router between sessions);
this note is the durable, synced record.
