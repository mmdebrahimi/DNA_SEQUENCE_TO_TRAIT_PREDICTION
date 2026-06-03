# Recommendation — minister run `2026-06-03-1309-ep4-v01-expec-recall`

## Stop: `blocked:user-only` (round 1) — interrogation receipts required before promotion

The minister GENERATED one candidate family to close gap `improve-expec-recall` and PARKED it on the
2-round interrogation gate. **You** must run `/interrogate-me` (user-only; Soraya cannot self-invoke) — twice
— against the candidate plan below, then re-invoke to promote + execute.

## Candidate family: `fam-per-gene-expec-scoring`

**Problem (current):** v0 resolver ExPEC recall = **9/12 = 0.75** on the H4 cohort. To clear ≥0.85 needs
**≥11/12**. Three ExPEC strains are missed because `SIDEROPHORES` and `CAPSULE_SERUM` are coarse single
clusters (present iff ANY member gene ≥0.80 coverage), so multi-gene extraintestinal burden is invisible:
- `JSMY00000000` → COMMENSAL (SIDEROPHORES + CAPSULE_SERUM, 0 strong adhesins)
- `JSPG00000000` → COMMENSAL (CAPSULE_SERUM only)
- `JSLG00000000` → AMBIGUOUS (P_FIMBRIAE only, 1 strong, 0 support)

**Proposal:** split `SIDEROPHORES` (iutA/iroN/fyuA/irp2/ireA/chuA/sitA) and `CAPSULE_SERUM`
(kpsMII/kpsM/iss/traT) into **per-gene presence**, derive a per-gene ExPEC score, and re-tune the resolver
so multi-gene iron-acquisition + serum-resistance burden lifts a strain over an ExPEC_COMPATIBLE
(LOW_CONFIDENCE) bar — **without** moving any non-CONFIDENT call into CONFIDENT (precision 1.0 invariant).
Re-tune against the cached coverage in `data/pathotype_cov_cache/` (instant reruns, no re-detection).

**Why not high-stakes:** local package code (`dna_decode/pathotype/{markers,resolve}.py`) + one new test +
one AC9 ledger row. No shipped skill, no global/synced state, no data migration, no auth → only the 2
interrogation receipts gate the promotion (no `/technical-plan` / `/brainstorm` required).

## Emit (run these, then tell me to resume)
```
/interrogate-me Candidate family fam-per-gene-expec-scoring: split SIDEROPHORES + CAPSULE_SERUM into per-gene presence and re-tune the v0 pathotype resolver to push ExPEC recall from 0.75 to >=0.85 on the 24-genome H4 cohort while keeping confident-supported-call precision at 1.0. Is per-gene support-marker scoring the right lever, or am I overfitting the 3 missed strains (JSMY/JSPG/JSLG) and risking EPEC false-positives?
```
```
/interrogate-me Round 2 on fam-per-gene-expec-scoring: pressure-test the precision-1.0 invariant. If per-gene ExPEC scoring promotes the 3 missed ExPEC strains, what stops it also promoting a borderline EPEC/commensal genome into a confident ExPEC call? Define the exact per-gene threshold + which calls stay LOW_CONFIDENCE so precision cannot regress.
```
