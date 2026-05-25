# DNA Decoder v0 closeout handoff - 2026-05-24

## Executive decision

- Recommended status: `close_v0_and_move_to_v0_1`
- Confidence: `HIGH`

Reason:

- current cipro v0 surface is implemented, leakage-safe, audit-aware, and documented
- remaining cipro interpretability work is explicitly scope-limited for v0
- next meaningful capability jump is no longer hardening; it is scope expansion

## What v0 achieved

Current v0 is the **cached-strain cipro predictor**.

It can:

- load the trained leakage-safe cipro model
- score a cached E. coli strain already present in the NT embedding cache
- emit JSON plus markdown sidecar
- propagate audit framing from the merge-gate sidecar
- emit provenance including CV strategy and primary CV AUROC

Current reference model state:

- drug scope: `ciprofloxacin` only
- CV strategy: `leave_one_accession_out`
- CV grouping: `assembly_accession`
- primary CV AUROC: `0.8697607851155182`

Canonical release packet:

- `reports/dna_decoder_v0_release_candidate_2026-05-24.md`

Canonical real example:

- `reports/dna_decoder_v0_release_candidate_example_2026-05-24.json`
- `reports/dna_decoder_v0_release_candidate_example_2026-05-24.md`

## What we can honestly claim now

We can honestly say:

> We can predict **ciprofloxacin resistance for cached E. coli strains** with good held-out performance, leakage-safe evaluation, and audit-aware output.

We should **not** yet say:

> We can identify whether a strain is resistant to an antibiotic

That broader statement is still too strong because:

- current v0 is cipro-only
- current v0 is cached-strain only
- current v0 does not accept arbitrary new genome input
- attribution is exploratory, not consistently mechanistic
- multi-drug support is not closed

## Why closing v0 here is reasonable

The remaining open issues are no longer v0 blockers.

### Predictive path

- green enough for v0
- leakage issue mitigated
- canonical predict path hardened

### Interpretability path

- useful as exploratory support
- not strong enough for a mechanism-first promise
- explicitly scope-limited in:
  - `reports/cipro_v0_scope_limit_decision_2026-05-23.md`

### Product contract

- relocked to the implemented surface in:
  - `wiki/decoder_v0_ux_and_success_criterion.md`

### Release / handoff surface

- now exists as a concrete gold-path artifact
- real audited example produced from the current reference model

## v0 closeout statement

Recommended internal closeout statement:

> v0 is complete as an audit-aware cached-strain cipro predictor. It is not a real genome-input decoder and not a multi-drug resistance classifier. Interpretability is included as exploratory support with explicit scope limits. Future work should move to v0.1 rather than continuing to broaden v0.

## Recommended v0.1 direction

v0.1 should be framed as **scope expansion**, not more v0 cleanup.

Primary v0.1 candidate themes:

1. **Real genome-input decode**
   - accept previously unseen genome input rather than only cached strain IDs
   - likely the most meaningful product capability jump

2. **Second drug follow-on**
   - current intended follow-on is cef
   - only after deciding whether genome-input support should come first

3. **Operational packaging**
   - make the decoder easier to hand off or run in a stable environment
   - lower priority than the first real capability expansion

## Recommended v0.1 first question

This is the first decision v0.1 should answer:

> Do we want the next milestone to be **real genome-input decode for cipro**, or **cached-strain follow-on expansion to cef**?

My recommendation:

- choose **real genome-input decode for cipro** first

Reason:

- it closes the biggest gap between current implementation and the original product intuition of a DNA decoder
- it improves product credibility more than adding a second cached-strain drug
- it creates a more reusable interface for later multi-drug expansion

## Suggested v0.1 starting plan

1. Define the accepted genome-input surface.
   - FASTA only vs FASTA + annotations
   - whether on-the-fly CDS extraction is in scope

2. Decide where embeddings are created for new inputs.
   - local
   - approved remote runtime
   - hybrid

3. Define the first end-to-end v0.1 success gate.
   - one previously unseen cipro genome
   - full JSON + markdown output
   - same audit-aware provenance block

4. Only after that, decide whether cef enters the same phase or the next one.

## Key files for the next session

- `wiki/decoder_v0_ux_and_success_criterion.md`
- `reports/dna_decoder_v0_release_candidate_2026-05-24.md`
- `reports/dna_decoder_v0_release_candidate_example_2026-05-24.md`
- `reports/cipro_v0_scope_limit_decision_2026-05-23.md`
- `README.md`

## Bottom line

- Yes: it is reasonable to treat this as the end of v0.
- Yes: the next phase should be v0.1.
- No: we should not yet make the broad claim that we can identify antibiotic resistance for arbitrary strains in general.
- Yes: we can make the narrower claim that we can predict **ciprofloxacin resistance for cached E. coli strains** with an audit-aware output.
