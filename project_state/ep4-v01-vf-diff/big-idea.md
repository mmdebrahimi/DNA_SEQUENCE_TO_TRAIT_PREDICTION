# Big Idea — EP-4 v0.1 VF side-by-side diff (completes the v0 spec)

<!-- mission-shape: in_place -->

> Bounded IN_PLACE minister mission. Scope: ship the LAST unfinished item of the v0 pathotype resolver
> spec — the **side-by-side comparison against CGE VirulenceFinder gene-call output** promised verbatim
> in the originating goal. This COMPLETES v0; after it lands, v0 is frozen as shipped+scope-capped and the
> project's frontier moves to the de-confounded-labels question (a later, separate epoch). One incumbent
> artifact (the v0 resolver CLI) is improved in place — NO family generation.
>
> Honest framing (the circularity caveat must ship in the output): the v0 resolver's own caller
> (`detect.py`, method `kmer_seed_coverage_v0`) already matches the SAME VirulenceFinder DB
> (`data/virulencefinder_db/virulence_ecoli.fsa`). So this diff is **lightweight-k-mer-seed-VF vs
> canonical-BLAST-VF over the same DB** — a method-vs-method AUDIT of the resolver's fast caller against
> canonical VirulenceFinder, NOT an independent-baseline comparison. `caller_is_independent_baseline:
> false` already records this; the diff section must carry the same honesty flag. This matches ledger
> Decision T1 (VF is both internal caller AND reference).

## What "done" means (concrete deliverable)
A `vf_diff` section in the CLI provenance JSON: for each VF DB gene, {canonical-VF call (present/%ID/cov),
resolver-k-mer-seed call (present/cov), agreement flag}, plus a per-cluster concordance summary and the
honesty flag. Canonical VF is run via real `virulencefinder` / `blastn` against the input FASTA.

## Mission Gaps

### gap: ship-vf-side-by-side-diff
<!-- strategy-budget: 2 -->
- test-exit-0 python -m pytest tests/test_pathotype_vf_diff.py
- project-state-row project_state/ecoli-pathotype-prediction-cli-2026-05-26.md:VF side-by-side diff shipped

## Endpoint Criteria
- test-exit-0 python -m pytest tests/test_pathotype_vf_diff.py
- project-state-row project_state/ecoli-pathotype-prediction-cli-2026-05-26.md:VF side-by-side diff shipped

## Known friction (resolve during the run)
- **BLAST+ (or KMA) is NOT installed** (`blastn` absent). Canonical VirulenceFinder shells out to one of
  them. The run must install BLAST+ (dep-install → minister `irreversible` gate, autonomous under money-only
  posture) OR route through the Docker VF path (`tools/docker_runner.py`) if a pinned VF image is preferred.
  Decide install-path during interrogation (draft-then-ratify).
- Pin VF software + DB version + DB checksum in the diff provenance (the project's standing version-pin
  discipline; the DB checksum is already emitted in the caller block).

## Test contract (what test_pathotype_vf_diff.py must assert)
1. The `vf_diff` section exists with the per-gene + per-cluster + honesty-flag shape.
2. On a known strain (e.g. an EHEC/STEC reference where stx is unambiguous), canonical VF and the resolver
   AGREE on the decisive cluster (stx present) — pins that the diff is wired correctly, not just structurally.
3. `caller_is_independent_baseline == false` is carried into the diff section (circularity honesty).
4. Offline-safe: if BLAST/VF is unavailable, the diff section degrades to `status: unavailable` with a
   reason (never silently drops it) — so the test runs in CI without the binary but asserts the contract.
