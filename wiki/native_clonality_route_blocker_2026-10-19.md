# The native clonality route's blocker was mis-identified — it is NOT distances

**Date:** 2026-10-19 · **Status:** `blocked:substrate` · **Artifact:** none (deliberately — see below)

## Correction (read first)

The Mission-Reconciliation plan treats the clonality route as blocked on **distances**:
Mash needs Docker/WSL, both unavailable on this host, so Steps 3-6 built and validated a
pip-installable `sourmash` producer instead (`dna_decode/eval/sourmash_distance.py`,
crosscheck `wiki/sourmash_vs_mash_crosscheck_2026-10-19.json`, 20/20 pairs, G8 bar amended).

That work is real and correct. **It does not unblock the lineage metrics.** A second
Docker dependency sits *upstream* of the distance step, and nothing had reached it before
today because the consumer was never run.

## Measured

`scripts/compute_lineage_metrics.py --distance-backend sourmash`, all ten provdisjoint
cohorts: **10/10 returned at the M4 reconcile gate with `recomputed=0`**, before the
distance producer was reached. Wall clock **3 s** — 580 genome sketches cannot occur in 3 s,
which is independent confirmation the backend never ran.

Cause, measured not inferred:

| input | state |
|---|---|
| `data/raw/<cohort>/selected.tsv` | **present** (10/10, recovered from server2 today) |
| `data/raw/<cohort>/amrfinder_runs/` | **absent, 0 of 10** |
| `data/raw/<cohort>/refseq/` | **absent, 0 of 10** (the 580 assemblies are on `E:`, reached via `--cohort-root`) |

`reconcile_raw_metrics` recomputes each isolate's prediction from the cached AMRFinder
output and refuses the cluster-weighted metric unless it reproduces the committed
confusion. With no runs, every recomputation is 0, and the gate correctly refuses.

**The committed artifacts cannot substitute.** Checked, not assumed: all
`wiki/provenance_disjoint_validation_*.json` carry aggregates only — `metrics`, `n_scored`,
`n_selected` — and **zero per-isolate keys**. Per-isolate predictions must be recomputed.

## The actual chain

lineage metrics → per-isolate predictions → `amrfinder_runs/` → AMRFinderPlus → Docker.

**Not claimed: that no non-Docker AMRFinder route exists here.** AMRFinderPlus is partly a
wrapper over blast, so a native route is conceivable and is **UNTESTED**. Six failed routes
would be evidence about six routes, not about the tool — the ECOFF lesson.

**CORRECTION (2026-10-20): this memo originally asserted BLAST+ 2.17.0 is installed natively
at `C:/Users/Farshad/ncbi-blast/bin` on this host. That is FALSE and I took it from CLAUDE.md
prose without checking** — exactly the error this repo's own orientation header warns about.
Measured through the DEPLOYED resolver: `vf_runner.find_blastn()` returns **None**; `$BLASTN_BIN`
unset, `blastn` not on PATH, and the file is absent under BOTH `~/ncbi-blast/bin` and
`C:/Users/Farshad/ncbi-blast/bin` (that profile exists; it has no `ncbi-blast`). So a native
BLAST+ is not available here as a starting point either — it would have to be installed
first, which is a separate decision.

## Two defects found in this run, both mine

1. **False provenance.** `main()` stamped `distance_backend` on every cell, including the
   ten that returned early. The artifact claimed `backend=sourmash` for a run sourmash
   contributed nothing to. Now stamped only on the path that reaches the producer;
   pinned by `test_an_EARLY_RETURN_cell_does_not_claim_a_backend_that_never_RAN`.
2. **Vacuous exit 0.** A run where no cohort reached the distance step printed
   `wrote ... (10 cells)` and exited 0 — indistinguishable from success. Now exits **4**
   with an explicit "exercised the backend ZERO times" line. This is the
   *output-that-has-only-ever-been-empty* pattern: the sourmash path in this consumer has
   still **never executed**, so it remains unvalidated end-to-end whatever the unit tests say.

No `wiki/provdisjoint_lineage_metrics_sourmash.json` was kept: an artifact named for a
backend that contributed no numbers is itself a trap.

## Namespace separation (the load-bearing safety property)

Each backend owns its own sidecar; `sidecar_path_for` is the only path selector and raises
on an unknown backend. Pointing one backend at the other's file **fails closed (exit 2)** —
the ten committed Mash cells share the `(organism, drug)` key with any sourmash
re-derivation, so a single file would overwrite them silently. Proven non-vacuous: the same
call shape on the owned path is accepted. The committed Mash artifact is byte-unchanged and
still holds its ten complete cells (pinned by test).

## What would unblock it

Docker on **server2** (its own blockers: session-0 launch + ~5.37 GB free RAM against a
~8 GB WSL2 default), or a measured answer to the untested native-AMRFinder question above.
