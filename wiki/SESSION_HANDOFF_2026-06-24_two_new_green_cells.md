# Session handoff — 2026-06-24 (overnight: two new GREEN cells built)

Built autonomously while you slept, per your "run as many of these such things as you want" permission.
Git clean + pushed (origin/main `0/0`). Project ledger `project_state/dna-decode-2026-05-11.md`.

## What shipped (both directive cells BUILT)
Two new deterministic typing decoders, mirroring `dna-ktype`/`dna-serotype` on the shared blastn engine:

| Cell | Console | Status | Real number |
|---|---|---|---|
| **S. pneumoniae capsular serotype** | `dna-pneumo-serotype` (`dna-decode pneumoserotype`) | SHIPPED + real DB built | **✅ 4/4 exact reference-control** (TIGR4→4, D39→2, Hungary19A-6→19A, Spain23F→23F; native blastn vs real PneumoCaT cps DB) |
| **Salmonella enterica serovar** | `dna-salmserovar` (`dna-decode salmserovar`) | SHIPPED (caller + synthetic control) | ⏳ data-engineering-gated (see below) |

Both: faithful-to-tool (`caller_is_independent_baseline=False`), offline-safe, registered in the dispatcher +
`pyproject` console scripts + `dna-decode list`, with namespace-separate report cards. **FROZEN AMR surface
byte-unchanged** (additive typing cells). 13 new tests (pure-logic + real-blastn controls), all green.

## The pneumococcus real number (the headline)
- Real cps DB built from PneumoCaT's Stage-1 reference (95 serotypes) — `scripts/build_pneumo_cps_db.py`
  (DB at `data/pneumoserotype_db/cps_references.fasta`, **gitignored**; rebuild from a PneumoCaT clone).
- Reference-control: **exact 4/4, serogroup 4/4** at 98.6–100% id / 100% cov. Strain names header-verified
  from ENA — verification **excluded 4 wrong-accession fetches** (returned *Beijerinckia*/*Korarchaeum*/
  *Leptothrix* — not pneumococcus). Artifact: `wiki/pneumo_serotype_reference_control_2026-06-24.json`.
- **Honest scope:** n=4 reference-control validates the real-DB + caller integration; the headline
  GREEN-VALIDATED number still needs the full **GPS Quellung cohort (11,810 isolates)** run.

## Why Salmonella has no real number yet (verified, not assumed)
Cloned + inspected SeqSero2's DB: it's `H_and_O_and_specific_genes.fasta` (H1=fliC / H2=fljB antigens,
extractable) **+ `antigens.pickle` = a k-mer DETECTION index, NOT a serovar formula table**. Unlike PneumoCaT
(flat per-serotype FASTA → pneumo DB built in minutes), SeqSero2 resolves serovar by an ALGORITHM + an
internal White-Kauffmann table. So the Salmonella real DB = genuine data-engineering (extract fliC/fljB +
assemble O-antigen alleles + source the formula→serovar TSV). Deferred — documented in
`wiki/salm_serovar_report_card.md`.

## To finish the GREEN-VALIDATED numbers (the runnable steps)
The full-cohort scorer is shipped: `scripts/serotype_cohort_validate.py --cell {pneumo,salm}` (native blastn,
no Docker, checkpointed). It needs:
1. **pneumo:** a `accession<tab>Quellung_serotype` TSV from the GPS metadata (the 11,810-isolate cohort) +
   genome fetching. The DB is already built. **Validate vs the WET-LAB Quellung label, NOT an in-silico tool**
   (circularity rail). Multi-hour cohort op — best when the host is free (the ktype run is still using it).
2. **salm:** build the real antigen DB first (above), then a `accession<tab>wet-lab-serovar` TSV (NARMS /
   PulseNet lab-serotyped, filtered to NOT-tool-predicted).
Note: no `datasets` CLI on this host — the runner fetches via it if present, else reads `--assembly-dir`;
ENA `https://www.ebi.ac.uk/ena/browser/api/fasta/<acc>` works for genome fetching.

## Other state
- **ktype full validation** still accruing in background (~204/447 at handoff; task may have rolled over —
  check `D:/dna_decode_cache/ktype_build/cohort_results.jsonl` line count + `COHORT_DONE`). Interim wzi-vs-
  serology ~0.55 (below the 0.745 ceiling). See the prior handoff `SESSION_HANDOFF_2026-06-24_pypi_published_ktype_running.md`.
- **PyPI:** `dna-decode 0.5.1` is live + public. The two new cells are NOT yet published (would be 0.5.2 —
  bump `pyproject` version + rebuild + `twine upload` when you want them on PyPI; revoke/reuse a real token).
- **Triage memo** (Round 2): `plans/Non_AMR_GREEN_Cell_Triage_Round2_2026-06-24.md` — pneumo + Salmonella were
  picks #1/#2; remaining slate (H. influenzae, E. coli K-antigen, virulence/toxin, catabolic) is GREEN-CALLER/
  NO_FREE_SOURCE.

## Standing discipline (unchanged)
FROZEN AMR surface sha-pinned; honesty tiers (validate vs wet-lab measured, never the in-silico tool);
money→hard-pause (everything tonight was free); no PyPI publish without you awake; commit-to-main = sync channel.
