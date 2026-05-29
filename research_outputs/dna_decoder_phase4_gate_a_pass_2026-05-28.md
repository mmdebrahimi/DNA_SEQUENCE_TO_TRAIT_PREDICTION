# DNA Decoder Phase 4 Gate A Pass — 2026-05-28

## Headline

- Gate A verdict: `PASS`
- caller runtime: `READY` in external mode
- local panel: `5 / 5` staged
- next gating item is **Gate B human validation**, not more local substrate work

## What completed

- installed `virulencefinder==3.2.0` into repo `.venv`
- cloned local VirulenceFinder DB under:
  - `data/processed/pathotype_gate_a/virulencefinder_db`
- downloaded 4 missing prototype genomes and froze the Gate A panel:
  - `ehec_edl933` -> `GCF_001997045.1`
  - `etec_h10407` -> `GCF_000210475.1`
  - `eaec_042` -> `GCF_000027125.1`
  - `upec_cft073` -> `GCF_000007445.1`
  - `commensal_mg1655` -> `GCF_000005845.2`
- installed local BLAST+ runtime under:
  - `data/processed/pathotype_gate_a/blast/ncbi-blast-2.17.0+/bin/blastn.exe`
- patched Gate A substrate to use the real runtime contract:
  - `python -m virulencefinder`
  - cloned VF DB
  - external `blastn`
- ran the full 5-strain Gate A call chain successfully
- populated manual verdict sheet and assembled final Gate A result

## Important runtime finding

- the original Gate A wrapper assumption was wrong
  - it looked for separate `virulencefinder` / `etecfinder` executables or Docker images
- actual working runtime on this machine is:
  - `python -m virulencefinder`
  - VirulenceFinder DB clone
  - local `blastn.exe`
- current Gate A `etecfinder` branch is an **alias** to the VF runtime scoped to `virulence_ecoli`
  - this was sufficient to make the decision table computable on the frozen 5-strain panel
  - future v0 spec work should still decide whether a stricter standalone ETECFinder DB source is required

## Final Gate A result

- final artifact:
  - `reports/pathotype_gate_a_result_2026-05-28.md`
  - `reports/pathotype_gate_a_result_2026-05-28.json`
- JSON headline:
  - `caller_runtime_ok = true`
  - `raw_outputs_usable = true`
  - `decision_table_computable = true`
  - `gate_a_verdict = PASS`

## Manual panel outcomes

- `commensal_mg1655` -> `COMMENSAL_LOW_MARKER_BURDEN`
- `ehec_edl933` -> `EHEC_COMPATIBLE`
- `etec_h10407` -> `ETEC_COMPATIBLE`
- `eaec_042` -> `EAEC_COMPATIBLE`
- `upec_cft073` -> `UPEC_COMPATIBLE`

## Key artifacts

- contract:
  - `reports/dna_decoder_phase4_gate_a_contract_2026-05-28.md`
- manifest:
  - `config/pathotype_gate_a_panel.json`
- preflight:
  - `data/processed/pathotype_gate_a/gate_a_preflight_external_2026-05-28.json`
- raw run:
  - `data/processed/pathotype_gate_a/run_summary.json`
- manual verdicts:
  - `data/processed/pathotype_gate_a/manual_verdicts/gate_a_manual_verdicts.csv`
- final result:
  - `reports/pathotype_gate_a_result_2026-05-28.md`

## Next step

- **Do not start the Horesh heavy build yet unless Gate B also passes.**
- best next action is now external / human:
  - complete Gate B user-pain validation
- if Gate B passes:
  - this machine becomes the workhorse for the Horesh H1-passing substrate build
