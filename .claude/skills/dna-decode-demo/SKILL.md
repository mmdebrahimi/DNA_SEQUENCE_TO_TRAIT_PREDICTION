---
name: dna-decode-demo
description: Run a dna-decode demo and summarize the AMR/typing output with its validation caveats. Research-use, never clinical.
---

# dna-decode demo

Give an agent user a clean one-command path to see what `dna-decode` does, with the scientific caveats
preserved.

## Instructions

1. Install in editable mode if `dna-decode` is not importable: `pip install -e ".[dev]"` (or `uv sync`).
2. Run `dna-decode list` and note the supported traits + each cell's validation tier.
3. Run one verified example (no Docker / no GPU needed):
   - zero-dependency: `dna-decode amr --drug efavirenz --observed RT:K103N`
   - from committed data: `dna-decode amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_000492655.1`
4. Summarize the result with these fields, verbatim from the tool output:
   - **input** (sample / observed mutation)
   - **command** run
   - **call** (R / S / abstain)
   - **driven-by** genes or mutations
   - **validation tier** (e.g. INDEPENDENT_MEASURED / INDEPENDENT_WETLAB / faithful-to-tool / NOT_CENSUSED) + provenance
   - **caveats** (blind spots the tool names; e.g. "an S call cannot rule out efflux/porin/regulatory resistance")
5. **Never** present the output as clinical guidance. **Never** claim organisms/drugs `dna-decode list`
   does not show. Preserve the validation tier — do not upgrade it.
