# Quickstart — install → run → read the output

`dna-decode` is research-use software for deterministic microbial genome→trait decoding. **Not a clinical tool.**

## 1. Install

```bash
pip install dna-decode
# or, from a clone, for development:
pip install -e ".[dev]"          # or: uv sync   (then prefix commands with `uv run`)
```

The shipped deterministic CLI imports no `torch`/GPU stack — a plain install is light and cross-platform.

## 2. See what it supports

```bash
dna-decode list
```

This is the **authoritative** supported-trait surface: every drug/organism/typing scheme the tool calls,
each with its honest validation tier and caveats. Prefer it over any static table.

## 3. Run a call (no Docker / no GPU)

```bash
# Zero-dependency: an observed mutation -> R/S + provenance
dna-decode amr --drug efavirenz --observed RT:K103N

# From a committed AMRFinder run (bacterial; no Docker for this path)
dna-decode amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_000492655.1
```

## 4. Read the output

Every call prints, in order:
- **CALL** — `R` / `S` (or an abstention) + a confidence tier + the number of determinants found,
- **driven by** — the exact genes / mutations behind the call (or "no curated determinants"),
- a one-line **mechanism note** + the **blind spots** the call cannot rule out,
- **validation:** the tier (e.g. `INDEPENDENT_MEASURED`, `INDEPENDENT_WETLAB`, faithful-to-tool,
  `NOT_CENSUSED`) + the provenance (which report card / cohort backs it).

`--json` / `--out result.json` emit the machine-readable record for agent pipelines.

## 5. Genome-FASTA input (optional, needs external tools)

```bash
dna-decode amr --drug ceftriaxone --genome-fasta X.fna --sample-id X     # bacterial: AMRFinder via Docker
```

The bacterial genome path runs AMRFinderPlus in a pinned container (Docker); the fungal/viral target-site
engines need NCBI BLAST+ on PATH. The `--amrfinder-run` and `--observed` paths above need neither.

## Next

- `docs/validation.md` — what each validation tier means + where the numbers come from.
- `AGENTS.md` — the agent-oriented setup/test/use guide + the scientific guardrails.
