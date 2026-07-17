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

## 3. Run a call (zero setup — no Docker, no BLAST, no downloads)

These four run on a bare `pip install`, offline, in seconds — the fastest way to see the tool work:

```bash
# a) resistance from an observed mutation -> R/S + the determinant that drove it
dna-decode amr --drug efavirenz --observed RT:K103N          # HIV-1 NNRTI -> CALL: R

# b) forward: a protein edit -> predicted molecular-effect (DMS-validated, Regime B)
dna-decode forward --mutation S2L --protein-seq MSIQHFRVALIPFFAAFCLPVFA

# c) inverse: a target percentile of damage -> proposed edits (ranks, never doses)
dna-decode inverse --protein-seq MSIQHFRVALIPFFAAFCLPVFA --target-percentile 0.05 --top-k 5

# d) a plant trait: Arabidopsis flowering habit from FRI/FLC allele calls
dna-decode flowering --fri Col --flc Col                     # -> SUMMER_ANNUAL_EARLY
```

Every one fails **loudly** on a bad input rather than guessing — e.g. `forward --mutation A6L` on a protein
whose 6th residue is `F` prints a WT-mismatch error and exits non-zero, never a silent wrong call.

Two more offline paths that need only committed data (still no Docker):

```bash
# from a committed AMRFinder run (bacterial)
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
