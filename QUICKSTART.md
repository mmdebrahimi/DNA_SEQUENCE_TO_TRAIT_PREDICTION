# Quickstart — `dna-decode` deterministic decoder

This quickstart uses **only** the default install — **no `[ml]` extra, no Docker, no external databases, no
network**. Every command below is executed end-to-end and asserted by `scripts/verify_quickstart.py` (and
`tests/test_verify_quickstart.py`), so it is *verified*, not aspirational.

> **Not a clinical tool.** Every call reports its own honest **validation tier** inline — read it.

## Install (deterministic core)

```bash
uv sync                 # or: pip install -e .
```

The default install is the **deterministic decoder** (AMR R/S + interpretable determinants + typing). It
does **not** pull `torch`/`transformers`/`xgboost` — those live in the optional `[ml]` extra and are only
needed for the foundation-model embedding track (a closed research negative, not the shipped product):

```bash
uv sync --extra ml      # only if you need the embedding/foundation-model experiments
```

External-tool workflows (genome-mode AMR via AMRFinder, the typing decoders via blastn, the genome-map
Bakta report) additionally need **Docker** and/or **BLAST+** and their databases — see the README Gotchas.
Without them those paths degrade to `unavailable` (never crash); the wheel-only paths below need none of it.

## Wheel-only decoder calls (no Docker, no DBs)

Antibiotic/antiviral/antifungal resistance from **observed substitutions** — pure-Python, instant:

```bash
# HIV-1 (RT NNRTI) -> R, validated against a FREE independent wet-lab label
uv run dna-amr --drug efavirenz   --observed RT:K103N    --sample-id demo
#   CALL: R   ...   validation: INDEPENDENT_WETLAB -- AUC 0.962 ...

# Candida auris (ERG11 azole target site) -> R
uv run dna-amr --drug fluconazole --observed ERG11:Y132F --sample-id demo
#   CALL: R   ...   validation: NO_FREE_PHENOTYPE_SOURCE ...

# SARS-CoV-2 (Mpro / nirmatrelvir) -> R
uv run dna-amr --drug nirmatrelvir --observed Mpro:E166V  --sample-id demo
#   CALL: R   ...   validation: IN_DISTRIBUTION (CoV-RDB), underpowered ...
```

Each call carries its **inline trust badge** (`validation:` line + a `validation` block in `--json-only`):
the cell's honest tier + headline metric + the standing report card it came from. A tier is **never
fabricated** and **never borrowed across organisms** — ask for a drug on the wrong organism and you get
`UNKNOWN (namespace_mismatch)`, not a borrowed number.

## Unified genome profile (cached AMR + offline typing)

Point at one genome assembly → every applicable decoder in one honest report. The AMR section reuses a
**cached** AMRFinder run (no Docker); the typing sections degrade to `unavailable` offline:

```bash
uv run dna-decode profile tests/fixtures/ecoli_mini/genome.fna \
  --amrfinder-run tests/fixtures/amr_mini
#   amr (Escherichia):  (organism ASSUMED E. coli -- pass --amr-organism for a non-E. coli genome)
#     ciprofloxacin  R  [2 det]  | INDEPENDENT_MEASURED acc 0.95 ...
#     ceftriaxone    R  [1 det]  | INDEPENDENT_MEASURED acc 0.919 ...
```

**Organism is self-disclosed:** omit `--amr-organism` and E. coli is *assumed* — the output stamps
`organism_assumed: true` (JSON) + a visible note, so it never masquerades as an explicit choice. For a
non-E. coli genome pass `--amr-organism Klebsiella` (etc.) so routing **and** the badge are correct.

Offline (no AMRFinder source) the AMR section degrades cleanly and the profile still exits 0:

```bash
uv run dna-decode profile tests/fixtures/ecoli_mini/genome.fna   # amr: [unavailable]
```

## Verify the quickstart yourself

```bash
uv run python scripts/verify_quickstart.py     # runs all of the above, asserts each; exit 0 = OK
uv run dna-decode list                         # what it decodes + per-trait routing
```
