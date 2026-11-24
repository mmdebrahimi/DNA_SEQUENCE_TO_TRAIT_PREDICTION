# dna-decode

`dna-decode` is a deterministic, interpretable genome-to-trait CLI.

Given a genome, a curated tool output, or an observed mutation, it reports:

- resistance or susceptibility calls for supported antibiotic, antiviral, and antifungal surfaces
- the exact genes or mutations driving the call
- typing/pathotype calls where supported
- trust badges, provenance, and explicit blind spots

It is mechanism-based, not an embedding model.

> Research use only. `dna-decode` is not a clinical diagnostic tool and must not be used to guide patient treatment.

## Install

```bash
pip install dna-decode
dna-decode list
```

Optional extras:

```bash
pip install "dna-decode[ml]"
pip install "dna-decode[quantize]"
```

## Quickstart

List the shipped surfaces:

```bash
dna-decode list
```

Two pure CLI examples that do not require a genome assembly:

```bash
dna-amr --drug efavirenz --observed RT:K103N --sample-id HIV_DEMO

dna-amr --drug fluconazole --observed ERG11:Y132F --sample-id CAURIS_DEMO
```

Representative assembly-driven commands:

```bash
dna-decode pathotype path/to/ecoli_assembly.fna --sample-id MY_STRAIN

dna-decode ktype path/to/klebsiella_assembly.fna --sample-id MY_STRAIN

dna-decode profile path/to/assembly.fna --sample-id MY_STRAIN
```

## What it includes

Current package surfaces include:

- AMR calling across bacteria, M. tuberculosis, fungi, HIV-1, SARS-CoV-2, and influenza
- E. coli pathotype calling
- plasmid, serotype, ResFinder, PointFinder, DisinFinder, MLST, and K-type typing surfaces
- Salmonella serovar and pneumococcal serotype typing
- human pharmacogenomics (`dna-pgx`)
- cross-decoder concordance, profile, and co-localization analyses

See the feature matrix:

- `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION/blob/main/docs/PYPI_FEATURE_MATRIX.md`

## Trust and validation

Every result should carry its own trust badge rather than borrowing confidence from another surface.

Common badges:

- `INDEPENDENT_MEASURED`
- `INDEPENDENT_WETLAB`
- `IN_DISTRIBUTION`
- `FAITHFUL_TO_TOOL`
- `ABSTAINS_BY_DESIGN`

Compact validation summary:

| Surface family | Primary trust shape | Notes |
|---|---|---|
| bacterial AMR | `INDEPENDENT_MEASURED` | measured AST on provenance-disjoint cohorts |
| M. tuberculosis | `INDEPENDENT_MEASURED` | deterministic WHO-catalog style rule on measured AST |
| HIV-1 | `INDEPENDENT_WETLAB` | isolate-level wet-lab fold-change validation |
| SARS-CoV-2 / influenza | mixed by cell | some surfaces are in-distribution or underpowered |
| fungal AMR | scope-limited | useful, but not a universal cross-lineage binary claim |
| typing / pathotype / PGx | varies by surface | often faithful-to-tool or phenotype-faithful rather than independent phenotype validation |

See:

- `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION/blob/main/docs/PYPI_VALIDATION_SUMMARY.md`

## Requirements by feature

Not every command is zero-setup after `pip install dna-decode`.

- Some commands are pure Python after install.
- Some need an organism-specific database directory.
- Some need BLAST+.
- Some AMR paths need AMRFinder outputs or Docker.

See:

- `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION/blob/main/docs/PYPI_FEATURE_MATRIX.md`

## Limitations

- not a clinical tool
- not a universal phenotype predictor
- trust tier varies by surface
- some surfaces are faithful-to-tool or in-distribution only
- some commands require external databases or helper tools
- the package prefers abstention over overclaiming when a surface is out of scope

## Documentation

- Quickstart: `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION/blob/main/docs/quickstart.md`
- Feature matrix: `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION/blob/main/docs/PYPI_FEATURE_MATRIX.md`
- Validation summary: `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION/blob/main/docs/PYPI_VALIDATION_SUMMARY.md`
- Repository: `https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`

## Development

For source-checkout work:

```bash
uv sync
uv run pytest tests/ -v
```
