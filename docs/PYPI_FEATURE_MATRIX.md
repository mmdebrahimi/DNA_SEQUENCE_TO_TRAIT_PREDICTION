# PyPI Feature Matrix

## Goal

Make first-use expectations explicit:

- which command exists
- what input it expects
- whether base install is enough
- what extra assets or tools it needs

## Matrix

| Command | Purpose | Typical input | Base install only? | Extra assets / tools | Notes |
|---|---|---|---|---|---|
| `dna-decode list` | enumerate shipped surfaces | none | yes | none | first smoke check |
| `dna-amr` | AMR trait call with explicit determinants | observed mutation, AMRFinder run, or genome FASTA | sometimes | genome mode may need AMRFinderPlus outputs or Docker; some modes need organism context | trust badge is load-bearing |
| `dna-pathotype` | E. coli pathotype compatibility | assembly FASTA | yes | none in base path | scope-limited to supported organism |
| `dna-plasmid` | plasmid Inc replicon typing | assembly FASTA | no | PlasmidFinder DB, BLAST+ | faithful-to-tool surface |
| `dna-serotype` | E. coli O:H serotype | assembly FASTA | no | SerotypeFinder-style DB, BLAST+ | unresolved loci should abstain |
| `dna-resfinder` | acquired AMR gene calling | assembly FASTA | no | ResFinder DB, BLAST+ | independent cross-tool check vs `amr` |
| `dna-pointfinder` | chromosomal AMR point mutations | assembly FASTA | no | PointFinder DB, BLAST+ | organism/locus scope matters |
| `dna-disinfinder` | disinfectant resistance genes | assembly FASTA | no | DisinFinder DB, BLAST+ | often paired with plasmid/co-localization |
| `dna-mlst` | MLST sequence type | assembly FASTA | no | PubMLST scheme, BLAST+ | incomplete profiles are not guessed |
| `dna-ktype` | Klebsiella capsule K-type | assembly FASTA | no | wzi DB, BLAST+ | faithful-to-tool, not a serology surrogate |
| `dna-salmserovar` | Salmonella serovar | assembly FASTA | no | antigen DB, BLAST+ | serovar only when formula resolves uniquely |
| `dna-pneumo-serotype` | pneumococcal capsular serotype | assembly FASTA | no | cps reference DB, BLAST+ | v0 is better at serogroup than exact within-serogroup calls |
| `dna-pgx` | human pharmacogenomics | phased VCF (GRCh38) | yes | none | pure Python deterministic caller |
| `dna-concordance` | compare AMR callers | decoder outputs / caller outputs | yes | depends on compared surfaces | second-opinion analysis |
| `dna-profile` | run-all summary report | assembly FASTA plus whatever each selected surface needs | no | may compose BLAST+, DBs, AMRFinder, Docker | highest setup-risk surface |
| `dna-coloc` | AMR × plasmid co-localization | AMR + plasmid evidence | no | contig-level evidence inputs | suggestive, not proof |

## Install language

Base install:

```bash
pip install dna-decode
```

Heavier extras:

```bash
pip install "dna-decode[ml]"
pip install "dna-decode[quantize]"
```

Meaning:

- base install: deterministic decoder CLI core
- `ml`: foundation-model / heavy research extras
- `quantize`: quantized-model extras where supported

## What the PyPI page should make explicit

- which commands are pure Python after install
- which commands require DB downloads
- which commands require BLAST+
- which commands require Docker
- which commands require precomputed AMRFinder or other external tool outputs

Without that, users will overestimate what `pip install dna-decode` alone enables.
