# examples

Tiny, **verified** runnable examples. All use data already committed to the repo — no Docker, no GPU, no
network, no large FASTA download.

## 1. Zero-dependency: observed mutation → R/S (wheel-only)

```bash
dna-decode amr --drug efavirenz --observed RT:K103N
```
Expected: `CALL: R` driven by `RT:K103N`, validation tier `INDEPENDENT_WETLAB` (Stanford HIVDB PhenoSense).
Try also `--drug fluconazole --observed ERG11:Y132F` (C. auris) or `--drug nirmatrelvir --observed Mpro:E166V`
(SARS-CoV-2). `--observed` is for the **target-site** engines (viral + fungal), where resistance is a known
substitution in one protein. **Bacterial** drugs are determinant-scan cells and take `--amrfinder-run` or
`--genome-fasta` instead — `--drug ciprofloxacin --observed gyrA:S83L` is rejected, by design, with
`--observed is fungal-only`. (Example 2 is the bacterial path.)

## 2. From a committed AMRFinder run (bacterial; no Docker)

```bash
dna-decode amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_000492655.1
```
Expected: a CALL + the driving determinants + a validation tier (`INDEPENDENT_MEASURED`, EBI AMR Portal).
Other committed accessions live under `data/amrfinder_runs/`; `--drug` accepts any drug shown by
`dna-decode list` for the bacterial engine.

## 3. The full supported surface

```bash
dna-decode list
```

## Machine-readable output

Add `--json` (or `--out result.json`) to any call to emit the structured record (call, determinants,
validation tier, provenance) for agent pipelines.

> These commands are exercised by `make example-observed` / `make example-amr`. Genome-FASTA input
> (`--genome-fasta`) additionally needs Docker (AMRFinder) or BLAST+ — see `docs/quickstart.md`.
