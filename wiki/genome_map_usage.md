# Genome-map — usage

The "Bakta honesty report": point at ONE microbial genome → an honest,
evidence-tiered per-feature map (4 tiers, phenotype only behind a validated
determinant wall, a DB-labelled unknown rate). Shipped 2026-06-18, spike
verdict **GO**. The achievable, label-free form of the north star — NOT a
learned phenotype predictor.

Package: `dna_decode/genome_map/`. Frozen AMR surface (`amr_rules.py` +
`calibrated_amr_rules.json`) is READ-only.

## Single genome (the primary surface)

```bash
# Live: annotate (Bakta db-light) + AMRFinder via Docker, then map.
MSYS_NO_PATHCONV=1 uv run python -m scripts.genome_map \
  --genome-fasta D:/dna_decode_cache/refseq/GCA_x/genome.fna \
  --organism Escherichia --sample-id GCA_x \
  --out-dir wiki/genome_map_GCA_x

# Hybrid: reuse a precomputed Bakta GFF (skips the slow annotation), run AMRFinder.
MSYS_NO_PATHCONV=1 uv run python -m scripts.genome_map \
  --genome-fasta X.fna --gff X.gff3 --organism Klebsiella_pneumoniae --sample-id X

# Offline / degraded: tiers from a provided GFF + the determinant cells only,
# no Docker (sets degraded_coverage). No phenotype overlay.
uv run python -m scripts.genome_map --gff X.gff3 --no-amrfinder --sample-id X
```

Outputs (per `--out-dir`): `genome_map_<id>.json` (per-feature map),
`genome_map_<id>_table.json` (flat table), `genome_map_<id>.md` (readable summary).

Flags: `--organism none` runs AMRFinder generic (no `-O`, for an organism it
doesn't curate); `--drugs cipro,ceftriaxone,...` scopes the overlay (default =
all `mic_tiers`-supported drugs).

## 3-genome prototype spike (the GO/NO-GO harness)

```bash
MSYS_NO_PATHCONV=1 uv run python -m scripts.genome_map_spike \
  --refseq-cache D:/dna_decode_cache/refseq
# -> wiki/genome_map_spike_verdict_<date>.{md,json} (per-genome 18MB maps gitignored)
```

## Tool-surface manifest (the feasibility gate — run first on a new host)

```bash
MSYS_NO_PATHCONV=1 uv run python -m scripts.genome_map_tool_surface \
  --accession GCA_x --fasta .../genome.fna --organism Escherichia
# -> wiki/genome_map_tool_surface_<date>.json (Bakta vocab + AMRFinder headers,
#    or a documented BAKTA_ANNOTATION_BLOCKED / AMRFINDER_BLOCKED)
```

## The 4 evidence tiers (precedence high → low)

| tier | fires when | phenotype claim? |
|---|---|---|
| `determinant-phenotype` | a HIGH-confidence AMRFinder determinant join (protein-id / coord) | **YES — only this tier** |
| `curated-molecular-function` | a named gene_symbol OR a specific product | no |
| `homology-only-hypothesis` | low-confidence wording (`putative`/`by similarity`/DUF/`domain-containing`) | no |
| `unknown` | `hypothetical protein` / empty product | no |

Headline metric: `unknown_under_bakta_db_light` (the db-light coverage caveat is
IN the field name — high unknown is tooling coverage, not biological unknown).

## Gotchas (real-data, see CLAUDE.md)

- AMRFinder `-n` uses ORIGINAL NCBI contig names; Bakta RENAMES contigs → the
  coord join reconciles by contig length. `Protein id` is `NA` in `-n` mode, so
  the coordinate join is the high-confidence path.
- Bakta emits a whole-contig `region` feature that the coord join excludes (else
  it swallows every determinant). A genome whose determinant joins are ALL
  symbol-fallback NO-GOs (the gene-symbol-trap guard).
- Bakta annotation is CPU-heavy (~10–40 min/genome on this host); a Docker wedge
  → `wsl --shutdown` to recover. Out of v1 scope: TB/fungal overlay (VCF-vs-GFF
  contract), hmmer/Pfam/eggNOG tiers, pathway/KEGG, a visual browser.
```
