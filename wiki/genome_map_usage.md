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
all `mic_tiers`-supported drugs); `--no-virulence` skips the VF virulence overlay
(see below).

## Virulence-determinant overlay (E. coli / Shigella, shipped 2026-06-21)

A 5th overlay tier (`dna_decode/genome_map/virulence_overlay.py`). Where a curated
VirulenceFinder (VF) allele is PRESENT in an E. coli/Shigella genome, it's surfaced
behind the SAME coordinate-join integrity gate + a presence-only wall as the AMR
`determinant-phenotype` tier (presence of a curated determinant, NEVER a learned
pathogenicity claim). The deterministic pathotype-resolver call is shown SEPARATELY
as a genome-level overlay (`genome_pathotype_call`, the virulence analog of the AMR
R/S call) — QC-gated, so a low-QC genome yields `AMBIGUOUS_LOW_QC` rather than a
confident commensal.

- VF caller: `vf_runner.run_canonical_vf(all_hits=True)` (NON-frozen) — raises
  `-max_target_seqs` so tandem/multi-copy alleles survive; adds a coord-retaining
  `per_hit` list + `db_sha`. `per_gene`/`per_cluster` stay best-hit byte-identical
  so `build_vf_diff` is unchanged.
- **no-`-parse_seqids` pin:** `makeblastdb` WITHOUT `-parse_seqids` keeps `sseqid` =
  the exact FASTA header first-token (the same token AMRFinder reports) so the shared
  `phenotype_overlay.build_contig_name_map` reconciles both overlays.
- All VF metrics live under SEPARATE keys (`virulence_join_quality`,
  `all_virulence_joins_symbol_fallback`, `virulence_determinant_feature_count`,
  `genome_pathotype_call`) — they never feed the AMR `all_joins_symbol_fallback`
  nor the AMR GO/NO-GO spike gate.
- `virulence_status ∈ {FULL, UNAVAILABLE_NO_BLASTN, SKIPPED_NON_ECOLI, SKIPPED_USER}`.
- DB: committed `data/virulencefinder_db/virulence_ecoli.fsa`. Offline-degrades to
  `UNAVAILABLE_NO_BLASTN` without `blastn`.
- Deferred: non-E.coli VF DBs; a virulence GO/NO-GO spike gate; folding virulence
  into the AMR tier. The 3-genome spike stays AMR-only in v1
  (`run_genome_map_for(virulence=True)` is opt-in; default False).

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

## The evidence tiers (precedence high → low)

| tier | fires when | phenotype claim? |
|---|---|---|
| `determinant-phenotype` | a HIGH-confidence AMRFinder determinant join (protein-id / coord) | presence-only (AMR R/S behind the wall) |
| `virulence-determinant` | a HIGH-confidence VF allele coord-join (E. coli/Shigella; shipped 2026-06-21) | presence-only (pathotype call is a SEPARATE genome-level overlay) |
| `curated-molecular-function` | a named gene_symbol OR a specific product | no |
| `homology-only-hypothesis` | low-confidence wording (`putative`/`by similarity`/DUF/`domain-containing`) | no |
| `unknown` | `hypothetical protein` / empty product | no |

Per-feature precedence: high-confidence AMR join → `determinant-phenotype` wins over
a high-confidence VF join → `virulence-determinant`; a symbol-fallback VF hit is
`secondary_evidence` only (excluded from the tier).

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
