# 1001 Genomes flowering time — FT10 / FT16 (the label side of the flowering cell)

Free, public, and **the exact phenotype Zhang & Jiménez-Gómez 2020 used** to associate FRI alleles with
flowering time. This is the label the flowering cell (`dna_decode/organism_rules/arabidopsis_flowering.py`)
should be scored against.

## Provenance

The paper's Methods name the source directly:

> *"Flowering time (from plants grown in long days at 16°C) and STRUCTURE group membership for each accession
> was obtained from the 1001 genomes project (1001 Genomes Consortium, 2016)
> (`http://1001genomes.org/tables/1001genomes-FT10-FT16_and_1001genomes-accessions.html`)"*

**That URL is dead (HTTP 404)** — 1001genomes.org was rebuilt as a SvelteKit site since 2020. The data lives
on at AraPheno, under the study *"1001genomes flowering time phenotypes"*:

```bash
# phenotype 261 = FT10, phenotype 262 = FT16 (long days at 16 C -- the paper's phenotype)
curl -sL https://arapheno.1001genomes.org/rest/phenotype/262/values.csv -o pheno_262.csv
# discover ids:  https://arapheno.1001genomes.org/rest/phenotype/list.json  (536 phenotypes)
```

| file | phenotype | n accessions | min | median | max |
|---|---|---|---|---|---|
| `pheno_261.csv` | FT10 (10 C) | 1,163 | 50.8 | 80.5 | 157.5 |
| `pheno_262.csv` | **FT16 (16 C)** | **1,123** | 32.0 | **70.3** | 145.0 |

Columns: `phenotype_name, accession_id, accession_name, accession_cs_number, accession_longitude,
accession_latitude, accession_country, phenotype_value, obs_unit_id`. **`accession_id` is the 1001 Genomes
id** — the join key to the FRI allele table.

**Both classes are well represented** (median split: FT16 560 late / 563 early), which matters: the earlier
AraPheno **DTF1** spot-check returned `DEGENERATE_NO_LABEL_VARIATION` — the cell could only call ~5
accessions and all were early, so a constant-`early` predictor beat it 5/5 vs 4/5. FT10/FT16 is both the
correct phenotype *and* a non-degenerate one. **Supersedes DTF1 (phenotype 703)** as the scoring label.

## Why this label is trustworthy — it reproduces the paper's biology unaided

Cross-referencing FT16 against the FRI status stated in the article body (no supplement needed):

| accession | FRI status per the paper | FT16 | vs median 70.3 |
|---|---|---|---|
| **Bil-7** | **functional** (their reference functional allele) | **117.0** | **LATE** |
| Col-0 | null — 16-bp frameshift, a044 | 38.0 | early |
| Ler-1 | null — promoter/start deletion, a062 | 40.3 | early |
| Wa-1 | **L294F** — full-length, non-functional | 41.8 | early |
| Wil-2 / Litva / Tottarp-2 / Est-1 | L294F carriers | 48.0 / 47.5 / 43.8 / 45.7 | early |
| Wil-1 | L294F carrier | 58.8 | early |
| Bå1-2 | **L276R** — full-length, non-functional | 45.3 | early |
| Cvi-0 | premature stop, a021 | 47.3 | early |
| NFA-10 | C-terminal truncation (last 115 aa), a019 | 42.5 | early |
| Van-0 | **functional allele, early accession** | 38.0 | early ← *the paper's own documented exception* |

The paper's central finding — that the single substitutions **L294F** and **L276R** disable a *full-length*
FRI protein — is visible directly: those carriers sit with the frameshift/stop nulls at 42–59 days while
functional Bil-7 sits at 117.

**Van-0 is the honest counterexample and it is in-scope for the cell**: a functional FRI allele in an early
accession ("may contain mutations in additional flowering time genes or in unsampled regulatory regions").
The cell caps FRI-route confidence at MEDIUM for exactly this class of case (cf. its Lz-0 anchor).

## STATUS UPDATE 2026-07-16 — Table S3 arrived; this label is now the BACKUP, not the primary

Table S3 landed (browser fetch) and it carries **its own `FT16_mean` column**, so the cell was scored
directly off S3 — this AraPheno copy was **not needed** for the scoring run. It stays committed as an
independent cross-check of the same phenotype (AraPheno FT16 n=1,123 vs S3's 854 phenotyped) and because it
is the *scriptable* route to the label; S3 is browser-only. See `../zhang2020/README.md` +
`wiki/flowering_tables3_score_2026-07-16.md`.

## Superseded: the allele side (now solved)

These ~14 article-text anchors are **not a scoring set** — only **one** (Bil-7) is late, so a constant-`early`
null scores 12/13 and cannot be beaten. Scoring needs the per-accession FRI allele for all 1,016 accessions:
**Table S3** ("Allele id and flowering time for each accession"), which is browser-only — see
`../zhang2020/README.md`. With it, the paper reports 245 accessions carrying loss-of-function alleles vs
~771 functional: both classes, real power.

**Possible full-independence alternative (unexplored):** call FRI functional status ourselves from the public
1001 Genomes variant data, making Table S3 a *validation* of our caller rather than an input — which is what
a decoder should do. **Non-trivial caveat:** the paper had to realign reads to a *modified* reference with
FRI-H51 substituted for FRI-Col-0, precisely because the Col-0 reference itself carries a null allele, and
they discarded >1/3 of raw calls as indel-alignment artefacts. Variants called against stock Col-0 are not
a drop-in substitute.
