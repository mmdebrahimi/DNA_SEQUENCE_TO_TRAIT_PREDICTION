# TB Portals → independent TB gold set — runbook (2026-06-22)

DUA **signed**. TB Portals is now the active path to a SCORED independent RIF/INH number (no author contact;
powered ~3,305 isolates with measured "Sequenced DST"). This runbook is the one place that says exactly what
YOU download and the single command chain that does the rest.

## Division of labor
- **You (one step, DUA-gated):** download the de-identified DST + accession export from TB Portals and drop it
  on disk. That's it — the labels are the only DUA-gated piece.
- **Me (autonomous after that):** the **genomes are PUBLIC on SRA**, so once the export is on disk I run
  schema-probe → candidate TSV → CRyPTIC leakage check → SRA fetch → variant-call → score, no further gating.

## Step 1 — Download from TB Portals (you)
1. Log in to the TB Portals data portal (post-DUA): **tbportals.niaid.nih.gov** → Data Download (or the
   Analytic API export).
2. Export a table that has, per isolate, ALL of:
   - an **isolate/condition/specimen id**,
   - **DST results** for **Rifampicin** and **Isoniazid** (the per-drug R/S calls),
   - the **DST method** (MGIT / LJ / proportion / Hain / Xpert / WGS …) — **REQUIRED** so we can keep only
     *phenotypic* (measured) calls,
   - a **genome accession** (SRA `SRR…` run, and/or `SAMN…` BioSample). This is what lets us fetch the genome
     and prove non-overlap with CRyPTIC.
3. Save it as CSV (or TSV) to the same place as the CRyPTIC dump, e.g.
   `D:/dna_decode_cache/data files donwload/tbportals_dst_export.csv`.

> If the export splits DST and genomic accessions into two files, download both — tell me the two paths and
> I'll join them on the isolate id.

## Step 2 — I run this (autonomous once the file is on disk)
```bash
# (a) PROBE first — pins the real schema + shows the method vocabulary (W0 discipline).
uv run python scripts/build_tbportals_candidates.py "D:/.../tbportals_dst_export.csv" --probe

# (b) BUILD the candidate TSV — keeps PHENOTYPIC only; DROPS molecular/genotypic DST (gate G1, circularity).
uv run python scripts/build_tbportals_candidates.py "D:/.../tbportals_dst_export.csv" \
    --build --out data/raw/tb_goldset/tbportals_candidates.tsv

# (c) VALIDATE schema + labels + class balance (target ~10 each SS / SR / RR).
uv run python scripts/validate_tb_goldset_candidates.py data/raw/tb_goldset/tbportals_candidates.tsv

# (d) LEAKAGE CHECK vs CRyPTIC (alias-aware: SRA SRR/SAMN vs CRyPTIC ENA ERR/SAMEA) -> per-drug manifest.
uv run python -m scripts.build_tb_goldset_manifest \
    --candidates data/raw/tb_goldset/tbportals_candidates.tsv

# (e) FETCH the public SRA reads + variant-call -> masked VCFs (Docker; targeted/whole-genome per locus need),
#     fill masked_vcf/regeno_vcf, re-run (d), then SCORE:
uv run python -m scripts.score_tb_independent_goldset --drug rifampicin
uv run python -m scripts.score_tb_independent_goldset --drug isoniazid
```

## The one honesty guard that matters here
TB Portals carries **both phenotypic and molecular/genotypic DST**. Scoring our WHO-catalogue determinant
rule against a *molecular* DST call is rule-vs-rule (**circular — the exact reason Thorpe 2024 was a NO-GO**).
`--build` keeps only phenotypic methods and **drops** molecular rows (logged); if the export has **no method
column** it **refuses** unless you explicitly assert phenotypic-only (`--no-method-column-ok`). The number
this produces is `INDEPENDENT_VALIDATION` precisely because the label is measured, not catalogue-derived.

## What "done" looks like
`score_tb_independent_goldset` emits the `INDEPENDENT_VALIDATION` arm (raw + lineage-collapsed sens/spec with
Wilson CI), separate from the in-distribution `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`. That is the
project's first non-circular, out-of-distribution TB number — the gold set's whole point.

## Provenance
Adapter `scripts/build_tbportals_candidates.py` (+ `tests/test_build_tbportals_candidates.py`, 7 tests).
Feeds the existing `validate_tb_goldset_candidates.py` → `build_tb_goldset_manifest.py`
(`tb_goldset.assert_independent_aliased`) → `score_tb_independent_goldset.py`. Shortlist
`research_outputs/tb_goldset_source_shortlist_2026-06-22.md`; Thorpe NO-GO context
`wiki/tb_goldset_thorpe2024_assessment_2026-06-22.md`.
