# How to build the TB "golden baseline" (independent gold set) — a playbook (2026-06-22)

Written for a software engineer (no TB-genomics background assumed). The current TB number is
*in-distribution* (WHO rule scored on CRyPTIC, which the rule was built from). A **golden / independent**
baseline scores the SAME frozen rule on TB isolates the rule has never seen → an honest, non-circular number.

## What you need, per isolate (the minimal contract)
1. **A measured RIF and/or INH DST result** — a lab phenotype `R`/`S` (NOT a software prediction).
2. **A genome** as a **masked VCF vs H37Rv `NC_000962.3`** (the variant calls our scorer reads). Optionally a
   regeno VCF for callability.
3. **An ENA accession** (run `ERR…` / sample `ERS…` / biosample `SAMEA…`) — used to PROVE it's not in CRyPTIC.

Independence rule (ratified): the isolate must **not** be in CRyPTIC. We check this by accession against the
CRyPTIC reuse table (`tb_goldset.cryptic_accessions` → `assert_independent`). Post-2022/2023 isolates are the
safest (CRyPTIC's June-2022 release + the WHO-2023 catalogue swept most earlier public TB WGS+DST).

## The hard part — and why it isn't an API fetch
TB genotype-phenotype data is concentrated in CRyPTIC + specialist resources; the general AST databases are
dry for TB (verified 2026-06-22: BV-BRC has **0** measured TB DST). So you must obtain a **specific** source.
Candidate source categories (verify specifics in a web-enabled session — my web search was filter-blocked):
- **Published post-2023 TB-WGS + DST studies** that deposit reads to ENA (national surveillance / clinical
  cohorts). Search scholar/ENA for 2024–2025 TB WGS resistance studies with an ENA project accession.
- **NIAID TB Portals** (depot of clinical TB WGS + DST), **ReSeqTB / the Relational Sequencing TB platform**.
- **National TB reference-lab WGS programs** (e.g. UK HSA, NL RIVM, ZA NICD) depositing to ENA with DST.
- A **small hand-curated set (~30 isolates)** is a legitimate first cut (the project's ratified
  "hand-curate-30-then-decision" option) — enough to get a first honest signal even if underpowered.

## Step-by-step (commands)
1. **Get the source data**: for each isolate, its sequencing reads (FASTQ) or a provided VCF, + its measured
   RIF/INH DST, + its ENA accession.
2. **Turn reads into a masked VCF vs H37Rv** (skip if the source already provides VCFs/variants):
   - Easiest: run **TBProfiler** (`tb-profiler profile -1 reads_1.fq.gz -2 reads_2.fq.gz`) — it aligns to
     H37Rv + emits per-isolate variants you can convert to the masked-VCF shape; OR
   - Standard: `minimap2 -ax sr H37Rv.fna reads.fq | samtools sort | bcftools call` → VCF (then apply the
     CRyPTIC mask). H37Rv = NC_000962.3.
   - (If the source provides a CRyPTIC-style `VARIANTS` parquet/table, reuse the genomic-nucleotide adapter
     in `scripts/score_tb_cryptic_parquet.py::parse_cryptic_variant` instead — no VCF needed.)
3. **Assemble a candidate TSV** (one row per isolate) with columns:
   `strain_id  ena_accession  masked_vcf  regeno_vcf  rif_label  inh_label`
4. **Leakage-check + build the manifest** (executor-run, one command):
   ```
   uv run python -m scripts.build_tb_goldset_manifest --candidates my_candidates.tsv
   ```
   This drops any CRyPTIC-present isolate (auditable `leaked_excluded.txt`) and writes
   `data/raw/tb_goldset/goldset_{rif,inh}.json`.
5. **Score the independent number** (the existing arm):
   ```
   uv run python -m scripts.score_tb_independent_goldset --drug rifampicin
   uv run python -m scripts.score_tb_independent_goldset --drug isoniazid
   ```
   → an `INDEPENDENT_VALIDATION` sens/spec (vs the in-distribution baseline), with the same lineage-collapse
   + Wilson-CI honesty rails.

## Honest feasibility + recommendation
- This is a **real acquisition + bioinformatics project**, not a quick fetch. The binding step is *obtaining a
  non-CRyPTIC source with measured DST + genomes* (an external/data wall), then variant-calling it.
- **Lowest-effort first cut:** a ~30-isolate hand-curated set from one post-2023 study that deposits reads +
  publishes DST. Underpowered but yields a first honest independent signal; scale later.
- **What's already built for you:** the leakage check (`tb_goldset.assert_independent`), the ingester
  (`build_tb_goldset_manifest`), and the scorer (`score_tb_independent_goldset`). Once you hand me a candidate
  TSV + the VCFs (or a source to fetch), the leakage→manifest→score path is a few commands and I can drive it.
- I can also, in a fresh web-enabled session, run the source survey (the part my filter blocked here) to hand
  you 2–3 concrete candidate datasets + their ENA accessions.
