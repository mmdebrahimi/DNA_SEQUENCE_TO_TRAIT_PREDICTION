# Variety-decoder roadmap — run as many as possible (laptop, while workhorse runs G2)

> 2026-06-08. Goal: maximize the genome→trait decoder VARIETY shipped on the laptop (no GPU, no money)
> while the workhorse runs the one GPU-heavy job (Path B / G2). Grounded by a 2026-06-08 feasibility probe.
> Current tool: 5 decoders (amr [bacterial+fungal], pathotype, plasmid, serotype, resfinder) on one shared
> curated-DB blastn engine (`dna_decode/typing/blast_caller.py`). North star: AI DNA decoder TOOL, breadth.

## Feasibility probe (2026-06-08) — what's actually buildable
| Candidate | Kind | DB / logic status | Verdict |
|---|---|---|---|
| amr⇄resfinder concordance | combined (no new DB) | both decoders exist; pure logic | **EASY / high VOI** |
| `dna-decode profile` (run-all) | combined (no new DB) | composes the 5 existing decoders | **EASY-MED / high VOI** |
| plasmid×(amr/resfinder) co-localization | combined + engine ext | engine outfmt lacks subject contig+coords (only `qseqid pident length qlen`) → needs a positions mode | **MED / high VOI** |
| PointFinder (chromosomal AMR point mutations) | new DB | `pointfinder_db/<species>/resistens-overview.txt` = HTTP 200; needs position-mapping logic (not presence) | **MED** |
| MLST (7-gene sequence type) | new DB | needs allele-profile→ST table; probed `mlst_db` raw path 404 → discover real layout first | **MED** |
| DisinFinder (biocide/disinfectant) | new DB | probed `disinfinder_db/disinfectant.fsa` 404 → filename differs; discover first | **LOW-MED (gated on discovery)** |
| S. aureus typing (spaTyper / SCCmecFinder) | new DB | not probed; spa/SCCmec are region/repeat typing (more than presence) | **MED (defer)** |

Insight: the **combined analyses** (compose the 5 decoders we already shipped) are simultaneously the
highest-VOI AND the most feasible — they turn separate callers into integrated analyses, which is where
product value compounds. New-DB decoders beyond simple presence (PointFinder/MLST) are real but heavier and
have DB-discovery friction. So the roadmap front-loads combined analyses.

## Wave 1 — combined analyses (no new DB; compose the existing 5) — ✅ DONE 2026-06-08 (commit 80d90fb)
Highest VOI/cost; zero new external dependency; each is pure logic + tests over shipped decoders.
1. **`amr`⇄`resfinder` concordance** — run both AMR callers on a genome, report per-gene/per-drug agreement
   (both-call / amr-only / resfinder-only), surfacing the independent cross-check resfinder was built for.
   New: `dna_decode/concordance/` (or a `dna-decode concordance` subcommand) + tests. ~1 unit.
2. **`dna-decode profile <genome>`** — one command runs every applicable decoder (amr-resfinder + plasmid +
   serotype + pathotype) on an assembly and emits a unified per-genome report (the "tell me everything"
   UX). Thin orchestration over existing CLIs; offline-safe per-decoder degrade. ~1 unit.

## Wave 2 — engine extension + co-localization (one small core change unlocks a high-VOI analysis)
3. **Engine positions mode** — extend `blast_caller.call_alleles` to optionally return subject contig +
   coords (`-outfmt "6 qseqid sseqid pident length qlen sstart send"`), behind a flag so existing callers
   are unchanged. ~0.5 unit + tests.
4. **`plasmid`×(`amr`/`resfinder`) co-localization** — link each acquired resistance gene to the plasmid
   replicon on the SAME contig (gene-on-replicon), answering "is THIS resistance plasmid-borne?" concretely
   rather than "are both present somewhere". Uses the positions mode. ~1 unit + tests.

## Wave 3 — new-DB decoders needing more than presence (DB-discovery first)
5. **PointFinder** — chromosomal AMR point-mutation decoder (`pointfinder_db/<species>/`), parsing
   `resistens-overview.txt` to map codon/position changes → resistance. Complements AMRFinder's POINT calls;
   per-species. First step: fetch + parse the overview format. ~1.5-2 units.
6. **MLST** — 7-gene sequence typing: best allele per locus → profile → ST via the species `profiles` table.
   First step: discover the real `mlst_db` layout (raw path 404). New logic (profile→ST lookup). ~1.5-2 units.

## Wave 4 — more presence-style curated-DB decoders (gated on DB discovery)
7. **DisinFinder** (biocide/disinfectant resistance) — gene presence like resfinder; discover the real DB
   filename first (probed name 404). ~0.5 unit once DB resolves.
8. **S. aureus typing** (spaTyper / SCCmecFinder) — repeat/region typing; heavier; defer unless requested.

## Cross-cutting (every wave)
- Each decoder/analysis: offline-safe degrade, uniform record schema, real-BLAST e2e test on synthetic
  fixtures + offline-degrade test, wired into `dna-decode` dispatcher + pyproject + the registry-contract
  test, DB downloaded on demand (not committed). Docs (README capability table + CHANGELOG) + ledger row
  per wave.

## Execution recommendation
Run **Wave 1 now** (2 combined analyses — highest VOI, zero dependency, fully testable offline). Then
Wave 2 (engine ext unlocks co-loc). Waves 3-4 are heavier + have DB-discovery friction — do after Waves 1-2
land, or on explicit request. Suggested scope per `/soraya` batch: one wave at a time, commit per item,
re-assess after each wave for the diminishing-returns plateau.
