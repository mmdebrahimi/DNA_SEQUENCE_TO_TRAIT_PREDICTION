> **⚠️ SUPERSEDED (2026-06-28).** This memo's premise — "the label wall blocks independence; acquisition
> (author emails / DUA) is the lever" — was already FALSE when written: the EBI AMR Portal (2026-06-23,
> `wiki/amr_portal_feasibility_result_2026-06-23.md`) is a FREE, public, isolate-level, MEASURED-AST source
> that broke the wall (74 powered cells; TB + E. coli/Salmonella/Klebsiella/Shigella/Campylobacter already
> SCORED_INDEPENDENT). The TB author emails are MOOT. Independence is now FREE + code-closable. Kept for
> audit trail. See `wiki/frontier_reassessment_2026-06-28.md`.

# Independent-label acquisition strategy — the best way to find the blocked data sources (2026-06-27)

The decoder's binding constraint is the independent LABEL (reproducibility freeze 2026-06-13). Two cells are
`blocked:external` on it — TB RIF/INH independent gold set, SARS-CoV-2 Mpro independent fold — plus the
general "next-epoch" need for a non-public wet-lab label. This is the ranked, verified acquisition plan: what
to pursue, in what order, by access mode. **The actual contact/DUA/spend is the USER's authority; this memo
makes that decision concrete and one-step-from-execution.**

## The meta-finding (load-bearing): more SEARCH won't help; OUTREACH will
Re-verified 2026-06-27 against 2 fresh 2025 TB candidates — the exhaustion pattern holds:
- **PLOS Glob. Public Health 2025** (`PRJEB68143`): per-isolate measured DST exists but n=17 (too small alone).
- **2025 rifampicin-heteroresistance study** (`PMC12310241`, 2917 genomes): measured pDST is aggregate/subset
  only — the same wall the prior 5× exhaustion (`wiki/tb_goldset_public_source_exhaustion_2026-06-22.md`) found.

So the free public per-isolate-measured path is exhausted across ≥7 checked sources. The remaining levers are
all **acquisition** (author-request → DUA/program-application → paid-last-resort), which is user-authority. The
ingestion pipelines are READY the instant any per-isolate measured label lands (TB: `scripts/score_tb_independent_goldset.py`
+ `organism_rules/tb_goldset.assert_independent`; SARS: `scripts/sarscov2_mpro_validate.py`; both leakage-/
operator-censoring-aware).

## Ranked acquisition paths — TB RIF/INH gold set

| # | Source | Access mode | Cleanness (leakage/circularity) | Effort | Likelihood | USER next action |
|---|---|---|---|---|---|---|
| 1 | **Author request — Clark/LSHTM** (Thorpe 2024 n=59 + Thawong 2025 n≈1826) | free email | CLEAN: genomes verified provdisjoint vs CRyPTIC (0 overlap); measured pDST exists, just not public per-isolate | minutes (email DRAFTED) | medium-high (genomes already public; ask is just the linked table) | **Send Email A** (`wiki/tb_goldset_author_emails_2026-06-22.md`) |
| 2 | **Author request — PLOS-GPH 2025** (`PRJEB68143`, n=17 per-isolate measured) | free email | CLEAN if post-CRyPTIC (accession-check); n=17 = under-powered first cut, combine with #1 | minutes (draft below) | high (per-isolate DST already in supplement; just confirm accession↔R/S linkage) | send the PLOS-GPH draft (below) |
| ~~3~~ | ~~ReSeqTB~~ **DEFUNCT — DROPPED (verified 2026-06-27)** | — | — | — | — | `www.reseqtb.org` now 301-redirects to a commercial QC vendor (domain lapsed); `platform.reseqtb.org` DNS dead. No longer a public download source. (It also carried circularity risk — WHO-catalogue family + CRyPTIC-era — so no loss.) |
| 4 | **NIAID TB Portals** (depot.tbportals.niaid.nih.gov) | DUA/application (genomic split open; clinical+phenotypic split access-gated) | CLEAN (clinical measured DST) once DUA granted | weeks (application) | medium | apply for the clinical/phenotypic data-use agreement |
| 5 | **National TB ref-lab WGS programs** (UK HSA / NL RIVM / ZA NICD) | DUA / collaboration | CLEAN | weeks-months | low-medium | institutional contact |

**Recommended TB order:** send #1 + #2 now (free, drafted/near-drafted, highest yield, CLEAN); verify #3
liveness in parallel (free but circularity-risky — treat as a bonus, not the headline); escalate to #4 DUA
only if #1/#2 go cold. A combined #1+#2 set (≈59+17, plus Thawong if granted) clears the ratified
"hand-curate-30-then-decision" first-cut bar.

## Ranked acquisition paths — SARS-CoV-2 Mpro independent fold

| # | Source | Access mode | Cleanness | Effort | Likelihood | USER next action |
|---|---|---|---|---|---|---|
| 1 | **Future CoV-RDB releases** (hivdb/covid-drdb-payload) | free (passive) | becomes independent ONLY for studies added AFTER the frozen catalog's source set | nil (watch releases) | medium over time | re-run the census periodically; score any NEW held-out study |
| 2 | **Author request — clinical nirmatrelvir-failure cohorts** (Mpro fold from treatment-failure isolates not in CoV-RDB) | free email | CLEAN (external clinical fold) | minutes-hours | low-medium (clinical Paxlovid resistance is rare → small) | contact authors of clinical nirmatrelvir-resistance reports for per-isolate Mpro mutation + measured fold |
| 3 | **prospective-lock the SARS cell** ✅ **BUILT 2026-06-27** | free, in-repo | CLEAN by time-construction | done (`scripts/sarscov2_prospective_lock.py` + manifest + 4 tests; reuses the AMR `is_prospective_eligible`) | high | re-run `sarscov2_mpro_validate.py` filtered to refs dated > lock when a new study lands |

**Recommended SARS order:** #1 is the zero-effort default (the v0 in-distribution number stands meanwhile);
#3 (prospective-lock by time) is the only path that produces an *independent-by-construction* number without
external acquisition — a future code move, not a data move. #2 is low-yield (clinical Paxlovid resistance is rare).
**Honest:** there is no free non-CoV-RDB Mpro-inhibitor fold corpus today (CoV-RDB IS the field aggregator).

## General / next-epoch independent label (the broadest lever)
Per `wiki/next_epoch_idea_anchor_prompt_2026-06-13.md`: the durable unlock is ACQUIRING a non-public
wet-lab/clinical label source that clears the 8 rejection gates by construction. The same access-mode ladder
applies (author-request → DUA → collaboration). The HIV PhenoSense precedent proves the model: one free
INDEPENDENT isolate-level wet-lab label corpus turned a whole cell from in-distribution to validated.

## Outreach kit (ready / near-ready)
- **TB Email A (Clark/LSHTM)** — DRAFTED, ready to send: `wiki/tb_goldset_author_emails_2026-06-22.md`.
- **TB PLOS-GPH 2025 draft** (new): *"Dear [authors of PRJEB68143] — I maintain an open, deterministic
  M. tuberculosis RIF/INH caller (WHO-2023-catalogue rule) and validate it ONLY against measured-phenotype
  benchmarks to stay non-circular. Your 17-isolate set links genomes to measured DST; would you share a
  per-isolate table of ENA accession ↔ measured RIF/INH R/S? Full credit + citation. — Farshad,
  mmdebrahimi@gmail.com"*
- **SARS template** (new): *"Dear [authors] — I maintain an open SARS-CoV-2 Mpro (nirmatrelvir) resistance
  caller and need an INDEPENDENT measured fold-change benchmark not derived from Stanford CoV-RDB. Would you
  share per-isolate Mpro substitutions + measured nirmatrelvir fold-change for your clinical isolates? — Farshad"*

## The decision (USER authority — this is the wall)
Pick which to pursue; everything upstream is done + ready:
Status after the 2026-06-27 sequential pass (moves 2 + 4 executed autonomously; 1 + 3 staged + parked at your authority):
1. **Send the free TB author emails** — STAGED + recipient-verified (5 ready in `wiki/tb_goldset_author_emails_2026-06-22.md`: A Clark, B Ethiopia, C India, D TB-Portals-followup, E Elton/PLOS-GPH). ← **your send** (outward to academics = your authority). Recommended: A + E first.
2. ~~Verify ReSeqTB~~ — **DONE: DEFUNCT, dropped** (domain lapsed).
3. **Apply for TB Portals DUA** — depot confirmed auth-gated (HTTP 403); the application is your authority. Email D is the follow-up once submitted.
4. **SARS prospective-lock-by-time** — **DONE: BUILT** (`scripts/sarscov2_prospective_lock.py`, manifest pinned, 4 tests). Accrues; v0 in-distribution number stands meanwhile.

Nothing here is further code-closable without a label landing. The pipelines are built; the ask is the user's
to send. See `wiki/tb_goldset_howto_2026-06-22.md` for the ingest-the-instant-it-lands runbook.

## Sources (verified 2026-06-27)
- ReSeqTB platform + open-access claim: [reseqtb.org](https://platform.reseqtb.org/), [C-Path WHO-adoption](https://c-path.org/global-health-partners-accelerate-uptake-of-genetic-sequencing-for-surveillance-and-diagnosis-of-drug-resistant-tuberculosis/), [PMC4583571](https://pmc.ncbi.nlm.nih.gov/articles/PMC4583571/)
- TB 2025 candidates: [PLOS GPH 2025 / PRJEB68143](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0004099), [PMC12310241 (heteroresistance)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12310241/)
- Prior exhaustion proof: `wiki/tb_goldset_public_source_exhaustion_2026-06-22.md`; SARS census: `wiki/sarscov2_mpro_v0.1_independence_infeasible_2026-06-27.md`
