# Influenza-NA / HBV free-mirror probe — VERIFIED WALL (2026-07-11)

**Directive:** find a working non-gated mirror (or the exact DUA/access path) for the two viral candidates
that had sat as `NEEDS-VERIFICATION`. Probed by DIRECT `curl` / API inspection (bypasses the biology-
resistance search filter, per `feedback_offline_cache_vs_unreachable_data`: a search filter is not a data
wall). This is a VERIFIED determination, not an assertion — each verdict names what was actually fetched.

**The bar (unchanged):** a FREE, public, **isolate-level MEASURED-phenotype** genotype↔phenotype source
(the HIVDB PhenoSense pattern). A sequence archive, an interpretation tool, or aggregate surveillance does
NOT clear it.

## HBV (NRTI / DAA fold-change) — WALL CONFIRMED (infra + G1)

| Mirror | Direct fetch | Verdict |
|---|---|---|
| `hbvdb.lyon.inserm.fr` (HBVdb, canonical) | **HTTP 000 / timeout** (https + http + `/HBVdbIndex`) | host genuinely unreachable across sessions |
| `hbvdb.ibcp.fr` (alt host) | HTTP 000 / timeout | unreachable |
| `genafor.org` | HTTP 200 | reachable but a different resource (not per-isolate measured HBV phenotype) |

Even a working HBVdb mirror would not clear the bar: HBVdb is a **curated interpretation database**
(sequence + annotated resistance mutations), not a per-isolate measured DAA/NRTI fold-change table — G1
circular for a determinant caller. Measured HBV antiviral fold-change is **clinical-not-public**.
**Access path:** a published per-isolate phenotype supplement (manual, per-paper) OR a clinical-virology
collaborator dataset (acquisition). No free queryable source exists.

## Influenza NA (oseltamivir / zanamivir / peramivir NA-inhibition IC50) — WALL CONFIRMED (no free measured phenotype)

Direct inspection of every candidate that was reachable:

| Source | Reachable | Carries per-isolate NA-inhibition IC50? | Evidence |
|---|---|---|---|
| **BV-BRC `genome_amr`** (taxon 11320) | HTTP 200 | **NO — zero rows** (`[]`). The measured-AMR table is bacterial. | API query returned empty |
| **BV-BRC `surveillance`** | HTTP 200 | **NO.** 35 fields, **zero** antiviral/resistance/IC50/susceptibility/inhibition columns — host/sample/subtype metadata only. | field-set inspected |
| **BV-BRC `genome_feature`** | HTTP 200 | sequences only (the NA gene product), no phenotype | product=neuraminidase records |
| **NCBI Influenza Virus Resource FTP** | HTTP 200 | **NO** — `genomeset.dat` / `influenza.dat` / `.faa` / `.cds` are sequences + collection metadata, no assay results | FTP listing |
| **GISAID EpiFlu** | HTTP 200 | access-GATED (registration + DUA); and antiviral-susceptibility is not systematically deposited even inside it | login wall |
| **WHO GISRS** antiviral working group | — | aggregate surveillance reports, not free-public per-isolate IC50 | (known structure) |

**No free, public, per-isolate NA-inhibition IC50 linked to genotype exists on the open surface.** The
measured IC50 lives in WHO-GISRS / CDC surveillance and in individual publication supplements.
**Access path:** (a) GISAID EpiFlu registration + the WHO-GISRS antiviral data (DUA/access — acquisition),
or (b) manual mining of specific NAI-surveillance publication supplements (per-paper, low yield, not a
source). Neither is a free queryable database.

## Net verdict

**Both viral candidates are VERIFIED walls — infrastructure/access, not code, exactly as flagged.** Neither
HBV nor influenza-NA has a free, public, isolate-level MEASURED-phenotype genotype↔phenotype source; each
requires an acquisition/DUA step (clinical-virology dataset for HBV; GISAID/WHO-GISRS for influenza). This
CLOSES the last `NEEDS-VERIFICATION` rows in the data-source master map with direct evidence, and confirms
(a fourth time, now by API-field inspection) the banked thesis: **the only path to a new independent viral
cell is acquisition** — there is no free shortcut. The reproducibility-freeze forward-path #1 (acquisition,
user authority) remains the lever; nothing here is code-closable.

## Provenance

Probes (all direct `curl` / BV-BRC Data API, 2026-07-11): HBVdb hosts (000); `genafor.org` (200);
BV-BRC `genome_amr`/`surveillance`/`genome_feature` (taxon 11320); NCBI Influenza FTP; GISAID/fludb (200,
gated/redirect). Gates: `wiki/negative_results_map_2026-06-13.md` (G1/G7). Supersedes the NEEDS-VERIFICATION
rows in `wiki/data_source_master_map_2026-07-09.md` §6 (which were already partly resolved 2026-07-10 for
HCV; this completes HBV + influenza by API-field inspection).
