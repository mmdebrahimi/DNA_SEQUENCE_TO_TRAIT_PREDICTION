# NCBI Pathogen Detection — host_disease Facet Audit (Tier-3 substrate-density audit, 2026-05-27)

> Audit memo. NOT a V1 13-column supported memo (no `<!-- memo-schema: ... -->` marker — this is a project-internal analysis sidecar). Sidecar to `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Rows 5+6 (which left NCBI Pathogen Detection's host-related metadata coverage as an honest gap).
>
> **Purpose:** estimate Tier-3 substrate density for COMMENSAL E. coli records via NCBI Pathogen Detection facets. Closes Action 4-new in `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md`.

## Pre-committed verdict bar (set BEFORE the query)

| Verdict | Trigger | Action |
|---|---|---|
| **PASS** | ≥75 commensal-indicating E. coli records confirmed accessible (via host_disease, isolation_source, or epi_type facet) | Include COMMENSAL as full v0 class; promote Tier-3 to substrate strategy as the commensal source |
| **PARTIAL** | 30-74 commensal-indicating records | Demote COMMENSAL to honesty-only label (`COMMENSAL_LOW_MARKER_BURDEN`); keep Tier-3 as supplementary source |
| **FAIL** | <30 records OR access blocked AND no other path | Drop COMMENSAL from v0 surface; emit `COMMENSAL_LOW_MARKER_BURDEN` only on absence-of-DEC-modules logic; document in handoff that v0 commensal class is honesty-only by construction |
| **HONEST-GAP** | Access blocked AND no count obtainable | Defer verdict; queue manual query path for user |

## Query mechanism discovery

NCBI Pathogen Detection Isolates Browser uses SOLR under the hood, exposed as a JavaScript-driven web UI. Three programmatic-access paths were evaluated:

| Path | Mechanism | Verdict for this audit |
|---|---|---|
| Web UI facets (interactive) | SOLR-backed filter panel at `https://www.ncbi.nlm.nih.gov/pathogens/isolates/#/search/taxgroup_name:%22E.coli%20and%20Shigella%22` | Returns SPA HTML shell via WebFetch — no rendered facet data accessible without a JS runtime. User-driven manual interaction is the only access route. |
| BigQuery (programmatic) | `ncbi-pathogen-detect.pdbrowser.isolates` table on Google Cloud Platform | Documented as the recommended programmatic path. **Requires GCP authentication.** Not accessible via WebFetch in this session. |
| Web download (TSV/CSV) | Download button on the Isolates Browser web UI | 100,000 row hard cap; would suffice for E. coli + host_disease facet retrieval. Requires browser interaction (download triggered by JS, not URL-fetchable). |

**Confirmed available metadata fields per NCBI Pathogen Detection help docs** (`https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/`):

| Field | Notes |
|---|---|
| `isolation_source` | Free-text (per substrate-survey Row 9 evidence; weakly structured) |
| `epi_type` | Isolation type — discrete (clinical / environmental/other / NULL) |
| `food_origin` | Food-source category |
| `taxgroup_name` | Organism group (e.g., `E.coli and Shigella`) |
| `geo_loc_name` | Geographic location |
| `collection_date` | ISO date |

**`host_disease` is NOT documented as a top-level facet field** on the Isolates Browser. The standard BioSample attribute `host_disease` exists upstream (NCBI BioSample), but the Pathogen Detection isolates index may not surface it as a facet. The help docs explicitly mention only the fields above. Hovering over column headers in the live web UI is the documented way to see the full field-name list — not accessible via WebFetch.

## Query attempted

| Attempt | URL | Result |
|---|---|---|
| 1 | `https://www.ncbi.nlm.nih.gov/pathogens/isolates/#/search/taxgroup_name:%22E.coli%20and%20Shigella%22` | Returned SPA shell only; no facet data extractable |
| 2 | `https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/` | Field list partial; no `host_disease` example queries surfaced |
| 3 | `https://www.ncbi.nlm.nih.gov/pathogens/docs/isolates_gcp/` | BigQuery schema list partial; `host_disease` column existence not confirmed in public docs |

**Access date:** 2026-05-27.

## Verdict

**HONEST-GAP.**

Programmatic facet access blocked for this discovery session. The COMMENSAL Tier-3 substrate-density question cannot be resolved from this machine without one of:

1. **User-driven manual query** via the live Isolates Browser web UI (recommended; ~10-30 min):
   - Navigate to `https://www.ncbi.nlm.nih.gov/pathogens/isolates/#/search/taxgroup_name:%22E.coli%20and%20Shigella%22`.
   - Hover over column headers to confirm the field-name for host-disease / health-state metadata (the live UI exposes the SOLR field name).
   - Apply the discovered field as a filter / facet.
   - Record top values + counts.
   - Cross-check with `isolation_source` filter (free-text; faceted into LexMapr-style top values).
   - Cross-check with `epi_type = environmental/other` or `NULL` as a commensal-proxy.
2. **BigQuery setup** (one-time cost ~30 min; programmatic thereafter):
   - GCP account + enabled BigQuery API.
   - Run `SELECT host_disease_or_equivalent_field, COUNT(*) FROM `ncbi-pathogen-detect.pdbrowser.isolates` WHERE taxgroup_name = 'E.coli and Shigella' GROUP BY 1 ORDER BY 2 DESC LIMIT 100`.
   - Free tier (1 TB/month) covers this query.
3. **Email pd-help@ncbi.nlm.nih.gov** with the field-name + facet question. Reply latency typically 1-2 business days.

## Workhorse-side recommendation (per scope guardrail)

The workhorse should:

- **NOT block on this verdict.** v0 design holds COMMENSAL as a class (per project ledger v3 11-class surface), but the threshold for activating it vs collapsing to `COMMENSAL_LOW_MARKER_BURDEN` honesty-only can be deferred until the manual query returns.
- **Default to honesty-only emission** if cohort assembly proceeds before the manual query returns: emit `COMMENSAL_LOW_MARKER_BURDEN` based purely on absence-of-DEC-modules + assembly-QC-PASS logic, with no positive commensal-class training.
- **Re-evaluate on count-return.** If manual query confirms ≥75 commensal-indicating E. coli records, retrofit COMMENSAL as a positive class in v0.1+.

## Reproducibility note

- Tool calls: 2 WebSearch + 3 WebFetch. All results captured in this memo.
- No JSON / response excerpts extractable (SPA shell only).
- No counts extractable from this session.
- Manual query path is the documented next step.

## Status (for the project ledger)

- **Action 4-new (NCBI Pathogen Detection `host_disease` facet query):** HONEST-GAP. Pending Decision row NOT yet retired in `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md`. Manual query path now appended to the Pending Decision's Notes column.
- **H1 status:** unchanged. Tier-3 substrate confirmation deferred to user-driven manual query.
- **H2 commensal floor:** still PENDING (Horesh N=2 + Tier-1b Whittam pending + Tier-3 unconfirmed).
