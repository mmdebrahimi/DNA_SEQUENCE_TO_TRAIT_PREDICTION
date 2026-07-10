"""Fetch the prospective-lock accrual cohort for the frozen AMR decoder (LA-3 accrual step).

Builds the `{biosample, first_public_date, gca, drug, label}` TSV that `scripts/prospective_lock_validate.py`
consumes — i.e. isolates that became PUBLIC strictly after the prospective-lock date (2026-06-13) and carry a
MEASURED (non-circular) AST phenotype + a downloadable assembly. Such isolates are leakage-free test cases by
construction (the frozen decoder could not have been tuned to data that did not yet exist).

DATA SOURCES (both free, public, no key — verified 2026-06-22):
  - BV-BRC `genome_amr` API → MEASURED AST: `resistant_phenotype` (Resistant/Susceptible) with
    `laboratory_typing_method != "Computational Prediction"` (excludes the circular ML-predicted rows), per
    organism (taxon_id) × drug; joined to BV-BRC `genome` for `assembly_accession`.
  - NCBI Datasets v2 `genome/accession/<GCA>/dataset_report` → the AUTHORITATIVE assembly `release_date`
    (the real "first public" date) + `biosample`. BV-BRC's own `date_inserted` LAGS NCBI, so it is used ONLY
    as a cheap NECESSARY pre-filter (release_date <= date_inserted always), never as the eligibility date.

FUNNEL (efficient + correct): BV-BRC `genome` where taxon + public + has-assembly + `date_inserted > lock`
(few recent rows) → BV-BRC `genome_amr` measured phenotype for those genomes × the cell's drugs → NCBI
Datasets release_date per GCA → `prospective_lock.is_prospective_eligible` (strictly-after lock, fail-closed).

HONEST: the prospective window is small + both NCBI assembly release and BV-BRC ingestion LAG, so the cohort
ACCRUES over time — `0 eligible` today is the expected, correct state, not a failure. This script is the
standing pipeline; re-run it periodically.

Output (per the user's directive) defaults to `D:/dna_decode_cache/data files donwload/`.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.eval.prospective_lock import SCORED_CELLS, is_prospective_eligible  # noqa: E402

DEFAULT_OUT_DIR = Path("D:/dna_decode_cache/data files donwload")
LOCK_DATE = "2026-06-13"
BVBRC = "https://www.bv-brc.org/api"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"

# SCORED-cell registry organism -> BV-BRC taxon_id (the dominant species/genus for the cell).
ORG_TAXON = {
    "Escherichia_coli_Shigella": 562,   # E. coli
    "Klebsiella": 573,                  # K. pneumoniae
    "Campylobacter": 194,               # genus Campylobacter
}


def cells_by_taxon() -> dict[int, set[str]]:
    """{taxon_id: {drug,...}} derived from the frozen SCORED grid (only mappable organisms)."""
    out: dict[int, set[str]] = {}
    for org, drug in SCORED_CELLS:
        t = ORG_TAXON.get(org)
        if t is not None:
            out.setdefault(t, set()).add(drug)
    return out


# ---- pure parsers / mappers (offline-tested) ----

def phenotype_to_label(resistant_phenotype: str | None) -> str | None:
    """BV-BRC `resistant_phenotype` -> 'R' / 'S' / None (skip Intermediate/ambiguous/blank)."""
    s = (resistant_phenotype or "").strip().lower()
    if s == "resistant":
        return "R"
    if s == "susceptible":
        return "S"
    return None  # intermediate / non-susceptible / blank -> not a clean R/S label


def parse_bvbrc_genomes(rows: list[dict]) -> dict[str, str]:
    """BV-BRC genome rows -> {genome_id: assembly_accession} (only public rows with a GCA/GCF)."""
    out: dict[str, str] = {}
    for r in rows:
        gid = str(r.get("genome_id", "")).strip()
        acc = (r.get("assembly_accession") or "").strip()
        if gid and acc.startswith(("GCA_", "GCF_")) and r.get("public", True):
            out[gid] = acc
    return out


def parse_bvbrc_amr(rows: list[dict], wanted_drugs: set[str]) -> list[dict]:
    """BV-BRC genome_amr rows -> [{genome_id, drug, label}] for measured rows with a clean R/S label."""
    out = []
    for r in rows:
        method = (r.get("laboratory_typing_method") or "").strip().lower()
        if method == "computational prediction":   # circular ML row -> exclude
            continue
        drug = (r.get("antibiotic") or "").strip().lower()
        if drug not in wanted_drugs:
            continue
        label = phenotype_to_label(r.get("resistant_phenotype"))
        gid = str(r.get("genome_id", "")).strip()
        if label and gid:
            out.append({"genome_id": gid, "drug": drug, "label": label,
                        "method": r.get("laboratory_typing_method", "")})
    return out


def parse_datasets_report(report_json: dict) -> dict:
    """NCBI Datasets dataset_report -> {release_date, biosample, status} (empty strings if absent)."""
    reports = report_json.get("reports") or [{}]
    ai = (reports[0] or {}).get("assembly_info", {}) or {}
    return {
        "release_date": (ai.get("release_date") or "")[:10],
        "biosample": ((ai.get("biosample") or {}).get("accession") or ""),
        "status": ai.get("assembly_status") or "",
    }


def build_cohort_rows(amr: list[dict], gid_to_gca: dict[str, str],
                      gca_release: dict[str, dict], lock_date: str) -> tuple[list[dict], dict]:
    """Assemble eligible cohort rows + a funnel-stats dict (PURE — release dates supplied, no network).

    A row is emitted IFF the isolate has a GCA, a resolved NCBI release_date+biosample, and is
    prospective-eligible (release_date strictly after lock_date)."""
    rows, stats = [], {"amr_records": len(amr), "with_gca": 0, "resolved": 0, "eligible": 0,
                       "excluded_pre_or_undatable": 0}
    for rec in amr:
        gca = gid_to_gca.get(rec["genome_id"])
        if not gca:
            continue
        stats["with_gca"] += 1
        rel = gca_release.get(gca)
        if not rel or not rel.get("release_date") or not rel.get("biosample"):
            continue
        stats["resolved"] += 1
        verdict = is_prospective_eligible(rel["release_date"], lock_date)
        if not verdict.eligible:
            stats["excluded_pre_or_undatable"] += 1
            continue
        stats["eligible"] += 1
        rows.append({"biosample": rel["biosample"], "first_public_date": rel["release_date"],
                     "gca": gca, "drug": rec["drug"], "label": rec["label"]})
    return rows, stats


def write_cohort_tsv(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["biosample\tfirst_public_date\tgca\tdrug\tlabel"]
    lines += [f"{r['biosample']}\t{r['first_public_date']}\t{r['gca']}\t{r['drug']}\t{r['label']}" for r in rows]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---- live HTTP seam (stdlib urllib — no dep-install on a disk-tight host) ----

class SourceUnavailable(RuntimeError):
    """The upstream API could not be queried (outage / timeout / error envelope).

    Load-bearing: an accrual pipeline's honest signal is "0 eligible isolates have accrued yet". A DEAD
    SOURCE must never be able to produce that same signal, or an outage silently reads as a real zero.
    Every fetch path raises this rather than returning an empty list.
    """


def _get_json(url: str, timeout: int = 60, retries: int = 3, backoff: float = 3.0) -> object:
    """GET JSON, retrying transient failures. Raises SourceUnavailable — never returns empty on error.

    BV-BRC answers an upstream outage with HTTP 200 carrying an ERROR ENVELOPE
    (`{"status":500,"message":"... 503 Service Unavailable ..."}`), not an HTTP error. Iterating that dict
    yields str keys, so the row parsers blow up with an obscure AttributeError deep in the funnel. Detect
    it here, at the boundary.
    """
    req = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": "dna_decode-prospective/1.0"})
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:   # noqa: S310 (trusted public APIs)
                payload = json.loads(resp.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001 — URLError / TimeoutError / JSONDecodeError all transient
            last = e
        else:
            # BV-BRC error envelope: a dict where a LIST of rows was expected.
            if isinstance(payload, dict) and int(payload.get("status") or 0) >= 400:
                last = SourceUnavailable(f"API error envelope: {str(payload.get('message'))[:120]}")
            else:
                return payload
        if attempt < retries - 1:
            time.sleep(backoff * (2 ** attempt))
    raise SourceUnavailable(f"{type(last).__name__ if last else 'unknown'}: {last}") from last


def _recent_genomes_query(clause: str, lock_date: str, limit: int) -> str:
    # BV-BRC RQL: separate &-clauses are ANDed; gt() on a date needs the FULL ISO timestamp with
    # URL-ENCODED colons (%3A) — date-only or raw-colon forms 400 (verified 2026-06-22).
    iso = f"{lock_date}T00:00:00.000Z".replace(":", "%3A")
    return (f"{clause}&eq(public,true)&gt(date_inserted,{iso})"
            f"&select(genome_id,assembly_accession,public,date_inserted)&limit({limit})")


def _bvbrc_recent_genomes(taxon: int, lock_date: str, limit: int) -> tuple[list[dict], str]:
    """Recent public genomes for `taxon`. Returns (rows, taxon_filter_used).

    `eq(taxon_id,N)` matches only genomes whose OWN taxid is N, so it MISSES strain-level taxids. Measured
    offline against the local BV-BRC genome export: only 95.4% of 85,099 E. coli genomes carry taxon_id 562;
    the remaining 4.6% sit under 2,422 strain taxids (244319 / 83334 = O157:H7 lineages — exactly the
    AST-tested pathogens this cohort wants). `eq(taxon_lineage_ids,N)` matches N and all descendants.

    The lineage field is NOT verified against the live API (it was down when this was written), so we TRY it
    and fall back to the narrow filter, stamping which one actually ran.
    """
    try:
        rows = _get_json(f"{BVBRC}/genome/?{_recent_genomes_query(f'eq(taxon_lineage_ids,{taxon})', lock_date, limit)}"
                         f"&http_accept=application/json", retries=2)
        if isinstance(rows, list):
            return rows, "taxon_lineage_ids"
    except SourceUnavailable:
        pass  # may be an unsupported-field 400 OR a real outage; the narrow query disambiguates
    rows = _get_json(f"{BVBRC}/genome/?{_recent_genomes_query(f'eq(taxon_id,{taxon})', lock_date, limit)}"
                     f"&http_accept=application/json")
    if not isinstance(rows, list):
        raise SourceUnavailable(f"unexpected payload type {type(rows).__name__}")
    return rows, "taxon_id_only"


def _bvbrc_amr_for_genomes(genome_ids: list[str], drugs: set[str], limit: int) -> list[dict]:
    if not genome_ids:
        return []
    gid_clause = ",".join(genome_ids)          # genome_ids are digits.dot -> no URL-encoding needed
    dr_clause = ",".join(sorted(drugs))        # drug names are lowercase words -> no encoding needed
    q = (f"in(genome_id,({gid_clause}))&in(antibiotic,({dr_clause}))"
         f"&select(genome_id,antibiotic,resistant_phenotype,laboratory_typing_method)&limit({limit})")
    return _get_json(f"{BVBRC}/genome_amr/?{q}&http_accept=application/json")


def _ncbi_release(gca: str) -> dict:
    return parse_datasets_report(_get_json(f"{DATASETS}/{urllib.parse.quote(gca)}/dataset_report"))


# ---------------------------------------------------------------------------------------------------
# SOURCE 2: NCBI Pathogen Detection (added 2026-07-10).
#
# WHY a second source: BV-BRC's API was fully down (HTTP 500 / a 503 error envelope) when the first
# accrual sweep ran, so no number could be earned. NCBI-PD is not merely a fallback — it is the SAME
# source the 10 frozen SCORED cells were built from, so an accrual cohort drawn from it is the most
# apples-to-apples prospective extension of that grid.
#
# Its `latest_snps/Metadata/PDG*.metadata.tsv` carries (verified live 2026-07-10, 67 cols):
#   asm_acc | AST_phenotypes ("ciprofloxacin=R;gentamicin=S;...") | biosample_acc
#   collection_date | sra_release_date | target_creation_date
# AST_phenotypes are SUBMITTED measured phenotypes (not the circular computational-prediction rows that
# BV-BRC mixes in and this script has to exclude by `laboratory_typing_method`).
# ---------------------------------------------------------------------------------------------------

PD_ORG_GROUPS: tuple[str, ...] = ("Campylobacter", "Escherichia_coli_Shigella", "Klebsiella")


def cells_by_group() -> dict[str, set[str]]:
    """{PD organism group: {drug,...}} from the frozen SCORED grid."""
    out: dict[str, set[str]] = {}
    for org, drug in SCORED_CELLS:
        if org in PD_ORG_GROUPS:
            out.setdefault(org, set()).add(drug)
    return out


def parse_ast_phenotypes(field: str | None, wanted_drugs: set[str]) -> dict[str, str]:
    """PD `AST_phenotypes` -> {drug: label} for wanted drugs, R/S only.

    REAL format, verified against the live PD metadata (2026-07-10):
        "ampicillin=ND,cefazolin=ND,ceftriaxone=R,ciprofloxacin=S,gentamicin=ND"
    i.e. COMMA-separated and wrapped in literal double quotes. Two traps, both hit in this repo:
      * the separator is a comma, NOT a semicolon (the `;` form appears only in the DERIVED
        `candidates.tsv`, not in the source);
      * the surrounding quotes ride on the FIRST and LAST tokens, so a naive `tok in ast.split(",")`
        silently never matches a drug that happens to sit at either end of the list.
    Values other than R/S (`ND` not-determined, `I`, `NS`) are DROPPED — the frozen decoder emits a
    binary R/S call, so those isolates carry no unambiguous ground truth.
    """
    raw = (field or "").strip()
    if not raw or raw == "NULL":
        return {}
    out: dict[str, str] = {}
    for part in raw.replace(";", ",").split(","):
        token = part.strip().strip('"').strip("'").strip()
        if "=" not in token:
            continue
        drug, _, val = token.partition("=")
        drug, val = drug.strip().lower(), val.strip().upper()
        if drug in wanted_drugs and val in ("R", "S"):
            out[drug] = val
    return out


def earliest_public_date(*dates: str | None) -> str | None:
    """The EARLIEST non-empty ISO date among the candidates (the conservative 'first public' bound).

    Prospective eligibility is fail-closed, so we take the earliest date any evidence supports. If none
    is available we return None, which `is_prospective_eligible` treats as INELIGIBLE.
    """
    clean = sorted(d[:10] for d in dates if d and len(d) >= 10 and d[:4].isdigit())
    return clean[0] if clean else None


def _pd_metadata_lines(group: str, timeout: int = 600):
    """Stream the PD metadata TSV line by line (these files are hundreds of MB — never buffer whole)."""
    from scripts.ncbi_pd_provenance_census import latest_metadata_url

    url = latest_metadata_url(group)
    try:
        resp = urllib.request.urlopen(url, timeout=timeout)  # noqa: S310 (trusted NCBI FTP-over-HTTPS)
    except Exception as e:  # noqa: BLE001
        raise SourceUnavailable(f"{type(e).__name__}: {e}") from e
    with resp:
        for raw in resp:
            yield raw.decode("utf-8", "replace").rstrip("\n")


def fetch_live_ncbi_pd(lock_date: str, sleep: float, row_cap: int | None = None,
                       ) -> tuple[list[dict], dict, dict[str, str], dict[str, str]]:
    """PD accrual funnel. Same contract as `fetch_live`: a failure is RECORDED, never a silent zero."""
    all_rows: list[dict] = []
    agg = {"pd_rows": 0, "with_ast_and_asm": 0, "prefilter_post_lock": 0,
           "resolved": 0, "eligible": 0, "excluded_pre_or_undatable": 0}
    group_status: dict[str, str] = {}
    group_filter: dict[str, str] = {}

    for group, drugs in sorted(cells_by_group().items()):
        try:
            candidates: list[dict] = []
            idx: dict[str, int] = {}
            n_rows = 0
            for i, line in enumerate(_pd_metadata_lines(group)):
                cols = line.split("\t")
                if i == 0:
                    idx = {c: j for j, c in enumerate(cols)}
                    missing = [c for c in ("asm_acc", "AST_phenotypes", "sra_release_date") if c not in idx]
                    if missing:
                        raise SourceUnavailable(f"PD metadata missing columns {missing}")
                    continue
                n_rows += 1
                if row_cap and n_rows > row_cap:
                    break

                def g(name: str) -> str:
                    j = idx.get(name, -1)
                    return cols[j].strip() if 0 <= j < len(cols) else ""

                asm = g("asm_acc")
                if not asm.startswith(("GCA_", "GCF_")):
                    continue
                labels = parse_ast_phenotypes(g("AST_phenotypes"), drugs)
                if not labels:
                    continue
                agg["with_ast_and_asm"] += 1
                sra_rel = g("sra_release_date")
                # cheap NECESSARY pre-filter: an isolate whose SRA was already public pre-lock can never
                # be eligible (its earliest public date is <= that). Fail-closed on a missing date.
                if not sra_rel or sra_rel[:10] <= lock_date:
                    continue
                agg["prefilter_post_lock"] += 1
                candidates.append({"asm": asm, "biosample": g("biosample_acc"),
                                   "sra_release_date": sra_rel[:10], "labels": labels})

            agg["pd_rows"] += n_rows
        except SourceUnavailable as e:
            group_status[group] = f"source_unavailable: {e}"
            print(f"  {group}: SOURCE UNAVAILABLE ({e}) — NOT counted as '0 accrued'")
            continue

        group_status[group] = "ok"
        group_filter[group] = "ncbi_pd_metadata"
        print(f"  {group}: {n_rows} PD rows -> {len(candidates)} post-lock candidates with measured AST + assembly")

        # confirm with the AUTHORITATIVE assembly release_date; take the EARLIEST public date.
        for c in candidates:
            try:
                rel = _ncbi_release(c["asm"])
            except Exception as e:  # noqa: BLE001 — one bad accession must not abort the sweep
                print(f"  WARN datasets fetch failed for {c['asm']}: {type(e).__name__}")
                agg["excluded_pre_or_undatable"] += 1
                continue
            time.sleep(sleep)
            agg["resolved"] += 1
            earliest = earliest_public_date(c["sra_release_date"], rel.get("release_date"))
            if not is_prospective_eligible(earliest, lock_date).eligible:
                agg["excluded_pre_or_undatable"] += 1
                continue
            for drug, label in c["labels"].items():
                all_rows.append({
                    "biosample": c["biosample"] or rel.get("biosample", ""),
                    "first_public_date": earliest,
                    "gca": c["asm"],
                    "drug": drug,
                    "label": label,
                })
            agg["eligible"] += 1

    return all_rows, agg, group_status, group_filter


def overall_status(taxon_status: dict) -> str:
    """OK only when EVERY taxon was successfully queried. Anything else must not read as a real zero."""
    if not taxon_status:
        return "SOURCE_UNAVAILABLE"
    oks = [t for t, s in taxon_status.items() if s == "ok"]
    if len(oks) == len(taxon_status):
        return "OK"
    return "SOURCE_UNAVAILABLE" if not oks else "PARTIAL_SOURCE_UNAVAILABLE"


def fetch_live(lock_date: str, per_taxon_limit: int, amr_limit: int,
               sleep: float) -> tuple[list[dict], dict, dict[int, str], dict[int, str]]:
    """Run the full funnel live. Returns (cohort_rows, stats, taxon_status, taxon_filter_used).

    A per-taxon failure is RECORDED and the sweep continues; it is never swallowed into "0 accrued".
    """
    all_rows, agg = [], {"amr_records": 0, "with_gca": 0, "resolved": 0, "eligible": 0,
                         "excluded_pre_or_undatable": 0, "recent_genomes": 0}
    taxon_status: dict[int, str] = {}
    taxon_filter: dict[int, str] = {}
    for taxon, drugs in sorted(cells_by_taxon().items()):
        try:
            genomes, filt = _bvbrc_recent_genomes(taxon, lock_date, per_taxon_limit)
            taxon_filter[taxon] = filt
        except SourceUnavailable as e:
            taxon_status[taxon] = f"source_unavailable: {e}"
            print(f"  taxon {taxon}: SOURCE UNAVAILABLE ({e}) — NOT counted as '0 accrued'")
            continue

        agg["recent_genomes"] += len(genomes)
        gid_to_gca = parse_bvbrc_genomes(genomes)
        if not gid_to_gca:
            taxon_status[taxon] = "ok"
            print(f"  taxon {taxon}: 0 recent public genomes with assembly since {lock_date} "
                  f"(genuine zero; filter={filt})")
            continue

        try:
            amr = parse_bvbrc_amr(_bvbrc_amr_for_genomes(list(gid_to_gca), drugs, amr_limit), drugs)
        except SourceUnavailable as e:
            taxon_status[taxon] = f"source_unavailable(amr): {e}"
            print(f"  taxon {taxon}: AMR query SOURCE UNAVAILABLE ({e})")
            continue

        gca_release: dict[str, dict] = {}
        for gca in {gid_to_gca[r["genome_id"]] for r in amr if r["genome_id"] in gid_to_gca}:
            try:
                gca_release[gca] = _ncbi_release(gca)
            except Exception as e:   # noqa: BLE001 — one bad GCA must not abort the sweep
                print(f"  WARN datasets fetch failed for {gca}: {type(e).__name__}")
            time.sleep(sleep)
        rows, stats = build_cohort_rows(amr, gid_to_gca, gca_release, lock_date)
        for k in stats:
            agg[k] = agg.get(k, 0) + stats[k]
        all_rows.extend(rows)
        taxon_status[taxon] = "ok"
        print(f"  taxon {taxon}: {len(genomes)} recent genomes -> {stats['amr_records']} measured AMR rows "
              f"-> {stats['eligible']} eligible post-lock (filter={filt})")
    return all_rows, agg, taxon_status, taxon_filter


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lock-date", default=LOCK_DATE)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--per-taxon-limit", type=int, default=500)
    ap.add_argument("--amr-limit", type=int, default=5000)
    ap.add_argument("--sleep", type=float, default=0.34, help="seconds between NCBI calls (rate politeness)")
    ap.add_argument("--offline", action="store_true", help="skip the live fetch (smoke the wiring only)")
    ap.add_argument("--source", choices=("ncbi_pd", "bvbrc"), default="ncbi_pd",
                    help="ncbi_pd (default): the SAME source the frozen SCORED cells came from. "
                         "bvbrc: the original API path (was fully down 2026-07-10).")
    ap.add_argument("--row-cap", type=int, default=None,
                    help="cap PD metadata rows per group (smoke runs only — a cap makes '0 eligible' "
                         "meaningless, so it is stamped into the status artifact)")
    args = ap.parse_args(argv)

    if args.source == "bvbrc":
        scope = sorted(cells_by_taxon())
    else:
        scope = sorted(cells_by_group())
    print(f"[prospective-cohort] lock_date={args.lock_date}  source={args.source}  scope={scope}  "
          f"out={args.out_dir}")
    if args.offline:
        print("[prospective-cohort] --offline: wiring OK, no fetch performed.")
        return 0

    if args.source == "bvbrc":
        rows, agg, taxon_status, taxon_filter = fetch_live(
            args.lock_date, args.per_taxon_limit, args.amr_limit, args.sleep)
    else:
        rows, agg, taxon_status, taxon_filter = fetch_live_ncbi_pd(
            args.lock_date, args.sleep, row_cap=args.row_cap)
    status = overall_status(taxon_status)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / "prospective_cohort.tsv"
    status_path = args.out_dir / "prospective_cohort_status.json"
    # A row-capped sweep did not see the whole source, so its "0 eligible" carries no ACCRUING meaning.
    if args.row_cap and status == "OK":
        status = "TRUNCATED_SMOKE_RUN"

    status_path.write_text(json.dumps({
        "_schema": "prospective-cohort-fetch-status-v1",
        "lock_date": args.lock_date,
        "source": args.source,
        "row_cap": args.row_cap,
        "overall_status": status,
        "taxon_status": {str(k): v for k, v in taxon_status.items()},
        "taxon_filter_used": {str(k): v for k, v in taxon_filter.items()},
        "funnel": agg,
        "n_eligible_rows": len(rows),
        "cohort_tsv_written": status == "OK",
        "coverage_caveat": (
            "taxon_id_only filters MISS strain-level taxids (measured offline: 95.4% of 85,099 BV-BRC "
            "E. coli genomes carry taxon_id 562; the other 4.6% sit under strain taxids incl. the O157:H7 "
            "lineages). Where taxon_filter_used == 'taxon_id_only' the cohort is an UNDER-COUNT."),
        "honesty": (
            "A cohort TSV is written ONLY when overall_status == OK. '0 eligible' is a meaningful ACCRUING "
            "signal exclusively under OK — a dead/degraded source must never be reportable as a real zero."),
    }, indent=2), encoding="utf-8")

    print(f"[prospective-cohort] funnel: {agg}")
    print(f"[prospective-cohort] overall_status={status}  -> {status_path}")

    if status == "TRUNCATED_SMOKE_RUN":
        print("[prospective-cohort] --row-cap set: wiring exercised, but a truncated sweep cannot claim "
              "ACCRUING. No cohort TSV written.")
        return 1

    if status != "OK":
        print(f"[prospective-cohort] REFUSING to write a cohort TSV: the source was degraded, so "
              f"'0 eligible' would be indistinguishable from a genuine zero. Per-scope: {taxon_status}")
        return 2 if status == "SOURCE_UNAVAILABLE" else 1

    write_cohort_tsv(rows, out_path)
    if rows:
        print(f"[prospective-cohort] wrote {len(rows)} eligible post-lock rows -> {out_path}")
    else:
        print(f"[prospective-cohort] ACCRUING: 0 eligible post-lock isolates yet (EXPECTED — small window + "
              f"NCBI/BV-BRC ingestion lag), and every taxon queried cleanly. Wrote header-only TSV -> "
              f"{out_path}. Re-run periodically as the cohort accrues.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
