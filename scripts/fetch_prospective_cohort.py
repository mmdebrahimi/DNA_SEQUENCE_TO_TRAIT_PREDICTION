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

def _get_json(url: str, timeout: int = 60) -> object:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "dna_decode-prospective/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:   # noqa: S310 (trusted public APIs)
        return json.loads(resp.read().decode("utf-8"))


def _bvbrc_recent_genomes(taxon: int, lock_date: str, limit: int) -> list[dict]:
    # BV-BRC RQL: separate &-clauses are ANDed; gt() on a date needs the FULL ISO timestamp with
    # URL-ENCODED colons (%3A) — date-only or raw-colon forms 400 (verified 2026-06-22).
    iso = f"{lock_date}T00:00:00.000Z".replace(":", "%3A")
    q = (f"eq(taxon_id,{taxon})&eq(public,true)&gt(date_inserted,{iso})"
         f"&select(genome_id,assembly_accession,public,date_inserted)&limit({limit})")
    return _get_json(f"{BVBRC}/genome/?{q}&http_accept=application/json")


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


def fetch_live(lock_date: str, per_taxon_limit: int, amr_limit: int, sleep: float) -> tuple[list[dict], dict]:
    """Run the full funnel live. Returns (cohort_rows, stats)."""
    all_rows, agg = [], {"amr_records": 0, "with_gca": 0, "resolved": 0, "eligible": 0,
                         "excluded_pre_or_undatable": 0, "recent_genomes": 0}
    for taxon, drugs in sorted(cells_by_taxon().items()):
        genomes = _bvbrc_recent_genomes(taxon, lock_date, per_taxon_limit)
        agg["recent_genomes"] += len(genomes)
        gid_to_gca = parse_bvbrc_genomes(genomes)
        if not gid_to_gca:
            print(f"  taxon {taxon}: 0 recent public genomes with assembly since {lock_date}")
            continue
        amr = parse_bvbrc_amr(_bvbrc_amr_for_genomes(list(gid_to_gca), drugs, amr_limit), drugs)
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
        print(f"  taxon {taxon}: {len(genomes)} recent genomes -> {stats['amr_records']} measured AMR rows "
              f"-> {stats['eligible']} eligible post-lock")
    return all_rows, agg


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lock-date", default=LOCK_DATE)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--per-taxon-limit", type=int, default=500)
    ap.add_argument("--amr-limit", type=int, default=5000)
    ap.add_argument("--sleep", type=float, default=0.34, help="seconds between NCBI calls (rate politeness)")
    ap.add_argument("--offline", action="store_true", help="skip the live fetch (smoke the wiring only)")
    args = ap.parse_args(argv)

    cells = cells_by_taxon()
    print(f"[prospective-cohort] lock_date={args.lock_date}  organisms/taxa={sorted(cells)}  "
          f"out={args.out_dir}")
    if args.offline:
        print("[prospective-cohort] --offline: wiring OK, no fetch performed.")
        return 0

    rows, agg = fetch_live(args.lock_date, args.per_taxon_limit, args.amr_limit, args.sleep)
    out_path = args.out_dir / "prospective_cohort.tsv"
    write_cohort_tsv(rows, out_path)
    print(f"[prospective-cohort] funnel: {agg}")
    if rows:
        print(f"[prospective-cohort] wrote {len(rows)} eligible post-lock rows -> {out_path}")
    else:
        print(f"[prospective-cohort] 0 eligible post-lock isolates yet (EXPECTED — small window + "
              f"NCBI/BV-BRC ingestion lag). Wrote header-only TSV -> {out_path}. Re-run periodically as the "
              f"cohort accrues.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
