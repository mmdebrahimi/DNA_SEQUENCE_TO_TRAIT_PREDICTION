"""Fetch AraPheno phenotypes and rank them by overlap with the 1001-Genomes genotyped accessions.

WHY: the F3 calibration's frozen bar needs MIN_TRAITS_FOR_SIGN_TEST (5) traits per control, and only 3
Arabidopsis phenotypes are on disk. The honest fix is to MEET the bar, not to lower it -- re-sizing a
frozen bar after seeing a result is exactly the move this project refused earlier the same week.

AraPheno lists 536 Arabidopsis phenotypes with a public REST API whose `values.json` shape matches the
phenotype CSVs already committed under `data/arabidopsis/` (verified: id 261 == FT10). Overlap with the
1,135 genotyped accessions is the binding quantity and cannot be known without fetching, so this caches
every phenotype once and ranks them.

SEQUENTIAL AND POLITE BY DESIGN: one request at a time with a delay. The recorded failure mode in this
project is fan-out -- 12 concurrent agents died on a rate limit with zero reports -- and this is a public
academic server, not a CDN.

Run:  uv run python scripts/fetch_arapheno_traits.py            # cache + rank
      uv run python scripts/fetch_arapheno_traits.py --limit 40 # bounded pass
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

BASE = "https://arapheno.1001genomes.org/rest"
CACHE = Path("D:/dna_decode_cache/arapheno")
KINSHIP_H5 = "D:/dna_decode_cache/arabidopsis_1001g/1001_SNP_MATRIX/kinship_ibs_mac5.hdf5"
MIN_OVERLAP = 250          # asserted: enough for 5-fold CV with room for grouped folds
DELAY_S = 0.35


def _get(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers={"User-Agent": "dna_decode-research/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return json.loads(fh.read().decode("utf-8", errors="replace"))


def genotyped_accessions() -> set[str]:
    """The 1,135 accessions the genotype matrix actually carries."""
    import h5py

    with h5py.File(KINSHIP_H5, "r") as f:
        return {s.decode() for s in f["accessions"][:]}


def fetch_list(force: bool = False) -> list[dict]:
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / "phenotype_list.json"
    if p.exists() and not force:
        return json.loads(p.read_text(encoding="utf-8"))
    d = _get(f"{BASE}/phenotype/list.json")
    p.write_text(json.dumps(d), encoding="utf-8")
    return d


def fetch_values(pid: int) -> list[dict] | None:
    """Cached; returns None on a fetch/parse failure (recorded, never silently treated as empty)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f"values_{pid}.json"
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:                             # noqa: BLE001 - corrupt cache -> refetch
            p.unlink(missing_ok=True)
    try:
        d = _get(f"{BASE}/phenotype/{pid}/values.json")
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError):
        return None
    except Exception:                                 # noqa: BLE001
        return None
    if not isinstance(d, list):
        return None
    p.write_text(json.dumps(d), encoding="utf-8")
    return d


def overlap_count(values: list[dict], geno: set[str]) -> tuple[int, int]:
    """(n_overlapping_accessions_with_a_finite_value, n_values_total)."""
    seen = set()
    for row in values:
        acc, val = row.get("accession_id"), row.get("phenotype_value")
        if acc is None or val is None:
            continue
        try:
            float(val)
        except (TypeError, ValueError):
            continue
        k = str(int(acc))
        if k in geno:
            seen.add(k)
    return len(seen), len(values)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=0, help="0 = all listed phenotypes")
    ap.add_argument("--min-overlap", type=int, default=MIN_OVERLAP)
    ap.add_argument("--out", type=Path,
                    default=REPO / "wiki" / f"arapheno_trait_overlap_{date.today().isoformat()}.json")
    a = ap.parse_args(argv)

    geno = genotyped_accessions()
    print(f"genotyped accessions: {len(geno)}", flush=True)
    lst = fetch_list()
    print(f"phenotypes listed: {len(lst)}", flush=True)
    ids = [p["phenotype_id"] for p in lst if p.get("phenotype_id") is not None]
    if a.limit:
        ids = ids[: a.limit]
    meta = {p["phenotype_id"]: p for p in lst}

    rows, failed = [], []
    for k, pid in enumerate(ids, 1):
        cached = (CACHE / f"values_{pid}.json").exists()
        v = fetch_values(pid)
        if not cached:
            time.sleep(DELAY_S)
        if v is None:
            failed.append(pid)
            continue
        n_ov, n_tot = overlap_count(v, geno)
        m = meta.get(pid, {})
        rows.append({"phenotype_id": pid, "name": str(m.get("name", ""))[:90],
                     "study": str(m.get("study", ""))[:120],
                     "n_overlap_genotyped": n_ov, "n_values": n_tot})
        if k % 50 == 0:
            print(f"  {k}/{len(ids)} fetched; usable so far="
                  f"{sum(1 for r in rows if r['n_overlap_genotyped'] >= a.min_overlap)}", flush=True)

    rows.sort(key=lambda r: -r["n_overlap_genotyped"])
    usable = [r for r in rows if r["n_overlap_genotyped"] >= a.min_overlap]
    out = {
        "schema": "arapheno-trait-overlap-v1", "generated": date.today().isoformat(),
        "n_genotyped_accessions": len(geno),
        "n_phenotypes_listed": len(lst), "n_fetched": len(rows), "n_failed": len(failed),
        "failed_ids": failed[:50],
        "min_overlap": a.min_overlap, "n_usable": len(usable),
        "usable": usable, "all_ranked": rows[:120],
        "honest_limits": [
            "Overlap counts accessions with a FINITE value that appear in the genotype matrix; it does "
            "NOT check the phenotype is quantitative-continuous, nor that it is independent of the "
            "others (many AraPheno traits come from one study and are correlated).",
            "A high overlap does not make a trait a good calibration target. Several of these are "
            "binary/ordinal scores, which a ridge will handle poorly.",
            "min_overlap is ASSERTED, not derived.",
            "Failures are RECORDED, never counted as zero-overlap -- a fetch failure and an empty "
            "phenotype must not look alike.",
        ],
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print()
    print(f"usable (overlap >= {a.min_overlap}): {len(usable)} of {len(rows)} fetched "
          f"({len(failed)} failed)")
    for r in usable[:20]:
        print(f"  {r['phenotype_id']:>5}  n={r['n_overlap_genotyped']:>4}  {r['name'][:60]}")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
