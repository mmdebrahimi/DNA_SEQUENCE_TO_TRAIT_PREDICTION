"""Collect + merge the 2 overnight full-manifest ESM2 shards -> the definitive large-N R2 number.

The overnight run (2026-07-21) sharded the 2569-assay held-out manifest across 2 Kaggle T4 kernels
(mavedb-full-esm2-sharda / -shardb). This pulls both (when terminal), merges by URN, and writes the packet:
overall median |Spearman| + shuffled control + n_scored + per-gene + human-only median.

  uv run --with kaggle python scripts/mavedb_full_esm2_collect.py            # pull + merge + write packet
  uv run python scripts/mavedb_full_esm2_collect.py --from-dir <dir>         # merge already-pulled JSONs
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import statistics as st
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
SHARDS = ["mavedb-full-esm2-sharda", "mavedb-full-esm2-shardb"]


def _pull(dest: Path) -> None:
    """Pull both shards' output via the Kaggle API IF both are terminal; else report + exit."""
    from kaggle.api.kaggle_api_extended import KaggleApi
    os.environ.setdefault("KAGGLE_CONFIG_DIR", os.path.expanduser("~/.kaggle"))
    api = KaggleApi(); api.authenticate()
    user = api.get_config_value("username") or "emanueleebrahimi"
    for slug in SHARDS:
        st_ = api.kernels_status(f"{user}/{slug}")
        status = str(getattr(st_, "status", st_))
        print(f"  {slug}: {status}")
        if "COMPLETE" not in status.upper():
            print(f"  -> {slug} not COMPLETE yet; re-run later.")
            raise SystemExit(2)
        d = dest / slug; d.mkdir(parents=True, exist_ok=True)
        api.kernels_output(f"{user}/{slug}", path=str(d))


def _load_results(dest: Path) -> list[dict]:
    out = []
    for f in glob.glob(str(dest / "**" / "mavedb_full_esm2_shard*_result.json"), recursive=True):
        out.extend(json.loads(Path(f).read_text(encoding="utf-8")).get("results", []))
    return out


def merge_and_report(results: list[dict], manifest: dict) -> dict:
    by_urn = {r["urn"]: r for r in results if r.get("rho") == r.get("rho")}  # drop NaN
    absr = sorted(abs(r["rho"]) for r in by_urn.values())
    shuf = [abs(r["rho_shuf"]) for r in by_urn.values() if r.get("rho_shuf") == r.get("rho_shuf")]
    bygene, human = {}, []
    for u, r in by_urn.items():
        g = manifest.get(u, {}).get("gene", "?").split()[0]
        bygene.setdefault(g, []).append(abs(r["rho"]))
        if manifest.get(u, {}).get("organism") == "Homo sapiens":
            human.append(abs(r["rho"]))
    return {
        "n_scored": len(by_urn),
        "median_abs_spearman": round(st.median(absr), 4) if absr else None,
        "median_abs_shuffled": round(st.median(shuf), 4) if shuf else None,
        "median_human": round(st.median(human), 4) if human else None,
        "n_human": len(human),
        "per_gene_top": {g: round(st.median(v), 3) for g, v in
                         sorted(bygene.items(), key=lambda kv: -st.median(kv[1]))[:15]},
        "n_genes": len(bygene),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from-dir", default="", help="dir with already-pulled shard result JSONs (skip pull)")
    a = ap.parse_args()
    dest = Path(a.from_dir) if a.from_dir else (ROOT / "data" / "raw" / "mavedb_full_esm2")
    if not a.from_dir:
        dest.mkdir(parents=True, exist_ok=True)
        _pull(dest)
    results = _load_results(dest)
    if not results:
        print(f"no shard result JSONs under {dest}"); return 1
    man = {r["urn"]: r for r in
           json.loads((ROOT / "wiki" / "mavedb_prospective_holdout_full_2026-07-21.json")
                      .read_text(encoding="utf-8"))["manifest"]}
    rep = merge_and_report(results, man)
    print(json.dumps(rep, indent=2))
    art = {"_schema": "mavedb-full-esm2-merged-v1", "date": _date.today().isoformat(),
           "shards": SHARDS, "model": "facebook/esm2_t33_650M_UR50D", **rep}
    (ROOT / "wiki" / f"mavedb_full_esm2_{_date.today().isoformat()}.json").write_text(
        json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nDEFINITIVE full-manifest ESM2: {rep['n_scored']} assays, median |Spearman| "
          f"{rep['median_abs_spearman']} (shuffled {rep['median_abs_shuffled']}; human {rep['median_human']} "
          f"on {rep['n_human']}). wrote wiki/mavedb_full_esm2_{_date.today().isoformat()}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
