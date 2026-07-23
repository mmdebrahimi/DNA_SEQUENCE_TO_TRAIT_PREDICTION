"""Build GEMME per-variant tables for the held-out hybrid proteins — the evolution modality for the 3-way.

The 2-way ESM2+ProSST hybrid ran on Kaggle. GEMME needs Docker (JET2/R/MUSCLE), which Kaggle can't run, so
GEMME is computed HERE (local Docker) and its per-variant tables are SHIPPED to Kaggle as a dataset — mirroring
how the ProSST structure tokens were pre-computed locally and shipped. The Kaggle notebook then combines
ESM2 + ProSST + the shipped GEMME table into the 3-way rank-average.

Per held-out protein (from wiki/holdout_hybrid_manifest.json): fetch the assay sequence, fetch its ColabFold
MSA (cached, single-query etiquette), run GEMME in Docker, keep the GEMME score for each of the assay's DMS
single-missense variants (ASSAY-LOCAL numbering, the same the manifest/Kaggle use). Checkpointed per protein
to wiki/gemme_holdout_tables.json (restartable). GEMME work dirs go to D: (C: is disk-tight).

  MSYS_NO_PATHCONV=1 uv run python scripts/build_gemme_holdout_tables.py --limit 30

Frozen AMR surface byte-unchanged (READ-only).
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward.gemme_scorer import run_gemme, GemmeUnavailable  # noqa: E402
from dna_decode.forward.msa_fetch import fetch_msa, MsaFetchError  # noqa: E402
from scripts.build_holdout_hybrid_manifest import assay_seq  # noqa: E402
from scripts.mavedb_prospective_holdout import parse_hgvs_pro  # noqa: E402

API = "https://api.mavedb.org/api/v1"
MANIFEST = Path("wiki/holdout_hybrid_manifest.json")
OUT = Path("wiki/gemme_holdout_tables.json")
GEMME_WORK = Path("D:/dna_decode_cache/gemme_work")


def dms_variants(urn: str) -> set[str]:
    """The assay's single-missense DMS variant keys 'wt{pos}alt' (assay-local numbering)."""
    with urllib.request.urlopen(f"{API}/score-sets/{urn}/scores", timeout=120) as r:
        text = r.read().decode("utf-8", "replace")
    out = set()
    for row in csv.DictReader(io.StringIO(text)):
        p = parse_hgvs_pro(row.get("hgvs_pro", ""))
        if p:
            out.add(f"{p[0]}{p[1]}{p[2]}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=30, help="max held-out proteins to GEMME-score")
    a = ap.parse_args()

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    tables = json.loads(OUT.read_text(encoding="utf-8")) if OUT.exists() else {}
    # process the largest-DMS proteins first (most informative); skip already-done
    order = sorted(manifest.items(), key=lambda kv: -(kv[1].get("seq_len") or 0))
    print(f"held-out proteins: {len(manifest)} | already GEMME-scored: {len(tables)}", flush=True)

    for urn, meta in order:
        if len([u for u in tables]) >= a.limit:
            break
        if urn in tables:
            continue
        g = meta.get("gene")
        try:
            seq = assay_seq(urn)
            if not seq:
                print(f"  {g:12s} SKIP no seq", flush=True)
                continue
            msa = fetch_msa(seq)                                   # ColabFold, cached to D:
            depth = sum(1 for l in msa.read_text(encoding='utf-8').splitlines() if l.startswith('>'))
            full = run_gemme(msa, seq, work_dir=GEMME_WORK / urn.replace(":", "_"))
            keep = dms_variants(urn)
            tbl = {m: full[m] for m in keep if m in full}
            tables[urn] = {"gene": g, "seq_len": len(seq), "msa_depth": depth,
                           "n_gemme": len(full), "n_kept": len(tbl), "table": tbl}
            OUT.write_text(json.dumps(tables), encoding="utf-8")   # checkpoint per protein (restartable)
            print(f"  {g:12s} OK depth={depth} gemme={len(full)} kept={len(tbl)}  [{len(tables)} done]", flush=True)
        except (GemmeUnavailable, MsaFetchError) as e:
            print(f"  {g:12s} SKIP {type(e).__name__}: {str(e)[:70]}", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"  {g:12s} ERR {type(e).__name__}: {str(e)[:70]}", flush=True)

    print(f"\ngemme tables: {OUT} ({len(tables)} proteins)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
