"""Capture the OpenGWAS study index (gwasinfo, ~50k GWAS datasets) — the auth-gated 5th summary-stat source.

OpenGWAS (MRC IEU) requires a free JWT token since 2024-05 (a USER credential — Soraya cannot generate it).
This is the one-step-ready capture: supply the token (env OPENGWAS_JWT, --token, or a --token-file) and it
queries /api/gwasinfo -> a committed compact study index (id, trait, year, sample_size, population, category,
ncase, ncontrol, nsnp). Same class as the other summary-stat indexes (data/summary_stat_sources/).

Get a token: log in at https://api.opengwas.io/profile/ -> generate token. Then:
    OPENGWAS_JWT=<token> uv run python scripts/capture_opengwas_index.py
The token is a password-equivalent; it is NEVER written to disk by this script (only used for the request).
"""
from __future__ import annotations

import json
import os
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "summary_stat_sources" / "opengwas_study_index.tsv"
_API = "https://api.opengwas.io/api/gwasinfo"
_KEEP = ["id", "trait", "year", "sample_size", "population", "category", "subcategory",
         "ncase", "ncontrol", "nsnp", "sex", "build"]


def _token(arg: str | None, token_file: str | None) -> str | None:
    if arg:
        return arg.strip()
    if token_file and Path(token_file).exists():
        return Path(token_file).read_text(encoding="utf-8").strip()
    return (os.environ.get("OPENGWAS_JWT") or "").strip() or None


def capture(token: str, out: Path) -> dict:
    req = urllib.request.Request(_API, headers={"Authorization": f"Bearer {token}",
                                                "Content-Type": "application/json",
                                                "User-Agent": "dna_decode/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = json.loads(r.read().decode("utf-8"))
    # gwasinfo returns either a dict keyed by id, or a list of study dicts
    studies = list(data.values()) if isinstance(data, dict) else data
    rows = [{k: str(s.get(k, "")).replace("\t", " ").replace("\n", " ")[:80] for k in _KEEP} for s in studies]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(_KEEP) + "\n")
        for r_ in rows:
            f.write("\t".join(r_[k] for k in _KEEP) + "\n")
    return {"source": "OpenGWAS", "n_studies": len(rows), "out": str(out),
            "bulk_recipe": "per-study harmonised VCF sumstats via the OpenGWAS API (associations/tophits) with the token"}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--token", default=None, help="OpenGWAS JWT (else OPENGWAS_JWT env or --token-file)")
    ap.add_argument("--token-file", default=None)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    tok = _token(a.token, a.token_file)
    if not tok:
        print("BLOCKED: no OpenGWAS JWT. Generate one (free) at https://api.opengwas.io/profile/ then:\n"
              "  OPENGWAS_JWT=<token> uv run python scripts/capture_opengwas_index.py", file=sys.stderr)
        return 3
    try:
        res = capture(tok, a.out)
    except Exception as e:
        print(f"ERROR querying OpenGWAS (token invalid/expired or API down?): {e}", file=sys.stderr)
        return 1
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
