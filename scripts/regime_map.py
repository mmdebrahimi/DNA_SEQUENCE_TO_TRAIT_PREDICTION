"""Emit the measured g->p regime map, and screen a proposal against it.

Read-only report; exit 0 when every regime's cited artifact resolves, 1 when one does not (a regime
whose evidence file is missing is not a regime -- it is a memory, which is the failure this module
exists to prevent).

    uv run python scripts/regime_map.py
    uv run python scripts/regime_map.py --screen natural organism zero_shot
    uv run python scripts/regime_map.py --screen constructed molecular supervised --catalog-exists

WHY. The boundary was mis-stated three times, always by compressing a SCOPED negative into a general
one. Prose went stale; this re-derives from the artifacts each run.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
WIKI = ROOT / "wiki"


def main() -> int:
    from dna_decode.eval.regime import REGIMES, screen_proposal

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--screen", nargs=3, metavar=("POPULATION", "ENDPOINT", "METHOD"),
                    help="screen one proposal instead of emitting the map")
    ap.add_argument("--catalog-exists", action="store_true",
                    help="a curated determinant catalog already covers this endpoint")
    args = ap.parse_args()

    if args.screen:
        res = screen_proposal(*args.screen, curated_catalog_exists=args.catalog_exists)
        print(json.dumps(res.as_dict(), indent=2))
        print(f"\n{'REFUSED' if res.refused else 'not refused'}: {res.verdict}")
        return 1 if res.refused else 0

    rows, missing = [], []
    for r in REGIMES:
        exists = (ROOT / r.artifact).exists()
        if not exists:
            missing.append(r.key)
        rows.append({**r.as_dict(), "artifact_exists": exists})

    print(f"\n{len(rows)} measured regimes\n")
    print(f"  {'population':12} {'endpoint':28} {'method':22} verdict")
    for r in rows:
        flag = "" if r["artifact_exists"] else "   [ARTIFACT MISSING]"
        print(f"  {r['population']:12} {r['endpoint']:28} {r['method']:22} {r['verdict']}{flag}")
    print("\nThe discriminating variable is POPULATION DESIGN, not organism complexity.")
    print("The natural-population negative is ZERO-SHOT-scoped; a supervised proposal gets conditions,")
    print("not a refusal -- compressing that scope has hidden a live direction three separate times.")

    out = {"schema": "learned-regime-map-v1", "generated": date.today().isoformat(),
           "discriminating_variable": "population design (constructed vs natural), NOT organism complexity",
           "scope_warning": ("the natural-population negative is ZERO-SHOT-scoped. A supervised "
                             "complement is the shipped architecture; refusing it would re-commit the "
                             "over-generalisation this map exists to prevent."),
           "regimes": rows, "artifacts_missing": missing}
    (WIKI / "learned_regime_map.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("\nwrote wiki/learned_regime_map.json")
    if missing:
        print(f"REFUSING to certify the map: cited artifacts missing for {missing}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
