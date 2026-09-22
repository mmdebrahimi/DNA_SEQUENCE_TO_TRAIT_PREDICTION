"""What this project IS, derived live. Run this before making any claim about project scope.

WHY THIS EXISTS. On 2026-08-29 a session spent a full day reporting "27 cells / 10 SCORED" as the
project's validated surface. That is the AMR provenance-disjoint arm alone; the registry holds 110 cells
with 28 independently measured. The number was ~3x understated, and it came from prose -- a written
figure in a wiki artifact that was true when written and is a subset now.

The fix is not a better summary. A summary is a written number, and a written number goes stale silently.
Every figure here is COMPUTED from the registries at run time, so it cannot drift from what ships:

  - version + entry points   <- pyproject.toml
  - CLI traits               <- dna_decode.cli.TRAITS
  - evidence tiers           <- dna_decode.data.cell_registry.cells()
  - AMR report card          <- wiki/decoder_validation_report_card.json (LABELLED as the AMR arm only)

READ-ONLY. No network, no Docker, no model load. Exit 0 always -- this is orientation, not a gate.

Run: uv run python scripts/project_status.py
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def version_and_entry_points() -> tuple[str, int]:
    """Parse pyproject WITHOUT a toml dependency -- section-scoped so a stray '=' can't inflate the count."""
    txt = (ROOT / "pyproject.toml").read_text(encoding="utf-8").splitlines()
    version, n, in_scripts = "?", 0, False
    for line in txt:
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            in_scripts = s == "[project.scripts]"
            continue
        if in_scripts and "=" in s and not s.startswith("#"):
            n += 1
        elif version == "?" and s.startswith("version"):
            version = s.split("=", 1)[1].strip().strip('"').strip("'")
    return version, n


def tier_counts() -> tuple[Counter, int]:
    from dna_decode.data.cell_registry import cells

    cs = cells()
    return Counter(c.evidence_tier.value for c in cs), len(cs)


def track_counts() -> Counter:
    from dna_decode.data.cell_registry import cells

    return Counter(c.track for c in cells())


def amr_report_card() -> dict:
    """The AMR arm's card. Returned WITH its scope, because reading it without one is the original bug."""
    f = ROOT / "wiki" / "decoder_validation_report_card.json"
    if not f.exists():
        return {"present": False}
    d = json.loads(f.read_text(encoding="utf-8"))
    rows = d.get("cells") or d.get("rows") or []
    return {"present": True, "n_rows": len(rows),
            "states": Counter(r.get("cell_state") or r.get("state") or "?" for r in rows),
            "generated": d.get("generated") or d.get("date") or "?"}


def n_traits() -> int:
    from dna_decode.cli import TRAITS

    return len(TRAITS)


def main() -> int:
    version, n_entry = version_and_entry_points()
    tiers, n_cells = tier_counts()
    tracks = track_counts()
    card = amr_report_card()

    print("=" * 78)
    print("  dna_decode -- LIVE project status (derived, not written down)")
    print("=" * 78)
    print(f"\n  A published multi-kingdom deterministic genotype->phenotype decoder.")
    print(f"  NOT an E. coli AMR research repo -- AMR is one track of {len(tracks)}.\n")
    print(f"  version              {version}")
    print(f"  console entrypoints  {n_entry}")
    print(f"  CLI traits           {n_traits()}")

    print(f"\n  EVIDENCE SURFACE -- {n_cells} registered cells, all tracks:")
    for tier, k in tiers.most_common():
        print(f"    {tier:24} {k:>4}")

    print(f"\n  BY TRACK:")
    for t, k in tracks.most_common():
        print(f"    {t:24} {k:>4}")

    print(f"\n  AMR PROVENANCE-DISJOINT REPORT CARD  <- ONE ARM, NOT THE WHOLE TOOL")
    if not card["present"]:
        print("    (absent -- rebuild: uv run python scripts/build_validation_report_card.py)")
    else:
        print(f"    rows {card['n_rows']}   generated {card['generated']}")
        for st, k in card["states"].most_common():
            print(f"      {st:30} {k:>4}")
        print(f"\n    This card covers bacterial AMR cells scored against NCBI-PD. Quoting its row")
        print(f"    count as the project's evidence surface understates it by ~{n_cells / max(card['n_rows'], 1):.1f}x.")

    # DERIVED from dna_decode.eval.regime.REGIMES, never restated. This block used to hardcode
    # "TEM-1 genome-edit, rho 0.761" -- and regime.py's own note for that same regime says, verbatim,
    # "do not quote 0.761 as the path's general strength" (the PEAR replication on a second
    # beta-lactamase came in at 0.352, so the honest range is 0.35-0.76). The one script every session
    # is told to run BEFORE making a scope claim was making the exact claim its data layer forbids.
    print()
    print("  THE REGIME MAP (what works, and the variable that decides it):")
    try:
        from dna_decode.eval.regime import REGIMES

        label = {"zero_shot": "zero-shot", "supervised": "supervised",
                 "deterministic_catalog": "catalog"}
        import textwrap

        for r in REGIMES:
            lhs = f"{r.population} + {r.endpoint} ({label.get(r.method, r.method)})"
            print(f"    {lhs:<56} {r.verdict}")
            # WRAPPED, never truncated. A hard [:150] cut this very block's warning mid-sentence at
            # "do not quote 0.761 as " -- silently amputating the load-bearing half of a caveat is
            # the same failure as omitting it.
            for line in textwrap.wrap(r.evidence, 104):
                print(f"        {line}")
            if r.note:
                for i, line in enumerate(textwrap.wrap(r.note, 98)):
                    print(f"        {'NOTE: ' if i == 0 else '      '}{line}")
    except Exception as exc:  # noqa: BLE001 - orientation must not die on an import
        # REFUSE rather than fall back to a remembered summary: a stale regime map is the specific
        # failure this script exists to prevent.
        print(f"    UNAVAILABLE - could not read dna_decode.eval.regime ({exc}).")
        print(f"    Read wiki/organism_gp_regime_correction_2026-08-29.md; do NOT quote from memory.")
    print(f"    The discriminating variable is POPULATION DESIGN, not organism complexity.")
    print(f"\n  BEFORE CLAIMING ANYTHING ABOUT SCOPE: re-derive it. Prose in CLAUDE.md and wiki/ was")
    print(f"  true when written; this script is true now.")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
