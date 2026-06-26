"""Certification capstone — ONE legible presentation of the whole validated decoder surface.

A THIN PRESENTATION LAYER over existing milestone evidence — NOT a new reducer gate. It reads the
already-shipped evidence (the Evidence-Contract Registry as the per-cell spine + the per-domain report
cards + the reproducibility freeze + the negative-results map) and lays them out in one document so a
reader sees the certified surface at a glance.

LOAD-BEARING HONESTY RAIL (the "boolean-verdict-on-unvalidated-model trap"): this capstone emits NO
aggregate pass/fail, NO "certified: true", NO averaged score. Reducing a surface of cells at DIFFERENT
honest tiers (independent_measured / near_independent / faithful_to_tool / knowledge_baseline /
not_censused / no_free_source) to one boolean would CERTIFY the weakest cell as strongly as the
strongest — the exact overclaim the registry exists to prevent. Each cell + each card keeps its own tier
VERBATIM; the capstone presents, it does not judge. `no_aggregate_verdict` is stamped True in the JSON.

It supersets `build_cross_kingdom_summary.py` (which predates the registry and covers only the 5 AMR-family
cards) by anchoring on the registry (67 cells / 5 tracks incl PGx/typing/finder + NOT_CENSUSED) and adding
the freeze + negative-map BOUNDARIES (what is closed, and the label wall).

Exit 0 always — a REPORT, not a gate.  ->  wiki/certification_capstone.{md,json}
"""
from __future__ import annotations

import datetime
import json
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WIKI = REPO / "wiki"
sys.path.insert(0, str(REPO))

from dna_decode.data.cell_registry import cells  # noqa: E402

# Per-domain report cards + the small set of headline keys to surface VERBATIM (no normalization — the
# heterogeneous shapes are intentional, each card's independence construction is distinct).
_CARDS: list[tuple[str, list[str]]] = [
    ("decoder_validation_report_card.json", ["honest_tier", "no_aggregate_headline", "state_counts"]),
    ("amr_portal_independent_report_card.json", ["n_cells", "n_scored_independent", "n_underpowered", "status_field"]),
    ("external_validation_report_card.json", ["honest_tier", "state_counts", "note"]),
    ("tb_report_card.json", ["headline_rule", "n_independent_drugs", "n_in_distribution_drugs"]),
    ("hiv_decoder_report_card.json", ["label_independence", "n_cells", "honest_caveats"]),
    ("pgx_report_card.json", ["note", "sources"]),
]

# Boundary documents — what is CLOSED + why (the label wall). Presented as references, not flattened.
_BOUNDARIES: list[str] = [
    "reproducibility_freeze_2026-06-13.md",
    "negative_results_map_2026-06-13.md",
]


def _load_json(name: str) -> dict | None:
    p = WIKI / name
    try:
        return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
    except Exception:
        return None


def _title_of(md_name: str) -> str | None:
    p = WIKI / md_name
    if not p.exists():
        return None
    first = p.read_text(encoding="utf-8").splitlines()[:1]
    return first[0].lstrip("# ").strip() if first else None


def _registry_census() -> dict:
    cs = list(cells())
    by_track = dict(Counter(c.track for c in cs))
    by_tier = dict(Counter(c.evidence_tier.value for c in cs))
    not_censused = sorted(c.cell_id for c in cs if c.evidence_tier.value == "not_censused")
    return {
        "total_cells": len(cs),
        "by_track": by_track,
        "by_evidence_tier": by_tier,
        "not_censused_cells": not_censused,  # surfaced explicitly — routable but unvalidated
    }


def _card_headlines() -> list[dict]:
    out = []
    for name, keys in _CARDS:
        d = _load_json(name)
        if d is None:
            out.append({"card": name, "status": "NOT_RUN"})
            continue
        out.append({"card": name, "status": "present",
                    "headline": {k: d.get(k) for k in keys if k in d}})
    return out


def _boundaries() -> list[dict]:
    return [{"document": b, "title": _title_of(b),
             "status": "present" if (WIKI / b).exists() else "MISSING"} for b in _BOUNDARIES]


def main() -> int:
    census = _registry_census()
    cap = {
        "schema": "certification-capstone-v0",
        "analysis_date": datetime.date.today().isoformat(),
        "what_this_is": ("A thin presentation layer over existing milestone evidence — NOT a clinical tool, "
                         "NOT a new gate, NOT an aggregate verdict. Each cell + card keeps its own honest tier."),
        "no_aggregate_verdict": True,  # LOAD-BEARING: this capstone never reduces the surface to a boolean
        "aggregate_verdict_disclaimer": ("Cells span 6 honest tiers; an aggregate pass/fail would certify the "
                                         "weakest cell as strongly as the strongest. No such field exists by design."),
        "registry_census": census,
        "domain_report_cards": _card_headlines(),
        "boundaries": _boundaries(),
        "sources": {"registry": "dna_decode/data/cell_registry.py",
                    "supersedes_partially": "scripts/build_cross_kingdom_summary.py (AMR-family only, pre-registry)"},
    }
    (WIKI / "certification_capstone.json").write_text(json.dumps(cap, indent=2), encoding="utf-8")

    L = [f"# Decoder certification capstone ({cap['analysis_date']})", "",
         f"> _{cap['what_this_is']}_", "",
         "**No aggregate verdict.** " + cap["aggregate_verdict_disclaimer"], "",
         "## Registry census (the per-cell spine)", "",
         f"- **Total cells:** {census['total_cells']}",
         f"- **By track:** " + ", ".join(f"{k}={v}" for k, v in sorted(census["by_track"].items())),
         f"- **By honest evidence tier:** " + ", ".join(f"{k}={v}" for k, v in sorted(census["by_evidence_tier"].items())),
         f"- **Routable-but-NOT_CENSUSED ({len(census['not_censused_cells'])}):** "
         + (", ".join(census["not_censused_cells"]) or "none"), "",
         "## Per-domain validated headlines (verbatim — NOT averaged)", "",
         "| card | status | headline |", "|---|---|---|"]
    for c in cap["domain_report_cards"]:
        hl = json.dumps(c.get("headline", {}), ensure_ascii=True) if c["status"] == "present" else "—"
        L.append(f"| {c['card']} | {c['status']} | {hl} |")
    L += ["", "## Boundaries (what is closed + the label wall)", ""]
    for b in cap["boundaries"]:
        L.append(f"- **{b['document']}** ({b['status']}): {b['title'] or '—'}")
    L += ["", "_A presentation of existing evidence; it certifies nothing the underlying cards do not. "
          "NOT a clinical tool._", ""]
    (WIKI / "certification_capstone.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print("[capstone -> wiki/certification_capstone.{md,json}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
