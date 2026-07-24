"""Build the standing HCMV decoder report card (read-only; exit 0 always — a report, not a gate).

Rolls up the shipped HCMV cell (`dna_decode/data/hcmv_amr.py`) into `wiki/hcmv_decoder_report_card.{md,json}`
+ the packaged copy the wheel carries, so the CLI trust surface renders a REAL tier for HCMV drugs instead
of UNKNOWN. The honest tier is **IN_DISTRIBUTION** (knowledge baseline): the catalog is curated FROM the
measured recombinant fold-changes it is scored against, so this is NOT independent validation — the SAME
posture as the SARS-CoV-2 CoV-RDB cell. An INDEPENDENT number is a CLOSED negative for free data: HCMV
phenotyping is per-mutation recombinant marker-transfer and the Chou compilations are its comprehensive
consensus, so no held-out per-isolate measured-phenotype resource disjoint from the catalog exists (research
2026-07-24, `wiki/free_label_acquisition_round2_2026-07-23.md`).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.hcmv_amr import (  # noqa: E402
    BENIGN_BY_GENE, GENES_FOR_DRUG, RESISTANCE_BY_GENE, all_supported_hcmv_drugs,
)

WIKI = REPO / "wiki"
PKG_CARDS = REPO / "dna_decode" / "report_cards"
TIER = "IN_DISTRIBUTION"
_POINT = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$")


def build() -> dict:
    genes = {}
    for g in ("UL97", "UL54", "UL56"):
        res = RESISTANCE_BY_GENE[g]
        genes[g] = {
            "n_resistance": len(res),
            "n_point": sum(1 for m in res if _POINT.match(m)),
            "n_indel": sum(1 for m in res if not _POINT.match(m)),
            "n_benign": len(BENIGN_BY_GENE[g]),
        }
    cells = []
    for drug in all_supported_hcmv_drugs():
        tgt = list(GENES_FOR_DRUG[drug])
        cells.append({
            "drug": drug, "organism": "HCMV", "genes": tgt, "tier": TIER,
            "state": "IN_DISTRIBUTION",
            "n_resistance": sum(genes[g]["n_resistance"] for g in tgt),
            "n_benign": sum(genes[g]["n_benign"] for g in tgt),
        })
    return {
        "card": "hcmv_decoder_report_card",
        "organism": "HCMV",
        "tier": TIER,
        "headline": "in-distribution (Chou recombinant fold-change); NOT independent",
        "n_cells": len(cells),
        "genes": genes,
        "provenance": ("catalog curated VERBATIM from Chou recombinant-phenotyping literature "
                       "(PMC3262590 / PMC5483911 / AAC 10.1128/aac.00922-18 / PMC9759347); AD169 numbering; "
                       "genome mode ref = Merlin NC_006273.2, integrity-gated"),
        "independence": ("CLOSED for free data (2026-07-24): HCMV phenotyping is per-mutation recombinant "
                         "marker-transfer + the Chou compilations are its comprehensive consensus, so no "
                         "held-out per-isolate measured-phenotype set disjoint from the catalog exists "
                         "(structurally in-distribution, like SARS-CoV-2 CoV-RDB)"),
        "cells": cells,
    }


def render_md(rc: dict, generated: str) -> str:
    lines = [f"# HCMV decoder report card ({generated})", "",
             f"**Tier: {rc['tier']}** — {rc['headline']}.", "",
             "Catalog per gene:", ""]
    lines.append("| gene | resistance (point / indel) | benign |")
    lines.append("|---|---|---|")
    for g, s in rc["genes"].items():
        lines.append(f"| {g} | {s['n_resistance']} ({s['n_point']} / {s['n_indel']}) | {s['n_benign']} |")
    lines += ["", "Cells (per drug):", "", "| drug | genes | tier |", "|---|---|---|"]
    for c in rc["cells"]:
        lines.append(f"| {c['drug']} | {'/'.join(c['genes'])} | {c['tier']} |")
    lines += ["", f"**Provenance:** {rc['provenance']}", "",
              f"**Independence:** {rc['independence']}", ""]
    return "\n".join(lines)


def main(argv=None) -> int:
    import datetime
    today = datetime.date.today().isoformat()
    rc = build()
    md = render_md(rc, today)
    for d in (WIKI, PKG_CARDS):
        d.mkdir(parents=True, exist_ok=True)
        (d / "hcmv_decoder_report_card.json").write_text(json.dumps(rc, indent=2), encoding="utf-8")
    (WIKI / "hcmv_decoder_report_card.md").write_text(md, encoding="utf-8")
    print(f"HCMV report card: tier={rc['tier']}  cells={rc['n_cells']}  -> wiki/ + dna_decode/report_cards/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
