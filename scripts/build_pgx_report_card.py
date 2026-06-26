"""Standing PGx trust-surface report card -- one honest roll-up of every shipped human-PGx cell.

Read-only consolidation (exit 0 always; a report, NOT a gate) -- the PGx analogue of the AMR
`decoder_validation_report_card`. Rows = the deployed PGx genes; columns = the validation evidence already
produced (GeT-RM consensus concordance on real 1000G, PharmCAT fixtures, independent functional-evidence
verdicts, trio Mendelian QC). NO aggregate headline; each cell's honest tier stands on its own. Reads the
committed `wiki/pgx_*` JSON sidecars; a missing sidecar renders as NOT_RUN, never a fabricated number.
"""
from __future__ import annotations

import datetime
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
W = REPO / "wiki"


def _load(name):
    p = W / name
    try:
        return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
    except Exception:
        return None


def main() -> int:
    getrm_c19 = _load("pgx_getrm_concordance_2026-06-25.json")
    getrm_c9 = _load("pgx_getrm_concordance_cyp2c9_2026-06-25.json")
    pharmcat = _load("pgx_cyp2c19_report_card.json")
    fe = _load("pgx_functional_evidence_2026-06-25.json")
    trio = _load("pgx_trio_mendelian_2026-06-25.json")

    fe_by_gene = {}
    for e in (fe or {}).get("evidence", []):
        fe_by_gene.setdefault(e["gene"], []).append(e["verdict"])

    def fe_summ(g):
        v = fe_by_gene.get(g)
        if not v:
            return "—"
        from collections import Counter
        c = Counter(v)
        return f"A{c['AGREE']}/D{c['DISAGREE']}/F{c['FLAG']}"

    def trio_summ(g):
        r = (trio or {}).get("genes", {}).get(g)
        return (f"{r['consistent']}/{r['n_callable']}" if r and r.get("status") == "ok" else "—")

    cells = [
        {"gene": "CYP2C19", "trait": "metabolizer phenotype (PM/IM/NM/RM/UM)",
         "getrm": getrm_c19 and getrm_c19.get("core_diplotype_hits"),
         "getrm_pct": getrm_c19 and getrm_c19.get("core_diplotype_concordance"),
         "pharmcat": pharmcat and pharmcat.get("core_diplotype_hits"),
         "functional_evidence": fe_summ("CYP2C19"), "trio_mendelian": trio_summ("CYP2C19"),
         "tier": "GeT-RM consensus (independent of consensus tools) + PharmCAT fixtures; phenotype faithful-to-CPIC",
         "residual": "non-core *4/*35 withheld (sentinel v0.1)"},
        {"gene": "CYP2C9", "trait": "metabolizer phenotype (activity-score)",
         "getrm": getrm_c9 and getrm_c9.get("core_diplotype_hits"),
         "getrm_pct": getrm_c9 and getrm_c9.get("core_diplotype_concordance"),
         "pharmcat": None,
         "functional_evidence": fe_summ("CYP2C9"), "trio_mendelian": trio_summ("CYP2C9"),
         "tier": "GeT-RM consensus; phenotype faithful-to-CPIC (activity-score)",
         "residual": "non-core *5/*8/*9/*11 withheld (sentinel v0.1); *6-indel/*61 residual"},
        {"gene": "VKORC1", "trait": "warfarin sensitivity (rs9923231)",
         "getrm": None, "getrm_pct": None, "pharmcat": None,
         "functional_evidence": fe_summ("VKORC1"), "trio_mendelian": "—",
         "tier": "single-SNP genotype->sensitivity (minus-strand encoded); not a star/diplotype system",
         "residual": "—"},
    ]

    rep = {
        "schema": "pgx-report-card-v0", "analysis_date": datetime.date.today().isoformat(),
        "note": ("Standing PGx trust surface -- a roll-up, NOT a gate (exit 0 always). No aggregate "
                 "headline; each cell's honest tier stands alone. CALLING is independently validatable vs "
                 "GeT-RM (free consensus panel); PHENOTYPE is faithful-to-CPIC (assigned, not measured)."),
        "cells": cells,
        "sources": {"getrm_cyp2c19": bool(getrm_c19), "getrm_cyp2c9": bool(getrm_c9),
                    "pharmcat_cyp2c19": bool(pharmcat), "functional_evidence": bool(fe),
                    "trio_mendelian": bool(trio)},
    }
    (W / "pgx_report_card.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")

    L = [f"# PGx decoder report card ({rep['analysis_date']})", "", f"_{rep['note']}_", "",
         "| gene | trait | GeT-RM core | PharmCAT | func-evidence (A/D/F) | trio Mendelian | residual |",
         "|---|---|---|---|---|---|---|"]
    for c in cells:
        g = f"{c['getrm']} ({c['getrm_pct']})" if c["getrm"] else "—"
        L.append(f"| {c['gene']} | {c['trait']} | {g} | {c['pharmcat'] or '—'} | "
                 f"{c['functional_evidence']} | {c['trio_mendelian']} | {c['residual']} |")
    L += ["", "## Honest tier per cell", ""]
    for c in cells:
        L.append(f"- **{c['gene']}:** {c['tier']}")
    L += ["", "_Validation axes: GeT-RM = consensus concordance on real 1000G (independent of the consensus "
          "tools); PharmCAT = reference-tool fixtures; func-evidence = non-CPIC cross-check of the function "
          "assignment (AGREE/DISAGREE/FLAG); trio = Mendelian calling-consistency on 1000G trios. "
          "NOT a clinical tool._", ""]
    (W / "pgx_report_card.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print("[report -> wiki/pgx_report_card.{md,json}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
