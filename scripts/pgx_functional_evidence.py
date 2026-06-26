"""Emit the PGx independent functional-evidence report card (Unit A).

Reads the curated `dna_decode.pgx.functional_evidence.EVIDENCE` (independent, non-CPIC signals per allele)
and writes `wiki/pgx_functional_evidence_2026-06-25.{md,json}`: the per-allele AGREE/DISAGREE/FLAG/NO_SIGNAL
cross-check vs CPIC's function. Offline (no fetch) -- the signals were grounded at curation time.
"""
from __future__ import annotations

import datetime
import json
from dataclasses import asdict
from pathlib import Path

from dna_decode.pgx.functional_evidence import EVIDENCE, SCHEMA, summary

REPO = Path(__file__).resolve().parent.parent


def main() -> int:
    s = summary()
    rep = {
        "schema": SCHEMA, "analysis_date": datetime.date.today().isoformat(),
        "purpose": ("Independent (non-CPIC) functional-evidence cross-check of each PGx allele's CPIC "
                    "function assignment -- the circularity-break. AGREE raises confidence; DISAGREE/FLAG "
                    "surface where 'faithful-to-CPIC' rests on clinical evidence sequence-signals miss."),
        "honesty": ("Independent signals are ORTHOGONAL to CPIC curation, NOT ground truth: missense -> "
                    "Ensembl VEP predictors (ML); stop/splice -> consequence class (fact); regulatory -> "
                    "documented expression effect (measured, primary literature). Small-N per-allele "
                    "annotation, NOT a concordance-%; GTEx-eQTL confirmation deferred (not asserted)."),
        "summary": s,
        "evidence": [asdict(e) for e in EVIDENCE],
    }
    out_json = REPO / "wiki" / "pgx_functional_evidence_2026-06-25.json"
    out_md = REPO / "wiki" / "pgx_functional_evidence_2026-06-25.md"
    out_json.write_text(json.dumps(rep, indent=2), encoding="utf-8")

    L = [f"# PGx independent functional-evidence cross-check ({rep['analysis_date']})", "",
         f"_{rep['purpose']}_", "",
         f"**Verdicts:** AGREE {s['AGREE']} / DISAGREE {s['DISAGREE']} / FLAG {s['FLAG']} / "
         f"NO_SIGNAL {s['NO_SIGNAL']}  (n={s['n']})", "",
         "| gene | allele | rsid | CPIC function | variant class | independent signal | verdict |",
         "|---|---|---|---|---|---|---|"]
    for e in EVIDENCE:
        L.append(f"| {e.gene} | {e.allele} | {e.rsid} | {e.cpic_function} | {e.variant_class} | "
                 f"{e.independent_signal} | **{e.verdict}** |")
    L += ["", "## Notes (the informative cases)", ""]
    for e in EVIDENCE:
        if e.verdict in ("DISAGREE", "FLAG"):
            L.append(f"- **{e.gene} {e.allele} ({e.verdict}):** {e.note}")
    L += ["", f"_{rep['honesty']}_", ""]
    out_md.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:10]))
    print(f"[report -> {out_md} + {out_json}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
