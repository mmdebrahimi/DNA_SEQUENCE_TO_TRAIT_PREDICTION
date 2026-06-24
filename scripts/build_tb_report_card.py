"""Standing report card for the M. tuberculosis cell — the durable, visible TB validation surface.

NAMESPACE-SEPARATE from the frozen NCBI-PD bacterial card, the AMR-Portal bacterial card, the HIV card, and
the external-cohort card (the shared-key silent-overwrite lesson: TB is the NON-frozen `organism_rules` cell,
its phenotype is BMD-MIC / genomic-VCF — not the bacterial card's CLSI-broth-microdilution honesty text — and
its independence is constructed differently; it must NOT be keyed into the frozen `canonical_cell_key` card).

Read-only roll-up (exit 0 always — a REPORT, not a gate) of the two TB validation tiers:
  IN-DISTRIBUTION    `wiki/tb_{rif,inh}_cryptic_parquet_baseline_*.json` (WHO catalogue scored on CRyPTIC;
                     the catalogue was built partly FROM CRyPTIC -> a knowledge baseline, NOT independent).
  INDEPENDENT        `wiki/tb_independent_amr_portal_scores.json` (+ `..._lineage_collapsed.json`) — the WHO
                     rule scored on the EBI AMR Portal provenance-disjoint, measured-AST cohort (free; no DUA).

HONEST TIER PER DRUG: the RAW per-isolate sens/spec is the headline for TB AMR (resistance is HOMOPLASIC —
acquired independently within sublineages — so a lineage-MAJORITY collapse measures the wrong question; the
lineage figure is a clonality DISCLOSURE, never the headline). Mirrors the project's bacterial
disclose-not-demote discipline. -> `wiki/tb_report_card.{md,json}`.
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WIKI = REPO / "wiki"
OUT_MD = WIKI / "tb_report_card.md"
OUT_JSON = WIKI / "tb_report_card.json"


def _latest(pattern: str) -> Path | None:
    fs = sorted(glob.glob(str(WIKI / pattern)))
    return Path(fs[-1]) if fs else None


def _indist_rows() -> list[dict]:
    out = []
    for drug, code in (("rifampicin", "rif"), ("isoniazid", "inh")):
        p = _latest(f"tb_{code}_cryptic_parquet_baseline_*.json")
        if not p:
            continue
        d = json.loads(p.read_text(encoding="utf-8"))
        raw, lc = d.get("raw", {}), d.get("lineage_collapsed", {})
        out.append({"tier": "IN_DISTRIBUTION_KNOWLEDGE_BASELINE", "drug": drug,
                    "n": d.get("n_isolates"), "raw_sens": raw.get("sens"), "raw_spec": raw.get("spec"),
                    "lineage_sens": lc.get("sens"), "lineage_spec": lc.get("spec"),
                    "source": "WHO catalogue on CRyPTIC (catalogue built partly from CRyPTIC -> in-distribution)"})
    return out


def _indep_rows() -> list[dict]:
    sp = WIKI / "tb_independent_amr_portal_scores.json"
    lcp = WIKI / "tb_independent_amr_portal_lineage_collapsed.json"
    if not sp.exists():
        return []
    scores = json.loads(sp.read_text(encoding="utf-8"))
    lin = json.loads(lcp.read_text(encoding="utf-8")) if lcp.exists() else {"drugs": {}}
    out = []
    for drug, m in scores.get("drugs", {}).items():
        lc = (lin.get("drugs", {}).get(drug, {}) or {}).get("lineage_collapsed", {})
        out.append({"tier": "PROVENANCE_DISJOINT_INDEPENDENT", "drug": drug,
                    "n": m.get("n_R", 0) + m.get("n_S", 0), "n_R": m.get("n_R"), "n_S": m.get("n_S"),
                    "raw_sens": m.get("sens"), "raw_spec": m.get("spec"), "raw_acc": m.get("accuracy"),
                    "raw_sens_ci95": m.get("sens_ci95"), "raw_spec_ci95": m.get("spec_ci95"),
                    "lineage_sens_disclosure": lc.get("sens"), "lineage_spec_disclosure": lc.get("spec"),
                    "source": "EBI AMR Portal (CABBAGE) provenance-disjoint, measured AST (free; no DUA)"})
    return out


def build() -> dict:
    indist, indep = _indist_rows(), _indep_rows()
    return {"schema": "tb-report-card-v1",
            "headline_rule": "RAW per-isolate is the TB-AMR headline; lineage is a clonality DISCLOSURE "
                             "(resistance is homoplasic). NAMESPACE-SEPARATE from the frozen bacterial cards.",
            "n_independent_drugs": len(indep), "n_in_distribution_drugs": len(indist),
            "in_distribution": indist, "independent": indep}


def _f(x):
    return f"{x:.3f}" if isinstance(x, (int, float)) else "—"


def _ci(c):
    return f" [{c[0]:.3f}, {c[1]:.3f}]" if isinstance(c, (list, tuple)) and len(c) == 2 else ""


def render_md(card: dict) -> str:
    L = ["# M. tuberculosis — validation report card",
         "",
         "Standing TB validation surface, **namespace-separate** from the frozen NCBI-PD / AMR-Portal / HIV / "
         "external cards (TB is the non-frozen `organism_rules` cell; its phenotype is BMD-MIC / genomic-VCF "
         "and its independence is constructed differently — it must never be keyed into the frozen "
         "`canonical_cell_key` card). **Headline = RAW per-isolate sens/spec** (TB resistance is HOMOPLASIC, so "
         "a lineage-majority collapse measures the wrong question; the lineage figure is a clonality "
         "DISCLOSURE, never the headline — mirrors the bacterial disclose-not-demote discipline).",
         ""]
    if card["independent"]:
        L += ["## INDEPENDENT — provenance-disjoint, measured AST (the gold-set saga's resolved holy grail)",
              "Free (no DUA, no author-contact). WHO-2023 catalogue rule applied UNCHANGED; isolates the rule "
              "was never tuned on.",
              "| Drug | n (R/S) | RAW sens (95% CI) | RAW spec (95% CI) | RAW acc | lineage disclosure (NOT headline) |",
              "|---|---|---|---|---|---|"]
        for r in card["independent"]:
            L.append(f"| {r['drug']} | {r.get('n_R')}/{r.get('n_S')} | {_f(r['raw_sens'])}{_ci(r.get('raw_sens_ci95'))} | "
                     f"{_f(r['raw_spec'])}{_ci(r.get('raw_spec_ci95'))} | {_f(r.get('raw_acc'))} | "
                     f"sens {_f(r.get('lineage_sens_disclosure'))} / spec {_f(r.get('lineage_spec_disclosure'))} |")
        L += ["", "Source: EBI AMR Portal (CABBAGE) provenance-disjoint cohort + assembly→H37Rv-VCF→WHO-rule "
              "pipeline (`scripts/run_tb_independent_amr_portal.py`). Memo "
              "`wiki/tb_independent_number_2026-06-23.md`; homoplasy/disclosure rationale "
              "`wiki/tb_independent_lineage_finding_2026-06-23.md`.", ""]
    if card["in_distribution"]:
        L += ["## IN-DISTRIBUTION — knowledge baseline (NOT independent)",
              "The WHO catalogue was built partly FROM CRyPTIC, so a CRyPTIC-scored number is in-distribution. "
              "Shown for comparison; the independent number above is the real external test.",
              "| Drug | n | RAW sens / spec | lineage (disclosure) |",
              "|---|---|---|---|"]
        for r in card["in_distribution"]:
            L.append(f"| {r['drug']} | {r.get('n')} | {_f(r['raw_sens'])} / {_f(r['raw_spec'])} | "
                     f"{_f(r['lineage_sens'])} / {_f(r['lineage_spec'])} |")
        L += ["", "Source: `wiki/tb_{rif,inh}_cryptic_parquet_baseline_*.json` "
              "(`wiki/tb_cryptic_parquet_baseline_2026-06-22.md`).", ""]
    L += ["## Honesty rails",
          "- **Independence is BioSample-resolution-CHECKED (upgraded 2026-06-23, `wiki/tb_independence_"
          "biosample_check.json`).** The ENA-side disjoint isolates (1,364 with an `ERS` accession) are "
          "already BioSample-grade — their `ERS` is string-matched DIRECTLY against CRyPTIC's `ENA_SAMPLE` "
          "(ERS), the same namespace. The NCBI-side (1,480 `SAMN`) are the only cross-archive risk vs the "
          "European CRyPTIC set; a bounded ENA-portal probe found **0/30** of their ENA-mirror accessions in "
          "CRyPTIC. The one irreducible residual is genomic RE-SUBMISSION (an isolate sequenced twice as "
          "distinct BioSample records) — which needs Mash genomic dedup, NOT accession resolution.",
          "- **Measured phenotype = non-circular** (BMD-MIC / measured DST); **WHO rule applied UNCHANGED.**",
          "- **RAW is the headline; lineage is disclosure** (homoplasy). The independent lineage figures "
          "(~0.44 / 0.32) match the in-distribution lineage figures (0.41 / 0.349), confirming the clonal "
          "structure of the R classes.",
          "- FROZEN bacterial AMR surface byte-unchanged; this card is READ-only. Rebuild: "
          "`uv run python scripts/build_tb_report_card.py`."]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    card = build()
    OUT_JSON.write_text(json.dumps(card, indent=2, default=str), encoding="utf-8")
    OUT_MD.write_text(render_md(card), encoding="utf-8")
    print(f"[tb-card] independent drugs: {card['n_independent_drugs']} | in-distribution: "
          f"{card['n_in_distribution_drugs']} -> {OUT_MD.name}, {OUT_JSON.name}")
    for r in card["independent"]:
        print(f"  INDEPENDENT {r['drug']}: raw sens={_f(r['raw_sens'])} spec={_f(r['raw_spec'])} "
              f"acc={_f(r.get('raw_acc'))} (n={r.get('n')})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
