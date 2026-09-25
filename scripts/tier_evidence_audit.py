"""Flag cells whose EVIDENCE TIER disagrees with the evidence their own contract describes.

WHY
---
`evidence_tier` records the CLASS of evidence behind a cell. This repo has repeatedly shipped cells whose
tier contradicted their own contract text, in BOTH directions, and every instance was found by reading the
contract against the tier rather than by any test:

  * `finder:any:forward` shipped NOT_CENSUSED / "nothing to score" while DMS-validated (Spearman 0.7315).
  * `pneumoserotype` sat FAITHFUL_TO_TOOL while its own card recorded independent Quellung validation.
  * `plasmid` / `disinfinder` were measured and left on the shared "never measured" default tuple.
  * `pathotype` and `mlst` both carried the never-measured default while holding enforced measured gates.

The repo's own rule is that **under-claiming is as much a trust-surface falsehood as over-claiming**, and
`tests/test_cell_registry.py` pins that rule for exactly ONE cell, by hand. Nothing checks it generally.

WHAT IT IS, HONESTLY
--------------------
A TRIAGE FUNNEL, not an oracle -- same posture as scripts/contract_number_audit.py and the staleness
auditor. Tier assignment is a JUDGEMENT about evidence class (is this comparator a wet-lab label, a
reference tool, or the rule's own source?), and a lexical reading of prose cannot settle it. So this
detects a narrow, checkable CONTRADICTION and emits it for adjudication; it never edits a tier.

The narrow contradiction, in two directions:
  * UNDER-CLAIM SUSPECT -- the tier asserts the cell was never measured, yet the contract reports a
    measured performance figure against a named comparator.
  * OVER-CLAIM SUSPECT -- the tier asserts an independent measurement, yet the contract reports no
    measured figure at all.

Deliberately NOT flagged: which of NEAR_INDEPENDENT vs FAITHFUL_TO_TOOL vs INDEPENDENT_MEASURED a measured
cell belongs in. That turns on whether the comparator is a tool or a wet-lab label -- exactly the judgement
this cannot make, and getting it wrong in the confident direction is how an over-claim ships.
"""
from __future__ import annotations

import datetime
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

PROSE_FIELDS = ("claim", "claim_status", "validation_slice", "label_provenance", "demotion_rule")

# Tiers that ASSERT the cell has not been measured against anything.
NEVER_MEASURED_TIERS = {"NOT_CENSUSED", "KNOWLEDGE_BASELINE", "FAITHFUL_TO_TOOL"}
# Tiers that ASSERT an independent measurement exists.
MEASURED_TIERS = {"INDEPENDENT_MEASURED", "NEAR_INDEPENDENT"}

# A measured PERFORMANCE figure. Requiring a metric NAME adjacent to a number keeps cohort sizes, allele
# frequencies and coordinates from reading as performance -- but the FIRST version of this regex required
# ONLY that, and flagged 46 of 115 cells as over-claiming. A 40% hit rate on a two-sided detector means the
# detector is miscalibrated, not that 40 cells over-claim: this corpus's dominant idiom is a BARE count
# ("85/85 core-comparable", "232/264", "64/64 exact"), which carries no metric word at all.
METRIC_RE = re.compile(
    r"\b(acc(?:uracy)?|sens(?:itivity)?|spec(?:ificity)?|ppv|npv|auroc|auc|f1|mcc|precision|recall|"
    r"spearman|pearson|concordance|agreement|comparable|balacc|r2|top-?\d|lift)\b[^.;|]{0,40}?"
    r"((?<![\w.])\d+\.\d+|\d+\s*/\s*\d+)", re.I)
# A bare count-over-count. In this corpus that is nearly always a concordance/hit count, which is exactly
# the evidence a measured tier rests on.
COUNT_RE = re.compile(r"(?<![\w./])\d{1,6}\s*/\s*\d{1,6}(?![\w./])")
# "r=+0.619" / "r = 0.42" -- a correlation reported without the word spearman/pearson.
CORR_RE = re.compile(r"\br\s*=\s*[+-]?\d*\.\d+", re.I)

# For a FAITHFUL_TO_TOOL cell, being measured against a REFERENCE TOOL is exactly what the tier means, so a
# measured figure is NOT a contradiction there. ktype (agreement 232/264 vs Kaptive) is the worked example:
# correctly FAITHFUL_TO_TOOL precisely because the comparator is a tool, not a wet-lab label.
TOOL_COMPARATOR = (
    "vs kaptive", "vs amrfinder", "vs seqsero2", "vs resfinder", "vs pharmcat", "reference tool",
    "faithful-to-tool", "faithful to tool", "tool comparator", "consensus", "in-distribution",
    "in distribution", "knowledge_baseline", "bounds agreement", "agreement-with-reference",
)

# Phrases that explicitly DENY a measurement exists. Their presence makes a never-measured tier coherent
# even when a number appears nearby (e.g. a null baseline, or a figure quoted from the literature).
DENIES_MEASUREMENT = (
    "no free", "validation_data_wall", "external wall", "never been scored", "no measured number",
    "not been measured", "no independent", "no measurement", "no cohort", "no phenotype source",
    "not censused", "nothing to score", "unmeasured", "cannot be tiered up", "no free validation",
    # An IN-DISTRIBUTION or INFERRED label is not an independent measurement, so a never-measured tier stays
    # coherent beside a figure. Adjudicated 2026-09-24: arabidopsis:flowering says outright that its
    # catalogue and its label trace to the same paper, and klebsiella:kleb's KL-types are prophage-host-LCA
    # INFERRED rather than observed. Both correctly hold KNOWLEDGE_BASELINE.
    "in-distribution", "in distribution", "lca-inferred", "lca_inferred", "inferred kl-type",
    "prophage-host-lca",
)


def _artifact_figures(blob: str) -> list[str]:
    """Measured figures reachable through an artifact the contract CITES.

    Load-bearing narrowing. The first over-claim rule was "no figure in the contract prose", which flagged
    39 of 115 cells -- because ~40 AMR cells are PROJECTED from shipped_decoder_surface and their evidence
    legitimately lives in the report card, not in the prose. A measured tier is defensible when the evidence
    is reachable EITHER in the prose OR through a cited artifact; only "reachable in neither" is a suspect.
    """
    from scripts.contract_number_audit import _cited_artifacts, _load_artifact_text
    out = []
    for rel in _cited_artifacts(blob):
        txt = _load_artifact_text(rel)
        if txt:
            out += _measured_figures(txt)[:3]
    return list(dict.fromkeys(out))[:6]


_CARD_INDEX: dict[tuple[str, str], dict] | None = None
# Matched as SUBSTRINGS of a card's field names, not as exact keys: card schemas differ per track and an
# exact allowlist silently missed all 25 HIV rows, whose metrics are named auc_call_separates_fold and
# catalog_balacc. Substring matching is what makes this work across cards written by different builders.
_METRIC_KEY_PARTS = ("acc", "sens", "spec", "auc", "auroc", "spearman", "pearson", "concordance",
                     "agreement", "balacc", "r2", "ppv", "npv", "f1", "mcc")


def _card_index() -> dict[tuple[str, str], dict]:
    """(organism, drug/target) -> the standing-report-card row, across every wiki/*report_card*.json.

    THE THIRD reachability path, and without it this audit is wrong about 39 cells. AMR and viral cells are
    PROJECTED from shipped_decoder_surface / _viral_contracts, so their contracts are deliberately short
    generic prose ("NCBI-PD provenance-disjoint stress test") with no figures and no citation -- their
    evidence lives in a STANDING card keyed by (organism, drug) that no individual cell names. Reading only
    prose + cited artifacts declared 39 of 115 cells over-claiming, which is a statement about the reader.
    """
    global _CARD_INDEX
    if _CARD_INDEX is None:
        idx: dict[tuple[str, str], dict] = {}
        for f in sorted((REPO / "wiki").glob("*report_card*.json")):
            try:
                doc = json.loads(f.read_text(encoding="utf-8", errors="replace"))
            except (OSError, json.JSONDecodeError):
                continue
            rows = doc.get("cells", doc) if isinstance(doc, dict) else doc
            if not isinstance(rows, list):
                continue
            for r in rows:
                if not isinstance(r, dict):
                    continue
                org = str(r.get("organism") or r.get("gene") or "").strip().lower()
                tgt = str(r.get("drug") or r.get("target") or r.get("trait") or "").strip().lower()
                if org or tgt:
                    idx[(org, tgt)] = r
        _CARD_INDEX = idx
    return _CARD_INDEX


def _card_figures(c) -> list[str]:
    """Metrics for this cell from a standing card, matched on (organism, target) case-insensitively.

    Also records, when a card row IS found but reports no metric, that the card DECLINES to score the cell
    (state UNDERPOWERED / ABSTAINS_BY_DESIGN) -- which is a different thing from no card row existing.
    """
    org = str(getattr(c, "organism", "") or "").strip().lower()
    tgt = str(getattr(c, "target", "") or "").strip().lower()
    idx = _card_index()
    # Gather EVERY plausible card row, then return the first that actually yields a metric. Taking the
    # first row that merely MATCHES was a bug: a loose ('human', '') fallback matched an unrelated row and
    # short-circuited before the gene lookup, so cyp2c19 read as over-claiming while its card row holds
    # getrm 72/72.
    candidates = [idx.get((org, tgt)), idx.get(("", tgt)), idx.get((tgt, "")), idx.get((org, ""))]
    candidates += [r for (o, _t), r in idx.items() if o == tgt]
    for row in candidates:
        if not row:
            continue
        out = []
        for k, v in row.items():
            if not any(part in k.lower() for part in _METRIC_KEY_PARTS):
                continue
            # Cards store a metric as a number OR as a string ("72/72", "n=65: sens 1.00"). A numeric-only
            # filter skipped every string-valued metric, which is how the PGx getrm field was missed.
            if isinstance(v, (int, float)):
                out.append(f"{k}={v}")
            elif isinstance(v, str) and re.search(r"\d", v):
                out.append(f"{k}={v[:40]}")
        if out:
            return out
    return []


# Suspects ADJUDICATED once as false positives, with the reason. An explicit table, like the one in
# contract_number_audit, because the alternative -- widening the metric vocabulary until nothing flags --
# is how a detector becomes vacuous.
ADJUDICATED_FALSE_POSITIVE: dict[str, str] = {
    "pgx:human:dpyd":
        "NOT an under-claim: DPYD is absent from the GeT-RM consolidated table entirely, so no free "
        "consensus truth exists for it (established 2026-09-22). Its contract reports deployment on 5 "
        "PGP-UK humans who were all *1/*1 -- a run with no true positives is not a measurement. The "
        "flagged figure was the fragment 'concordance = v0.1', a version string, not a metric.",
    "pgx:human:cyp2c19":
        "NOT an over-claim: wiki/pgx_report_card.json records getrm '72/72' for this gene. The detector "
        "misses it because that card names its metric fields by DOMAIN (getrm / pharmcat / trio_mendelian) "
        "rather than by metric (acc / sens / auc), and matching every card's private vocabulary is "
        "unbounded. This names the limitation instead of widening the vocabulary until nothing flags.",
}


def _tier_name(c) -> str:
    return str(getattr(c, "evidence_tier", "")).split(".")[-1]


def _measured_figures(blob: str) -> list[str]:
    out = [m.group(0).strip() for m in METRIC_RE.finditer(blob)]
    out += [m.group(0).strip() for m in CORR_RE.finditer(blob)]
    out += [m.group(0).strip() for m in COUNT_RE.finditer(blob)]
    return list(dict.fromkeys(out))


def audit() -> dict:
    sys.path.insert(0, str(REPO))
    from dna_decode.data.cell_registry import cells

    rows, suspects = [], []
    for c in cells():
        blob = " ".join(str(getattr(c, f, "") or "") for f in PROSE_FIELDS)
        low = blob.lower()
        tier = _tier_name(c)
        figs = _measured_figures(blob)
        denies = [p for p in DENIES_MEASUREMENT if p in low]

        status, why = "CONSISTENT", ""
        tool_cmp = [p for p in TOOL_COMPARATOR if p in low]
        # FAITHFUL_TO_TOOL + a tool comparator is COHERENT, not a contradiction: that is what the tier says.
        if tier == "FAITHFUL_TO_TOOL" and tool_cmp:
            denies = denies + [f"tool-comparator: {tool_cmp[0]}"]
        if tier in NEVER_MEASURED_TIERS and figs and not denies:
            status = "UNDER_CLAIM_SUSPECT"
            why = (f"tier {tier} asserts the cell was not measured, but the contract reports "
                   f"{len(figs)} measured figure(s) and carries no phrase denying a measurement")
        elif tier in MEASURED_TIERS and not figs:
            art_figs = _artifact_figures(blob) or _card_figures(c)
            if not art_figs:
                status = "OVER_CLAIM_SUSPECT"
                why = (f"tier {tier} asserts a measurement, but no measured figure is reachable either in "
                       f"the contract prose or through any artifact it cites")
            else:
                why = f"measured figure reachable via cited artifact or standing card: {art_figs[0]}"

        if status != "CONSISTENT" and c.cell_id in ADJUDICATED_FALSE_POSITIVE:
            why = "ADJUDICATED FALSE POSITIVE -- " + ADJUDICATED_FALSE_POSITIVE[c.cell_id]
            status = "CONSISTENT"
        row = {"cell_id": c.cell_id, "tier": tier, "status": status,
               "measured_figures": figs[:6], "denial_phrases": denies[:4],
               "tool_comparator_phrases": tool_cmp[:3], "why": why}
        rows.append(row)
        if status != "CONSISTENT":
            suspects.append(row)

    return {
        "schema": "tier-evidence-audit-v1",
        "analysis_date": datetime.date.today().isoformat(),
        "posture": "TRIAGE FUNNEL, not an oracle. Tier assignment is a JUDGEMENT about evidence class "
                   "(wet-lab label vs reference tool vs the rule's own source) and lexical reading cannot "
                   "settle it. This detects a narrow CONTRADICTION and emits it for adjudication; it never "
                   "edits a tier.",
        "not_flagged_on_purpose": "which measured tier a cell belongs in (INDEPENDENT_MEASURED vs "
                                  "NEAR_INDEPENDENT vs FAITHFUL_TO_TOOL) -- that turns on whether the "
                                  "comparator is a tool or a wet-lab label, the judgement this cannot make",
        "n_cells": len(rows),
        "n_under_claim_suspect": sum(1 for r in rows if r["status"] == "UNDER_CLAIM_SUSPECT"),
        "n_over_claim_suspect": sum(1 for r in rows if r["status"] == "OVER_CLAIM_SUSPECT"),
        "verdict": "NOTHING_TO_ADJUDICATE" if not suspects else "ADJUDICATION_REQUIRED",
        "suspects": suspects,
        "cells": rows,
    }


def main(argv=None) -> int:
    rep = audit()
    if rep["n_cells"] == 0:
        print("REFUSED: no cells audited -- a clean report over zero cells is not a result.",
              file=sys.stderr)
        return 3
    stamp = rep["analysis_date"]
    (REPO / "wiki" / f"tier_evidence_audit_{stamp}.json").write_text(json.dumps(rep, indent=2),
                                                                   encoding="utf-8")
    for r in rep["suspects"]:
        print(f"{r['status']:20} {r['cell_id']:42} tier={r['tier']}")
        print(f"    {r['why']}")
        if r["measured_figures"]:
            print(f"    figures: {'; '.join(r['measured_figures'][:4])}")
    print(f"\n{rep['verdict']}  under-claim {rep['n_under_claim_suspect']} / over-claim "
          f"{rep['n_over_claim_suspect']} across {rep['n_cells']} cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
