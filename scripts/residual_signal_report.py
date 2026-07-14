"""CLI: residual-signal report from a cross-axis lineage-deconfound artifact (2026-07-13).

The g8-residual-detector product surface. Reads a committed `crossaxis_lineage_deconfound_*.json` and
emits a per-feature residual-signal report (md + json): which genome features carry mechanism signal
that SURVIVES leave-one-clade-out de-confounding (GENERALIZES) vs which are clade-mediated / possibly
clonal (LINEAGE_MEDIATED). Read-only over committed artifacts; frozen decoder surface untouched.

  uv run python scripts/residual_signal_report.py \
    wiki/crossaxis_lineage_deconfound_determinant_2026-07-12.json \
    --out wiki/residual_signal_ecoli_determinant.md
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.eval.residual_detector import GENERALIZES, LINEAGE_MEDIATED, UNTESTED, load_report  # noqa: E402


def _render_md(rep: dict) -> str:
    m = rep["meta"]
    tc = rep["tier_counts"]
    lines = [
        f"# Residual-signal report — {m.get('organism')} / {m.get('axis_label')}",
        "",
        f"generated {_date.today().isoformat()} · dna_decode.eval.residual_detector (g8) · "
        "signal-provenance, NOT a phenotype prediction",
        "",
        "**What this is.** For each feature, whether its genotype→axis association SURVIVES "
        "leave-one-clade-out (Mash-clade GroupKFold) de-confounding. "
        f"`{GENERALIZES}` = real signal beyond lineage; `{LINEAGE_MEDIATED}` = clade-concentrated "
        f"(may be clonal structure); `{UNTESTED}` = no cross-axis entry.",
        "",
        f"- source verdict: `{m.get('source_verdict')}`",
        f"- median AUC naive {m.get('median_auc_naive')} → clade-grouped {m.get('median_auc_clade_grouped')}",
        f"- Mash partition: threshold {m.get('mash_threshold')}, {m.get('n_clades')} clades, "
        f"largest-clade fraction {m.get('largest_clade_frac')}",
        f"- **tiers:** {GENERALIZES}={tc.get(GENERALIZES,0)} · {LINEAGE_MEDIATED}={tc.get(LINEAGE_MEDIATED,0)} · "
        f"{UNTESTED}={tc.get(UNTESTED,0)} ({m.get('n_features')} features)",
        "",
        "## Per-feature (residual signal — GENERALIZES first, strongest de-confounded AUC on top)",
        "",
        "| feature | tier | family | n | AUC naive | AUC clade-grouped | drop |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for r in rep["per_feature"]:
        lines.append(f"| `{r['feature_id']}` | {r['tier']} | {r['family']} | {r['n_present']} | "
                     f"{r['auc_naive']} | {r['auc_clade_grouped']} | {r['drop']} |")
    lines += ["", "## Gene-family rollup (residual signal by family)", "",
              "| family | generalizes | lineage-mediated | untested | total |",
              "|---|---:|---:|---:|---:|"]
    for fam, c in sorted(rep["family_rollup"].items(), key=lambda kv: -kv[1][GENERALIZES]):
        lines.append(f"| {fam} | {c[GENERALIZES]} | {c[LINEAGE_MEDIATED]} | {c[UNTESTED]} | {c['total']} |")
    lines += ["", "## Honest caveats", ""] + [f"- {c}" for c in rep["honest_caveats"]]
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("artifact", help="crossaxis_lineage_deconfound_*.json")
    ap.add_argument("--out", type=Path, required=True, help="output .md (a .json sidecar is written too)")
    a = ap.parse_args(argv)

    rep = load_report(a.artifact)
    a.out.write_text(_render_md(rep), encoding="utf-8")
    json_out = a.out.with_suffix(".json")
    json_out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
    tc = rep["tier_counts"]
    print(f"[residual] {rep['meta']['organism']}/{rep['meta']['axis_label']}: "
          f"{tc[GENERALIZES]} generalizes / {tc[LINEAGE_MEDIATED]} lineage-mediated / {tc[UNTESTED]} untested "
          f"({rep['meta']['n_features']} features) -> {a.out} (+ {json_out.name})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
