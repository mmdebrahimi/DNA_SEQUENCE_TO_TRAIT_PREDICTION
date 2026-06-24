"""Top-level CROSS-KINGDOM validation summary — one legible view of the whole validated decoder surface.

The project now spans bacteria + M. tuberculosis + viral cells, each with its OWN standing report card and
its OWN (different) independence construction. This is a READ-ONLY roll-up of those cards into one view, so a
reader can see the validated surface at a glance — WITHOUT flattening the distinct honesty tiers into a
misleading aggregate (the project's "no aggregate headline" discipline; each card's independence claim is
preserved verbatim, NOT averaged).

Reads (each namespace-separate, by design — the shared-key-overwrite lesson):
  - `wiki/decoder_validation_report_card.json`   — bacterial NCBI-PD provenance-disjoint cells (frozen surface)
  - `wiki/amr_portal_independent_report_card.json`— bacterial EBI AMR Portal independent measured-AST cells
  - `wiki/tb_report_card.json`                    — M. tuberculosis (independent AMR Portal + in-distribution CRyPTIC)
  - `wiki/hiv_decoder_report_card.json`           — HIV (Stanford HIVDB PhenoSense, free independent label)
  - `wiki/external_validation_report_card.json`   — external-cohort revalidation arm

Exit 0 always — a REPORT, not a gate. -> `wiki/cross_kingdom_validation_summary.{md,json}`.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WIKI = REPO / "wiki"


def _load(name: str) -> dict | None:
    p = WIKI / name
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


def collect() -> list[dict]:
    """One summary row per validation surface, preserving each card's distinct independence tier."""
    rows = []

    # 1. Bacteria — NCBI-PD provenance-disjoint (the frozen deployed surface)
    d = _load("decoder_validation_report_card.json")
    if d:
        sc = d.get("state_counts", {})
        rows.append({"kingdom": "Bacteria", "surface": "NCBI-PD provenance-disjoint (frozen)",
                     "independence": "isolate-level provenance-disjoint (different submitter/lab/country); "
                     "NOT methodology-independent",
                     "summary": f"{sc.get('SCORED', 0)} SCORED / {len(d.get('cells', []))} cells",
                     "card": "wiki/decoder_validation_report_card.md"})

    # 2. Bacteria — EBI AMR Portal independent measured-AST
    d = _load("amr_portal_independent_report_card.json")
    if d:
        rows.append({"kingdom": "Bacteria", "surface": "EBI AMR Portal independent (measured AST)",
                     "independence": "BioSample/GCA disjoint vs CRyPTIC + tuning cohorts; measured wet-lab AST "
                     "(non-circular); E.coli/Salmonella/Klebsiella/Shigella 0.83–0.995 acc",
                     "summary": f"{d.get('n_scored_independent', 0)} SCORED_INDEPENDENT / {d.get('n_cells', 0)} cells",
                     "card": "wiki/amr_portal_independent_report_card.md"})

    # 3. M. tuberculosis — independent (AMR Portal) + in-distribution (CRyPTIC)
    d = _load("tb_report_card.json")
    if d:
        ind = d.get("independent", [])
        nums = "; ".join(f"{r['drug'][:3].upper()} acc {r.get('raw_acc'):.3f}" for r in ind if r.get("raw_acc"))
        rows.append({"kingdom": "Bacteria (M. tuberculosis)", "surface": "EBI AMR Portal INDEPENDENT (WHO rule)",
                     "independence": "BioSample-resolution-checked (ENA-side biosample-grade; NCBI-side 0/30 "
                     "cross-archive); RAW headline (homoplasy -> lineage is disclosure); measured AST",
                     "summary": f"{nums or 'n/a'} (N~2,845, full cohort)",
                     "card": "wiki/tb_report_card.md"})

    # 4. HIV — Stanford HIVDB PhenoSense (free independent label)
    d = _load("hiv_decoder_report_card.json")
    if d:
        rows.append({"kingdom": "Virus (HIV)", "surface": "Stanford HIVDB PhenoSense",
                     "independence": d.get("label_independence", "free independent isolate-level wet-lab "
                     "fold-change (non-circular vs HIVDB's own interpretation)"),
                     "summary": f"{d.get('n_cells', '?')} cells (NNRTI/NRTI/PI/INSTI/CAI)",
                     "card": "wiki/hiv_decoder_report_card.md"})

    # 5. External-cohort revalidation arm
    d = _load("external_validation_report_card.json")
    if d:
        rows.append({"kingdom": "Bacteria (external cohort)", "surface": "external-cohort revalidation arm",
                     "independence": d.get("note", "frozen v0.5.0 decoder re-validated on an INDEPENDENT "
                     "measured-MIC cohort; BioSample-level leakage preflight"),
                     "summary": f"{d.get('n_cells', 0)} cells",
                     "card": "wiki/external_validation_report_card.md"})
    return rows


def build() -> dict:
    rows = collect()
    return {"schema": "cross-kingdom-validation-summary-v1",
            "no_aggregate_headline": True,
            "principle": "Each surface has a DIFFERENT independence construction + honesty tier; they are "
            "PRESERVED, never averaged into a single accuracy. The binding constraint of the project (a free "
            "independent measured-phenotype label) is broken across bacteria AND TB; SARS-CoV-2/influenza/"
            "fungal cells are no-free-source or in-distribution; learned-embedding expansion is a closed "
            "0-for-4 negative (see wiki/negative_results_map_2026-06-13.md).",
            "n_surfaces": len(rows), "surfaces": rows}


def render_md(card: dict) -> str:
    L = ["# Cross-kingdom validation summary — the whole validated decoder surface",
         "",
         "One legible view of every validation surface. **No aggregate headline** — each surface has a "
         "DIFFERENT independence construction + honesty tier, preserved verbatim (never averaged). See each "
         "card for per-cell detail.",
         "",
         "| Kingdom | Surface | Headline | Independence tier (honest) | Card |",
         "|---|---|---|---|---|"]
    for r in card["surfaces"]:
        L.append(f"| {r['kingdom']} | {r['surface']} | {r['summary']} | {r['independence']} | "
                 f"`{r['card']}` |")
    L += ["",
          "## The arc (what this surface represents)",
          "- **The binding constraint of the project — a FREE, independent, measured-phenotype label — is "
          "broken across bacteria AND M. tuberculosis.** The deterministic decoder is independently validated "
          "on E. coli / Salmonella / Klebsiella / Shigella (EBI AMR Portal, measured AST) **and** "
          "M. tuberculosis (WHO rule, N~2,845), all FREE (no DUA, no author-contact).",
          "- **HIV** is validated against a free independent wet-lab fold-change (Stanford HIVDB PhenoSense) — "
          "the project's first independent-label win.",
          "- **SARS-CoV-2 / influenza / fungal** cells are in-distribution or no-free-phenotype-source (their "
          "honesty tiers say so on their own cards).",
          "- **Learned-embedding expansion is a CLOSED 0-for-4 negative** (cipro within-lineage, pathotype, "
          "Arabidopsis ×2) — `wiki/negative_results_map_2026-06-13.md`. The validated shippable artifact is "
          "the DETERMINISTIC decoder suite, not the foundation-model embedding bet.",
          "",
          "## Honesty discipline (preserved, not flattened)",
          "Each surface's independence is a DIFFERENT construction: NCBI-PD provenance-disjoint (submitter/lab/"
          "country) ≠ EBI AMR Portal accession-disjoint measured-AST ≠ TB BioSample-resolution-checked ≠ HIV "
          "free wet-lab fold-change ≠ external-cohort measured-MIC. This summary lists them side by side; it "
          "does NOT compute a single cross-kingdom accuracy (that would be a category error).",
          "",
          "Rebuild: `uv run python scripts/build_cross_kingdom_summary.py`."]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    card = build()
    (WIKI / "cross_kingdom_validation_summary.json").write_text(
        json.dumps(card, indent=2, default=str), encoding="utf-8")
    (WIKI / "cross_kingdom_validation_summary.md").write_text(render_md(card), encoding="utf-8")
    print(f"[cross-kingdom] {card['n_surfaces']} validation surfaces -> "
          "cross_kingdom_validation_summary.{md,json}")
    for r in card["surfaces"]:
        print(f"  {r['kingdom']:<28} {r['surface']:<42} {r['summary']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
