"""Essentiality decoder report card + the E4 cross-organism transfer eval (reproducible).

Read-only roll-up (exit 0 always -- a report, NOT a gate), mirroring the AMR/forward report cards:
per-organism honest tier, NO aggregate headline. The E. coli cell is validated by size+composition
(labels walled); the human cell is the cross-organism TRANSFER AUROC on the citable BAGEL CEG2/NEG
reference. Needs the D: data (gene_info + BAGEL sets); degrades to NOT_EVALUATED if absent.
"""
from __future__ import annotations
import gzip, json
from pathlib import Path
import numpy as np
from dna_decode.essentiality.core_decoder import score_gene

TGT = Path("D:/dna_decode_cache/essentiality")
W = Path("wiki")


def human_transfer_auroc():
    gi = TGT / "Homo_sapiens.gene_info.gz"
    if not (gi.exists() and (TGT/"CEGv2.txt").exists() and (TGT/"NEGv1.txt").exists()):
        return None
    desc = {}
    with gzip.open(gi, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); gid=h.index("GeneID"); sym=h.index("Symbol")
        dsc=h.index("description"); ty=h.index("type_of_gene")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[ty] == "protein-coding": desc[p[gid]] = (p[sym], p[dsc])
    def load(fn):
        return [desc[c.split("\t")[2]] for c in (TGT/fn).read_text().splitlines()[1:]
                if len(c.split("\t")) >= 3 and c.split("\t")[2] in desc]
    ceg, neg = load("CEGv2.txt"), load("NEGv1.txt")
    se = np.array([score_gene(g, d).core_score for g, d in ceg])
    sn = np.array([score_gene(g, d).core_score for g, d in neg])
    # SHARED AUROC definition, imported from the ladder runner, so the card and the ladder cannot report
    # two different AUROCs for the same rung. Previously computed inline here with mannwhitneyu.
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    from essentiality_transfer_ladder import auroc as shared_auroc
    auroc = shared_auroc(se, sn)
    return {"n_essential": len(ceg), "n_nonessential": len(neg), "auroc": round(float(auroc), 4),
            "sens": round(float((se >= 2).mean()), 4), "spec": round(float((sn < 2).mean()), 4)}


LADDER_GLOB = "essentiality_transfer_ladder_*.json"


def load_ladder():
    """Newest transfer-ladder artifact, or None. Read-only; the card never triggers a ladder run."""
    hits = sorted(W.glob(LADDER_GLOB))
    if not hits:
        return None
    return json.loads(hits[-1].read_text(encoding="utf-8"))


def ladder_rows(art):
    """Ladder rungs as card rows. AUGMENT-ONLY: these are APPENDED; no existing row is touched.

    Each rung carries its OWN honest tier -- COVERAGE_SCORED for a rung that scored, WALL_<reason> for
    one that could not -- so a walled rung renders as a wall rather than vanishing from the table. There
    is deliberately no aggregate over them: the card's standing invariant is per-organism tiers only.
    """
    if not art:
        return []
    rows = []
    for r in sorted(art["rungs"], key=lambda x: -x["depth"]):
        if r.get("scored"):
            rows.append({
                "organism": r["organism"],
                "cell": "transfer ladder rung (depth %d from E. coli)" % r["depth"],
                "tier": "COVERAGE_SCORED",
                "metric": ("coverage_lift %.4f (null p95 %.4f / MAX %.4f); coverage(ess) %.4f, "
                           "coverage(non) %.4f; phrasing-adjusted %.4f; AUROC %.4f is SECONDARY and "
                           "NOT cross-rung comparable"
                           % (r["coverage_lift"], r["null_p95"], r["null_max"],
                              r["coverage_essential"], r["coverage_nonessential"],
                              r["coverage_lift_adjusted"],
                              r["auroc_SECONDARY_not_cross_rung_comparable"])),
                "validation": ("decoder applied UNCHANGED; coverage_lift is class-conditioned so it is "
                               "base-rate robust and cross-rung comparable, which AUROC is not. Label "
                               "technology %s; class sourcing %s"
                               % (r["technology"], r["class_sourcing_mode"])),
            })
        else:
            rows.append({
                "organism": r["organism"],
                "cell": "transfer ladder rung (depth %d from E. coli)" % r["depth"],
                "tier": r["wall"],
                "metric": "none -- a walled rung carries no coverage number by construction",
                "validation": "route tried: %s" % r.get("route_tried", "(not recorded)"),
            })
    return rows


def main():
    ht = human_transfer_auroc()
    ec = None
    ecp = W / "essentiality_ecoli_v0_1_auroc_2026-07-28.json"
    if ecp.exists():
        ec = json.loads(ecp.read_text())
    ecoli_row = (
        {"organism": "Escherichia coli K-12", "cell": "conserved-core v0.1", "tier": "AUROC_SCORED",
         "metric": f"AUROC {ec['auroc']} vs null 0.5 (Goodall-TraDIS gold-standard, n={ec['n']}, "
                   f"{ec['n_essential']} ess/{ec['n_nonessential']} non, base rate {ec['base_rate']}); "
                   f"sens {ec['sens']} spec {ec['spec']} prec {ec['precision']}",
         "validation": "real per-gene AUROC vs the Goodall 2018 mBio Table S1 gold-standard (CC-BY); "
                       "high-precision moderate-recall -- catches the universal core, misses the E. coli-"
                       "specific essential tail (the E3 learned-complement target)"}
        if ec else
        {"organism": "Escherichia coli K-12", "cell": "conserved-core v0", "tier": "COMPOSITION_VALIDATED",
         "metric": "208/4318 predicted essential (known essentialome ~300)",
         "validation": "size + composition match the known essentialome; per-gene AUROC pending labels (walled)"})
    rows = [
        ecoli_row,
        {"organism": "Homo sapiens", "cell": "cross-organism transfer (E4)",
         "tier": "TRANSFER_SCORED" if ht else "NOT_EVALUATED",
         "metric": (f"AUROC {ht['auroc']} vs null 0.50 (BAGEL CEG2 n={ht['n_essential']} / NEG n={ht['n_nonessential']}); "
                    f"sens {ht['sens']} spec {ht['spec']}") if ht else "D: data absent",
         "validation": "universal core (ribosome/tRNA-synth/translation/polymerase) transfers cross-kingdom at "
                       "high precision; human-specific core (proteasome 0/53, spliceosome 0/49) MISSED -> "
                       "per-organism catalogue extension is the follow-on" if ht else "-"},
    ]
    ladder = load_ladder()
    lrows = ladder_rows(ladder)
    rows = rows + lrows          # APPEND only -- the two pre-existing rows above are untouched

    card = {"schema": "essentiality-report-card-v1", "generated": "2026-07-28",
            "note": "single-gene KO -> essential/non-essential; conserved-core R1 decoder; per-organism honest "
                    "tier, NO aggregate headline; E. coli composition-validated, human transfer-AUROC (BAGEL)",
            "organisms": rows}
    (W/"essentiality_report_card.json").write_text(json.dumps(card, indent=2), encoding="utf-8")
    md = ["# Essentiality decoder report card (standing trust surface)", "",
          "Single-gene KO -> essential/non-essential, via the conserved-core R1 decoder. Per-organism honest",
          "tier; **no aggregate headline**. E. coli validated by composition (labels walled); human = the",
          "cross-organism TRANSFER AUROC on the citable BAGEL CEG2/NEG reference.", "",
          "| organism | cell | tier | metric | validation |", "|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['organism']} | {r['cell']} | {r['tier']} | {r['metric']} | {r['validation']} |")
    if lrows:
        v = ladder["verdict"]
        md += ["", "## Cross-organism transfer ladder", "",
               "Primary metric is **`coverage_lift`** = coverage(essential) - coverage(non-essential),",
               "which is class-conditioned and therefore base-rate robust and cross-rung comparable.",
               "AUROC is a SECONDARY and is **not** compared across rungs (different sampling frames:",
               "BAGEL is two curated extremes at base rate 0.431 vs a genome-wide 0.0928).", "",
               "Verdict: **`%s`** - %s" % (v["verdict"], v["reason"]),
               "Scored %d / walled %d. Frozen bar %s - a DESCRIPTIVE-CONSISTENCY LOCK, not a"
               % (ladder["n_scored"], ladder["n_walled"], json.dumps(ladder["frozen_thresholds"])),
               "frozen-before-the-numbers endpoint test (it is derived from the two rungs that",
               "motivated the question). Full detail + the five honest limits:",
               "`wiki/essentiality_transfer_ladder_%s.md`." % ladder["date"], "",
               "| organism | rung | tier | coverage_lift |", "|---|---|---|---|"]
        for r in ladder["rungs"]:
            md.append("| %s | depth %d | %s | %s |"
                      % (r["organism"], r["depth"],
                         "COVERAGE_SCORED" if r.get("scored") else r["wall"],
                         ("%.4f" % r["coverage_lift"]) if r.get("scored") else "-"))
    md += ["", "## Honest scope",
           "- The conserved-core decoder is the R1 PRIOR: high-precision, conservative-recall; captures the",
           "  UNIVERSAL essential core, misses lineage-specific core (the R2/per-organism-catalogue target).",
           "- E. coli per-gene AUROC + a learned E3 complement are gated on gold-standard labels (see",
           "  `wiki/essentiality_label_wall_2026-07-28.md`); human labels (BAGEL CEG2/NEG) ARE available.",
           "- Regenerate: `scripts/build_essentiality_report_card.py` (needs D: gene_info + BAGEL sets)."]
    (W/"essentiality_report_card.md").write_text("\n".join(md), encoding="utf-8")
    print("wrote essentiality_report_card.{md,json}")
    for r in rows: print(f"  [{r['tier']}] {r['organism']}: {r['metric']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
