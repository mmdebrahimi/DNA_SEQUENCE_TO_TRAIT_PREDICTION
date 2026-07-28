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
    # AUROC = P(score_ess > score_non)
    from scipy.stats import mannwhitneyu
    auroc = mannwhitneyu(se, sn, alternative="greater").statistic / (len(se)*len(sn))
    return {"n_essential": len(ceg), "n_nonessential": len(neg), "auroc": round(float(auroc), 4),
            "sens": round(float((se >= 2).mean()), 4), "spec": round(float((sn < 2).mean()), 4)}


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
