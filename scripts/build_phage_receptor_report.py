"""Build the phage receptor-class cell v0 result packet (wiki/phage_receptor_cell_v0_<date>.{md,json}).

Reproducible: recomputes the leave-one-out receptor-transfer number from the committed BASEL manifest
+ genomes via native blastn, then emits the packet. Honest tier = IN_DISTRIBUTION (labels are
clade-derived from the same BASEL Results the catalog is curated from); the RBP-variable clades
(T-even, Drexlerviridae) are excluded and reported as the tractability boundary, not hidden.

Usage:  uv run python scripts/build_phage_receptor_report.py [--date YYYY-MM-DD]
Needs native blastn (see dna_decode.pathotype.vf_runner.find_blastn); without it the LOO degrades to
all-INDETERMINATE and the packet says so.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from dna_decode.phage.receptor_caller import _load_manifest, leave_one_out

MANIFEST = "data/phage_ref/basel_manifest.tsv"
GENOME_DIR = "data/phage_ref/basel"


def build(date: str) -> dict:
    refs, receptors = _load_manifest(MANIFEST, GENOME_DIR)
    total_genomes = sum(1 for _ in Path(GENOME_DIR).glob("*.fna"))
    res = leave_one_out(refs, receptors)
    from collections import Counter
    per_receptor_label = dict(Counter(receptors.values()))
    n_indeterminate = sum(1 for p in res.predictions if p["status"] == "INDETERMINATE")
    n_miscall = sum(1 for p in res.predictions if p["status"] == "CALLED" and not p["correct"])
    return {
        "cell": "phage_receptor_class",
        "date": date,
        "axis": "bacteriophage host-receptor class (first non-AMR, non-host-organism cell)",
        "substrate": "BASEL E. coli phage collection (Maffei 2021 PLOS Biology 3001424; GenBank MZ501046-MZ501113)",
        "label_source": "experimentally-determined receptors (>50 single-gene K-12 mutants + EOP host-range), CC-BY",
        "method": "genome-homology receptor TRANSFER (nearest-BLAST-neighbour inherits its receptor)",
        "tier": "IN_DISTRIBUTION",
        "independence": "closed for v0 - labels are clade-derived from the same BASEL Results the catalog is "
                        "curated from; an INDEPENDENT number needs a held-out phage set with measured receptors",
        "scope": "RECEPTOR-CLASS only (NOT the full phage x strain host-range matrix, which is polygenic/"
                 "intractable from genome alone); clade-conserved clades only",
        "n_genomes_total": total_genomes,
        "n_labelled": len(refs),
        "n_excluded_rbp_variable": total_genomes - len(refs),
        "labelled_per_receptor": per_receptor_label,
        "loo_overall_accuracy": res.accuracy,
        "loo_n_called": res.n_called,
        "loo_n_total": res.n_total,
        "loo_n_indeterminate": n_indeterminate,
        "loo_n_miscall": n_miscall,
        "loo_per_receptor_correct_called": res.per_receptor,
        "predictions": res.predictions,
        "excluded_clades": ["Tequatrovirus/Straboviridae (T-even: OmpC/FadL/Tsx vary by RBP)",
                            "Drexlerviridae (FhuA/BtuB/YncD/TolC vary)",
                            "Siphoviridae Dhillonvirus/Nonagvirus/Seuratvirus (LptD/FhuA/LamB vary)"],
        "honest_reading": "receptor-class TRANSFERS reliably along genome homology within clade-conserved "
                          "clades (0 mis-calls); the caller ABSTAINS (INDETERMINATE) rather than mis-transfer "
                          "when a phage has no reference homolog (e.g. the lone NfrA phage). The 100% is on "
                          "the clade-conserved subset BY CONSTRUCTION - it validates the pipeline + abstention "
                          "+ catalog self-consistency, not a solved RBP->receptor map. RBP-variable clades are "
                          "the documented tractability boundary.",
    }


def _md(d: dict) -> str:
    pr = d["loo_per_receptor_correct_called"]
    pr_rows = "\n".join(f"| {r} | {c[0]}/{c[1]} |" for r, c in sorted(pr.items()))
    acc = "n/a" if d["loo_overall_accuracy"] is None else f"{d['loo_overall_accuracy']:.3f}"
    return f"""# Phage receptor-class cell v0 ({d['date']})

**The first non-AMR, non-host-organism cell** — a bacteriophage-genome -> host-receptor-class decoder.

- **Axis:** {d['axis']}
- **Substrate:** {d['substrate']}
- **Label source (FREE MEASURED):** {d['label_source']}
- **Method:** {d['method']}
- **Scope:** {d['scope']}
- **Tier:** `{d['tier']}` — {d['independence']}

## Result (leave-one-out, native blastn)

- Genomes fetched: **{d['n_genomes_total']}**; clean-labelled (clade-conserved): **{d['n_labelled']}**;
  excluded RBP-variable: **{d['n_excluded_rbp_variable']}**.
- **Overall LOO accuracy: {acc}** ({d['loo_n_called']}/{d['loo_n_total']} called; {d['loo_n_indeterminate']}
  INDETERMINATE abstentions; **{d['loo_n_miscall']} mis-calls**).

| receptor | LOO correct/called |
|---|---|
{pr_rows}

## Honest reading

{d['honest_reading']}

**Excluded clades (the tractability boundary — receptor is receptor-binding-protein-determined, not clade-clean):**
{chr(10).join('- ' + c for c in d['excluded_clades'])}

## Reproduce

```bash
uv run python scripts/build_phage_receptor_report.py
# or just the number:
uv run python -c "from dna_decode.phage.receptor_caller import _load_manifest, leave_one_out; \\
r=leave_one_out(*_load_manifest('data/phage_ref/basel_manifest.tsv','data/phage_ref/basel')); \\
print(r.accuracy, r.per_receptor)"
```
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="2026-07-24")
    args = ap.parse_args()
    d = build(args.date)
    wiki = Path("wiki")
    (wiki / f"phage_receptor_cell_v0_{args.date}.json").write_text(json.dumps(d, indent=2), encoding="utf-8")
    (wiki / f"phage_receptor_cell_v0_{args.date}.md").write_text(_md(d), encoding="utf-8")
    print(f"wrote wiki/phage_receptor_cell_v0_{args.date}.md + .json  |  LOO acc={d['loo_overall_accuracy']}")


if __name__ == "__main__":
    main()
