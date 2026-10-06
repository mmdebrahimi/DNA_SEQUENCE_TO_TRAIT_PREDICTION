"""`dna-decode essentiality` / `dna-essentiality` — single-gene KO -> essential / non-essential.

    dna-decode essentiality --gene ftsZ --product "cell division protein FtsZ"
    dna-decode essentiality --feature-table GCF_..._feature_table.txt.gz --json

The DETERMINISTIC conserved-core decoder: predicts essentiality from gene FUNCTION (product text) via a
curated universal-core function catalogue (translation / replication / transcription / cell-envelope /
division). Label-INDEPENDENT (reads function, not a label). Offline, no deps.

Tier: KNOWLEDGE_BASELINE. Validated vs gold-standard: E. coli AUROC 0.695 (Goodall 2018 TraDIS, genome-wide);
composition matches the known essentialome (208/4318, translation/envelope/replication-dominated). The
conserved-core is HIGH-PRECISION, conservative-recall (captures the universal core, misses the lineage-specific
tail). The learned E3 complement is a SEPARATE WITHIN-organism 5-fold CV — E. coli 0.795 and human 0.911, each
against its OWN within-organism baseline — and does NOT lift the 0.580 CROSS-ORGANISM transfer number; reading
it as transfer evidence is a documented trap. NOT a clinical tool. See
wiki/essentiality_{decoder_v0,ecoli_v0_1_auroc,e4_transfer,e3_learned,e3_human}_2026-07-28.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path
from dna_decode.data.cell_evidence_line import evidence_one_line

_SCOPE = ("scope: conserved-core deterministic prior; high-precision, conservative-recall; validated vs "
          "gold-standard (E. coli AUROC 0.695 genome-wide). Captures the universal essential core, misses the "
          "lineage-specific tail. NOT clinical.")


def _iter_feature_table(path: str):
    """Thin alias — the single definition now lives in `annotation_join.iter_feature_table`.

    Promoted so the transfer-ladder arm and this CLI cannot drift to two readers of the same file
    format. Kept as a module-level name because it is referenced below and in tests.
    """
    from dna_decode.essentiality.annotation_join import iter_feature_table

    yield from iter_feature_table(path)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode essentiality",
        description="Single-gene KO -> essential/non-essential (deterministic conserved-core decoder).",
        epilog=_SCOPE)
    ap.add_argument("--gene", help="gene symbol (e.g. ftsZ)")
    ap.add_argument("--product", help="gene product / function text (e.g. 'cell division protein FtsZ')")
    ap.add_argument("--feature-table", help="NCBI feature_table.txt(.gz): score every CDS gene genome-wide")
    ap.add_argument("--threshold", type=float, default=2.0, help="core-score essential threshold (default 2.0)")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.essentiality.core_decoder import score_gene, decode_genome

    if args.feature_table:
        if not Path(args.feature_table).exists():
            print(f"error: cannot read --feature-table: {args.feature_table}", file=sys.stderr); return 2
        genes = list(_iter_feature_table(args.feature_table))
        if not genes:
            print("error: no CDS rows with a 'symbol'/'name' column found in the feature table", file=sys.stderr)
            return 2
        calls = decode_genome(genes, threshold=args.threshold)
        ess = [c for c in calls if c.prediction == "essential"]
        out = {"n_genes": len(calls), "n_predicted_essential": len(ess),
               "essential_genes": [c.gene for c in ess],
               "scope": _SCOPE, "tier": "KNOWLEDGE_BASELINE"}
        if args.json:
            print(json.dumps(out, indent=2))
        else:
            print(f"genes: {out['n_genes']} | predicted essential: {out['n_predicted_essential']}")
            print("essential:", ", ".join(out["essential_genes"][:40]) + (" ..." if len(ess) > 40 else ""))
            print(_SCOPE)
        return 0

    if not args.product and not args.gene:
        print("error: give --gene/--product (single) or --feature-table (genome-wide)", file=sys.stderr)
        return 2
    c = score_gene(args.gene or "", args.product or "", threshold=args.threshold)
    d = c.as_dict(); d["scope"] = _SCOPE; d["tier"] = "KNOWLEDGE_BASELINE"
    if args.json:
        print(json.dumps(d, indent=2))
    else:
        print(f"{c.gene or '(no symbol)'}: {c.prediction.upper()}  (core-score {c.core_score:.1f}"
              f"{'; matched ' + ', '.join(c.matched) if c.matched else ''})")
        print(_SCOPE)
        # The cell's MEASURED evidence, from its committed registry contract. Appended BESIDE the
        # scope/caveat line above, never replacing it: a caveat states the method's limits, this
        # states what was actually measured. Evidence carried only in the registry is not a disclosure.
        _ev = evidence_one_line("dna-essentiality")
        if _ev:
            print(f"  {_ev}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
