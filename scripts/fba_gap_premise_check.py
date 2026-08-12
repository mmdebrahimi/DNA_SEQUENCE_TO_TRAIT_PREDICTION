"""Does MODEL-GAP structure explain FBA's essentiality false negatives? (Track C's premise, tested)

Track C of the design epoch (`wiki/design_epoch_plan_2026-08-07.md`) proposes protein-function prediction
over dark-matter genes -> candidate missing reactions -> better FBA accuracy. That is a large build, and it
rests on one premise worth testing FIRST:

    the model's FALSE NEGATIVES -- genes essential in vivo that FBA calls dispensable -- are there because
    the reactions they control are BLOCKED or touch a DEAD-END metabolite (i.e. the model is incomplete),
    rather than because of regulation / kinetics / moonlighting / medium (which gap-filling cannot reach).

FALSIFIER, pre-registered before running: FN gap-adjacency must be materially ABOVE the TN rate among the
genes the model calls dispensable (Fisher two-sided p < 0.05). Equal rates falsify the premise as written.

Also separates the CHEAPER lever: reactions blocked under the model's default medium but NOT when every
exchange is opened are a MEDIUM mis-specification (a config fix), not a structural gap needing new
biochemistry. If most of the blocked set were medium-caused, Track C would be the wrong tool to reach for.

    uv run python scripts/fba_gap_premise_check.py --organism yeast

Exit 0 = premise SUPPORTED; 1 = FALSIFIED; 2 = not scorable (no gold standard / no model).
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from datetime import date
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.essentiality_labels import ESSENTIALITY_LABEL_SOURCES, parse_essential  # noqa: E402
from dna_decode.fba.gapfill import model_dead_ends  # noqa: E402
from dna_decode.fba.model import gene_essentiality, load_model, wildtype_growth  # noqa: E402


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact on a 2x2 via hypergeometric tails. No scipy dependency.

    Two-sided is the sum of every table at least as extreme (probability <= the observed table's), which
    is the standard definition and does NOT assume symmetry -- relevant here because the margins are very
    unequal (92 FN vs 703 TN).
    """
    n = a + b + c + d
    if n == 0 or (a + b) == 0 or (a + c) == 0:
        return float("nan")
    row1, col1 = a + b, a + c

    def p(k: int) -> float:
        return comb(row1, k) * comb(n - row1, col1 - k) / comb(n, col1)

    obs = p(a)
    lo, hi = max(0, col1 - (n - row1)), min(row1, col1)
    return float(sum(p(k) for k in range(lo, hi + 1) if p(k) <= obs * (1 + 1e-9)))


def structural_vs_medium_blocked(model) -> dict:
    """Split the blocked set into STRUCTURAL (blocked even with every exchange open) and MEDIUM-caused.

    `processes=1` is deliberate: cobra's FVA uses multiprocessing, which on Windows re-imports the calling
    module in each child and crashes any caller lacking an `if __name__ == "__main__"` guard.
    """
    from cobra.flux_analysis import find_blocked_reactions  # noqa: PLC0415

    default_blocked = set(find_blocked_reactions(model, processes=1))
    with model:
        for r in model.exchanges:
            r.lower_bound = -1000.0
        structural = set(find_blocked_reactions(model, processes=1))
    freed = default_blocked - structural
    return {
        "n_reactions": len(model.reactions),
        "n_blocked_default_medium": len(default_blocked),
        "n_blocked_structural": len(structural),
        "n_unblocked_by_opening_medium": len(freed),
        "structural_fraction_of_blocked": round(len(structural) / max(1, len(default_blocked)), 4),
        "blocked_ids_default": sorted(default_blocked),
        "blocked_ids_structural": sorted(structural),
    }


def gap_adjacency(model, blocked: set[str], dead_end_metabolites: set[str]):
    """A gene is gap-adjacent if ANY reaction it controls is blocked or touches a dead-end metabolite."""
    def is_adj(gene_id: str) -> bool:
        g = model.genes.get_by_id(gene_id)
        return any(r.id in blocked or any(m.id in dead_end_metabolites for m in r.metabolites)
                   for r in g.reactions)
    return is_adj


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", default="yeast")
    ap.add_argument("--label-file", default=None, help="local gold standard (skip the network fetch)")
    ap.add_argument("--frac", type=float, default=0.01, help="essentiality growth threshold (of wild type)")
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    if a.organism not in ESSENTIALITY_LABEL_SOURCES:
        print(f"no essentiality gold standard registered for '{a.organism}'", file=sys.stderr)
        return 2
    kind, url = ESSENTIALITY_LABEL_SOURCES[a.organism]
    text = (Path(a.label_file).read_text(encoding="utf-8", errors="replace") if a.label_file
            else urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "dna-decode"}),
                                        timeout=90).read().decode("utf-8", "replace"))
    ess = parse_essential(kind, text)

    model = load_model(organism=a.organism)
    wt = wildtype_growth(model)
    dead_ids = {d.metabolite for d in model_dead_ends(model)}
    blk = structural_vs_medium_blocked(model)
    blocked = set(blk["blocked_ids_default"])
    print(f"{model.id}: {len(model.genes)} genes, {len(model.reactions)} reactions, wt {wt:.4f}")
    print(f"dead-end metabolites {len(dead_ids)}/{len(model.metabolites)} | "
          f"blocked {blk['n_blocked_default_medium']} (default) / {blk['n_blocked_structural']} (structural)")

    fba = gene_essentiality(model, frac=a.frac)
    pred = {g: v[1] for g, v in fba.items()}
    is_adj = gap_adjacency(model, blocked, dead_ids)

    cells: dict[str, list[str]] = {"TP": [], "FP": [], "FN": [], "TN": []}
    for gid, p in pred.items():
        e = gid in ess
        cells["TP" if (e and p) else "FN" if e else "FP" if p else "TN"].append(gid)
    adj_counts = {k: sum(is_adj(g) for g in v) for k, v in cells.items()}
    rates = {k: (round(adj_counts[k] / len(v), 4) if v else None) for k, v in cells.items()}

    fn_adj, fn_not = adj_counts["FN"], len(cells["FN"]) - adj_counts["FN"]
    tn_adj, tn_not = adj_counts["TN"], len(cells["TN"]) - adj_counts["TN"]
    p_value = fisher_exact_two_sided(fn_adj, fn_not, tn_adj, tn_not)
    supported = bool(rates["FN"] and rates["TN"] and rates["FN"] > rates["TN"] and p_value < 0.05)

    for k in ("TP", "FP", "FN", "TN"):
        print(f"   {k}: {len(cells[k]):4d} genes | gap-adjacent {adj_counts[k]:4d} "
              f"({100 * (rates[k] or 0):5.1f}%)")
    print(f"\nFN vs TN gap-adjacency: {100 * rates['FN']:.1f}% vs {100 * rates['TN']:.1f}% "
          f"| Fisher two-sided p = {p_value:.4g}")
    print(f"PREMISE {'SUPPORTED' if supported else 'FALSIFIED'}")

    result = {
        "record": "fba-gap-premise-check-v1",
        "date": a.date,
        "organism": a.organism,
        "model": model.id,
        "label_source": f"{kind}: {url}",
        "n_model_genes_scored": len(pred),
        "confusion_counts": {k: len(v) for k, v in cells.items()},
        "gap_adjacent_counts": adj_counts,
        "gap_adjacent_rates": rates,
        "fn_vs_tn_fisher_p": None if p_value != p_value else round(p_value, 8),
        "premise_supported": supported,
        "blocked_reactions": {k: v for k, v in blk.items() if not k.endswith("_ids_default")
                              and not k.endswith("_ids_structural")},
        "verdict": ("Gap-adjacency predicts WHICH essential genes FBA misses; the gap-filling premise is "
                    "worth building on." if supported else
                    "FN and TN are equally gap-adjacent -- model gaps do NOT explain the false negatives, "
                    "so gap-filling is the wrong lever for this error set."),
        "caveats": [
            "CORRELATION, not proof of repair: this shows gaps mark the FNs, NOT that filling a specific "
            "gap flips a specific gene. The decisive test is a measured MCC delta after gap-filling.",
            "The mechanism is near-tautological in one direction -- a gene whose reactions are all blocked "
            "carries zero flux and MUST be called dispensable. The informative part is that this holds "
            "MORE for truly-essential genes (FN) than for truly-dispensable ones (TN).",
            "Gap-adjacency is a coarse proxy: 'controls >=1 blocked or dead-end-touching reaction'.",
            "Essentiality is medium-dependent; the model's DEFAULT medium is used throughout.",
        ],
    }
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"fba_gap_premise_{a.organism}_{a.date}.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    print(f"wrote {outdir / f'fba_gap_premise_{a.organism}_{a.date}.json'}")
    return 0 if supported else 1


if __name__ == "__main__":
    raise SystemExit(main())
