"""External-validity of the v0 cross-axis ExPEC rule on the Horesh 2021 independent-label set.

Closes the `markers.py` scope-limit "K=1 cross-axis tested in-sample on N=24". The shipped rule
(`meets_cross_axis_support`: >=1 iron-acquisition gene AND >=1 capsule/serum gene) is applied to the
recorded virulence-gene calls (`Other_Vir_genes`) of the Horesh curated, isolation-source-labelled
strains:
  - POSITIVES: ExPEC (blood/urine isolation source -> label independent of the resolver's marker rules).
  - NEGATIVES (specificity proxy): EPEC + ETEC (intestinal/feces) and Not-determined/commensal.

SCOPE (honest): this validates the RULE LOGIC against an independent-LABEL population, holding
gene-detection constant (Horesh's ariba/VirulenceFinder calls are the same gene-detection FAMILY the
resolver uses). It is NOT a from-FASTA re-run and does NOT re-validate gene detection. For EPEC, the rule
firing is harmless in the live resolver because the DEC-module gate resolves LEE (eae) ABOVE the ExPEC
branch; the pure specificity signal is the commensal/Not-determined negatives.

Run: uv run python scripts/expec_rule_external_validity.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.pathotype.expec_score import iron_capsule_counts, meets_cross_axis_support
from dna_decode.pathotype.markers import CLUSTER_MARKERS, EXPEC_STRONG

F1 = Path("data/horesh_2021/F1_genome_metadata.csv")

# strong-marker prefixes (>=2 distinct clusters -> confident UPEC in the live resolver)
_STRONG_PREFIXES = {c: tuple(CLUSTER_MARKERS[c]) for c in EXPEC_STRONG}


def gene_symbols(cell) -> list[str]:
    """Parse 'iroN_6_AF449498;iss_9_CP002167' -> ['iron','iss'] (lowercased gene symbol = token[0])."""
    if not isinstance(cell, str) or not cell.strip():
        return []
    return [tok.split("_", 1)[0].lower() for tok in cell.split(";") if tok.strip()]


def pergene_cov_from_symbols(symbols: list[str]) -> dict[str, float]:
    """Build the {support-prefix -> 1.0} dict the resolver's cross-axis rule consumes. A support prefix
    is 'present' (cov 1.0) iff some gene symbol startswith-matches it (same match as detect.build_vf_index)."""
    from dna_decode.pathotype.markers import EXPEC_SUPPORT_GENE_PREFIXES
    cov = {}
    for p in EXPEC_SUPPORT_GENE_PREFIXES:
        if any(s.startswith(p) for s in symbols):
            cov[p] = 1.0
    return cov


def _strong_count(symbols: list[str]) -> int:
    return sum(1 for prefixes in _STRONG_PREFIXES.values()
               if any(s.startswith(p) for s in symbols for p in prefixes))


def evaluate(df: pd.DataFrame, label: str) -> dict:
    n = len(df)
    fired = 0
    confident_strong = 0
    has_genes = 0
    for cell in df["Other_Vir_genes"]:
        syms = gene_symbols(cell)
        if syms:
            has_genes += 1
        cov = pergene_cov_from_symbols(syms)
        if meets_cross_axis_support(cov):
            fired += 1
        if _strong_count(syms) >= 2:
            confident_strong += 1
    return {"label": label, "n": n, "has_gene_calls": has_genes,
            "cross_axis_fired": fired, "cross_axis_rate": round(fired / n, 3) if n else None,
            "confident_strong_ge2": confident_strong,
            "confident_strong_rate": round(confident_strong / n, 3) if n else None}


def main() -> int:
    if not F1.exists():
        print(f"ERROR: {F1} not found (download Figshare 13270073 file 25552514).", file=sys.stderr)
        return 2
    df = pd.read_csv(F1)
    groups = {
        "ExPEC (independent positive)": df[df["Pathotype"] == "ExPEC"],
        "EPEC (intestinal; DEC-gated negative)": df[df["Pathotype"] == "EPEC"],
        "ETEC (intestinal; DEC-gated negative)": df[df["Pathotype"] == "ETEC"],
        "Not determined (commensal-ish negative)": df[df["Pathotype"] == "Not determined"],
    }
    rows = [evaluate(g, name) for name, g in groups.items()]
    print(f"{'group':<42} {'N':>5} {'genes':>6} {'x-axis':>7} {'rate':>6} {'strong>=2':>9} {'rate':>6}")
    for r in rows:
        print(f"{r['label']:<42} {r['n']:>5} {r['has_gene_calls']:>6} {r['cross_axis_fired']:>7} "
              f"{r['cross_axis_rate']!s:>6} {r['confident_strong_ge2']:>9} {r['confident_strong_rate']!s:>6}")
    # recall is measured ONLY on strains that have gene calls (no-call rows are a detection gap, not a rule miss)
    ex = groups["ExPEC (independent positive)"]
    ex_called = ex[ex["Other_Vir_genes"].apply(lambda c: bool(gene_symbols(c)))]
    ex_fired_called = sum(1 for c in ex_called["Other_Vir_genes"]
                          if meets_cross_axis_support(pergene_cov_from_symbols(gene_symbols(c))))
    print()
    print(f"ExPEC cross-axis RECALL on gene-called strains: {ex_fired_called}/{len(ex_called)} "
          f"= {ex_fired_called/len(ex_called):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
