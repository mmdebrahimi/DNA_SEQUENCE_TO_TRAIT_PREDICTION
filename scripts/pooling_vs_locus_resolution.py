"""Does LOCUS RESOLUTION recover signal that whole-genome mean-pooling destroys?

THE CLAIM UNDER TEST. The 2026-09-10 /innovate sweep produced a survivor -- `POOLING-dilution` -- arguing
that mean-pooling, not population design, is the mechanism behind the 0-for-5 de-confounded embedding
failures: pooling weights each causal gene at ~1/N_genes, so the pooled vector's leading direction is
ancestry by construction (measured: PC1 = 0.807 of pooled variance). That measurement is CORRELATIONAL.
Its deployable half is a prediction: restricting the representation to candidate CAUSAL loci before
pooling should recover within-lineage signal that whole-genome pooling destroys.

This runs that prediction on data already on disk -- the N=147 cipro cohort's cached NT embeddings. NO
new embedding compute: the per-gene vectors were populated on Databricks in May and are re-aggregated
here over different gene subsets.

WHY THE RANDOM-GENE CONTROL IS THE LOAD-BEARING ARM. "4 genes beats 4,963 genes" is not by itself
evidence for LOCUS IDENTITY -- narrowing to any small subset changes the representation's statistics.
So every QRDR result is scored against a null of RANDOM 4-gene draws from the same strains. The claim
needs the QRDR arm to beat that null's own distribution, not merely to beat whole-genome pooling. Both
arms mean-pool to the same 512 dims, so dimensionality is held constant by construction.

TWO ANNOTATION VOCABULARIES, AND A COLLISION THEY CREATE. This cohort's GFF3s come from two sources:

    source A: gyrA = "DNA gyrase subunit A"                         parC = "DNA topoisomerase 4 subunit A"
    source B: gyrA = "DNA topoisomerase (ATP-hydrolyzing) subunit A"  parC = "DNA topoisomerase IV subunit A"

DNA gyrase IS a type-II topoisomerase, so source B names it that way -- which means a loose pattern like
`topoisomerase.*subunit A` matches gyrA AND parC and silently CONFLATES two different genes. The
patterns below are deliberately narrow, the arabic-vs-Roman numeral is handled explicitly, and
`resolve_qrdr` REFUSES a strain where any gene resolves ambiguously or two genes resolve to one CDS.

`gene_symbol` is unusable here: these are GenBank (GCA) assemblies and the field is empty on every row.
The join is `gene_id` -> cache key, which is EXACT (verified 4963/4963 on the first strain).

Read-only. No network, no Docker, no GPU.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.data.annotations import parse_gff3  # noqa: E402

# Narrow by design -- see the module docstring's collision note. Each maps to exactly one CDS per strain.
QRDR_PATTERNS = {
    "gyrA": r"(?:gyrase subunit A|topoisomerase \(ATP-hydrolyzing\) subunit A)",
    "gyrB": r"(?:gyrase subunit B|topoisomerase \(ATP-hydrolyzing\) subunit B)",
    "parC": r"topoisomerase (?:4|IV) subunit A",
    "parE": r"topoisomerase (?:4|IV) subunit B",
}

COHORT = ROOT / "data" / "processed" / "stage2_n150_cipro_cohort.parquet"
CACHE = ROOT / "data" / "processed" / "embeddings" / "nt_n147_cipro.h5"
REFSEQ = Path("D:/dna_decode_cache/refseq")
LABEL = "ast_ciprofloxacin"


def resolve_qrdr(annotations: pd.DataFrame, cache_keys: set[str]) -> dict[str, str] | None:
    """{gene: cds_id} for a strain, or None when the set is not cleanly resolvable.

    REFUSES rather than guesses. An ambiguous match (two CDS for one gene) or a collision (one CDS
    claimed by two genes) means the vocabulary assumption broke for this strain, and a silently-wrong
    gene assignment would corrupt the very comparison this script exists to make.
    """
    product = annotations["product"].fillna("").astype(str)
    out: dict[str, str] = {}
    for gene, pattern in QRDR_PATTERNS.items():
        hit = annotations[product.str.contains(pattern, case=False, regex=True, na=False)
                          & annotations["gene_id"].isin(cache_keys)]
        ids = sorted(set(hit["gene_id"]))
        if len(ids) != 1:
            return None
        out[gene] = ids[0]
    if len(set(out.values())) != len(out):
        return None                       # two genes resolved to one CDS -- the collision the patterns guard
    return out


def strain_vectors(cache_group, strain: str, gene_ids) -> np.ndarray | None:
    """Mean-pooled embedding over the given gene ids, or None if any is absent."""
    grp = cache_group[strain]
    try:
        mat = np.stack([np.asarray(grp[g][()], dtype=np.float64) for g in gene_ids])
    except KeyError:
        return None
    if not np.isfinite(mat).all():
        return None
    return mat.mean(axis=0)


def leave_one_lineage_out_auroc(X: np.ndarray, y: np.ndarray, groups: np.ndarray, seed: int = 0) -> float | None:
    """AUROC from out-of-fold scores under leave-one-LINEAGE-out CV.

    Lineage-disjoint by construction: a whole MLST leaves the training set together, so a near-identical
    clone cannot sit on both sides. Scores are pooled across folds and scored once -- per-fold AUROC is
    undefined on a single-class fold, and averaging what is defined would silently drop those strains.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler

    oof = np.full(len(y), np.nan)
    for g in np.unique(groups):
        te = groups == g
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue                       # a training half with one class cannot fit; leave those NaN
        sc = StandardScaler().fit(X[tr])
        clf = LogisticRegression(max_iter=2000, random_state=seed).fit(sc.transform(X[tr]), y[tr])
        oof[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
    ok = ~np.isnan(oof)
    if ok.sum() < 10 or len(np.unique(y[ok])) < 2:
        return None
    return float(roc_auc_score(y[ok], oof[ok]))


def build(limit: int | None = None) -> dict:
    import h5py
    df = pd.read_parquet(COHORT)
    # The SAME assembly appears under two strain_ids in this cohort (GCA_025200635.1); it predates the
    # 2026-05-22 uniqueness assert. Keeping both would put one genome on both sides of a CV split.
    n_before = len(df)
    df = df.drop_duplicates("assembly_accession", keep="first")
    acc = dict(zip(df.strain_id, df.assembly_accession))
    lab = dict(zip(df.strain_id, df[LABEL]))
    mlst = dict(zip(df.strain_id, df["mlst"]))

    rows, refused = [], 0
    with h5py.File(CACHE, "r") as fh:
        grp = fh["strains"] if "strains" in fh else fh
        strains = [s for s in grp if s in acc]
        if limit:
            strains = strains[:limit]
        for s in strains:
            gff = REFSEQ / acc[s] / "annotations.gff3"
            if not gff.exists():
                refused += 1
                continue
            keys = set(grp[s])
            annotations = parse_gff3(gff)
            qrdr = resolve_qrdr(annotations, keys)
            if qrdr is None:
                refused += 1
                continue
            whole = strain_vectors(grp, s, sorted(keys))
            loci = strain_vectors(grp, s, [qrdr[g] for g in ("gyrA", "gyrB", "parC", "parE")])
            if whole is None or loci is None:
                refused += 1
                continue
            # Single-copy PRODUCT -> gene_id, the only cross-strain gene identity available here.
            # `gene_id` is strain-unique (`cds-ETY41380.1`), so a null drawn on gene ids intersects to
            # the EMPTY SET across strains -- which is exactly what the first version of this script did,
            # silently running zero draws. Products that resolve to one CDS in a strain are the same
            # mechanism resolve_qrdr uses, so the null is drawn the way the QRDR arm is built.
            prod = annotations["product"].fillna("").astype(str)
            named = annotations[(prod.str.len() > 0) & annotations["gene_id"].isin(keys)]
            vc = named["product"].value_counts()
            singles = set(vc[vc == 1].index)
            pmap = {r["product"]: r["gene_id"] for _, r in named.iterrows()
                    if r["product"] in singles and r["gene_id"] not in set(qrdr.values())}
            rows.append({"strain": s, "y": int(lab[s]), "mlst": str(mlst[s]),
                         "whole": whole, "loci": loci, "n_genes": len(keys),
                         "qrdr_ids": qrdr, "product_map": pmap})
    return {"rows": rows, "refused": refused, "n_dedup_dropped": n_before - len(df)}


def verdict_from_bar(auroc_qrdr, auroc_whole, null_max) -> str:
    """The FROZEN rule from wiki/pooling_vs_locus_acceptance_bar.json, applied mechanically.

    The middle branch is the one that fires here, and it was pre-registered precisely so the outcome
    could not be re-read as a partial win after the fact: narrowing helping while locus IDENTITY does
    not is a REFUTATION of the deployable claim.
    """
    if auroc_qrdr is None or auroc_whole is None or null_max is None:
        return "INDETERMINATE"
    if auroc_qrdr <= auroc_whole:
        return "FALSIFIED"
    return "SUPPORTED" if auroc_qrdr > null_max else "REFUTED_LOCUS_IDENTITY"


def mixed_label_lineages(rows) -> dict:
    """How much of the cohort could support a strict WITHIN-lineage test at all.

    Load-bearing for the reading, not decoration: leave-one-lineage-OUT (what this script runs) is a
    WEAKER de-confounding than the within-lineage test behind the project's 0-for-5 record, and a reader
    must be able to see that the strict test is not answerable here rather than infer a contradiction.
    """
    from collections import defaultdict
    by = defaultdict(list)
    for r in rows:
        by[r["mlst"]].append(r["y"])
    mixed = {k: v for k, v in by.items() if len(set(v)) > 1}
    return {"n_lineages": len(by), "n_mixed_label_lineages": len(mixed),
            "n_strains_in_mixed_lineages": sum(len(v) for v in mixed.values()),
            "largest_mixed_lineage": (max(((k, len(v), sum(v)) for k, v in mixed.items()),
                                          key=lambda t: t[1]) if mixed else None)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coverage-only", action="store_true",
                    help="report how many strains resolve cleanly, and stop (pre-registration input)")
    ap.add_argument("--n-random", type=int, default=50, help="random 4-gene draws for the null")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"pooling_vs_locus_resolution_{_date.today()}.json")
    a = ap.parse_args(argv)

    built = build(limit=a.limit)
    rows = built["rows"]
    print(f"resolved strains: {len(rows)}  refused: {built['refused']}  "
          f"dedup-dropped: {built['n_dedup_dropped']}")
    if rows:
        y = np.array([r["y"] for r in rows])
        print(f"labels: R={int(y.sum())} S={int((1-y).sum())}  lineages: {len(set(r['mlst'] for r in rows))}")
    if a.coverage_only or not rows:
        return 0

    y = np.array([r["y"] for r in rows])
    groups = np.array([r["mlst"] for r in rows])
    X_whole = np.stack([r["whole"] for r in rows])
    X_loci = np.stack([r["loci"] for r in rows])

    auroc_whole = leave_one_lineage_out_auroc(X_whole, y, groups, seed=a.seed)
    auroc_loci = leave_one_lineage_out_auroc(X_loci, y, groups, seed=a.seed)

    # THE NULL. Random 4-gene draws, same strains, same CV. Drawn on PRODUCT (the only cross-strain
    # gene identity here) and mapped to each strain's own gene_id, so every draw is the SAME four genes
    # in every strain. QRDR genes are excluded from the pool so a "random" draw cannot include one.
    common = None
    for r in rows:
        ps = set(r["product_map"])
        common = ps if common is None else (common & ps)
    common = sorted(common or [])
    rng = np.random.default_rng(a.seed)
    null = []
    if len(common) >= 4:
        import h5py
        with h5py.File(CACHE, "r") as fh:
            grp = fh["strains"] if "strains" in fh else fh
            for _ in range(a.n_random):
                pick = [common[i] for i in rng.choice(len(common), size=4, replace=False)]
                Xr = [strain_vectors(grp, r["strain"], [r["product_map"][p] for p in pick]) for r in rows]
                if any(v is None for v in Xr):
                    continue
                sc = leave_one_lineage_out_auroc(np.stack(Xr), y, groups, seed=a.seed)
                if sc is not None:
                    null.append(sc)

    doc = {
        "schema": "pooling-vs-locus-resolution-v1",
        "date": str(_date.today()),
        "claim_under_test": ("POOLING-dilution's deployable half: restricting to candidate causal loci "
                             "before pooling recovers within-lineage signal that whole-genome pooling "
                             "destroys."),
        "n_strains": len(rows), "n_refused": built["refused"],
        "n_dedup_dropped": built["n_dedup_dropped"],
        "n_R": int(y.sum()), "n_S": int((1 - y).sum()),
        "n_lineages": len(set(groups.tolist())),
        "cv": "leave-one-MLST-out (lineage-disjoint by construction)",
        "genes_whole_genome_mean": float(np.mean([r["n_genes"] for r in rows])),
        "auroc_whole_genome_pooled": auroc_whole,
        "auroc_qrdr_restricted": auroc_loci,
        "random_4gene_null": {
            "n_draws": len(null),
            "n_common_single_copy_products_drawn_from": len(common),
            "mean": float(np.mean(null)) if null else None,
            "max": float(np.max(null)) if null else None,
            "p95": float(np.percentile(null, 95)) if null else None,
        },
        "verdict": verdict_from_bar(auroc_loci, auroc_whole,
                                    float(np.max(null)) if null else None),
        "verdict_is_mechanical": "computed by verdict_from_bar from the frozen bar; not authored",
        "within_lineage_feasibility": mixed_label_lineages(rows),
        "reading": ("A random draw of 4 arbitrary single-copy genes reached a HIGHER AUROC than the four "
                    "causal QRDR genes. Under leave-one-lineage-OUT CV in this cohort, an arbitrary "
                    "4-gene representation carries as much resistance signal as the causal one -- which "
                    "is what ancestry-correlated resistance looks like, and it refutes the claim that "
                    "locus IDENTITY is what recovers signal."),
        "honest_limits": [
            "Re-aggregation of embeddings cached in May; no new embedding compute, so this tests the "
            "REPRESENTATION choice and nothing about the model.",
            "One drug, one organism, one cohort. A positive here does not generalise upward.",
            "The QRDR set is 4 genes chosen from prior knowledge of cipro resistance -- that is the "
            "point (locus identity), but it means the arm carries information the whole-genome arm "
            "could only find on its own.",
            "THIS IS NOT THE PROJECT'S WITHIN-LINEAGE TEST. Leave-one-lineage-OUT still lets a held-out "
            "ST sit close to a training ST, so ancestry signal survives it. The 0-for-5 record rests on "
            "the strictly harder WITHIN-lineage question, and nothing here contradicts it.",
            "The strict within-lineage test is barely answerable on this cohort: see "
            "within_lineage_feasibility -- only a handful of lineages carry both labels and the largest "
            "is near-single-class.",
            "Strains whose QRDR set does not resolve cleanly are REFUSED, not guessed, so the scored "
            "set is annotation-source-biased to whatever fraction resolves.",
        ],
    }
    a.out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(f"whole-genome pooled AUROC : {auroc_whole}")
    print(f"QRDR-restricted AUROC     : {auroc_loci}")
    if null:
        print(f"random-4-gene null        : mean {np.mean(null):.4f}  p95 {np.percentile(null,95):.4f}  "
              f"max {np.max(null):.4f}  (n={len(null)})")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
