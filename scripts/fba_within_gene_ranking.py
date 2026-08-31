"""Does the model's growth ratio rank a gene's OWN conditions correctly? The de-confounded switch test.

THE LEVER THIS CLOSES. `conditional_essentiality.continuous_readout` scores the raw knockout growth ratio
as a ranking POOLED over every gene x condition cell and reports AUROC ~0.59 -- above chance, which reads
as "the 1% cutoff is discarding real signal". But a pooled ranking is dominated by the GENE MAIN EFFECT: a
gene that is essential in all four media has a low ratio in all four, and contributes four correctly-ranked
positives without the model ever having switched. A pooled 0.59 is achievable with EXACTLY ZERO
within-gene signal.

The conditional-switch question is strictly within-gene: for one gene, is the ratio LOWER in the conditions
where that gene is actually essential? Conditioning on the gene removes its main effect by construction --
the same de-confounding idiom this project already uses for lineage, clonality and ancestry, applied to the
axis that was never conditioned on.

PRE-REGISTERED BEFORE THE FIRST RUN:
  * PRIMARY: mean within-gene AUROC over NON-FLAT conditionally-essential genes.
  * PASS  : > 0.60 AND permutation-null p < 0.05  -> the cutoff discards real switch signal; a per-gene
            RELATIVE rule (rank a gene's own conditions) is worth building.
  * FAIL  : ~0.5 -> the model's cross-condition variation is uninformative even where it exists; the
            bottleneck is the model/data, not the readout, and this lever closes.
  * MUST-HOLD: the flat fraction reproduces the 2026-08-12 artifact's ~64%, proving this runs on the same
            subset. If it does not, the numbers are not comparable and nothing here is interpretable.

FLAT GENES ARE THE POINT, NOT AN EXCLUSION. A gene whose ratio is identical across all four media has
within-gene AUROC exactly 0.5 by ties -- including it would drag the mean toward chance and hide whatever
the varying genes do. Both are reported: `mean_auroc_all` (honest headline over every gene) and
`mean_auroc_nonflat` (the pre-registered primary, which asks whether variation, WHERE IT EXISTS, points
the right way).

DETERMINISM FIRST. Degenerate LP optima shift mid-range ratios between processes, which is why the
committed artifact quotes AUROC as ~0.60 rather than to four decimals. This runs single-process and
`--repeat 2` re-solves everything independently and reports whether the primary agrees to 3 decimals. A
bar cleared by a run whose variance exceeds the effect is not a result.

Run: uv run python scripts/fba_within_gene_ranking.py [--repeat 2]
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

FLAT_EPS = 1e-9
DEPLOYED_FRAC = 0.01   # the shipped essentiality cutoff: growth below 1% of wildtype = essential
PRIMARY_BAR = 0.60
N_PERM = 2000


# Where each axis's DEPLOYED exact-set number lives. Read from the artifact, never hardcoded: quoting the
# 4-media 3/67 beside a 25-condition carbon result would compare a ceiling against the wrong baseline and
# make the ranking lever look ~4x more valuable than it is on the axis where it is best measured.
DEPLOYED_SOURCE = {
    "media4": ("fba_conditional_essentiality_ecoli_*.json", ("model_scored", "switch")),
    "carbon": ("fba_conditional_carbon_*.json", ("switch",)),
}


def deployed_exact_set(axis: str) -> dict | None:
    """Pull this axis's own deployed exact-set from its newest committed artifact."""
    import glob

    pattern, path = DEPLOYED_SOURCE.get(axis, (None, None))
    if not pattern:
        return None
    files = sorted(glob.glob(str(ROOT / "wiki" / pattern)))
    if not files:
        return None
    node = json.loads(pathlib_read(files[-1]))
    for key in path:
        node = node.get(key, {})
    # Schema drift between two generations of the same producer: the carbon artifact carries
    # `n_scored_exact_set` (added with `exclude_cells`), the older 4-media one only
    # `n_conditionally_essential`. Accept either rather than silently reporting "unknown".
    m = node.get("exact_set_match")
    n = node.get("n_scored_exact_set", node.get("n_conditionally_essential"))
    if m is None or n is None:
        return None
    return {"match": m, "n": n, "rate": round(m / n, 4), "source": Path(files[-1]).name}


def pathlib_read(f: str) -> str:
    return Path(f).read_text(encoding="utf-8")


def within_gene_auroc(essential: list[bool], ratio: list[float],
                      tol: float = FLAT_EPS) -> float | None:
    """AUROC over ONE gene's conditions: does a lower ratio mark the essential ones? PURE.

    Returns None when the gene is not two-sided (all-essential or all-dispensable) -- undefined, not 0.5.

    THE TIE TOLERANCE IS LOAD-BEARING, and the first version got it wrong. Flatness was judged with a
    1e-9 tolerance while the tie test used exact float equality, so a gene whose four ratios differ in the
    15th decimal was called flat AND scored as if its ordering meant something: 36 of 41 flat genes
    returned an AUROC other than 0.5, purely from LP noise. The two tests must share one tolerance. The
    arithmetic that exposed it: 41 genes at exactly 0.5 plus 26 at 0.718 cannot average 0.6045.
    """
    pos = [r for r, e in zip(ratio, essential) if e]
    neg = [r for r, e in zip(ratio, essential) if not e]
    if not pos or not neg:
        return None
    # P(ratio_essential < ratio_dispensable), ties counted a half. Lower ratio = more essential.
    wins = sum(0.5 if abs(p - n) <= tol else (1.0 if p < n else 0.0) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def permutation_p(genes: list[tuple[list[bool], list[float]]], observed: float,
                  n_perm: int = N_PERM, seed: int = 0) -> float:
    """Shuffle the essential/dispensable labels WITHIN each gene; how often does chance beat observed?

    Within-gene shuffling holds each gene's ratio profile AND its number of essential conditions fixed, so
    the null isolates exactly one thing: whether the pairing of condition to ratio is informative.
    """
    import random

    rng = random.Random(seed)
    ge = 0
    for _ in range(n_perm):
        vals = []
        for ess, rat in genes:
            shuffled = ess[:]
            rng.shuffle(shuffled)
            a = within_gene_auroc(shuffled, rat)
            if a is not None:
                vals.append(a)
        if vals and sum(vals) / len(vals) >= observed:
            ge += 1
    return (ge + 1) / (n_perm + 1)


def compute_ratios_carbon(db: str | None, threshold: float) -> tuple[dict, list, list[str]]:
    """The 25-source Fitness Browser Keio CARBON axis -- an INDEPENDENT, far better-powered replication.

    The 4-media Orth panel gives each gene only 4 conditions, so its within-gene AUROC can take just a
    handful of discrete values and rests on 26 varying genes. 25 conditions over ~200 genes is a
    materially finer measurement of the same quantity, on a different substrate with a different label
    source (Keio transposon fitness, not Orth's curated E/N calls).

    Loading mirrors `scripts/fba_conditional_carbon_validate.py` exactly so the two are comparable --
    same `carbon_conditions` mapping, same `load_records` gene filter and threshold, same
    `apply_carbon_condition(all_carbon=...)` which closes every OTHER carbon exchange (without that the
    conditions are not actually distinct media).
    """
    from cobra.flux_analysis import single_gene_deletion

    from dna_decode.fba.conditional_essentiality import conditionally_essential_genes
    from dna_decode.fba.fitness_browser import (apply_carbon_condition, carbon_conditions, load_records,
                                                open_db)
    from dna_decode.fba.model import load_model, wildtype_growth

    conn = open_db(db)
    model = load_model()
    conds = carbon_conditions(conn, model)
    keys = sorted(conds)
    records = load_records(conn, conds, gene_filter={g.id for g in model.genes}, threshold=threshold)
    subset = conditionally_essential_genes(records)
    print(f"  {len(conds)} mappable carbon sources | {len(records)} complete rows | "
          f"{len(subset)} conditionally essential", flush=True)
    if not subset:
        return {}, [], keys

    genes = [r.gene_id for r in subset]
    all_ex = tuple(conds.values())
    ratios: dict[str, dict[str, float]] = {}
    for n, cond in enumerate(keys, 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            wt = wildtype_growth(model)
            rat: dict[str, float] = {}
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
                for idx, row in res.iterrows():
                    key = row["ids"] if "ids" in res.columns else idx
                    gid = next(iter(key)) if not isinstance(key, str) else key
                    g = row["growth"]
                    rat[gid] = 0.0 if g != g else g / wt
            ratios[cond] = rat
        if n % 5 == 0 or n == len(keys):
            print(f"    [{n:2d}/{len(keys)}] deletions done", flush=True)
    return ratios, subset, keys


def compute_ratios_nitrogen(db: str | None, threshold: float) -> tuple[dict, list, list[str]]:
    """The Keio NITROGEN axis -- a THIRD substrate, run to test a PREDICTION rather than to repeat a result.

    PRE-REGISTERED PREDICTION: flatness rose 61.2% (4 media) -> 68.2% (25 carbon sources). If a gene's
    ratio is flat because the AXIS itself carries little dynamic range, nitrogen should be flatter still --
    six of its thirteen conditions give an IDENTICAL wildtype growth (0.92593), so the model literally
    cannot distinguish those media. A confirmation turns flatness from an observation into a mechanism
    with a predictor; a refutation means flatness is a property of the genes, not of the axis.

    Mirrors `dna_decode/fba/nitrogen.py`'s own loaders exactly (`nitrogen_conditions` +
    `load_nitrogen_records` + `apply_nitrogen_condition(all_nitrogen=...)`, which closes every OTHER
    nitrogen exchange while holding glucose as the carbon source).
    """
    from cobra.flux_analysis import single_gene_deletion

    from dna_decode.fba.conditional_essentiality import conditionally_essential_genes
    from dna_decode.fba.fitness_browser import open_db
    from dna_decode.fba.model import load_model, wildtype_growth
    from dna_decode.fba.nitrogen import (apply_nitrogen_condition, load_nitrogen_records,
                                         nitrogen_conditions)

    conn = open_db(db)
    model = load_model()
    conds = nitrogen_conditions(conn, model)
    keys = sorted(conds)
    records = load_nitrogen_records(conn, conds, gene_filter={g.id for g in model.genes},
                                    threshold=threshold)
    subset = conditionally_essential_genes(records)
    print(f"  {len(conds)} mappable nitrogen sources | {len(records)} complete rows | "
          f"{len(subset)} conditionally essential", flush=True)
    if not subset:
        return {}, [], keys

    genes = [r.gene_id for r in subset]
    all_ex = tuple(conds.values())
    ratios: dict[str, dict[str, float]] = {}
    wt_seen: list[float] = []
    for n, cond in enumerate(keys, 1):
        with model:
            apply_nitrogen_condition(model, conds[cond], all_nitrogen=all_ex)
            wt = wildtype_growth(model)
            wt_seen.append(round(wt, 5))
            rat: dict[str, float] = {}
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
                for idx, row in res.iterrows():
                    key = row["ids"] if "ids" in res.columns else idx
                    gid = next(iter(key)) if not isinstance(key, str) else key
                    g = row["growth"]
                    rat[gid] = 0.0 if g != g else g / wt
            ratios[cond] = rat
        if n % 4 == 0 or n == len(keys):
            print(f"    [{n:2d}/{len(keys)}] deletions done", flush=True)
    print(f"  distinct wildtype growths across the axis: {len(set(wt_seen))}/{len(wt_seen)} "
          f"<- the axis's own dynamic range", flush=True)
    return ratios, subset, keys


def compute_ratios(organism: str = "ecoli") -> tuple[dict, list, dict]:
    """Re-solve the deletion panel. Single-process by design (see the determinism note in the docstring)."""
    from cobra.flux_analysis import single_gene_deletion

    from dna_decode.fba.conditional_essentiality import (CONDITIONS, apply_condition,
                                                         conditionally_essential_genes, load_labels)
    from dna_decode.fba.model import load_model, wildtype_growth

    records = load_labels()
    model = load_model(organism=organism)
    ids = {g.id for g in model.genes}
    scored = [r for r in records if r.gene_id in ids]

    ratios: dict[str, dict[str, float]] = {}
    wt_by_cond: dict[str, float] = {}
    for c in CONDITIONS:
        with model:
            apply_condition(model, c)
            wt = wildtype_growth(model)
            wt_by_cond[c] = round(wt, 6)
            rat: dict[str, float] = {}
            if wt > 1e-6:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(r.gene_id) for r in scored], processes=1)
                for idx, row in res.iterrows():
                    key = row["ids"] if "ids" in res.columns else idx
                    gid = next(iter(key)) if not isinstance(key, str) else key
                    g = row["growth"]
                    # NaN is a genuinely infeasible LP -> total loss of growth. Documented at length in
                    # wiki/fba_infeasibility_finding_2026-08-13.md: do NOT "fix" this to a skip.
                    rat[gid] = 0.0 if g != g else g / wt
            ratios[c] = rat
    return ratios, conditionally_essential_genes(scored), wt_by_cond


def score(ratios: dict, cond_genes: list, keys: list[str] | None = None) -> dict:
    if keys is None:
        from dna_decode.fba.conditional_essentiality import CONDITIONS
        keys = sorted(CONDITIONS)
    per_gene: list[dict] = []
    for r in cond_genes:
        if not all(r.gene_id in ratios.get(c, {}) for c in keys):
            continue
        ess = [bool(r.experimental[c]) for c in keys]
        rat = [ratios[c][r.gene_id] for c in keys]
        a = within_gene_auroc(ess, rat)
        if a is None:                       # not two-sided -> the switch question is undefined for it
            continue
        per_gene.append({"gene_id": r.gene_id, "auroc": a,
                         "flat": (max(rat) - min(rat)) < FLAT_EPS,
                         "spread": round(max(rat) - min(rat), 6), "_ess": ess, "_rat": rat})

    # What could a RELATIVE rule buy in the metric this project tracks (exact-set match)? For each gene
    # take k = its TRUE number of essential conditions and call the k lowest-ratio conditions essential.
    # It is handed k, which a deployed rule would have to infer, so this is a strict CEILING on any
    # rank-based per-gene rule and NEVER a deployable number -- the same rail `continuous_readout`
    # carries on its oracle threshold.
    def oracle_hit(e: list[bool], rt: list[float]) -> bool:
        """A FLAT gene has no ordering, so the oracle is UNDEFINED for it -- never a hit.

        The first version omitted this and reported 23/67. On a flat gene every ratio is equal, so
        `sorted` falls back to stable index order and "the k lowest conditions" means "the first k
        conditions"; 12 of the 41 flat genes matched their true pattern purely because of how the
        condition list happens to be ordered. That is tie-breaking luck presented as a ceiling. Flat
        genes are structurally unreachable by ANY change of readout -- the model emits one number for
        all four media -- so they belong in the denominator and never in the numerator.
        """
        if (max(rt) - min(rt)) < FLAT_EPS:
            return False
        k = sum(e)
        order = sorted(range(len(rt)), key=lambda i: rt[i])
        # An AMBIGUOUS top-k boundary is not a hit either -- the partial version of the flat-gene trap.
        # If the k-th and (k+1)-th ratios tie, which condition lands inside the selection is decided by
        # stable-sort index order, not by the model. Two genes at AUROC 0.833 were "hitting" that way.
        # With this guard the ceiling equals the count of AUROC==1.0 genes exactly, as it must: top-k
        # selection is right precisely when every essential condition ranks below every dispensable one.
        if k < len(rt) and abs(rt[order[k - 1]] - rt[order[k]]) <= FLAT_EPS:
            return False
        pred = [False] * len(rt)
        for i in order[:k]:
            pred[i] = True
        return pred == e

    hits = [oracle_hit(g["_ess"], g["_rat"]) for g in per_gene]

    # THE DENOMINATOR QUESTION, and it is not cosmetic. A gene whose ratio is flat gets a CONSTANT
    # threshold call, and a constant call can never match a two-sided truth -- every conditionally
    # essential gene is essential in some conditions and not others by definition. Scoring the model
    # against all 217 genes therefore scores it against a target it structurally cannot hit, and the
    # published 23/217 = 10.6% is that number. Verified on the committed carbon artifact before relying
    # on it: `commit_strata.predicted_constant` = 184 genes, 0 exact-set matches.
    #
    # Three nested strata, each answering a different question:
    #   flat            -- the model emits ONE ratio for every condition. It cannot distinguish them.
    #   varies_subthr   -- the ratio moves but never crosses the cutoff, so the CALL is still constant.
    #                      This is exactly the stratum a ranking rule could rescue.
    #   commits         -- the call itself varies. The only stratum where an exact-set hit is possible.
    for g in per_gene:
        rat = g["_rat"]
        call = [r < DEPLOYED_FRAC for r in rat]
        g["deployed_hit"] = (call == g["_ess"])
        g["stratum"] = ("flat" if g["flat"]
                        else "commits" if len(set(call)) > 1
                        else "varies_subthr")
    # Is "70% right when it commits" impressive? Unanchored it is just a number. The null: if the model
    # kept the SAME number of essential conditions it actually predicted but placed them at random, how
    # many exact matches would it get? A gene whose predicted count differs from the truth's count can
    # never match at all, so its chance contribution is 0 -- which is itself part of the difficulty.
    from math import comb

    chance_hits = 0.0
    for g in per_gene:
        call = [r < DEPLOYED_FRAC for r in g["_rat"]]
        if len(set(call)) <= 1:
            continue                       # a constant call cannot match a two-sided truth, chance 0
        k_pred, k_true, n = sum(call), sum(g["_ess"]), len(call)
        if k_pred == k_true and 0 < k_true < n:
            chance_hits += 1.0 / comb(n, k_true)

    # SECOND NULL, and a strictly harder one. The comb-based null treats conditions as INTERCHANGEABLE:
    # it asks only "how many placements of k essential conditions are there". If true essentiality is
    # concentrated in a few substrates AND the model tends to break on those same substrates, a model
    # could beat that null by learning the marginal shape rather than the per-gene placement.
    # This one shuffles the TRUTH matrix preserving BOTH margins -- every gene keeps its number of
    # essential conditions and every condition keeps its number of essential genes -- while the model's
    # predictions stay fixed. It reuses the repo's tested curveball implementation, which raises if a
    # margin ever breaks rather than silently returning an invalid null.
    from dna_decode.fba.nulls import margin_preserving_null

    gene_ids = [g["gene_id"] for g in per_gene]
    conds = tuple(keys)
    # `nulls` consumes {condition: {gene: bool}} -- the orientation every FBA caller uses. Building it
    # gene-keyed raised KeyError on the first shuffle, which is the right way for this to fail.
    truth_calls = {c: {g["gene_id"]: g["_ess"][i] for g in per_gene} for i, c in enumerate(keys)}
    model_call = {g["gene_id"]: [r < DEPLOYED_FRAC for r in g["_rat"]] for g in per_gene}
    committing = {g["gene_id"] for g in per_gene if len(set(model_call[g["gene_id"]])) > 1}

    def _exact_hits(shuffled_truth: dict) -> float:
        return float(sum(1 for gid in committing
                         if model_call[gid] == [bool(shuffled_truth[c][gid]) for c in keys]))

    marginal = margin_preserving_null(gene_ids, conds, truth_calls, _exact_hits,
                                      n_draws=200, seed0=0) if committing else None

    strata = {}
    for name in ("flat", "varies_subthr", "commits"):
        members = [g for g in per_gene if g["stratum"] == name]
        strata[name] = {
            "n_genes": len(members),
            "deployed_exact": sum(1 for g in members if g["deployed_hit"]),
            "oracle_exact": sum(1 for g, h in zip(per_gene, hits)
                                if g["stratum"] == name and h),
        }
    nonflat = [(g["_ess"], g["_rat"]) for g in per_gene if not g["flat"]]
    aurocs_all = [g["auroc"] for g in per_gene]
    aurocs_nf = [g["auroc"] for g in per_gene if not g["flat"]]
    n_flat = sum(1 for g in per_gene if g["flat"])

    return {"n_genes_scored": len(per_gene), "n_flat": n_flat, "n_nonflat": len(nonflat),
            "flat_fraction": round(n_flat / len(per_gene), 4) if per_gene else None,
            "mean_auroc_all": round(sum(aurocs_all) / len(aurocs_all), 4) if aurocs_all else None,
            "mean_auroc_nonflat": round(sum(aurocs_nf) / len(aurocs_nf), 4) if aurocs_nf else None,
            "strata": strata,
            "deployed_exact_set_recomputed": sum(1 for g in per_gene if g["deployed_hit"]),
            "chance_exact_hits_among_committing": round(chance_hits, 4),
            "marginal_preserving_null": marginal,
            "oracle_relative_exact_set": sum(hits),
            "oracle_relative_exact_set_nonflat": sum(
                h for h, g in zip(hits, per_gene) if not g["flat"]),
            "per_gene": [{k: v for k, v in g.items() if not k.startswith("_")} for g in per_gene],
            "_nonflat_pairs": nonflat}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repeat", type=int, default=2, help="independent re-solves; determinism check")
    ap.add_argument("--organism", default="ecoli")
    ap.add_argument("--axis", default="media4", choices=("media4", "carbon", "nitrogen"),
                    help="media4 = the 4-media Orth panel; carbon = the 25-source Keio axis")
    ap.add_argument("--db", default=None, help="feba.db path (carbon axis only)")
    ap.add_argument("--threshold", type=float, default=-2.0)
    args = ap.parse_args()

    runs = []
    for i in range(args.repeat):
        print(f"solving deletion panel, run {i + 1}/{args.repeat} (single-process) ...", flush=True)
        if args.axis == "carbon":
            ratios, cond_genes, keys = compute_ratios_carbon(args.db, args.threshold)
        elif args.axis == "nitrogen":
            ratios, cond_genes, keys = compute_ratios_nitrogen(args.db, args.threshold)
        else:
            ratios, cond_genes, _ = compute_ratios(args.organism)
            keys = None
        s = score(ratios, cond_genes, keys)
        runs.append(s)
        print(f"  genes {s['n_genes_scored']}  flat {s['n_flat']} ({s['flat_fraction']})  "
              f"mean AUROC all {s['mean_auroc_all']}  nonflat {s['mean_auroc_nonflat']}")

    primary = runs[0]
    # Determinism is a TOLERANCE question, not an equality one: degenerate LP optima shift mid-range
    # ratios between processes, so bit-identical repeats are not achievable and demanding them would
    # reject a real effect. The honest test is whether the between-run SPREAD is small against the
    # distance from the bar -- a bar cleared by a run whose variance exceeds the effect is not a result.
    vals = [r["mean_auroc_nonflat"] for r in runs if r["mean_auroc_nonflat"] is not None]
    spread = round(max(vals) - min(vals), 4) if len(vals) > 1 else None
    margin = round(min(vals) - PRIMARY_BAR, 4) if vals else None
    agree = (spread is not None and margin is not None and spread < margin) if len(vals) > 1 else None

    print("\npermutation null (within-gene label shuffle) ...", flush=True)
    p = permutation_p(primary["_nonflat_pairs"], primary["mean_auroc_nonflat"] or 0.5)

    obs = primary["mean_auroc_nonflat"]
    # The WORST repeat must clear the bar, not the first one -- otherwise the verdict is a coin flip on
    # which run happened to be reported.
    passed = bool(vals) and min(vals) > PRIMARY_BAR and p < 0.05 and (agree is not False)
    # The ~64% flat-fraction must-hold pins the 4-media panel to the 2026-08-12 artifact. The carbon axis
    # is a DIFFERENT substrate with 25 conditions, so it has no prior to reproduce -- asserting one there
    # would be inventing a bar. Reported as None rather than silently passing or failing.
    must_hold = (primary["flat_fraction"] is not None and 0.55 <= primary["flat_fraction"] <= 0.75
                 ) if args.axis == "media4" else None

    out = {"schema": "fba-within-gene-ranking-v1", "generated": date.today().isoformat(),
           "organism": args.organism, "axis": args.axis, "pre_registered": {"bar_mean_auroc_nonflat": PRIMARY_BAR,
                                                         "bar_perm_p": 0.05,
                                                         "must_hold_flat_fraction": "0.55-0.75 (~0.64)"},
           "n_genes_scored": primary["n_genes_scored"], "n_flat": primary["n_flat"],
           "n_nonflat": primary["n_nonflat"], "flat_fraction": primary["flat_fraction"],
           "mean_auroc_all": primary["mean_auroc_all"],
           "strata": primary["strata"],
           "deployed_exact_set_recomputed": primary["deployed_exact_set_recomputed"],
           "chance_exact_hits_among_committing": primary["chance_exact_hits_among_committing"],
           "marginal_preserving_null": primary["marginal_preserving_null"],
           "oracle_relative_exact_set": primary["oracle_relative_exact_set"],
           "oracle_relative_exact_set_nonflat": primary["oracle_relative_exact_set_nonflat"],
           "deployed_exact_set_this_axis": deployed_exact_set(args.axis),
           "oracle_note": ("the relative-rule exact-set is handed each gene's TRUE essential-condition "
                           "count k, which a deployed rule must infer. It is a CEILING, not a number "
                           "that could ship."),
           "mean_auroc_nonflat": obs, "permutation_p": round(p, 5), "n_permutations": N_PERM,
           "repeats": args.repeat, "repeat_values": vals,
           "between_run_spread": spread, "margin_over_bar_worst_run": margin,
           "spread_smaller_than_margin": agree,
           "must_hold_flat_fraction_reproduced": must_hold,
           "verdict": ("PASS" if passed else "FAIL"),
           "interpretation": (
               "Within-gene AUROC conditions on the gene, removing the gene main effect that the pooled "
               "readout's ~0.59 cannot separate from a real switch. PASS would mean the 1% cutoff discards "
               "genuine conditional signal and a per-gene RELATIVE rule is worth building. FAIL means the "
               "model's cross-condition variation is uninformative even where it exists, so the readout is "
               "not the bottleneck."),
           "per_gene": primary["per_gene"]}

    tag = "" if args.axis == "media4" else f"_{args.axis}"
    for f in (ROOT / "wiki" / f"fba_within_gene_ranking{tag}_{out['generated']}.json",):
        f.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print(f"\n  pooled comparison ....... continuous_readout AUROC ~0.59 (gene main effect NOT removed)")
    print(f"  within-gene, all genes .. {out['mean_auroc_all']}   (n={out['n_genes_scored']})")
    print(f"  within-gene, non-flat ... {obs}   (n={out['n_nonflat']})  <- PRE-REGISTERED PRIMARY")
    print(f"  permutation p ........... {out['permutation_p']}")
    dep = out["deployed_exact_set_this_axis"]
    dep_s = f"{dep['match']}/{dep['n']}" if dep else "unknown"
    head = (out["oracle_relative_exact_set"] - dep["match"]) if dep else None
    print(f"  oracle relative exact-set {out['oracle_relative_exact_set']}/{out['n_genes_scored']} "
          f"vs THIS AXIS's deployed {dep_s}"
          + (f"  -> headroom {head:+d} genes" if head is not None else "")
          + "   <- CEILING, handed true k")
    st = out["strata"]
    print("")
    print("  ANATOMY -- three nested strata (a constant call can NEVER match a two-sided truth):")
    for name, label in (("flat", "flat: one ratio for all conds"),
                        ("varies_subthr", "varies, never crosses cutoff"),
                        ("commits", "call varies -- MODEL COMMITS")):
        r = st[name]
        rate = f"{r['deployed_exact'] / r['n_genes']:.0%}" if r["n_genes"] else "  -"
        print(f"    {label:32} n={r['n_genes']:>4}  deployed {r['deployed_exact']:>3} ({rate:>4})  "
              f"oracle {r['oracle_exact']:>3}")
    c = st["commits"]
    if c["n_genes"]:
        print(f"  WHEN THE MODEL COMMITS it is exactly right "
              f"{c['deployed_exact']}/{c['n_genes']} = {c['deployed_exact'] / c['n_genes']:.0%}, "
              f"and it commits on {c['n_genes'] / out['n_genes_scored']:.0%} of genes.")
        print(f"  CHANCE baseline (conditions interchangeable): {out['chance_exact_hits_among_committing']}"
              f" expected exact hits.")
        mp = out["marginal_preserving_null"]
        if mp and mp.get("mean") is not None:
            print(f"  CHANCE baseline (BOTH margins preserved, {mp['n_draws']} curveball draws): "
                  f"mean {mp['mean']}  sd {mp.get('sd')}  p95 {mp.get('p95')}  max {mp['max']}")
    print(f"  flat fraction ........... {out['flat_fraction']}  must-hold reproduced: {must_hold}")
    print(f"  repeats ................. {vals}  spread {spread} vs margin-over-bar {margin} "
          f"-> spread<margin: {agree}")
    print(f"\n  VERDICT: {out['verdict']} against the pre-registered bar "
          f"(>{PRIMARY_BAR} and p<0.05)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
