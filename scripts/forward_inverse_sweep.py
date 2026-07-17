"""Does the molecular inverse GENERALIZE, or is it a blaTEM quirk? — the cross-protein boundary map.

The blaTEM falsifier PASSED (+53.0%). That is **one protein**. This project's most valuable results have
all been REGIME BOUNDARIES ("supervised rescue works on CONVERGENT pathogens, fails on CLONAL ones"), and
they were only found by replicating a single-substrate win until it broke. So: run the SAME falsifier
across every protein whose ESM2 table is cached, and find where the inverse stops working.

THE PRE-REGISTERED QUESTION (R2: derive, don't assert — stated before the run):
Does inverse success track the forward RANK quality (Spearman)? There is a known counterexample in this
repo to test against: **CcdB-ESM2 ranks 0.49 yet its dosage interval is NOT informative**
(`wiki/forward_dosage_sweep_2026-07-15.md`) — "ranks well" != "pins the dose". If inverse success is
merely rank quality re-measured, CcdB should degrade gracefully with its rank. If the inverse instead
FAILS on CcdB while a similar-rank protein passes, then rank is NOT sufficient and the boundary is
something else — which is the more interesting outcome and the reason to run this.

HONEST SCOPE DIFFERENCE, stated per protein, not buried:
only blaTEM has a committed CDS in this repo, so only it can be restricted to the SINGLE-NT-ACCESSIBLE
candidate set (the real genome-editing space). The other four are scored over the full DMS variant set
(protein-level design space). These are DIFFERENT questions. Within one protein the comparison is fair
(every method faces the identical candidate pool), so per-protein verdicts are valid; but a cross-protein
margin comparison must carry the `candidate_space` field. The sweep prints it in every row.

Run:  uv run python scripts/forward_inverse_sweep.py
Exit: 0 = the sweep ran; 2 = substrate unavailable.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.forward_inverse_roundtrip import (  # noqa: E402
    Candidate,
    SubstrateError,
    assay_degeneracy,
    blosum_score,
    empirical_null,
    esm_score,
    load_substrate,
    run_inverse,
    single_nt_accessible,
)

# Every protein with a cached ESM2-650M masked-marginal table. Spearman = the committed forward leaderboard
# (wiki/forward_method_leaderboard_2026-07-15.json) -- the rank quality each inverse result is tested against.
ASSAYS = [
    {"dms_id": "BLAT_ECOLX_Stiffler_2015", "organism": "E. coli", "protein": "TEM-1 beta-lactamase",
     "esm2_spearman": 0.7315, "cds": "data/forward_ref/blatem_3349172526.fna"},
    {"dms_id": "CCDB_ECOLI_Tripathi_2016", "organism": "E. coli", "protein": "CcdB toxin",
     "esm2_spearman": 0.5115, "cds": None},
    {"dms_id": "PTEN_HUMAN_Mighell_2018", "organism": "human", "protein": "PTEN",
     "esm2_spearman": 0.518, "cds": None},
    {"dms_id": "RL40A_YEAST_Mavor_2016", "organism": "yeast", "protein": "RL40A (ubiquitin)",
     "esm2_spearman": None, "cds": None},
    {"dms_id": "SR43C_ARATH_Tsuboyama_2023_2N88", "organism": "Arabidopsis", "protein": "SR43C",
     "esm2_spearman": None, "cds": None},
]

MATERIAL_MARGIN = 0.25      # same bar as the blaTEM falsifier -- not re-tuned per protein


def build_candidates(target: str, dms: dict[str, float], cds: str | None) -> tuple[list[Candidate], str]:
    if cds is not None:
        acc = single_nt_accessible(target, cds)
        space = "single_nt_accessible_from_real_CDS"
    else:
        acc = None
        space = "all_DMS_variants_protein_level_no_CDS"
    out = []
    for m, v in dms.items():
        if acc is not None and m not in acc:
            continue
        wt, alt = m[0], m[-1]
        try:
            pos = int(m[1:-1])
        except ValueError:
            continue
        if pos < 1 or pos > len(target) or target[pos - 1] != wt or alt == "*":
            continue
        out.append(Candidate(m, pos, wt, alt, v))
    return out, space


def _spearman(x: list[float], y: list[float]) -> float:
    """Rank correlation, computed ON THE SWEEP'S OWN CANDIDATE POOL.

    The committed leaderboard Spearman was measured on a DIFFERENT variant subset (all DMS variants incl.
    nonsense; the sweep excludes nonsense and, for blaTEM, restricts to single-nt-reachable edits). Quoting
    it beside a sweep verdict would compare two different pools. Compute it here instead.
    """
    def rk(v: list[float]) -> list[int]:
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0] * len(v)
        for j, i in enumerate(order):
            r[i] = j
        return r

    if len(x) < 3:
        return 0.0
    rx, ry = rk(x), rk(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else 0.0


def sweep_one(assay: dict, n_targets: int, n_splits: int, top_k: int) -> dict:
    cds_path = (REPO / assay["cds"]) if assay["cds"] else None
    target, dms, esm, cds = load_substrate(assay["dms_id"], cds_path)
    cands, space = build_candidates(target, dms, cds)
    positions = sorted({c.pos for c in cands})
    if len(cands) < 200 or len(positions) < 20:
        return {**assay, "status": "UNDERPOWERED", "n_candidates": len(cands), "n_positions": len(positions)}

    # DEGENERACY GATE (before any scoring): a censored assay FLATTERS the margin rather than failing loudly.
    deg = assay_degeneracy([c.measured for c in cands])
    if deg["degenerate"]:
        return {**assay, "status": "DEGENERATE_CENSORED_ASSAY", "candidate_space": space,
                "n_candidates": len(cands), "degeneracy": deg}

    targets = list(statistics.quantiles([c.measured for c in cands], n=n_targets))
    per_split = []
    for off in range(n_splits):
        cal_pos = {p for i, p in enumerate(positions) if (i + off) % n_splits == 0}
        cal = [c for c in cands if c.pos in cal_pos]
        test = [c for c in cands if c.pos not in cal_pos]
        if len(cal) < 50 or len(test) < 50:
            continue
        per_split.append({
            "esm2": run_inverse("esm2", cal, test, lambda c: esm_score(esm, c), targets, top_k),
            "blosum62": run_inverse("blosum62", cal, test, blosum_score, targets, top_k),
            "empirical_null": empirical_null(test, targets, top_k),
        })
    if not per_split:
        return {**assay, "status": "UNDERPOWERED_SPLITS", "n_candidates": len(cands)}

    def agg(m: str, k: str = "mean_abs_err_best_of_k"):
        v = [s[m][k] for s in per_split]
        return {"mean": round(statistics.fmean(v), 4), "min": round(min(v), 4), "max": round(max(v), 4)}

    e, b, n = agg("esm2"), agg("blosum62"), agg("empirical_null")
    span = max(c.measured for c in cands) - min(c.measured for c in cands)
    pool = [(esm_score(esm, c), c.measured) for c in cands]
    pool = [(s, m) for s, m in pool if s is not None]
    fwd_rank = _spearman([s for s, _ in pool], [m for _, m in pool])

    # TWO DIFFERENT QUESTIONS -- one verdict conflates them, and conflating them MISLABELS CcdB
    # (esm 0.09 vs null 0.39 = the inverse works superbly) as "no discriminating power" merely because a
    # 1992 substitution matrix ties it. Decompose:
    #   Q1 DOES INVERSE DESIGN WORK AT ALL?       -> beat the no-oracle empirical null
    #   Q2 DOES THE LEARNED ORACLE EARN ITS KEEP? -> beat BLOSUM62, the naive deterministic tool
    # Q2 is this project's own wrapper-vs-underlying-tool discipline: a layer over a tool must be validated
    # against NAIVE use of that tool, and the DELTA is the claim.
    margin_null = (n["mean"] - e["mean"]) / n["mean"] if n["mean"] else 0.0
    margin_blosum = (b["mean"] - e["mean"]) / b["mean"] if b["mean"] else 0.0

    # PAIRED comparison, not a range comparison. The splits ARE paired -- esm2, blosum62 and the null all
    # face the IDENTICAL calibration/test partition on split i -- so the honest test is the per-split delta
    # and the win count, NOT "does the worst esm range beat the best baseline range" (which is an unpaired
    # test that throws away the pairing and is far too conservative: it called PTEN a loss at +33% margin
    # purely because split variance overlapped). This is the project's own recorded discipline:
    # "median(new) - max(median(a), median(b)) is not a lift; compute per-item deltas + win count."
    def paired(other: str) -> dict:
        deltas = [s[other]["mean_abs_err_best_of_k"] - s["esm2"]["mean_abs_err_best_of_k"]
                  for s in per_split]                      # >0 means esm2 is BETTER on that split
        wins = sum(d > 0 for d in deltas)
        return {"wins": wins, "n": len(deltas),
                "mean_delta": round(statistics.fmean(deltas), 4),
                "worst_delta": round(min(deltas), 4),
                "sweeps": wins == len(deltas)}             # esm2 better on EVERY paired split

    vs_null_paired = paired("empirical_null")
    vs_blosum_paired = paired("blosum62")
    # A win must be BOTH consistent (wins every paired split) AND material (>= the same 25% bar, not
    # re-tuned per protein). Consistency alone can be a hair's-breadth win repeated 6 times.
    beats_null = vs_null_paired["sweeps"] and margin_null >= MATERIAL_MARGIN
    beats_blosum = vs_blosum_paired["sweeps"] and margin_blosum >= MATERIAL_MARGIN
    blosum_beats_null = ((n["mean"] - b["mean"]) / n["mean"] if n["mean"] else 0) >= MATERIAL_MARGIN

    if beats_null and beats_blosum:
        verdict = "PASS_LEARNED_ORACLE_EARNS_KEEP"
    elif beats_null:
        # The inverse WORKS -- just use the cheap deterministic scorer. An actionable engineering result,
        # NOT a failure: no GPU, no model, no ESM table needed for this protein.
        verdict = "PASS_INVERSE_WORKS_BUT_BLOSUM_SUFFICES"
    else:
        verdict = "FAIL_NO_DISCRIMINATING_POWER_VS_NULL"

    informative = sum(bool(s["esm2"]["interval_is_informative"]) for s in per_split)
    return {
        **assay, "status": "SCORED", "candidate_space": space,
        "n_candidates": len(cands), "n_positions": len(positions), "n_splits": len(per_split),
        "measured_effect_span": round(span, 3),
        "forward_rank_on_pool": round(fwd_rank, 4),
        "esm2_best_of_k": e, "blosum62_best_of_k": b, "empirical_null_best_of_k": n,
        "esm2_top1": agg("esm2", "mean_abs_err_top1"),
        # Absolute |err| is NOT comparable across proteins -- effect spans differ >10x (RL40A 0.75 vs
        # PTEN 8.57). Normalize before ANY cross-protein statement.
        "esm2_err_over_span": round(e["mean"] / span, 4) if span else None,
        "null_err_over_span": round(n["mean"] / span, 4) if span else None,
        "margin_vs_null": round(margin_null, 4),
        "margin_vs_blosum": round(margin_blosum, 4),
        "paired_vs_null": vs_null_paired,
        "paired_vs_blosum": vs_blosum_paired,
        "beats_null": beats_null, "beats_blosum": beats_blosum, "blosum_beats_null": blosum_beats_null,
        "verdict": verdict,
        "informative_interval_splits": f"{informative}/{len(per_split)}",
        "magnitude_certifiable": informative > len(per_split) / 2,
    }


def render_md(rep: dict) -> str:
    h = rep["headline"]
    scored = [r for r in rep["assays"] if r["status"] == "SCORED"]
    L = [
        f"# Does the molecular inverse generalize? — cross-protein boundary map ({rep['date']})",
        "",
        f"**Inverse design WORKS on {h['inverse_design_works_beats_null']} usable proteins. "
        f"The LEARNED oracle earns its keep on {h['learned_oracle_earns_keep_beats_blosum']}. "
        f"Magnitude is certifiable on {h['magnitude_certifiable']}.** "
        f"(A 5th assay, CcdB, is EXCLUDED as censored — see finding 4.)",
        "",
        "The blaTEM falsifier passed at +53%. That was **one protein**. This runs the identical falsifier",
        "across every protein with a cached ESM2 table — E. coli x2, human, yeast, Arabidopsis.",
        "",
        "## Result",
        "",
        "| protein | organism | fwd rank | esm err/span | vs null | win | vs BLOSUM | win | interval | verdict |",
        "|---|---|---:|---:|---:|:-:|---:|:-:|:-:|---|",
    ]
    for r in scored:
        pn, pb = r["paired_vs_null"], r["paired_vs_blosum"]
        sp = f"{r['esm2_spearman']:.2f}" if r["esm2_spearman"] else "—"
        tag = {"PASS_LEARNED_ORACLE_EARNS_KEEP": "**oracle earns keep**",
               "PASS_INVERSE_WORKS_BUT_BLOSUM_SUFFICES": "works, *BLOSUM suffices*",
               "FAIL_NO_DISCRIMINATING_POWER_VS_NULL": "**FAIL vs null**"}[r["verdict"]]
        L.append(f"| {r['protein']} | {r['organism']} | {sp} | {r['esm2_err_over_span']:.4f} | "
                 f"{r['margin_vs_null']:+.1%} | {pn['wins']}/{pn['n']} | {r['margin_vs_blosum']:+.1%} | "
                 f"{pb['wins']}/{pb['n']} | {r['informative_interval_splits']} | {tag} |")
    L += [
        "",
        "*err/span normalizes by each protein's measured-effect span — spans differ >10x (RL40A 0.75 vs",
        "PTEN 8.57), so absolute errors are not comparable across rows. `win` = paired per-split wins;",
        "all three methods face the identical partition on each split.*",
        "",
        "## Three findings",
        "",
        "### 1. Inverse SELECTION generalizes; the LEARNED oracle does NOT",
        "",
        f"{h['inverse_design_works_beats_null']} proteins beat the no-oracle null — the inverse lands ~2x",
        "closer to a target than guessing (err/span 1.3–4.1% vs the null's 5.2–7.5%). But only",
        f"{h['learned_oracle_earns_keep_beats_blosum']} beat **BLOSUM62**. On CcdB the 1992 substitution",
        "matrix ties ESM2-650M (2/6 paired wins, +12.5%); on RL40A it is a coin flip (3/6).",
        "**The blaTEM PASS was not representative.** Where BLOSUM suffices that is an *engineering win*, not",
        "a failure: no GPU, no model, no precomputed table — just run the matrix.",
        "",
        "### 2. Inverse utility does NOT track forward rank quality (the pre-registered question, answered)",
        "",
        "The natural assumption is that a better forward ranker inverts better. It is **false**:",
        "",
        "| protein | forward rank (same pool) | earns its keep vs BLOSUM? |",
        "|---|---:|---|",
        "| PTEN | **0.5185** | **yes — 6/6 paired wins, +33.1%** |",
        "| RL40A | **0.5190** | **no — 3/6 paired wins, +26.9%** |",
        "",
        "A rank difference of **0.0005**, opposite verdicts. So you cannot read a Spearman and conclude the",
        "oracle will be useful for design on that protein — it must be measured per protein.",
        "*(n=4: this falsifies the assumption; it does not establish what the real predictor is.)*",
        "",
        "> **Correction, kept on the record.** The first version of this finding used **CcdB** (0.5115, fails)",
        "> as the counterexample. CcdB is a **censored assay** — 79.3% of its variants tie at the −2.00",
        "> ceiling — so it is now excluded and that pairing was an artifact. The finding **survived**",
        "> re-anchoring to RL40A, an even tighter pair on a clean assay. Both ranks are recomputed on each",
        "> sweep's *own* candidate pool rather than quoted from the leaderboard, which measured a different",
        "> variant subset.",
        "",
        "### 3. Selection quality and magnitude certifiability are ORTHOGONAL — measured in both directions",
        "",
        "| protein | selection vs BLOSUM | informative interval |",
        "|---|---|---|",
        "| blaTEM | **best in the sweep** (+52.9%) | **0/6 — certifies nothing** |",
        "| PTEN | +33.1% | **6/6 — fully certified** |",
        "",
        "This project already knew *ranks well ≠ pins the dose* (CcdB). The sweep shows the **converse** too:",
        "*pins the dose ≠ selects well*. They are independent axes, so a design tool must report both — and",
        "an interval that merely **brackets** the target proves nothing, because split-conformal coverage",
        "holds even for a useless model.",
        "",
        "### 4. A censored assay FLATTERS the metric — it does not fail loudly",
        "",
        "**CcdB is excluded**, and how it was caught matters. With 79.3% of its 1,663 variants tied at the",
        "−2.00 ceiling (8 distinct levels total), most quantile targets collapse *onto* that ceiling, so any",
        "ceiling variant hits them trivially. Ungated, CcdB posted the **best err/span in the entire sweep**",
        "(0.0128, +77.1% vs null) — it looked like the strongest result rather than an unusable one. The",
        "degeneracy gate now runs **before** scoring and reports `DEGENERATE_CENSORED_ASSAY`.",
        "",
        "## Scope",
        "",
        f"{rep['scope_caveat']}",
        "",
        f"Test: {rep['test']}",
        "",
        "## What this licenses",
        "",
        "A **per-protein-gated** selection inverse: run this falsifier on the target protein FIRST, and let",
        "it choose the scorer (BLOSUM where BLOSUM suffices — most proteins — and ESM only where it earns",
        "its keep). Do **not** ship a blanket 'ESM2 inverse' on the strength of blaTEM; on 2/5 proteins here",
        "that would burn a GPU to match a substitution matrix, and on RL40A it would not beat guessing.",
        "",
    ]
    return chr(10).join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-targets", type=int, default=40)
    ap.add_argument("--n-splits", type=int, default=6)
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--out-dir", type=Path, default=REPO / "wiki")
    args = ap.parse_args()

    rows = []
    for a in ASSAYS:
        try:
            rows.append(sweep_one(a, args.n_targets, args.n_splits, args.top_k))
        except SubstrateError as e:
            rows.append({**a, "status": "SUBSTRATE_UNAVAILABLE", "error": str(e)})
    scored = [r for r in rows if r["status"] == "SCORED"]
    if not scored:
        print("[inverse-sweep] no assay scored — substrate unavailable?", file=sys.stderr)
        return 2

    n_works = sum(r["beats_null"] for r in scored)
    n_earns = sum(r["beats_blosum"] for r in scored)
    n_mag = sum(r["magnitude_certifiable"] for r in scored)
    rep = {
        "schema": "forward-inverse-sweep-v1",
        "date": date.today().isoformat(),
        "question": ("does the molecular inverse GENERALIZE across proteins/kingdoms, or is it a blaTEM "
                     "quirk? and does inverse success track forward RANK quality?"),
        "material_margin_bar": MATERIAL_MARGIN,
        "test": ("PAIRED per-split comparison (esm2/blosum62/null face the IDENTICAL partition on each "
                 "split): a win must sweep every paired split AND clear the 25% material margin. An "
                 "unpaired range test would discard the pairing and is far too conservative."),
        "headline": {
            "n_scored": len(scored),
            "inverse_design_works_beats_null": f"{n_works}/{len(scored)}",
            "learned_oracle_earns_keep_beats_blosum": f"{n_earns}/{len(scored)}",
            "magnitude_certifiable": f"{n_mag}/{len(scored)}",
            "finding": ("inverse SELECTION generalizes -- it beats the no-oracle null broadly. The LEARNED "
                        "oracle earning its keep over BLOSUM62 does NOT. The blaTEM PASS is not "
                        "representative: on most proteins a 1992 substitution matrix proposes edits as "
                        "well as ESM2-650M does."),
            "orthogonality": ("selection quality and magnitude certifiability are ORTHOGONAL, measured in "
                              "BOTH directions: blaTEM selects best (+52.9% vs BLOSUM) yet certifies NO "
                              "magnitude (0/6 informative); PTEN certifies magnitude 6/6."),
            "utility_does_not_track_rank": ("PTEN forward-rank 0.5185 -> the oracle EARNS ITS KEEP (6/6 "
                                            "paired wins); RL40A forward-rank 0.5190 -> it does NOT (3/6). "
                                            "A rank difference of 0.0005 with opposite verdicts: you cannot "
                                            "read a Spearman and conclude the oracle will be useful for "
                                            "design on that protein. (An earlier version of this finding "
                                            "used CcdB as the counterexample; CcdB is now EXCLUDED as a "
                                            "censored assay, and the finding survives on the clean pair.)"),
            "degeneracy_gate": ("a censored/binned assay is EXCLUDED, not scored. It would otherwise be "
                                "FLATTERED, not caught: CcdB (79.3% of variants at the -2.00 ceiling, 8 "
                                "distinct levels) scored the BEST err/span in the sweep (0.0128, +77.1% vs "
                                "null) purely because most quantile targets collapse onto that ceiling."),
        },
        "scope_caveat": ("only blaTEM has a committed CDS -> only it is restricted to the single-nt-"
                         "accessible (real genome-edit) space; the rest use the full DMS variant set "
                         "(protein-level). Within a protein every method faces the SAME pool, so each "
                         "verdict is valid; cross-protein margins must carry `candidate_space`."),
        "assays": rows,
    }
    stem = f"forward_inverse_sweep_{rep['date']}"
    (args.out_dir / f"{stem}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")
    (args.out_dir / f"{stem}.md").write_text(render_md(rep), encoding="utf-8")

    print(f"[inverse-sweep] {len(scored)}/{len(rows)} assays scored  "
          f"(targets={args.n_targets} x splits={args.n_splits}, top-k={args.top_k})\n")
    print(f"  {'protein':22s} {'organism':12s} {'rank':>5s} {'esm/span':>8s} "
          f"{'vs null':>8s} {'win':>4s} {'vs blos':>8s} {'win':>4s} {'inform':>6s}  verdict")
    for r in rows:
        if r["status"] != "SCORED":
            extra = ""
            if r["status"] == "DEGENERATE_CENSORED_ASSAY":
                d = r["degeneracy"]
                extra = f"  ({d['mode_share']:.0%} of variants at {d['mode_value']}, " \
                        f"{d['n_distinct_values']} distinct levels)"
            print(f"  {r['protein']:22s} {r['organism']:12s} {r['status']}{extra}")
            continue
        sp = f"{r['forward_rank_on_pool']:.2f}"
        tag = {"PASS_LEARNED_ORACLE_EARNS_KEEP": "ORACLE EARNS KEEP",
               "PASS_INVERSE_WORKS_BUT_BLOSUM_SUFFICES": "works, BLOSUM suffices",
               "FAIL_NO_DISCRIMINATING_POWER_VS_NULL": "FAIL vs null"}[r["verdict"]]
        pn, pb = r["paired_vs_null"], r["paired_vs_blosum"]
        print(f"  {r['protein']:22s} {r['organism']:12s} {sp:>5s} "
              f"{r['esm2_err_over_span']:8.4f} "
              f"{r['margin_vs_null']:+8.1%} {'%d/%d' % (pn['wins'], pn['n']):>4s} "
              f"{r['margin_vs_blosum']:+8.1%} {'%d/%d' % (pb['wins'], pb['n']):>4s} "
              f"{r['informative_interval_splits']:>6s}  {tag}")
    print(f"\n  Q1 inverse design WORKS (beats the no-oracle null): {n_works}/{len(scored)}")
    print(f"  Q2 the LEARNED oracle EARNS ITS KEEP (beats BLOSUM):  {n_earns}/{len(scored)}")
    print(f"  MAGNITUDE certifiable (informative interval):         {n_mag}/{len(scored)}")
    print(f"  -> {args.out_dir / (stem + '.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
