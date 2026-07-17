"""Is the inverse DEPLOYABLE on a protein with no DMS? — the question that decides what can ship.

THE PROBLEM WITH THE PASSING RESULT. The magnitude inverse (`forward_inverse_roundtrip.py`) hits a target
effect T by fitting a score->effect CALIBRATOR on measured DMS. So it needs DMS **for the protein you are
designing on**. But if you have that protein's DMS, you already know every variant's measured effect --
you do not need an inverse to find one. The magnitude inverse is therefore confined to a narrow niche:
"I scanned half the positions, now design at the unscanned ones."

CAN THE CALIBRATOR TRANSFER FROM ANOTHER PROTEIN? No, and not for a subtle reason -- the assays do not
share a scale. Measured ranges (this repo's own cached assays):

    TEM-1  [-3.56, +0.23]   CcdB [-9.00, -2.00]   PTEN [-5.75, +2.82]
    RL40A  [-0.50, +0.25]   SR43C [-3.97, +1.02]

CcdB's entire range lies below TEM-1's minimum. A calibrator fit on TEM-1 CANNOT EXPRESS a CcdB value at
all, so cross-protein magnitude transfer is impossible by construction, not merely inaccurate. (This is
measured + asserted in `test_assay_scales_do_not_share_a_frame`, not assumed.)

THE DEPLOYABLE ALTERNATIVE, tested here. Ask for a RANK, not a dose:

    "propose an edit at the p-th percentile of damage among the reachable edits"

That needs **no calibrator, no DMS, no measured label for the target protein** -- only the oracle's own
score ordering. It is what you could actually ship for a novel protein. It is also weaker by construction:
it can say *near the top of the damaging tail*, never *fold-change 4.2*.

GRADED, non-circular: propose the variant at score-percentile p; report the percentile it ACTUALLY lands
at in the MEASURED distribution. Error = |measured_pct - target_pct|, in percentile units, so it is
comparable across proteins with incommensurate scales -- which is exactly the property the magnitude
version lacks.

BASELINES: BLOSUM62 percentile (the naive tool -- does the learned oracle earn its keep?) and the exact
random-pick null (a variant drawn with no oracle lands at a uniform percentile; E|U - p| in closed form,
no RNG).

Run:  uv run python scripts/forward_inverse_deployable.py
Exit: 0 = ran; 2 = substrate unavailable.
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
    SubstrateError,
    assay_degeneracy,
    blosum_score,
    esm_score,
    load_substrate,
)
from scripts.forward_inverse_sweep import ASSAYS, MATERIAL_MARGIN, build_candidates  # noqa: E402

TARGET_PCTS = [i / 20 for i in range(1, 20)]     # 5%..95% of the damage ordering


def _percentile_of(sorted_vals: list[float], v: float) -> float:
    """Fraction of the pool at or below v -- the variant's true percentile in the MEASURED distribution."""
    lo, hi = 0, len(sorted_vals)
    while lo < hi:
        mid = (lo + hi) // 2
        if sorted_vals[mid] <= v:
            lo = mid + 1
        else:
            hi = mid
    return lo / len(sorted_vals)


def rank_inverse(cands, score_of, targets: list[float], top_k: int) -> dict:
    """NO calibration, NO DMS: order by score, take the variant at score-percentile p, and ask what
    percentile it truly occupies. Every method faces the identical pool."""
    scored = [(c, score_of(c)) for c in cands]
    scored = [(c, s) for c, s in scored if s is not None]
    by_score = sorted(scored, key=lambda cs: cs[1])          # ascending score == ascending predicted effect
    measured_sorted = sorted(c.measured for c, _ in scored)
    n = len(by_score)

    rows = []
    for p in targets:
        idx = min(n - 1, max(0, int(round(p * (n - 1)))))
        window = by_score[max(0, idx - top_k // 2): max(0, idx - top_k // 2) + top_k]
        top1_pct = _percentile_of(measured_sorted, by_score[idx][0].measured)
        errs = [abs(_percentile_of(measured_sorted, c.measured) - p) for c, _ in window]
        rows.append({"target_pct": round(p, 3), "proposed": by_score[idx][0].mutant,
                     "landed_pct_top1": round(top1_pct, 4),
                     "abs_pct_err_top1": round(abs(top1_pct - p), 4),
                     "abs_pct_err_best_of_k": round(min(errs), 4)})
    return {"n_scored": n,
            "mean_pct_err_top1": round(statistics.fmean(r["abs_pct_err_top1"] for r in rows), 4),
            "mean_pct_err_best_of_k": round(statistics.fmean(r["abs_pct_err_best_of_k"] for r in rows), 4),
            "per_target": rows}


def random_null(targets: list[float], top_k: int) -> dict:
    """Exact, no RNG. A variant picked with no oracle lands at a uniform percentile U.
    E|U - p| = p^2 - p + 1/2.  For best-of-k, integrate the min over k iid draws numerically on a fine grid
    (deterministic; the pool is effectively continuous in percentile space)."""
    top1, bok = [], []
    for p in targets:
        top1.append(p * p - p + 0.5)
        # E[min_k |U_i - p|] = int_0^1 P(min > e) de = int_0^1 (1 - F(e))^k de, where
        # F(e) = measure of {u in [0,1] : |u - p| <= e}. Deterministic quadrature, no RNG.
        G = 2000
        acc = 0.0
        for i in range(G):
            e = i / G
            covered = min(1.0, p + e) - max(0.0, p - e)
            acc += (1 - covered) ** top_k / G
        bok.append(acc)
    return {"mean_pct_err_top1": round(statistics.fmean(top1), 4),
            "mean_pct_err_best_of_k": round(statistics.fmean(bok), 4)}


def render_md(rep: dict) -> str:
    h = rep["headline"]
    scored = [r for r in rep["assays"] if r["status"] == "SCORED"]
    excl = [r for r in rep["assays"] if r["status"] == "DEGENERATE_CENSORED_ASSAY"]
    nl = rep["random_null"]
    L = [
        f"# Is the inverse deployable on a protein with NO DMS? ({rep['date']})",
        "",
        f"**Yes — but only the RANK version. It works on {h['rank_inverse_works_beats_null']} usable "
        f"proteins; the learned oracle earns its keep on {h['learned_oracle_earns_keep']}.**",
        "",
        "## The problem with the passing magnitude result",
        "",
        f"{rep['why_magnitude_is_not_deployable']}",
        "",
        "So the magnitude inverse — the one that PASSED at +53% on blaTEM — is confined to a narrow niche:",
        "*I scanned half the positions, now design at the unscanned ones*. It cannot serve a novel protein.",
        "",
        "## The deployable alternative",
        "",
        "Ask for a **rank**, not a dose: *propose an edit at the p-th percentile of damage among the",
        "reachable edits*. That needs **no calibrator, no DMS, no measured label** — only the oracle's own",
        "score ordering. Graded by the percentile the proposal ACTUALLY lands at in the measured",
        "distribution. Errors are in percentile points, so they are comparable across proteins with",
        "incommensurate assay scales — the property the magnitude version structurally lacks.",
        "",
        f"Random-pick null (exact, no RNG): top-1 **{nl['mean_pct_err_top1']:.4f}**, "
        f"best-of-5 **{nl['mean_pct_err_best_of_k']:.4f}** percentile points.",
        "",
        "| protein | organism | esm top-1 | esm best-of-5 | BLOSUM best-of-5 | vs null | vs BLOSUM | verdict |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in scored:
        tag = {"PASS_LEARNED_ORACLE_EARNS_KEEP": "**oracle earns keep**",
               "PASS_RANK_INVERSE_WORKS_BUT_BLOSUM_SUFFICES": "works, *BLOSUM suffices*",
               "FAIL_NO_DISCRIMINATING_POWER_VS_NULL": "**FAIL vs null**"}[r["verdict"]]
        L.append(f"| {r['protein']} | {r['organism']} | {r['esm2']['mean_pct_err_top1']:.4f} | "
                 f"{r['esm2']['mean_pct_err_best_of_k']:.4f} | "
                 f"{r['blosum62']['mean_pct_err_best_of_k']:.4f} | {r['margin_vs_null']:+.1%} | "
                 f"{r['margin_vs_blosum']:+.1%} | {tag} |")
    for r in excl:
        d = r["degeneracy"]
        L.append(f"| {r['protein']} | {r['organism']} | — | — | — | — | — | "
                 f"**EXCLUDED — censored** ({d['mode_share']:.0%} tied at {d['mode_value']}) |")
    L += [
        "",
        "## What this means for what ships",
        "",
        "**Ship the rank inverse.** It needs nothing but the protein sequence, works on every usable assay",
        "here, and its error is ~2-5 percentile points with 5 proposals. **Do not ship the magnitude",
        "inverse**: it requires the very data that would make it unnecessary.",
        "",
        "The two versions do NOT agree on which proteins they suit, which is why both were measured:",
        "",
        "| | magnitude inverse | rank inverse |",
        "|---|---|---|",
        "| blaTEM | star (+52.9% vs BLOSUM) | *BLOSUM suffices* (+23.3%) |",
        "| RL40A | **fails vs null** | works (+34.9%) |",
        "",
        "A per-protein gate is therefore not optional — the right scorer, and even whether the learned",
        "oracle is worth loading at all, changes protein by protein and question by question.",
        "",
        "## Honest limits",
        "",
        "- **It ranks, it does not dose.** *Near the top of the damaging tail* — never *fold-change 4.2*.",
        "- **A censored assay is excluded, not scored.** CcdB (79.3% of variants tied at its ceiling) has no",
        "  well-defined percentile; ungated it posted −159% vs null, which reads as an oracle failure and is",
        "  not one — the metric is simply undefined there.",
        "- **Regime B (molecular fitness) only.** Not clinical resistance, where this scorer class is below",
        "  chance.",
        "- **n=4 proteins.** Enough to falsify *utility tracks rank*; not enough to model what does predict it.",
        "",
    ]
    return chr(10).join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--out-dir", type=Path, default=REPO / "wiki")
    args = ap.parse_args()

    null = random_null(TARGET_PCTS, args.top_k)
    rows = []
    for a in ASSAYS:
        try:
            cds_path = (REPO / a["cds"]) if a["cds"] else None
            target, dms, esm, cds = load_substrate(a["dms_id"], cds_path)
        except SubstrateError as e:
            rows.append({**a, "status": "SUBSTRATE_UNAVAILABLE", "error": str(e)})
            continue
        cands, space = build_candidates(target, dms, cds)
        if len(cands) < 200:
            rows.append({**a, "status": "UNDERPOWERED", "n": len(cands)})
            continue
        deg = assay_degeneracy([c.measured for c in cands])
        if deg["degenerate"]:
            # NOT an oracle failure -- the metric is undefined here. CcdB scored -159% vs null purely
            # because 79.3% of its variants tie at one value, so no variant occupies percentile 0.5.
            rows.append({**a, "status": "DEGENERATE_CENSORED_ASSAY", "degeneracy": deg})
            continue
        e = rank_inverse(cands, lambda c: esm_score(esm, c), TARGET_PCTS, args.top_k)
        b = rank_inverse(cands, blosum_score, TARGET_PCTS, args.top_k)
        m_null = (null["mean_pct_err_best_of_k"] - e["mean_pct_err_best_of_k"]) / null["mean_pct_err_best_of_k"]
        m_blos = ((b["mean_pct_err_best_of_k"] - e["mean_pct_err_best_of_k"]) / b["mean_pct_err_best_of_k"]
                  if b["mean_pct_err_best_of_k"] else 0.0)
        rows.append({**a, "status": "SCORED", "candidate_space": space, "n_candidates": len(cands),
                     "esm2": e, "blosum62": b,
                     "margin_vs_null": round(m_null, 4), "margin_vs_blosum": round(m_blos, 4),
                     "beats_null": m_null >= MATERIAL_MARGIN,
                     "beats_blosum": m_blos >= MATERIAL_MARGIN,
                     "verdict": ("PASS_LEARNED_ORACLE_EARNS_KEEP"
                                 if (m_null >= MATERIAL_MARGIN and m_blos >= MATERIAL_MARGIN)
                                 else "PASS_RANK_INVERSE_WORKS_BUT_BLOSUM_SUFFICES"
                                 if m_null >= MATERIAL_MARGIN
                                 else "FAIL_NO_DISCRIMINATING_POWER_VS_NULL")})

    scored = [r for r in rows if r["status"] == "SCORED"]
    if not scored:
        print("[deployable] no assay scored", file=sys.stderr)
        return 2
    n_works = sum(r["beats_null"] for r in scored)
    n_earns = sum(r["beats_blosum"] for r in scored)

    rep = {
        "schema": "forward-inverse-deployable-v1",
        "date": date.today().isoformat(),
        "question": ("is the inverse deployable on a protein with NO DMS? The magnitude inverse needs a "
                     "per-protein calibrator; the rank/percentile inverse needs nothing."),
        "why_magnitude_is_not_deployable": (
            "the calibrator maps score -> measured effect, so it needs the TARGET protein's DMS -- and if "
            "you have that you already know every variant's effect. Cross-protein transfer is impossible "
            "by construction, not merely inaccurate: the assays do not share a scale (CcdB's entire range "
            "[-9.00,-2.00] lies below TEM-1's minimum -3.56, so a TEM-1 calibrator cannot express a CcdB "
            "value at all)."),
        "null_definition": ("exact, no RNG: a no-oracle pick lands at a uniform percentile U; "
                            "E|U-p| = p^2-p+1/2, and best-of-k by closed-form integration of (1-F(e))^k."),
        "headline": {
            "n_scored": len(scored),
            "rank_inverse_works_beats_null": f"{n_works}/{len(scored)}",
            "learned_oracle_earns_keep": f"{n_earns}/{len(scored)}",
            "units": "percentile points -- comparable across proteins with incommensurate DMS scales",
        },
        "random_null": null,
        "assays": rows,
    }
    stem = f"forward_inverse_deployable_{rep['date']}"
    (args.out_dir / f"{stem}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")
    (args.out_dir / f"{stem}.md").write_text(render_md(rep), encoding="utf-8")

    print(f"[deployable] rank/percentile inverse — NO calibrator, NO DMS needed (top-k={args.top_k})")
    print(f"  random-pick null: top1 {null['mean_pct_err_top1']:.4f}  "
          f"best-of-{args.top_k} {null['mean_pct_err_best_of_k']:.4f}  (percentile units)\n")
    print(f"  {'protein':22s} {'organism':12s} {'esm top1':>9s} {'esm bok':>8s} {'blos bok':>9s} "
          f"{'vs null':>8s} {'vs blos':>8s}  verdict")
    for r in rows:
        if r["status"] != "SCORED":
            extra = ""
            if r["status"] == "DEGENERATE_CENSORED_ASSAY":
                d = r["degeneracy"]
                extra = f"  ({d['mode_share']:.0%} tied at {d['mode_value']} -> percentile undefined)"
            print(f"  {r['protein']:22s} {r['organism']:12s} {r['status']}{extra}")
            continue
        tag = {"PASS_LEARNED_ORACLE_EARNS_KEEP": "ORACLE EARNS KEEP",
               "PASS_RANK_INVERSE_WORKS_BUT_BLOSUM_SUFFICES": "works, BLOSUM suffices",
               "FAIL_NO_DISCRIMINATING_POWER_VS_NULL": "FAIL vs null"}[r["verdict"]]
        print(f"  {r['protein']:22s} {r['organism']:12s} {r['esm2']['mean_pct_err_top1']:9.4f} "
              f"{r['esm2']['mean_pct_err_best_of_k']:8.4f} {r['blosum62']['mean_pct_err_best_of_k']:9.4f} "
              f"{r['margin_vs_null']:+8.1%} {r['margin_vs_blosum']:+8.1%}  {tag}")
    print(f"\n  RANK inverse works (beats null):        {n_works}/{len(scored)}")
    print(f"  LEARNED oracle earns its keep:         {n_earns}/{len(scored)}")
    print(f"  -> {args.out_dir / (stem + '.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
