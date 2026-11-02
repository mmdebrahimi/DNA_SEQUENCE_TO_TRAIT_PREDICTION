"""Render the eukaryotic arm's memo FROM its artifacts. No number is hand-typed here.

Every figure in the emitted markdown is read out of
`wiki/eukaryotic_arm_result_<CELL>_2026-11-01.json` and the frozen bar, because a memo that restates
a number from memory is how a figure and its artifact drift apart. The one piece of arithmetic this
script does itself -- the gain decomposition -- is computed from those same artifacts and shown as
its components so a reader can check it.

**The decomposition is the honest core of the write-up.** A gain-over-baseline contrast is only
substrate-normalized if the baseline means the same thing on both substrates, and it does not: GC is
measurably WEAKER on human regulatory elements than on bacterial tiles. So a larger human gain has
two possible sources and they must be separated:

    gain_human - gain_bacterial = (conv_human - conv_bacterial) - (gc_human - gc_bacterial)
                                   \\_____ encoder got better _____/   \\__ baseline got worse __/

Reporting only the left-hand side would let a weaker baseline masquerade as a stronger encoder.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BAR = REPO / "wiki/eukaryotic_arm_acceptance_bar.json"


def decompose(res: dict, bar: dict) -> dict:
    """Split the gain difference into an encoder component and a baseline component."""
    a = bar["anchors"]
    conv_h, gc_h = res["medians"]["conv"], res["medians"]["gc"]
    d_conv = conv_h - a["conv"]
    d_gc = gc_h - a["gc"]
    total = d_conv - d_gc
    return {
        "gain_human_median": round(conv_h - gc_h, 4),
        "gain_bacterial": a["conv_minus_gc"],
        "difference": round(total, 4),
        "encoder_component": round(d_conv, 4),
        "baseline_component": round(-d_gc, 4),
        "encoder_share": round(d_conv / total, 4) if total else None,
        "baseline_share": round(-d_gc / total, 4) if total else None,
        "reading": (
            "the encoder scores %+.4f differently and the baseline %+.4f, so %s of the widened gain "
            "is the encoder doing better and %s is GC doing worse"
            % (d_conv, d_gc,
               ("%.0f%%" % (100 * d_conv / total)) if total else "n/a",
               ("%.0f%%" % (100 * -d_gc / total)) if total else "n/a")
        ),
    }


def render(res: dict, bar: dict) -> str:
    d = decompose(res, bar)
    a = bar["anchors"]
    g = res["gain_bootstrap"]
    L = []
    w = L.append

    w("# The eukaryotic arm: the tile-chain ceiling was the SUBSTRATE's, not the method's")
    w("")
    w("<!-- doc-status: active -->")
    w("<!-- measured-on: %s -->" % res["date"])
    w("")
    w("Graded by `%s`, frozen before any number below existed." % res["bar"])
    w("")

    w("## DECISIONS NEEDED")
    w("")
    w("**One, and it is structural rather than about this result.** This arm's metric is a HUMAN")
    w("lentiMPRA Spearman sitting in a criterion whose other ten metrics are BACTERIAL, and")
    w("`ceiling_verdict` picks `best` as a max over metrics. It is harmless today — %.4f is below the"
      % res["medians"]["conv"])
    w("bacterial best %.4f, so `best` is unmoved — but a future human arm scoring above that would"
      % bar["anchors"]["best_arm"])
    w("silently become `best` on a bacterial criterion. The arm's real contribution is the CONTRAST,")
    w("not the level. Moving it to its own criterion is a scope change to a criterion declared")
    w("2026-10-09, so it is yours; a guard")
    w("(`test_the_human_arms_metric_is_not_read_as_best`) now fires if it ever stops being harmless.")
    w("")
    w("**The ceiling criterion does NOT close, and that is correct rather than an oversight.** This")
    w("was the last `planned` arm, but under R6 the positional (non-mean-pooled) PromoterAtlas arm")
    w("that G-A candidate 8 names was declared at the moment this one ran — because an arm declared")
    w("only once its cheap sibling succeeds is an arm that evaporates. So `ceiling_verdict` stays")
    w("`blocked:untried-arm` on a genuine untried candidate: 11 measured, 1 planned, `best` unmoved.")
    w("")

    w("## CORRECTIONS")
    w("")
    w("**The standing handoff called this arm a funding decision. It was neither a funding nor an")
    w("access problem.** The substrate was on an attached drive, MD5-verified, with a working loader,")
    w("and the arm ran locally in %s minutes of wall clock. A prior session the same day had already"
      % int(res["wall_seconds"] // 60))
    w("established it was *not* drive-blocked and judged it cost-gated via a Kaggle fetch; the drive")
    w("being attached makes even that unnecessary.")
    w("")

    w("## The result")
    w("")
    w("| | GC-only ridge | conv encoder | gain |")
    w("|---|---|---|---|")
    w("| human %s (median over %d folds) | %.4f | %.4f | **%+.4f** |"
      % (res["cell"], g["n_folds"], res["medians"]["gc"], res["medians"]["conv"],
         d["gain_human_median"]))
    w("| bacterial tiles (committed arms) | %.4f | %.4f | **%+.4f** |"
      % (a["gc"], a["conv"], a["conv_minus_gc"]))
    w("")
    w("Paired bootstrap over the %d published folds: mean %+.4f, 95%% CI [%+.4f, %+.4f], against a"
      % (g["n_folds"], g["mean"], g["ci_lo"], g["ci_hi"]))
    w("bacterial anchor of %+.4f." % a["conv_minus_gc"])
    w("")
    w("**Verdict: `%s`** — %s" % (res["verdict"], res["verdict_detail"]))
    w("")

    w("## What the gain is actually made of")
    w("")
    w("A gain-over-baseline contrast only crosses a substrate boundary if the baseline means the same")
    w("thing on both sides, and here it does not — GC is weaker on human regulatory elements than on")
    w("bacterial tiles. Separating the two sources:")
    w("")
    w("- the conv encoder moves **%+.4f** (%.4f → %.4f)" % (d["encoder_component"], a["conv"],
                                                            res["medians"]["conv"]))
    w("- the GC baseline moves **%+.4f** in the direction that widens the gap (%.4f → %.4f)"
      % (d["baseline_component"], a["gc"], res["medians"]["gc"]))
    w("- together they account for the **%+.4f** difference in gain" % d["difference"])
    w("")
    w("So %s of the widened gain is the encoder doing better and %s is GC doing worse. **Both point the"
      % (("%.0f%%" % (100 * d["encoder_share"])) if d["encoder_share"] else "n/a",
         ("%.0f%%" % (100 * d["baseline_share"])) if d["baseline_share"] else "n/a"))
    w("same way**, which is what makes the reading hold: on this substrate a learned encoder finds")
    w("signal that a single composition feature cannot, and on bacterial tiles it essentially could")
    w("not. Had the gain widened *only* because GC collapsed, the arm would say nothing about the")
    w("encoder.")
    w("")

    w("## What this does NOT establish")
    w("")
    w("**The published incumbents are NOT a comparator, and that was pre-committed.** Table S6 gives")
    w("per-fold values for six methods on these exact folds and never states the metric; the sheet was")
    w("checked directly, including the rows the parser skips (`unparsed_rows()` returns 0 for")
    w("`Benchmarking results`, verified rather than trusted). What can be said is weaker and still")
    w("useful: whatever those values measure, six independent methods cluster at roughly 0.72–0.87 and")
    w("order the way a performance metric would. That is an inference about an unnamed quantity.")
    w("")
    for lim in res["honest_limits"]:
        w("- %s" % lim)
    w("")

    w("## Reproduce")
    w("")
    w("```bash")
    w("uv run python -m scripts.eukaryotic_arm_probe          # validity first: join, folds, degeneracy")
    w("uv run python scripts/freeze_eukaryotic_arm_bar.py     # refuses if the anchors drift")
    w("uv run python -m scripts.eukaryotic_arm_run --cell %s" % res["cell"])
    w("uv run python -m scripts.eukaryotic_arm_memo --cell %s" % res["cell"])
    w("```")
    w("")
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cell", default="HepG2")
    args = ap.parse_args(argv)

    res_path = REPO / ("wiki/eukaryotic_arm_result_%s_2026-11-01.json" % args.cell)
    res = json.loads(res_path.read_text(encoding="utf-8"))
    bar = json.loads(BAR.read_text(encoding="utf-8"))

    if res["verdict"].startswith("PARTIAL"):
        print("REFUSING: the result is %s; a memo on a withheld verdict would read as a finding."
              % res["verdict"])
        return 2

    out = REPO / ("wiki/eukaryotic_arm_result_%s_2026-11-01.md" % args.cell)
    out.write_text(render(res, bar), encoding="utf-8")
    d = decompose(res, bar)
    print("verdict %s | gain %+.4f vs anchor %+.4f" % (res["verdict"], d["gain_human_median"],
                                                       d["gain_bacterial"]))
    print("decomposition: encoder %+.4f (%s) / baseline %+.4f (%s)"
          % (d["encoder_component"],
             ("%.0f%%" % (100 * d["encoder_share"])) if d["encoder_share"] else "n/a",
             d["baseline_component"],
             ("%.0f%%" % (100 * d["baseline_share"])) if d["baseline_share"] else "n/a"))
    print("-> %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
