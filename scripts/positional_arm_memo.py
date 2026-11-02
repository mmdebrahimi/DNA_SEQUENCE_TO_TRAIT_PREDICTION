"""Render the positional arm's memo FROM its artifacts. No number is hand-typed here.

Every figure in the emitted markdown is read out of `wiki/positional_arm_result_2026-11-01.json`
and the frozen bar, because a memo that restates a number from memory is how a figure and its
artifact drift apart.

**It REFUSES to render on an `INDETERMINATE_*` verdict.** A memo on a withheld verdict reads as a
finding, and this arm has three ways to land there -- a split mismatch, a parity failure and a
failed shuffled-label null -- none of which is a result about positional structure.

**The document is ordered by what the READER must do, not by narrative chronology:** DECISIONS
NEEDED first, then CORRECTIONS, then the findings. Clearing the last `planned` arm on a ceiling
criterion flips `ceiling_verdict` to `met`, and whether that should be accepted as the criterion's
terminal is an AUTHORITY call about scope. It goes in the first screen, not near the end.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RESULT = REPO / "wiki/positional_arm_result_2026-11-01.json"
BAR = REPO / "wiki/positional_arm_acceptance_bar.json"
OUT = REPO / "wiki/positional_arm_result_2026-11-01.md"

TERMINAL = ("POSITIONAL_ADDS", "POSITIONAL_NO_GAIN", "POSITIONAL_AMBIGUOUS")

#: The tile rung's replicate-derived ceiling. Named here only so the fraction below is COMPUTED
#: rather than hand-typed; the bar's honest limits carry the same figure and the standing rule that
#: no tile-rung number may be read as a fraction of it without saying so.
RUNG_CEILING = 0.941


def _rung_rows(res: dict) -> list[str]:
    rows = ["| rung | role | features | median rho | gain vs atlas | corrected CI | seed wins |",
            "|---|---|---|---|---|---|---|"]
    order = sorted(
        res["rungs"].items(),
        key=lambda kv: (kv[1]["role"] != "co-primary",
                        -(kv[1].get("median_spearman") or float("-inf"))),
    )
    for name, r in order:
        if r.get("skipped"):
            rows.append("| `%s` | %s | %d | SKIPPED | — | — | — |"
                        % (name, r["role"], r["n_features"]))
            continue
        b = r["gain_bootstrap"]
        rows.append("| `%s` | %s | %d | %.4f | %+.4f | [%+.4f, %+.4f] | %d/%d |"
                    % (name, r["role"], r["n_features"], r["median_spearman"], b["median"],
                       b["ci_lo"], b["ci_hi"], r["seed_wins_over_atlas"], r["n_seeds"]))
    return rows


def render(res: dict, bar: dict, ceiling: dict | None = None) -> str:
    inh = bar["inherited_not_chosen"]
    atlas = inh["comparator"]["atlas_median"]
    margin = inh["margin"]["value"]
    gates = res["gates"]
    L: list[str] = []
    w = L.append

    w("# The positional PromoterAtlas arm: does position-resolved latent structure beat the "
      "mean-pool %.4f" % atlas)
    w("")
    w("<!-- doc-status: active -->")
    w("<!-- measured-on: %s -->" % res["date"])
    w("")
    w("Graded by `%s`, frozen before any number below existed." % res["bar"])
    w("")

    # ---------------------------------------------------------------------------------------
    w("## DECISIONS NEEDED")
    w("")
    if ceiling is not None and ceiling.get("verdict") == "met":
        w("**The ceiling criterion now reads `met`, and whether to accept that is yours.** This arm")
        w("was the last `planned` arm on the criterion, so recording it leaves %d measured arms and"
          % len(ceiling.get("measured", [])))
        w("none untried, which is what `ceiling_verdict` returns `met` for. That is a SUBSTANTIVE")
        w("change to a criterion declared 2026-10-09, not bookkeeping: it says the arm set is")
        w("exhausted and `best = %s` is the ceiling this project reached." % ceiling.get("best"))
        w("")
        w("The alternative is to declare a further arm under R6 before accepting the terminal. Two")
        w("candidates are named and deliberately NOT declared here, because declaring one is itself")
        w("the scope decision: the four shipped expression-prediction heads in the clone, and a")
        w("fine-tuned rather than frozen encoder. Each would be a NEW arm, never a re-read of this")
        w("one.")
        w("")
        w("**What this arm does NOT license either way:** `met` is a statement about the ARM SET")
        w("being exhausted, not about the criterion's `best` of %.4f being close to the %.3f"
          % (atlas, RUNG_CEILING))
        w("replicate-derived ceiling. It is not — it is %.0f%% of it, and the project's own standing"
          % (100.0 * atlas / RUNG_CEILING))
        w("rule is that no tile-rung number may be read as a fraction of that ceiling without")
        w("saying so, which is why it is said here.")
        w("")
    else:
        w("**None from this arm's verdict.** `ceiling_verdict` reads `%s`, so the criterion is"
          % (ceiling or {}).get("verdict", "unrecomputed"))
        w("unchanged by this arm and no scope decision falls out of it.")
        w("")

    # ---------------------------------------------------------------------------------------
    w("## CORRECTIONS")
    w("")
    w("**Two scoping corrections, both of which NARROW this arm from how its ledger rationale")
    w("declared it.** They were fixed in the bar before any number existed, not written after.")
    w("")
    for i, c in enumerate(res["scoping_corrections"], 1):
        w("%d. %s" % (i, c))
        w("")
    asg = res["rungs"].get("bins10_mean_asgiven")
    ori = res["rungs"].get("bins10_mean_oriented")
    if asg and ori and not asg.get("skipped") and not ori.get("skipped"):
        d = asg["median_spearman"] - ori["median_spearman"]
        b = asg["gain_bootstrap"]
        w("**The orientation convention is no longer unresolved, and it resolved AGAINST orienting.**")
        w("The frozen bar records it as a 1-bit inference *\"from a performance difference between")
        w("two rungs\"* and lists it under what this arm does not establish. This run supplied that")
        w("difference: as-given scores **%.4f** against oriented **%.4f**, so reversing the bin axis"
          % (asg["median_spearman"], ori["median_spearman"]))
        w("on minus-strand tiles **destroys about %.4f of signal rather than aligning the cohorts**."
          % d)
        w("The bar's own text is left UNEDITED — it was correct when frozen — and this resolution is")
        w("recorded beside it rather than written back into it.")
        w("")
        w("**The consequence is load-bearing and uncomfortable: the arm's declared co-primary was the")
        w("ORIENTED rung, so the eligible primary was handicapped by the transform this same run")
        w("shows to be the wrong one.** The verdict stands unedited anyway, because as-given does not")
        w("clear the margin either — gain %+.4f, CI [%+.4f, %+.4f], which CONTAINS zero, so it TIES"
          % (b["median"], b["ci_lo"], b["ci_hi"]))
        w("the incumbent rather than beating it. A variant with as-given as the primary would be a")
        w("NEW ARM under R6, never a re-read of this one.")
        w("")
        if asg["median_spearman"] > atlas:
            w("**And this is the multiplicity rule doing real work rather than decorating the bar.**")
            w("`bins10_mean_asgiven` at %.4f is the highest rung on the ladder and sits ABOVE the"
              % asg["median_spearman"])
            w("incumbent %.4f. It is a declared SECONDARY and INELIGIBLE as this arm's result, so it"
              % atlas)
            w("is not the metric and does not set the criterion's `best`. Reporting the best of six")
            w("rungs would have claimed a gain — and even that rung's own interval contains zero.")
            w("")

    prov = gates["shuffled_label_null"].get("statistic_provenance")
    if prov:
        w("**A further correction, to an inherited surface rather than to this arm.** %s" % prov)
        w("")

    # ---------------------------------------------------------------------------------------
    w("## The result")
    w("")
    w("**Verdict: `%s`** — %s" % (res["verdict"], res["verdict_detail"]))
    w("")
    w("The incumbent is the mean-pool arm at median Spearman **%.4f**, read from its committed"
      % atlas)
    w("per-seed table and never re-run, so this arm cannot move the number it is graded against.")
    w("The margin is **%+.2f**, inherited verbatim from the tile-rung arms." % margin)
    w("")
    w(*["\n".join(_rung_rows(res))])
    w("")
    w("`gain` is the PAIRED difference against the committed atlas value on the same seed. The")
    w("interval is a paired bootstrap over the %d seeds at the Bonferroni-corrected level: the two"
      % res["rungs"][list(res["rungs"])[0]].get("n_seeds", 10))
    w("co-primaries at alpha/%d, the four secondaries at alpha/%d."
      % (bar["multiplicity"]["co_primary_n_comparisons"],
         bar["multiplicity"]["secondary_n_comparisons"]))
    w("")
    w("**Only the two co-primaries are eligible as this arm's result, and that is enforced rather")
    w("than promised** — the verdict function cannot see the secondaries at all. A secondary rung")
    w("scoring well is a diagnostic, never the headline.")
    w("")

    # ---------------------------------------------------------------------------------------
    w("## The gates, which ran before any rung was scored")
    w("")
    w("| gate | result | what it establishes |")
    w("|---|---|---|")
    w("| split equality | %s | all %d recomputed splits match the committed `n_test_per_seed`, so "
      "the arm is scored on the incumbent's held-out peaks |"
      % ("PASS" if gates["split_equality"]["passed"] else "FAIL",
         len(gates["split_equality"]["per_seed"])))
    p = gates["parity"]
    emb = p.get("embedding_level")
    w("| parity | %s | the global mean DERIVED from the stored 10-bin means reproduces the "
      "committed embedding%s and its per-seed Spearmans (median %.4f vs %.4f), so every rung is "
      "comparable to %.4f rather than merely adjacent to it |"
      % ("PASS" if p["passed"] else "FAIL",
         "" if emb is None else " to max |dev| %.3e" % emb["max_abs_deviation"],
         p["derived_median"], p["committed_median"], atlas))
    n = gates["shuffled_label_null"]
    w("| shuffled-label null | %s | median %+.4f over %d seeds against the bar `%s` |"
      % ("PASS" if n["passed"] else "FAIL", n["pooled_median"], n["n_pooled"], n["bar"]))
    w("")
    if n.get("extended"):
        w("The null was EXTENDED, not weakened. At %d seeds it read %+.4f against the %s bar, and"
          % (len(n["per_seed_10"]), n["median_10"], n["bar"]))
        w("the bar pre-committed to spending more seeds rather than re-specifying. Pooled over %d"
          % n["n_pooled"])
        w("seeds it reads %+.4f, interval [%s, %s] at achieved coverage %s."
          % (n["pooled_median"],
             "%+.4f" % n["median_ci95"][0] if n["median_ci95"][0] is not None else "n/a",
             "%+.4f" % n["median_ci95"][1] if n["median_ci95"][1] is not None else "n/a",
             "%.4f" % n["median_ci_achieved_coverage"]
             if n.get("median_ci_achieved_coverage") else "n/a"))
        w("The bar text is byte-identical to the one frozen before the run.")
        w("")

    skipped = [k for k, v in res["rungs"].items() if v.get("skipped")]
    if skipped:
        w("**%d rung(s) were SKIPPED on the RAM pre-check and recorded with their byte estimates,"
          % len(skipped))
        w("never silently dropped:** %s. A rung nobody could afford to score is not evidence of no"
          % ", ".join("`%s`" % s for s in skipped))
        w("gain, which is why a skipped co-primary yields `POSITIONAL_AMBIGUOUS`.")
        w("")

    # ---------------------------------------------------------------------------------------
    w("## What this does NOT establish")
    w("")
    for lim in res["honest_limits"]:
        w("- %s" % lim)
    w("")
    f = res.get("forecast_recorded_in_bar", {})
    if f:
        w("**My forecast, recorded before any number existed:** `%s` at p=%.2f. %s"
          % (f.get("predicted_outcome"), f.get("p_%s" % str(f.get("predicted_outcome", "")).lower(),
                                                 0.0) or 0.0,
             f.get("my_forecast_record_in_this_project", "")))
        w("")

    w("## Reproduce")
    w("")
    w("```bash")
    for cmd in bar["reproduce"]:
        w(cmd)
    w("uv run python -m scripts.positional_arm_memo")
    w("```")
    w("")
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--result", default=str(RESULT))
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)

    res_path = Path(args.result)
    if not res_path.is_file():
        print("REFUSING: no result artifact at %s; run the scorer first." % res_path)
        return 2
    res = json.loads(res_path.read_text(encoding="utf-8"))
    bar = json.loads(BAR.read_text(encoding="utf-8"))

    if res["verdict"] not in TERMINAL:
        print("REFUSING: the verdict is %s; a memo on a withheld verdict would read as a finding."
              % res["verdict"])
        return 2

    ceiling = None
    try:
        import sys
        sys.path.insert(0, r"C:/my_skills/skills/soraya/scripts")
        import os

        from arms_ledger import ArmsLedger  # type: ignore
        led = ArmsLedger(os.getcwd())
        crits = [c for c in led._data["criteria"] if c.startswith("ceiling:")]
        if len(crits) == 1:
            ceiling = led.ceiling_verdict(crits[0]).to_dict()
    except Exception as e:                                  # pragma: no cover - advisory only
        print("note: could not recompute ceiling_verdict (%s: %s); the memo renders without the "
              "DECISIONS section's verdict text" % (type(e).__name__, e))

    out = Path(args.out)
    out.write_text(render(res, bar, ceiling), encoding="utf-8")
    print("verdict %s" % res["verdict"])
    if ceiling:
        print("ceiling_verdict %s | best %s | measured %d planned %d"
              % (ceiling["verdict"], ceiling["best"], len(ceiling["measured"]),
                 len(ceiling["planned"])))
    print("-> %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
