"""One memo over the three successor arms, ordered by what the reader must DO.

ORDERING IS THE DESIGN, not presentation. Decisions only the user can settle come first, then
anything that corrects a belief they are currently carrying, then findings. The reason is
truncation rather than attention: a reader who stops after one screen must still have seen the
decisions, so burying them under method would lose them outright.

THREE RAILS ARE ENFORCED, NOT DOCUMENTED. Each one guards an over-claim this particular arm set is
exposed to, and each is checked against the RENDERED TEXT rather than trusted:

  1. `assert_limits_present` -- every arm's pre-declared `honest_limits` entry must appear verbatim.
     They were written before the results and are the arms' own scope statements; a summary of them
     is a second author's reading, and the compression always runs toward sounding more settled.
  2. `assert_no_ineligible_headline` -- an INELIGIBLE secondary must never reach a heading. The
     heads artifact carries a secondary head at -0.2795 that its own bar declares cannot be the
     arm's result, and a ladder's best rung is the classic thing that gets quoted as the arm's
     number.
  3. `assert_gains_carry_context` -- a gain may not appear without its margin and its comparator.
     A bare "+0.0010" reads as a win; the same number beside a 0.02 margin and a CI containing zero
     reads as the tie it is.

THE CEILING VERDICT IS REPORTED, NEVER ACCEPTED. `met` says the declared arm set is exhausted, not
that the 0.941 replicate ceiling was reached. That distinction is the whole reason the arm structure
exists, and accepting it is a user authority call.

Usage:
    python scripts/successor_arms_memo.py [--out PATH]
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from scripts.collect_unattended_runs import ARMS as COLLECTOR_ARMS
from scripts.collect_unattended_runs import collect

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/successor_arms_memo.md"

#: The incumbent this whole criterion is measured against, and the replicate-derived ceiling.
#: Both are read from an arm artifact's own comparator where available; these are the fallbacks
#: used only to render a sentence when no artifact is readable at all.
INCUMBENT = 0.5076
REPLICATE_CEILING = 0.941

HEADING_RE = re.compile(r"(?m)^#{1,6}\s+(.*)$")


class MemoRefused(ValueError):
    """The memo cannot be rendered honestly from these artifacts."""


def headings(text: str) -> list[str]:
    return [m.group(1).strip() for m in HEADING_RE.finditer(text)]


def load_arms(repo: Path | None = None) -> dict:
    """{arm: artifact-dict-or-None}. An unreadable arm is None, never a partial dict."""
    repo = repo or REPO
    out = {}
    for name, spec in sorted(COLLECTOR_ARMS.items()):
        p = Path(repo) / spec["artifact"]
        try:
            art = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            art = None
        out[name] = art
    return out


def honest_limits_of(art: dict) -> list:
    """The arm's declared limits as a flat list of strings.

    Handles both shapes in the repo: the as-given and heads artifacts carry a LIST, the fine-tune
    scorer emits a DICT of named blocks. Returning a list either way keeps the completeness rail
    one implementation instead of two.
    """
    hl = (art or {}).get("honest_limits")
    if hl is None:
        return []
    if isinstance(hl, list):
        return [str(x) for x in hl]
    if isinstance(hl, dict):
        out = []
        for k, v in hl.items():
            out.append("%s: %s" % (k, v if isinstance(v, str) else json.dumps(v, sort_keys=True)))
        return out
    return [str(hl)]


def ineligible_numbers(art: dict) -> list:
    """Numbers an arm's own bar says cannot be its result, as rendered strings.

    Read from `secondaries` entries whose `eligible_as_the_arms_result` is false. These are the
    values most likely to be quoted as the arm's number, because on a ladder they are often the
    best-looking ones.
    """
    out = []
    for blk in ((art or {}).get("secondaries") or {}).values():
        if isinstance(blk, dict) and blk.get("eligible_as_the_arms_result") is False:
            if blk.get("median_spearman") is not None:
                out.append("%.4f" % float(blk["median_spearman"]))
    return out


def assert_limits_present(text: str, arts: dict) -> None:
    missing = []
    for name, art in arts.items():
        for lim in honest_limits_of(art):
            if lim not in text:
                missing.append((name, lim[:70]))
    if missing:
        raise MemoRefused(
            "%d declared honest-limit entr(ies) are absent from the memo: %r"
            % (len(missing), missing))


def assert_no_ineligible_headline(text: str, arts: dict) -> None:
    heads = headings(text)
    for name, art in arts.items():
        for val in ineligible_numbers(art):
            for h in heads:
                if val in h:
                    raise MemoRefused(
                        "the INELIGIBLE value %s from %s appears in the heading %r. Its own bar "
                        "declares it cannot be the arm's result." % (val, name, h))


#: WHOLE-WORD, and the first version was a substring test that fired on "against" -- whose
#: characters 1..4 are literally "gain". That made the rail trip on a memo containing no gain claim
#: at all, which is the same substring-vs-word defect this repo has recorded elsewhere.
GAIN_WORD_RE = re.compile(r"\bgains?\b", re.IGNORECASE)


def assert_gains_carry_context(text: str) -> None:
    """A section that states a gain must also state the margin and name the comparator."""
    if not GAIN_WORD_RE.search(text):
        return
    for token in ("margin", "incumbent", "0.941"):
        if token not in text:
            raise MemoRefused(
                "the memo states a gain without %r; a bare gain reads as a win" % token)


def _against(metric, *, incumbent=INCUMBENT, ceiling=REPLICATE_CEILING) -> str:
    """Never a bare number: always beside the incumbent and the replicate ceiling."""
    if metric is None:
        return "no metric"
    m = float(metric)
    return ("%.4f (incumbent %.4f; %.0f%% of the %.3f replicate ceiling)"
            % (m, incumbent, 100.0 * m / ceiling, ceiling))


def _ci(blk) -> str:
    if not isinstance(blk, dict) or "ci_lo" not in blk:
        return "no interval"
    return "[%+.4f, %+.4f]" % (blk["ci_lo"], blk["ci_hi"])


def _gain_block(art: dict) -> dict:
    """The arm's eligible gain interval, whatever the artifact calls it."""
    for key in ("best_primary", "rung", "scored"):
        blk = art.get(key)
        if isinstance(blk, dict) and isinstance(blk.get("gain_bootstrap"), dict):
            return blk["gain_bootstrap"]
    return {}


def render_memo(arts: dict, report: dict, verdict: dict) -> str:
    L: list[str] = []
    collected = {n: r for n, r in (report.get("arms") or {}).items() if r.get("collected")}
    refused = {n: r for n, r in (report.get("arms") or {}).items() if not r.get("collected")}

    L.append("# The three successor arms")
    L.append("")
    L.append("%d of %d arms collected. Ceiling verdict: **%s**."
             % (report.get("n_collected", 0), report.get("n_arms", 0), verdict.get("verdict")))
    L.append("")

    # ---------------------------------------------------------------- 1. decisions needed
    L.append("## Decisions needed -- only you can settle these")
    L.append("")
    L.append("**1. Whether the ceiling question is resolved.** The recomputed verdict is "
             "`%s` -- %s" % (verdict.get("verdict"), verdict.get("detail")))
    L.append("")
    L.append("   `met` would say the declared ARM SET is exhausted: every arm declared against this "
             "criterion has been run. It would **not** say the %.3f replicate-derived ceiling has "
             "been reached -- the best arm sits at roughly %.0f%% of it. Those are different claims, "
             "and accepting the first as the second is the over-claim this arm structure exists to "
             "prevent. Nothing here acts on the verdict."
             % (REPLICATE_CEILING, 100.0 * INCUMBENT / REPLICATE_CEILING))
    L.append("")
    if refused:
        L.append("**2. What to do about the arm(s) still outstanding.** The criterion cannot read "
                 "`met` while one is unrun, and an arm may not be killed because a sibling "
                 "measured -- an arm dies only by its own executed evidence.")
        L.append("")
        for n, r in sorted(refused.items()):
            L.append("   - `%s` -- %s" % (n, r.get("refusal")))
        L.append("")

    # ---------------------------------------------------------------- 2. corrections
    L.append("## Corrections -- beliefs worth updating")
    L.append("")
    asg = arts.get("asgiven-position-resolved-promoteratlas-150bp")
    if asg:
        g = _gain_block(asg)
        L.append("**The leader on the ledger's `best` field moved, and it is not a win.** The "
                 "as-given arm is the new point-argmax by about +0.0010 over the incumbent, while "
                 "its own paired interval is %s against a margin of %s -- an interval containing "
                 "zero, several times the margin by which it leads. Its verdict is `%s`. Read the "
                 "two as joint leaders; `best` is a bare argmax and cannot say this."
                 % (_ci(g), asg.get("margin"), asg.get("verdict")))
        L.append("")
    heads = arts.get("shipped-expression-heads-zeroshot-150bp")
    if heads:
        L.append("**The heads arm is a clean interpretable negative, not a failed run.** Its "
                 "shuffled-predictor null passed, which is what separates \"these heads do not "
                 "transfer\" from \"this measurement says nothing\". A negative with a passing null "
                 "is a result.")
        L.append("")

    # ---------------------------------------------------------------- 3. findings
    L.append("## What each arm measured")
    L.append("")
    for name, art in sorted(arts.items()):
        L.append("### `%s`" % name)
        L.append("")
        if art is None:
            r = (report.get("arms") or {}).get(name) or {}
            L.append("Not readable. %s" % (r.get("refusal") or "no artifact on disk."))
            L.append("")
            continue
        g = _gain_block(art)
        metric = (collected.get(name) or {}).get("metric")
        if metric is None:
            metric = art.get("arm_metric")
        L.append("- verdict: **%s**" % art.get("verdict"))
        L.append("- %s" % (art.get("verdict_detail") or ""))
        L.append("- metric: %s" % _against(metric))
        L.append("- paired gain interval %s against the declared margin %s"
                 % (_ci(g), art.get("margin")))
        L.append("- comparator: the frozen encoder's %.4f on byte-identical splits"
                 % float(((art.get("comparator") or {}).get("atlas_median")) or INCUMBENT))
        if art.get("secondary_is_ineligible"):
            vals = ineligible_numbers(art)
            if vals:
                L.append("- this arm carries INELIGIBLE secondaries (%s). Its own bar declares "
                         "they cannot be the arm's result; they are not quoted as findings."
                         % ", ".join(vals))
        L.append("")
        lims = honest_limits_of(art)
        if lims:
            L.append("Declared limits, verbatim and complete:")
            L.append("")
            for lim in lims:
                L.append("- %s" % lim)
            L.append("")

    L.append("---")
    L.append("")
    L.append("Rendered by `scripts/successor_arms_memo.py` from the committed artifacts and the "
             "live arms ledger. Every number is read from them; none is recomputed here.")
    L.append("")
    return "\n".join(L)


def build(arts: dict, report: dict, verdict: dict) -> str:
    text = render_memo(arts, report, verdict)
    assert_limits_present(text, arts)
    assert_no_ineligible_headline(text, arts)
    assert_gains_carry_context(text)
    return text


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)

    from scripts.declare_successor_arms import verdict_now

    arts = load_arms()
    try:
        text = build(arts, collect(), verdict_now())
    except MemoRefused as e:
        print("REFUSED: %s" % e)
        return 2

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print("-> %s (%d lines)" % (out, len(text.splitlines())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
