"""Render the K562 eukaryotic-arm result to a markdown memo.

WHY THIS NEEDS A GUARD AND NOT A STYLE NOTE. The artifact carries `verdict` and
`provisional_verdict_if_complete` SIDE BY SIDE. A provisional verdict is what the arm WOULD have
concluded had the substrate been complete; presenting it as the verdict is the single over-claim
this artifact is unusually exposed to, because both fields are strings of the same shape sitting in
the same object. So the rail is ENFORCED: `assert_provisional_not_headlined` re-reads the rendered
text and REFUSES if the provisional value appears in any heading. A comment asking a future editor
to be careful would not survive the edit it exists to prevent.

THE REFUSAL AND THE LIMITS ARE COPIED VERBATIM, NEVER SUMMARISED. Both exist because someone
decided a comparison could not honestly be made; a paraphrase is a second author's reading of that
decision, and the compression always runs the same direction -- toward sounding more settled. The
`honest_limits` list is additionally CHECKED for completeness after rendering, so a future edit
cannot quietly drop one.

ORDERED BY WHAT THE READER MUST DO. The verdict first, then what the verdict does NOT license,
then the limits, and only then the numbers and the method. A reader who stops after one screen
should have the scope, not the method.

Usage:
    python scripts/eukaryotic_arm_memo.py [--artifact PATH] [--out PATH]
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ARTIFACT = REPO / "wiki/eukaryotic_arm_result_K562_2026-11-01.json"
OUT = REPO / "wiki/eukaryotic_arm_result_K562_2026-11-01.md"

SCHEMA = "eukaryotic-arm-result-v1"

#: Any markdown heading line. Used by the enforced rail, so it must match every heading level.
HEADING_RE = re.compile(r"(?m)^#{1,6}\s+(.*)$")


class MemoRefused(ValueError):
    """The memo cannot be rendered honestly from this artifact."""


def headings(text: str) -> list[str]:
    return [m.group(1).strip() for m in HEADING_RE.finditer(text)]


def assert_provisional_not_headlined(text: str, art: dict) -> None:
    """REFUSE if the provisional verdict reached a heading.

    `None` is the common case and the check is then vacuous ON THIS ARTIFACT -- which is exactly
    why the test suite drives it with a synthetic artifact that carries a value. A guard whose only
    exercise is against a null field has never been tested.
    """
    prov = art.get("provisional_verdict_if_complete")
    if prov is None:
        return
    prov = str(prov).strip()
    if not prov:
        return
    for h in headings(text):
        if prov in h:
            raise MemoRefused(
                "the provisional verdict %r appears in the heading %r. A provisional verdict is "
                "what the arm WOULD have concluded on a complete substrate; headlining it states "
                "a result the arm did not reach." % (prov, h))


def assert_limits_complete(text: str, art: dict) -> None:
    """Every declared limit must be present verbatim, so none can be dropped by a later edit."""
    missing = [lim for lim in (art.get("honest_limits") or []) if lim not in text]
    if missing:
        raise MemoRefused(
            "%d of %d honest limits are absent from the rendered memo: %r"
            % (len(missing), len(art.get("honest_limits") or []), missing))


def _fold_table(per_fold: list) -> list[str]:
    rows = ["| fold | conv | GC | gain | n_train | n_test |",
            "|---|---|---|---|---|---|"]
    for f in per_fold:
        rows.append("| %s | %.4f | %.4f | %+.4f | %s | %s |"
                    % (f.get("fold"), f.get("conv_spearman"), f.get("gc_spearman"),
                       f.get("gain"), f.get("n_train"), f.get("n_test")))
    return rows


def render_memo(art: dict) -> str:
    """The memo text. Pure: no file access, so it is testable on a synthetic artifact."""
    if art.get("schema") != SCHEMA:
        raise MemoRefused("artifact schema is %r, expected %r" % (art.get("schema"), SCHEMA))
    for required in ("verdict", "verdict_detail", "honest_limits"):
        if art.get(required) in (None, ""):
            raise MemoRefused("artifact carries no %r; the memo would state a result it does not "
                              "have" % required)

    verdict = str(art["verdict"])
    cell = art.get("cell", "?")
    b = art.get("gain_bootstrap") or {}
    med = art.get("medians") or {}
    anc = art.get("bacterial_anchor") or {}
    sub = art.get("substrate") or {}
    tr = art.get("target_reparameterization") or {}
    prov = art.get("provisional_verdict_if_complete")

    L: list[str] = []
    L.append("# Eukaryotic arm, %s cell line -- %s" % (cell, verdict))
    L.append("")
    L.append("Artifact: `wiki/%s`  ·  bar: `%s`  ·  date: %s"
             % (ARTIFACT.name, art.get("bar"), art.get("date")))
    L.append("")

    # ---- 1. the verdict
    L.append("## Verdict: %s" % verdict)
    L.append("")
    L.append("%s." % str(art["verdict_detail"]).rstrip("."))
    L.append("")
    L.append("Criterion: %s" % art.get("criterion"))
    L.append("")
    L.append("- Human conv median Spearman **%s**, GC-only baseline **%s**"
             % (med.get("conv"), med.get("gc")))
    L.append("- Gain CI over %s folds: **[%+.4f, %+.4f]**, mean %+.4f"
             % (b.get("n_folds"), b.get("ci_lo", float("nan")), b.get("ci_hi", float("nan")),
                b.get("mean", float("nan"))))
    L.append("- Bacterial anchor gain (the quantity it is compared against): **%s**"
             % anc.get("conv_minus_gc"))
    L.append("")

    # ---- 2. what it does not say
    L.append("## What this verdict does NOT say")
    L.append("")
    if prov is None:
        L.append("- The artifact records **no** provisional verdict, so there is no "
                 "would-have-been result to read here.")
    else:
        L.append("- A provisional verdict is recorded in the artifact's "
                 "`provisional_verdict_if_complete` field. It is **not** this arm's result: it is "
                 "what the arm would have concluded on a complete substrate, and it is "
                 "deliberately not headlined.")
    L.append("- The comparison that crosses the substrate boundary is the **within-substrate "
             "gain**, not raw Spearman. %s" % anc.get("why_this_is_the_anchor", ""))
    L.append("")

    # ---- 3. the refusals, verbatim
    ref = art.get("refused") or {}
    if ref:
        L.append("## Refused, and what was done instead")
        L.append("")
        for name, blk in ref.items():
            L.append("### `%s`" % name)
            L.append("")
            L.append("**Why refused:** %s" % blk.get("why"))
            L.append("")
            L.append("**Done instead:** %s" % blk.get("what_is_done_instead"))
            L.append("")

    # ---- 4. the limits, verbatim and complete
    L.append("## Honest limits")
    L.append("")
    for lim in art["honest_limits"]:
        L.append("- %s" % lim)
    L.append("")

    # ---- 5. substrate and method
    L.append("## Substrate")
    L.append("")
    L.append("- %s rows measured of %s designed (join rate %s)"
             % (sub.get("n_measured"), sub.get("n_designed"), sub.get("join_rate")))
    L.append("- sequence: %s bp -- %s" % (sub.get("sequence_length"), sub.get("sequence_column")))
    L.append("- published folds: %s" % sub.get("n_folds_published"))
    L.append("")
    if tr:
        L.append("## Target reparameterization")
        L.append("")
        L.append("- applied: `%s`" % tr.get("applied"))
        L.append("- computed over: %s" % tr.get("computed_over"))
        L.append("- why: %s" % tr.get("why"))
        L.append("- evaluation unaffected: %s" % tr.get("evaluation_unaffected"))
        L.append("")

    L.append("## Bacterial anchor")
    L.append("")
    L.append("- best arm %s, conv %s, GC %s, replicate ceiling %s"
             % (anc.get("best_arm"), anc.get("conv"), anc.get("gc"), anc.get("replicate_ceiling")))
    L.append("- fraction of ceiling captured by best: %s"
             % anc.get("fraction_of_ceiling_captured_by_best"))
    L.append("- headroom (best minus GC): %s" % anc.get("headroom_best_minus_gc"))
    L.append("- corroborating artifact: `%s`" % anc.get("corroborating_artifact"))
    L.append("")

    if art.get("per_fold"):
        L.append("## Per fold")
        L.append("")
        L.extend(_fold_table(art["per_fold"]))
        L.append("")

    L.append("---")
    L.append("")
    L.append("Wall time %s s. Rendered by `scripts/eukaryotic_arm_memo.py` from the artifact "
             "named above; every number here is read from it, none recomputed."
             % art.get("wall_seconds"))
    L.append("")
    return "\n".join(L)


def build(art: dict) -> str:
    """Render and then CHECK the rendered text against the artifact's own rails."""
    text = render_memo(art)
    assert_provisional_not_headlined(text, art)
    assert_limits_complete(text, art)
    return text


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--artifact", default=str(ARTIFACT))
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)

    p = Path(a.artifact)
    if not p.is_file():
        print("REFUSED: artifact absent: %s" % p)
        return 2
    try:
        text = build(json.loads(p.read_text(encoding="utf-8")))
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
