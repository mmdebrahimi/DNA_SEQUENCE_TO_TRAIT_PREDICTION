"""Does combining the LEARNED encoder with the BIOPHYSICAL scalar beat the scalar alone?

    uv run python scripts/glm_tile_hybrid_arm.py [--seeds 10]

**The one untried combination on this rung.** The Promoter Calculator was measured AGAINST the conv
encoder (+0.0401 over it, 4/4 seeds) and the 55-feature decomposition was measured against the scalar
(`NO_GAIN`, paired -0.0040, CI95 includes zero). Nobody combined the learned and the biophysical
scorer. That gap is worth one run because this repo has a MEASURED rule about when such a combine pays:

> a naive rank-average of ORTHOGONAL modalities beats the strongest single modality on 84-90% of
> ProteinGym proteins -- `ESM2+GEMME` lifts **even though GEMME alone LOSES to ESM2**.

So "the conv loses to max_tx" is not an argument against trying. **Orthogonality is.**

**And the precondition is measured FIRST, not assumed, because the failure mode is also measured.** On
CTX-M-14 the same rank-average LOST -- ESM2+ProSST 0.204 against ESM2's 0.352 -- because ProSST alone
sat at chance (-0.040), so averaging it in halved the signal. The conv sits at 0.3234 against GC's
0.3112 with a paired CI that includes zero. A negative here is therefore a legitimate, expected-possible
outcome and is recorded as a measured negative (R6/R7), never as a deleted arm.

**No fitting, no calibrator, no second-level model.** Both prediction vectors are already out-of-sample
on each seed's held-out peaks, and `rank_average_hybrid` has no parameters, so the combine introduces
no leakage. It RANKS; it does not dose.

**The control is what makes a positive interpretable.** `max_tx` rank-averaged with a SHUFFLED conv
vector: if that lifts too, the gain came from rank-averaging two vectors rather than from the conv's
content, and the headline is an artifact.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from glm_tile_predictions import (  # noqa: E402
    DEFAULT_OUT,
    load_table,
    reconcile_problems,
    require_real,
)
from promcalc_decomposition_arm import paired  # noqa: E402

from dna_decode.forward.variant_effect import rank_average_hybrid  # noqa: E402
from dna_decode.glm.genomewide import spearman  # noqa: E402

#: FROZEN before any hybrid number was read, MIRRORING the decomposition arm's bar so the two
#: "does more information help?" results are directly comparable.
PREREGISTERED = {
    "question": "does rank-averaging the learned conv encoder with the biophysical max_tx scalar beat "
                "the scalar alone on the tile rung?",
    "primary_metric": "spearman",
    "comparator": "max_tx",
    "margin": 0.02,
    "min_seeds_winning": 8,
    "n_seeds_expected": 10,
    "mirrors": "scripts/promcalc_decomposition_arm.py PREREGISTERED (margin 0.02, 8/10 seeds)",
    "precondition": {
        "name": "orthogonality",
        "rule": "a naive rank-average lifts only when the modalities are ORTHOGONAL; a near-chance "
                "second modality REDUCES the signal (measured: ESM2+ProSST 0.204 vs ESM2 0.352 on "
                "CTX-M-14, because ProSST alone was at chance)",
        "measured_first": True,
        "not_a_gate": "the orthogonality number is REPORTED, not used to skip the run -- a prediction "
                      "that a combine will fail is not a substitute for measuring it",
    },
    "control": "max_tx rank-averaged with a SHUFFLED conv vector; if this lifts, the gain is an "
               "artifact of rank-averaging rather than of the conv's content",
    "verdicts": {
        "HYBRID_ADDS": "paired median >= +0.02 over max_tx AND >= 8/10 seeds win AND control clean",
        "WEAK": "gain > 0 but below the margin or the seed count",
        "NO_GAIN": "no positive paired gain -- the learned encoder carries nothing orthogonal here",
        "INDETERMINATE_CONTROL_ALSO_LIFTS": "the shuffled-conv control lifts too; the combine, not "
                                            "the conv, is producing the number",
        "INDETERMINATE_RECONCILE_FAILED": "the prediction table's published anchors drift",
        "INDETERMINATE_MISSING_SCORER": "conv or max_tx absent from the table",
        "INDETERMINATE_SHORT_ARM": "an arm scored fewer seeds than were loaded, so its median is "
                                   "taken over a different seed set than its comparator's and the "
                                   "two are not comparable",
    },
    #: ADDITIVE second question, frozen BEFORE the first stacking number was read. It is a NEW key,
    #: never an edit to `margin` / `min_seeds_winning` / `comparator` above: those three are pinned by
    #: `test_frozen_bar_mirrors_the_decomposition_arm` and describe the question already ANSWERED
    #: (HYBRID_ADDS, +0.0339). Re-fitting an answered question's bar is the failure being avoided.
    "stacking": {
        "question": "does a THIRD channel add over the MEASURED 2-channel hybrid, rather than over "
                    "the scalar alone?",
        "primary_metric": "spearman",
        "comparator": "hybrid_conv_maxtx",
        "margin": 0.02,
        "min_seeds_winning": 8,
        "mirrors": "the parent bar above (margin 0.02, 8 of 10 seeds), itself mirroring "
                   "scripts/promcalc_decomposition_arm.py PREREGISTERED. Reused VERBATIM rather than "
                   "derived, so the two 'does more information help?' results stay directly "
                   "comparable; the inheritance is recorded as an inheritance, not as a derivation.",
        #: WRITTEN BEFORE THE RUN. An unrecorded prediction is what lets a negative read as a surprise
        #: and a positive read as retro-rationalised.
        "predicted_outcome": "decomposition-as-third is predicted to LOSE, and to lose by MORE than "
                             "gc did. Measured pairwise orthogonality (10-seed median) puts "
                             "decomposition <-> max_tx at +0.7812, against gc <-> max_tx's +0.5508, "
                             "so decomposition is the MORE redundant of the two available third "
                             "channels -- and gc-as-third already lost (0.4217 against the "
                             "2-channel's 0.4254). The expected verdict is therefore NO_STACK for "
                             "both new combos.",
        #: The scope a negative DOES and DOES NOT close. Stated before the result so it cannot be
        #: narrowed afterwards to fit one.
        "negative_scope": "a negative closes THESE FOUR SCORERS on THIS rung, not stacking in "
                          "general. conv and gc are both composition-family; max_tx and "
                          "decomposition are both readings of the same Promoter Calculator output, "
                          "max_tx being a scalar reduction of it. So the four span only TWO modality "
                          "families, and a channel from a NEW family is untested here. A negative is "
                          "evidence about the available channels, never about the mechanism.",
        #: NARROWED 2026-10-26, AFTER the result was known -- which needs saying plainly rather than
        #: hiding in a commit message. The original `WEAK_STACK` text read "gain > 0 but below the
        #: margin or the seed count", and the code implemented it faithfully: `median > 0` alone. So a
        #: gain reproducing on 10 of 10 seeds with a CI excluding zero and positive NOISE on 4 of 10
        #: with a CI spanning zero received the SAME machine-readable verdict. That is the identical
        #: collapse `INDETERMINATE_STACKING_NOT_MEASURED` exists to prevent two keys below -- two
        #: materially different findings sharing one string -- and it was guarded there and not here,
        #: in the same function.
        #:
        #: The narrowing introduces NO NEW THRESHOLD: it reuses `min_seeds_winning`, already frozen at
        #: 8 and already load-bearing for `STACKS`. Inventing a second seed count after seeing the
        #: number would be the post-hoc parameter this whole arm is built to avoid.
        #:
        #: It is legitimate only because it is OUTCOME-NEUTRAL on the measured statistics, and that is
        #: CHECKED, not argued: see `test_the_NARROWED_bar_moved_no_measured_verdict`. If a verdict had
        #: moved, the correct response was to revert this edit, not to defend it.
        "verdicts": {
            "STACKS": "paired median >= +0.02 over hybrid_conv_maxtx AND >= 8/10 seeds win",
            "WEAK_STACK": "gain > 0 AND >= min_seeds_winning seeds win, but below the margin -- a "
                          "REPRODUCIBLE gain that fell short of the effect-size bar",
            "INDETERMINATE_GAIN_NOT_REPRODUCIBLE": "gain > 0 but FEWER than min_seeds_winning seeds "
                                                   "win -- a positive gain that did NOT reproduce, "
                                                   "which is neither a measured stack nor a measured "
                                                   "negative and must not be read as either",
            "NO_STACK": "no positive paired gain -- the third channel carries nothing the measured "
                        "2-channel hybrid lacks",
            "INDETERMINATE_STACKING_NOT_MEASURED": "no superset combo produced a paired statistic -- "
                                                   "a MISSING measurement, which must not be read as "
                                                   "a negative one",
        },
    },
}

COMBOS = {
    #: The two ORIGINAL entries, byte-identical. The published 0.4254 and 0.4217 must stay reproducible.
    "hybrid_conv_maxtx": ("conv", "max_tx"),
    "hybrid_conv_maxtx_gc": ("conv", "max_tx", "gc"),
    #: The 3-channel the orthogonality argument was used to DECLINE. Added because declining a run on
    #: a prediction is what this file's own `precondition.not_a_gate` forbids.
    "hybrid_conv_maxtx_decomp": ("conv", "max_tx", "decomposition"),
    #: The 4-channel. Included because with all four channels on disk it is free, and because a
    #: 3-channel negative ALONE cannot distinguish "a third channel never helps here" from "this
    #: particular third channel does not". Its cost is a second number against one bar -- recorded in
    #: `honest_limits` as a multiple-comparison exposure rather than left implicit.
    "hybrid_conv_maxtx_gc_decomp": ("conv", "max_tx", "gc", "decomposition"),
}

#: The two paired keys that PUBLISHED artifacts and tests already read. Every other combo gets a
#: generated `_paired_<combo>_vs_comparator`. Renaming these two for symmetry would break backward
#: compatibility for a cosmetic gain, so they are pinned here instead.
LEGACY_PAIRED_KEY = {
    "hybrid_conv_maxtx": "_paired_hybrid_vs_comparator",
    "hybrid_conv_maxtx_gc": "_paired_hybrid3_vs_comparator",
}
#: Same reasoning for the human-facing labels in `render_md`.
LEGACY_PAIRED_LABEL = {
    "hybrid_conv_maxtx": "hybrid",
    "hybrid_conv_maxtx_gc": "hybrid+gc",
}


def stack_base() -> tuple[str, ...]:
    """The channels of the combo the stacking bar compares AGAINST — derived, never restated."""
    return COMBOS[PREREGISTERED["stacking"]["comparator"]]


def required_channels() -> tuple[str, ...]:
    """Channels without which the arm's PRIMARY (already-answered) question cannot be asked.

    Deliberately NOT the union over `COMBOS`. A table carrying conv+max_tx can still answer the
    published question, and refusing the whole arm because a LATER-added optional combo wants a
    channel it lacks would make the answered result unreproducible from any older substrate. The
    hazard that extension was reaching for -- a combo silently scored on a seed SUBSET -- is closed
    instead by `skipped_combos` (uniform skip, recorded) plus the `_short_arm` refusal.
    """
    return stack_base()


def all_channels() -> tuple[str, ...]:
    """Every channel named by any COMBOS entry, in first-appearance order."""
    seen: list[str] = []
    for names in COMBOS.values():
        for n in names:
            if n not in seen:
                seen.append(n)
    return tuple(seen)


def channel_pairs() -> list[tuple[str, str]]:
    """Every unordered pair of channels any combo uses — the orthogonality matrix's row set."""
    ch = all_channels()
    return [(a, b) for i, a in enumerate(ch) for b in ch[i + 1:]]


def superset_combos() -> list[str]:
    """Combos that are a STRICT superset of the stacking comparator — the ones that stack."""
    base = set(stack_base())
    return [k for k, names in COMBOS.items() if set(names) > base]


def paired_key(combo: str) -> str:
    return LEGACY_PAIRED_KEY.get(combo, f"_paired_{combo}_vs_comparator")


def stacking_pair_key(combo: str) -> str:
    """The key carrying `combo` measured against the 2-channel hybrid rather than against the scalar.

    This is the statistic the stacking question needs and which did not previously exist in any form.
    """
    return f"_paired_{combo}_vs_hybrid2"


def paired_label(combo: str, pre: dict = PREREGISTERED) -> str:
    """Human-facing label; legacy combos keep their published wording verbatim."""
    if combo in LEGACY_PAIRED_LABEL:
        return LEGACY_PAIRED_LABEL[combo]
    extra = [c for c in COMBOS[combo] if c not in set(stack_base())]
    return "hybrid+" + "+".join(extra) if extra else combo


def skipped_combos(tables: list[dict]) -> dict[str, list[str]]:
    """Combos whose channels are not present in EVERY table, with the missing channels named.

    A combo is run on all tables or on none. The alternative -- scoring whichever seeds happen to
    carry the channel -- produces a median over a seed subset that a reader comparing medians would
    not notice, which is exactly the hazard this guard exists for.
    """
    out: dict[str, list[str]] = {}
    for name, names in COMBOS.items():
        gaps = sorted({n for n in names if any(n not in t["predictions"] for t in tables)})
        if gaps:
            out[name] = gaps
    return out


def as_table(pred: list[float]) -> dict[str, float]:
    """A prediction vector as the `{key: score}` table `rank_average_hybrid` consumes.

    Keys are the POSITION in this seed's test vector as a string. Position, not tile id: every vector
    in one seed's table is already aligned to the same `tile_indices`, and using a positional key makes
    a length mismatch impossible to paper over -- `rank_average_hybrid` raises on an empty intersection.
    """
    return {str(i): float(v) for i, v in enumerate(pred)}


def combine(preds: dict[str, list[float]], names: tuple[str, ...]) -> list[float]:
    """Rank-average `names` into one score vector, in the original index order."""
    tables = [as_table(preds[n]) for n in names]
    combined = rank_average_hybrid(tables)
    n = len(preds[names[0]])
    return [combined[str(i)] for i in range(n)]


def shuffled(vec: list[float], seed: int) -> list[float]:
    """A permutation of `vec` -- the control's second modality. Same distribution, no content."""
    out = list(vec)
    random.Random(1000 + seed).shuffle(out)
    return out


def orthogonality(a: list[float], b: list[float]) -> float:
    """Spearman between two scorers' predictions. Near 1 = redundant; near 0 = orthogonal."""
    return spearman(a, b)


def verdict(arms: dict, pre: dict = PREREGISTERED) -> tuple[str, str]:
    """Mechanical application of the frozen bar, ordered so a refusal wins over a headline."""
    if arms.get("_missing_scorer"):
        return ("INDETERMINATE_MISSING_SCORER",
                f"the table lacks {arms['_missing_scorer']}; nothing to combine")
    if arms.get("_reconcile_problems"):
        return ("INDETERMINATE_RECONCILE_FAILED",
                "; ".join(arms["_reconcile_problems"]))
    if arms.get("_short_arm"):
        return ("INDETERMINATE_SHORT_ARM",
                "; ".join(f"{k} scored {v['scored']} of {v['expected']} seeds"
                          for k, v in sorted(arms["_short_arm"].items()))
                + " -- medians over different seed sets are not comparable")
    p = arms[paired_key(pre["stacking"]["comparator"])]
    ctrl = arms.get("_paired_control_vs_comparator") or {}
    ctrl_lifts = (ctrl.get("median") or 0.0) >= pre["margin"] and \
                 (ctrl.get("wins") or 0) >= pre["min_seeds_winning"]
    if ctrl_lifts:
        return ("INDETERMINATE_CONTROL_ALSO_LIFTS",
                f"the shuffled-conv control gains {ctrl['median']:+.4f} ({ctrl['wins']}/{ctrl['n']}) "
                f"-- indistinguishable from the real combine's {p['median']:+.4f}; the number is an "
                f"artifact of rank-averaging, not of the conv")
    if p["median"] >= pre["margin"] and p["wins"] >= pre["min_seeds_winning"]:
        return ("HYBRID_ADDS",
                f"hybrid {p['median']:+.4f} over {pre['comparator']}, {p['wins']}/{p['n']} seeds, "
                f"CI95 {p['ci95']}")
    if p["median"] > 0:
        return ("WEAK",
                f"hybrid {p['median']:+.4f} over {pre['comparator']}, {p['wins']}/{p['n']} seeds, "
                f"CI95 {p['ci95']} -- positive but short of the {pre['margin']:+.2f} margin or the "
                f"{pre['min_seeds_winning']}-seed bar")
    return ("NO_GAIN",
            f"hybrid {p['median']:+.4f} over {pre['comparator']}, {p['wins']}/{p['n']} seeds, "
            f"CI95 {p['ci95']} -- the learned encoder carries nothing orthogonal the scalar lacks")


def stacking_verdict_for(p: dict, pre: dict = PREREGISTERED) -> tuple[str, str]:
    """Apply the stacking bar to ONE paired statistic (a combo measured against the 2-channel)."""
    s = pre["stacking"]
    if p["median"] >= s["margin"] and p["wins"] >= s["min_seeds_winning"]:
        return ("STACKS", f"{p['median']:+.4f} over {s['comparator']}, {p['wins']}/{p['n']} seeds, "
                          f"CI95 {p['ci95']}")
    # REPRODUCIBLE but under the effect-size bar. The seed condition is what separates this from
    # positive noise; without it both land here and the verdict field stops carrying information.
    if p["median"] > 0 and p["wins"] >= s["min_seeds_winning"]:
        return ("WEAK_STACK",
                f"{p['median']:+.4f} over {s['comparator']}, {p['wins']}/{p['n']} seeds, "
                f"CI95 {p['ci95']} -- REPRODUCIBLE ({p['wins']}/{p['n']} clears the "
                f"{s['min_seeds_winning']}-seed bar) but short of the {s['margin']:+.2f} margin")
    # Positive and NOT reproducible. Neither a stack nor a negative, so it says so rather than
    # borrowing the credibility of a result that won on every seed.
    if p["median"] > 0:
        return ("INDETERMINATE_GAIN_NOT_REPRODUCIBLE",
                f"{p['median']:+.4f} over {s['comparator']}, only {p['wins']}/{p['n']} seeds, "
                f"CI95 {p['ci95']} -- the gain is positive but did NOT reproduce: "
                f"{p['wins']} of {p['n']} seeds win, against the {s['min_seeds_winning']} the bar "
                f"requires; this is neither a measured stack nor a measured negative and must not be "
                f"read as either")
    return ("NO_STACK",
            f"{p['median']:+.4f} over {s['comparator']}, {p['wins']}/{p['n']} seeds, "
            f"CI95 {p['ci95']} -- this third channel carries nothing the measured 2-channel lacks")


def stacking_verdict(arms: dict, pre: dict = PREREGISTERED) -> tuple[str, str, dict]:
    """The stacking question, per superset combo plus one aggregate.

    SEPARATE from `verdict()` on purpose: the top-level field describes the question already ANSWERED
    (does conv add over max_tx) and is read by published artifacts and by `test_verdict_*`. Folding a
    second question into it would change what that field means.

    An ABSENT statistic returns `INDETERMINATE_STACKING_NOT_MEASURED`, never `NO_STACK` -- a missing
    measurement and a measured negative are different findings and this is the seam where they would
    otherwise collapse.
    """
    per: dict[str, dict] = {}
    for combo in superset_combos():
        pst = arms.get(stacking_pair_key(combo))
        if not pst:
            continue
        v, detail = stacking_verdict_for(pst, pre)
        per[combo] = {"verdict": v, "detail": detail}
    if not per:
        return ("INDETERMINATE_STACKING_NOT_MEASURED",
                "no superset combo produced a paired statistic against the 2-channel hybrid; the "
                "stacking question was not measured and must not be read as answered", per)
    # Ordered, best-first. The new key MUST be in this tuple: without it an arm whose every combo is
    # positive-but-unreproducible falls through to the `NO_STACK` return below -- reporting a measured
    # NEGATIVE for a set of positive gains, which re-creates the very defect this change removes, one
    # level up at the aggregate.
    # EVERY combo is named, with its own verdict, on every return path. The aggregate stays best-of
    # by design, but a detail string listing only the combos that cleared presents best-of as if it
    # were all-of -- and `honest_limits` tells the reader to "read the per-combo block, not only the
    # aggregate", which is exactly what the old string defeated. Routed through one helper so the hit
    # path and the fallback cannot drift apart again: the defect was that the hit path joined over
    # `hits` while the fallback below already joined correctly over `sorted(per)`.
    def _all_named() -> str:
        return "; ".join(f"{c} [{per[c]['verdict']}]: {per[c]['detail']}" for c in sorted(per))

    for want in ("STACKS", "WEAK_STACK", "INDETERMINATE_GAIN_NOT_REPRODUCIBLE"):
        if any(r["verdict"] == want for r in per.values()):
            return (want, _all_named(), per)
    return ("NO_STACK", _all_named(), per)


def render_md(art: dict) -> str:
    a = art["arms"]
    L = [f"# Tile-rung hybrid arm (conv ⊕ max_tx) — {art['date']}", "",
         f"**VERDICT: `{art['verdict']}`**", "", art["verdict_detail"], "",
         "## The precondition, measured first", "",
         f"Spearman between the conv and `max_tx` prediction vectors, per seed: "
         f"median **{art['orthogonality']['median']:.4f}** "
         f"(min {art['orthogonality']['min']:.4f}, max {art['orthogonality']['max']:.4f}).", "",
         "Near 1 would mean the two scorers are redundant and a rank-average cannot add; near 0 means "
         "they are orthogonal, which is the condition under which the combine is measured to pay.", "",
         "## Arms", "", "| arm | median spearman | per-seed |", "|---|---|---|"]
    for name in sorted(k for k in a if not k.startswith("_")):
        arm = a[name]
        L.append(f"| `{name}` | {arm['median']:+.4f} | "
                 + ", ".join(f"{v:+.3f}" for v in arm["per_seed"]) + " |")
    L += ["", "## Paired margins (same splits, per seed)", ""]
    cmp_name = PREREGISTERED["comparator"]
    lines = [(paired_key(c), f"{paired_label(c)} vs {cmp_name}") for c in COMBOS]
    lines.append(("_paired_control_vs_comparator",
                  f"CONTROL (shuffled conv) vs {cmp_name}"))
    for key, label in lines:
        pr = a.get(key)
        if pr:
            L.append(f"- **{label}**: median {pr['median']:+.4f}, {pr['wins']}/{pr['n']} seeds, "
                     f"CI95 {pr['ci95']}")

    # The SECOND question. The top-level `verdict` field describes the 2-channel question only, so a
    # reader who reads it alone learns nothing about stacking -- which makes this section
    # load-bearing rather than decorative.
    if art.get("stacking_verdict"):
        base_combo = PREREGISTERED["stacking"]["comparator"]
        L += ["", f"## Does a THIRD channel stack on top of `{base_combo}`?", "",
              f"**STACKING VERDICT: `{art['stacking_verdict']}`**", "",
              art.get("stacking_verdict_detail", ""), ""]
        per = art.get("stacking_per_combo") or {}
        if per:
            L += ["| combo | verdict | margin over the 2-channel |", "|---|---|---|"]
            for combo in sorted(per):
                pr = a.get(stacking_pair_key(combo)) or {}
                marg = (f"{pr['median']:+.4f}, {pr['wins']}/{pr['n']} seeds, CI95 {pr['ci95']}"
                        if pr else "—")
                L.append(f"| `{combo}` | `{per[combo]['verdict']}` | {marg} |")
            L.append("")
        s = art.get("preregistered", {}).get("stacking", {})
        if s:
            L += ["**Predicted BEFORE the run:** " + s.get("predicted_outcome", ""), "",
                  "**What a negative does and does not close:** " + s.get("negative_scope", ""), ""]
        if a.get("_skipped_combos"):
            L += ["**Combos SKIPPED (channel absent from at least one table) — a missing measurement, "
                  "not a negative one:**", ""]
            L += [f"- `{k}`: missing {', '.join(v)}" for k, v in sorted(a["_skipped_combos"].items())]
            L.append("")

    # The measurement that refuted "redundant by construction". Most-redundant first, so the reader
    # meets the stacking obstacle before the pair that worked.
    mx = art.get("orthogonality_matrix") or {}
    if mx:
        L += ["", "## Pairwise orthogonality — every channel pair any combo uses", "",
              "Near 1 = redundant, so a rank-average cannot add. Near 0 = orthogonal, the condition "
              "under which the combine is measured to pay. Ordered most-redundant first.", "",
              "| pair | median | min | max |", "|---|---|---|---|"]
        for key, m in mx.items():
            ca, cb = key.split("__", 1)
            L.append(f"| `{ca}` ↔ `{cb}` | {m['median']:+.4f} | {m['min']:+.4f} | "
                     f"{m['max']:+.4f} |")
        L.append("")

    L += ["", "## Honest limits", ""] + [f"- {x}" for x in art["honest_limits"]]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", default=str(DEFAULT_OUT))
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    in_dir = Path(a.in_dir)
    paths = sorted(in_dir.glob("seed-*.json.gz"))[: a.seeds]
    if not paths:
        raise SystemExit(f"no seed tables in {in_dir}; run scripts/glm_tile_predictions.py first")
    tables = [load_table(p) for p in paths]
    require_real(tables)

    arms: dict = {}
    recon: list[str] = []
    for t in tables:
        recon += reconcile_problems(t["seed"], t["spearman_recomputed"])
    if recon:
        arms["_reconcile_problems"] = recon

    missing = [s for s in required_channels()
               if any(s not in t["predictions"] for t in tables)]
    if missing:
        arms["_missing_scorer"] = missing
        v, detail = verdict(arms)
        print(f"VERDICT: {v}\n  {detail}")
        return 1

    # A combo whose channels are not in EVERY table is skipped UNIFORMLY and RECORDED, rather than
    # scored on whichever seeds happen to carry them. Recorded, because an absent combo and a combo
    # that lost must not produce the same empty output.
    skipped = skipped_combos(tables)
    if skipped:
        arms["_skipped_combos"] = dict(sorted(skipped.items()))
        for k, gaps in sorted(skipped.items()):
            print(f"  SKIPPING combo {k}: absent from at least one table -> {', '.join(gaps)}",
                  flush=True)

    runs: dict[str, list[float]] = {k: [] for k in
                                    (*COMBOS, "max_tx", "conv", "control_shuffled_conv")}
    orth: list[float] = []
    orth_matrix: dict[tuple[str, str], list[float]] = {pr: [] for pr in channel_pairs()}
    base = stack_base()
    for t in tables:
        preds, meas = t["predictions"], t["measured"]
        orth.append(orthogonality(preds["conv"], preds["max_tx"]))
        for ca, cb in channel_pairs():
            if ca in preds and cb in preds:
                orth_matrix[(ca, cb)].append(orthogonality(preds[ca], preds[cb]))
        runs["max_tx"].append(spearman(preds["max_tx"], meas))
        runs["conv"].append(spearman(preds["conv"], meas))
        for name, names in COMBOS.items():
            # `all(n in preds ...)` is belt-and-braces BEHIND `skipped`, not a substitute for it: if
            # `skipped_combos` ever under-reports, this skips the table instead of raising KeyError
            # and the `_short_arm` guard below then REFUSES. Three layers, each catching the previous
            # one's failure -- and it is what makes `_short_arm` reachable by mutation at all, since
            # on real data every table carries every channel.
            if name not in skipped and all(n in preds for n in names):
                runs[name].append(spearman(combine(preds, names), meas))
        ctrl = dict(preds)
        ctrl["conv"] = shuffled(preds["conv"], t["seed"])
        runs["control_shuffled_conv"].append(spearman(combine(ctrl, base), meas))
        bits = " | ".join(f"{n} {runs[n][-1]:+.4f}" for n in COMBOS if runs[n])
        print(f"  seed {t['seed']}: orth {orth[-1]:+.4f} | max_tx {runs['max_tx'][-1]:+.4f} | "
              f"conv {runs['conv'][-1]:+.4f} | {bits} | "
              f"ctrl {runs['control_shuffled_conv'][-1]:+.4f}", flush=True)

    # Any arm that scored a DIFFERENT number of seeds than were loaded. Previously `if vals` accepted
    # 3 of 10 and recorded `"n": 3`, which a reader comparing medians would not notice.
    short = {k: {"scored": len(v), "expected": len(tables)}
             for k, v in sorted(runs.items()) if v and len(v) != len(tables)}
    if short:
        arms["_short_arm"] = short
        v, detail = verdict(arms)
        print(f"VERDICT: {v}\n  {detail}")
        return 1

    for k, vals in runs.items():
        if vals:
            arms[k] = {"median": round(statistics.median(vals), 4), "n": len(vals),
                       "per_seed": [round(v, 4) for v in vals]}
    # One paired statistic per combo, so a combo can never land as a median with no CI that no
    # verdict judges. The two legacy key names are preserved by `paired_key`.
    for name in COMBOS:
        if runs[name]:
            arms[paired_key(name)] = paired(runs[name], runs["max_tx"])
    # The statistic the STACKING question needs: each superset combo against the 2-channel hybrid
    # rather than against the scalar. This did not previously exist in any form.
    base_combo = PREREGISTERED["stacking"]["comparator"]
    for name in superset_combos():
        if runs[name] and runs[base_combo]:
            arms[stacking_pair_key(name)] = paired(runs[name], runs[base_combo])
    arms["_paired_control_vs_comparator"] = paired(runs["control_shuffled_conv"], runs["max_tx"])
    arms["_paired_conv_vs_comparator"] = paired(runs["conv"], runs["max_tx"])

    v, detail = verdict(arms)
    sv, sv_detail, sv_per = stacking_verdict(arms)
    art = {
        "schema": "glm-tile-hybrid-arm-v1",
        "date": str(date.today()),
        "substrate": "GEO GSE144621 peak tiles, leave_peak_out, ceiling 0.941",
        "arm_name": "promcalc-plus-conv-rank-average-hybrid-150bp",
        "n_seeds": len(tables),
        "preregistered": PREREGISTERED,
        "orthogonality": {"per_seed": [round(o, 4) for o in orth],
                          "median": round(statistics.median(orth), 4),
                          "min": round(min(orth), 4), "max": round(max(orth), 4)},
        "orthogonality_matrix": {
            f"{ca}__{cb}": {"per_seed": [round(x, 4) for x in vals],
                            "median": round(statistics.median(vals), 4),
                            "min": round(min(vals), 4), "max": round(max(vals), 4)}
            for (ca, cb), vals in sorted(
                (kv for kv in orth_matrix.items() if kv[1]),
                key=lambda kv: -statistics.median(kv[1]))
        },
        "arms": arms,
        "verdict": v,
        "verdict_detail": detail,
        "stacking_verdict": sv,
        "stacking_verdict_detail": sv_detail,
        "stacking_per_combo": sv_per,
        "honest_limits": [
            "The combine is label-free and parameter-free (a rank-average), so it RANKS and does not "
            "dose -- the same deployability class as the inverse cell.",
            "Both inputs are out-of-sample on each seed's held-out peaks, so no second-level leakage "
            "is introduced; but the pool is MEASURED peak tiles, not generator output.",
            "One rung, one experiment, one combine rule. A negative does not show that no learned "
            "representation can add here -- only that THIS encoder, rank-averaged with THIS scalar, "
            "does not.",
            "`max_tx` and `conv` are both fitted by ridge on the same split, so their errors may be "
            "correlated through the split rather than through the biology.",
            "The STACKING question puts more than one superset combo against ONE pre-registered bar "
            "and the aggregate stacking verdict is the best of them. That is a multiple-comparison "
            "exposure: with k combos the chance one clears the margin by noise exceeds the bar's "
            "nominal rate, so read the per-combo block, not only the aggregate.",
            "Orthogonality here is Spearman between PREDICTION vectors. Two scorers can be "
            "rank-correlated and still carry complementary error structure; a rank-average is only "
            "sensitive to the rank relationship, which is why this is the measure the arm uses -- it "
            "is not a general claim that the channels share information.",
            "SCOPE on the direction finding, and it bounds what this arm showed: `max_tx` and "
            "`decomposition` are NOT two modalities, they are two readings of ONE Promoter "
            "Calculator output -- scripts/glm_tile_predictions.py passes the same decomposition dict "
            "to both `tile_features()` and `max_tx_scalar()`, the latter being a scalar REDUCTION of "
            "it (`log1p(max Tx_rate)`). So their +0.7812 redundancy is mechanically expected, and "
            "the pair that is most redundant while still stacking is a REDUCTION/SUPERSET pair. "
            "Pairwise prediction-vector Spearman is refuted as a screen ACROSS THAT PAIR; it is NOT "
            "shown invalid between genuinely different modality families, which remain untested "
            "here. Read the result as 'finer-grained access to a channel already present pays', not "
            "as 'a new kind of information pays'.",
            "The CI is a normal approximation on the MEAN with z=1.96, while every point estimate "
            "reported beside it is the MEDIAN -- two different estimands -- and at n=10 the correct "
            "critical value is t(9)=2.262, so these intervals are ~15% narrower than their nominal "
            "level. RE-DERIVED under t and the conclusion SURVIVES on all three combos "
            "(decomposition [0.0103, 0.0177] and gc+decomposition [0.0030, 0.0165] still exclude "
            "zero; gc still includes it), with sign tests corroborating independently (p=0.00098 / "
            "0.0107 / 0.828). Nothing published here is wrong. It is NOT fixed because `paired()` "
            "lives in promcalc_decomposition_arm and is SHARED with the published parent arm, so "
            "correcting it moves parent numbers and needs its own pre-registration.",
        ],
    }
    out = Path(a.out) if a.out else Path(f"wiki/glm_tile_hybrid_arm_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(art, indent=1), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(art), encoding="utf-8")
    print(f"\nVERDICT: {v}\n  {detail}")
    print(f"STACKING VERDICT: {sv}\n  {sv_detail}")
    print(f"wrote {out} + {out.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
