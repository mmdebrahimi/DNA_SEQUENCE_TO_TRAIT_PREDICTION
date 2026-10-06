"""The cross-organism transfer-ladder runner — coverage-primary, reconcile-first, every rung reported.

Scores the UNCHANGED conserved-core essentiality decoder at increasing phylogenetic distance from the
organism it was tuned on, to separate (A) continuous decay with distance from (B) a cliff at the eukaryote
boundary. Design + metric live in `dna_decode/eval/transfer_ladder.py`; schema gating in
`dna_decode/essentiality/label_sources.py`; the join in `annotation_join`; the phrasing subtraction in
`phrasing_floor`.

THREE RAILS, each because something already went wrong without it:

  RECONCILE FIRST, AND BLOCKING. Before any new rung is touched, the two rungs that already have
  published numbers are recomputed and must match. This is the gate that caught a silent Goodall
  mis-join (3936/372 against the true 351/3432) during planning; an attribution built on an unreconciled
  prediction set is how a plausible wrong story gets published.

  EVERY RUNG IS REPORTED. A rung that cannot be scored emits a NAMED wall plus the route tried -- never
  dropped from the table. A ladder silently missing its hard rungs would read as a clean short ladder.

  THRESHOLDS FROZEN HERE, FROM THE TWO EXISTING RUNGS, AND HONESTLY LABELLED. They are derived from
  measured permutation noise (below) rather than asserted at design time -- but they are derived from the
  two numbers that MOTIVATED the question, so this is a DESCRIPTIVE-CONSISTENCY LOCK, not a
  frozen-before-the-numbers test of the endpoint hypothesis. Saying otherwise would be the
  pre-registration theatre this repo has twice recorded.

    uv run python scripts/essentiality_transfer_ladder.py --self-check   # offline, no D:, no network
    uv run python scripts/essentiality_transfer_ladder.py                # reconcile + score what exists
    uv run python scripts/essentiality_transfer_ladder.py --emit         # + wiki artifact

Exit 0 = a verdict was produced · 3 = INDETERMINATE (too few rungs, or a null not cleared) · 2 = the
reconcile gate failed, which means nothing downstream is trustworthy.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from dna_decode.essentiality import annotation_join as aj          # noqa: E402
from dna_decode.essentiality import label_sources as ls            # noqa: E402
from dna_decode.essentiality import phrasing_floor as pf           # noqa: E402
from dna_decode.essentiality.core_decoder import score_gene        # noqa: E402
from dna_decode.eval import transfer_ladder as tl                  # noqa: E402

CACHE = Path("D:/dna_decode_cache/essentiality")
WIKI = ROOT / "wiki"
NULL_DRAWS = 400
NULL_SEED = 20261006

# ---------------------------------------------------------------------------------------------------
# RECONCILE TARGETS — the committed numbers. Perturbing any of these must make the gate raise.
# ---------------------------------------------------------------------------------------------------
RECONCILE = {
    "ecoli": {"coverage_lift": 0.3125, "auroc": 0.6952, "n_essential": 351, "n_nonessential": 3432},
    "human": {"coverage_lift": 0.1633, "auroc": 0.5805, "n_essential": 681, "n_nonessential": 899},
}


class ReconcileMismatch(RuntimeError):
    """A committed number did not reproduce. Nothing downstream may be trusted."""


# ---------------------------------------------------------------------------------------------------
# FROZEN THRESHOLDS + the arithmetic that produced them, kept together so the derivation ships with the
# number. Measured 2026-10-06 with NULL_DRAWS=400 at NULL_SEED:
#
#   rung      observed_lift   null_mean   null_p95   null_MAX
#   E. coli      0.3125        -0.0011     0.0298     0.0581
#   human        0.1633        +0.0001     0.0239     0.0446
#
#   plateau_tol = ceil(max null_MAX) = ceil(0.0581) -> 0.06
#       Two rungs closer together than the label-permutation null can resolve are indistinguishable. The
#       null MAX is used rather than p95 deliberately: this repo's recorded practice is to clear the
#       null's MAXIMUM, not its 95th percentile.
#   cliff_drop  = 2 * plateau_tol -> 0.12
#       A drop must be twice the noise floor before it is called a cliff rather than scatter.
#   min_rungs   = 3
#       Two rungs are the existing observation, not a shape.
#
# NOT VACUOUS, checked: the observed E. coli->human spread is 0.1492, which exceeds cliff_drop (so a
# cliff IS reachable) and exceeds plateau_tol (so a plateau is NOT trivially satisfied).
# ---------------------------------------------------------------------------------------------------
FROZEN_THRESHOLDS = {"plateau_tol": 0.06, "cliff_drop": 0.12, "min_rungs": 3}
THRESHOLD_DERIVATION = {
    "measured_nulls": {"ecoli": {"null_p95": 0.0298, "null_max": 0.0581},
                       "human": {"null_p95": 0.0239, "null_max": 0.0446}},
    "plateau_tol": "ceil(max null_MAX across the two existing rungs) = ceil(0.0581) -> 0.06; the null "
                   "MAX rather than p95, matching this repo's practice of clearing the maximum",
    "cliff_drop": "2 x plateau_tol -> 0.12; a drop must be twice the noise floor to be a cliff",
    "min_rungs": "3; two rungs are the motivating observation, not a shape",
    "observed_two_rung_spread": 0.1492,
    "non_vacuity": "0.1492 > cliff_drop 0.12 (a cliff is reachable) and > plateau_tol 0.06 (a plateau is "
                   "not trivially satisfied)",
    "HONEST_LABEL": "DESCRIPTIVE-CONSISTENCY LOCK, not a frozen-before-the-numbers endpoint test: these "
                    "are derived from the two rungs that motivated the question, so a bar built from "
                    "them cannot also test them.",
}


# ---------------------------------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------------------------------
def auroc(scores_pos, scores_neg) -> float:
    """SHARED definition. Mann-Whitney U / (n*m) = P(score_pos > score_neg), ties at 0.5.

    The report card imports THIS so the card and the ladder cannot report two different AUROCs for one
    rung. Reported as the ladder's SECONDARY metric only -- it is not compared across rungs.
    """
    from scipy.stats import mannwhitneyu
    import numpy as np
    a, b = np.asarray(list(scores_pos), float), np.asarray(list(scores_neg), float)
    return float(mannwhitneyu(a, b, alternative="greater").statistic / (len(a) * len(b)))


def permutation_null(all_scores, n_essential, draws=NULL_DRAWS, seed=NULL_SEED) -> dict:
    """Label-shuffle null for coverage_lift. Expected 0 by construction; reports mean, p95 and MAX."""
    rnd = random.Random(seed)
    out = []
    pool = list(all_scores)
    for _ in range(draws):
        rnd.shuffle(pool)
        out.append(tl.coverage_lift(pool[:n_essential], pool[n_essential:]))
    out.sort()
    return {"draws": draws, "seed": seed,
            "null_mean": round(sum(out) / len(out), 4),
            "null_p95": round(out[int(0.95 * len(out))], 4),
            "null_max": round(out[-1], 4)}


# ---------------------------------------------------------------------------------------------------
# rung loaders for the two rungs that EXIST. The three new rungs route through label_sources, which
# refuses until a W0 record verifies them -- so they emit a wall rather than a number.
# ---------------------------------------------------------------------------------------------------
def _have(*names) -> bool:
    return all((CACHE / n).exists() for n in names)


def load_ecoli_rows():
    from essentiality_e3_learned import load_labels
    lab = load_labels()
    ann = aj.load_annotation(str(CACHE / "ecoli_k12_feature_table.txt.gz"), "ncbi_feature_table")
    ess = [k for k, v in lab.items() if v == 1]
    non = [k for k, v in lab.items() if v == 0]
    return ann, ess, non, "gene_symbol"


def load_human_rows():
    from essentiality_missed_vocabulary import human_rows
    ceg, neg = human_rows()
    ann = {g: "%s %s" % (g, d) for g, d in (ceg + neg)}
    return ann, [g for g, _ in ceg], [g for g, _ in neg], "gene_symbol"


def score_rung(rung, ann, ess_ids, non_ids, key_space) -> dict:
    """The per-rung pipeline. Returns a record; NEVER raises on a scoreless rung."""
    rec = {"key": rung.key, "organism": rung.scientific_name, "taxid": rung.taxid,
           "depth": tl.depth_from_reference(rung), "technology": rung.technology,
           "label_source_id": rung.label_source_id, "source_url": rung.source_url,
           "class_sourcing_mode": ls.RUNG_LABEL_PLAN[rung.key]["mode"], "scored": False}
    jr = aj.join_labels(ann, ess_ids, non_ids, key_space=key_space)
    rec["join"] = jr.as_dict()
    if jr.wall:
        rec["wall"] = jr.wall
        return rec

    ess = [(g, t) for g, t, is_e in jr.rows if is_e]
    non = [(g, t) for g, t, is_e in jr.rows if not is_e]
    s_ess = [score_gene(g, t).core_score for g, t in ess]
    s_non = [score_gene(g, t).core_score for g, t in non]

    rec["n_essential"] = len(s_ess)
    rec["n_nonessential"] = len(s_non)
    rec["base_rate"] = round(len(s_ess) / (len(s_ess) + len(s_non)), 4)
    rec["coverage_essential"] = round(tl.coverage(s_ess), 4)
    rec["coverage_nonessential"] = round(tl.coverage(s_non), 4)
    rec["coverage_lift"] = round(tl.coverage_lift(s_ess, s_non), 4)
    fired_e = sum(1 for s in s_ess if s != 0)
    fired_n = sum(1 for s in s_non if s != 0)
    rec["precision_where_fires"] = (round(fired_e / (fired_e + fired_n), 4)
                                    if (fired_e + fired_n) else None)
    rec.update(permutation_null(s_ess + s_non, len(s_ess)))
    raw, adj, detail = pf.coverage_lift_adjusted(ess, non,
                                                 lambda g, t: score_gene(g, t).core_score)
    rec["coverage_lift_adjusted"] = adj
    rec["phrasing_floor"] = {k: detail[k] for k in
                             ("n_phrasing_floor", "phrasing_floor_fraction_of_missed",
                              "coverage_essential_adjusted", "floor_is_a_FLOOR_not_an_estimate")}
    rec["auroc_SECONDARY_not_cross_rung_comparable"] = round(auroc(s_ess, s_non), 4)
    rec["scored"] = True
    return rec


def reconcile_or_raise(records) -> dict:
    """BLOCKING. Every committed number must reproduce before anything else is trusted."""
    checked = {}
    for key, want in RECONCILE.items():
        got = next((r for r in records if r["key"] == key and r.get("scored")), None)
        if got is None:
            raise ReconcileMismatch(
                "rung %r did not score, so its committed numbers cannot be reconciled -- refusing to "
                "proceed on an unreconciled prediction set" % key)
        for field, dp in (("coverage_lift", 4), ("auroc_SECONDARY_not_cross_rung_comparable", 3)):
            w = want["auroc"] if field.startswith("auroc") else want["coverage_lift"]
            g = got[field]
            if round(g, dp) != round(w, dp):
                raise ReconcileMismatch(
                    "%s %s: recomputed %.4f but the committed value is %.4f" % (key, field, g, w))
        for field in ("n_essential", "n_nonessential"):
            if got[field] != want[field]:
                raise ReconcileMismatch(
                    "%s %s: recomputed %d but committed %d" % (key, field, got[field], want[field]))
        checked[key] = {"coverage_lift": got["coverage_lift"],
                        "auroc": got["auroc_SECONDARY_not_cross_rung_comparable"]}
    return checked


def walled_record(rung, reason: str, wall: str) -> dict:
    return {"key": rung.key, "organism": rung.scientific_name, "taxid": rung.taxid,
            "depth": tl.depth_from_reference(rung), "technology": rung.technology,
            "label_source_id": rung.label_source_id, "source_url": rung.source_url,
            "class_sourcing_mode": ls.RUNG_LABEL_PLAN[rung.key]["mode"],
            "scored": False, "wall": wall, "route_tried": reason}


def build_records() -> list[dict]:
    records = []
    for rung in tl.RUNGS:
        if rung.key == "ecoli":
            if not _have("goodall_TableS1_essential.xlsx", "ecoli_k12_feature_table.txt.gz"):
                records.append(walled_record(rung, "Goodall workbook or E. coli feature table absent "
                                             "from %s" % CACHE, aj.WALL_ANNOTATION_UNREACHABLE))
                continue
            records.append(score_rung(rung, *load_ecoli_rows()))
        elif rung.key == "human":
            if not _have("CEGv2.txt", "NEGv1.txt", "Homo_sapiens.gene_info.gz"):
                records.append(walled_record(rung, "BAGEL sets or gene_info absent from %s" % CACHE,
                                             aj.WALL_ANNOTATION_UNREACHABLE))
                continue
            records.append(score_rung(rung, *load_human_rows()))
        else:
            # Routed through the gate on purpose: it refuses until a W0 record verifies the schema, so an
            # unfetched rung produces a NAMED WALL rather than a number from a guessed parser.
            ok, why = ls.w0_status(rung.label_source_id, ls.latest_w0_record())
            records.append(walled_record(rung, why, "WALL_SCHEMA_UNVERIFIED"))
    return records


def _self_check() -> int:
    """Offline: the reconcile gate is live, and a walled rung cannot carry a coverage number."""
    fake = [{"key": "ecoli", "scored": True, "coverage_lift": 0.3125,
             "auroc_SECONDARY_not_cross_rung_comparable": 0.6952,
             "n_essential": 351, "n_nonessential": 3432},
            {"key": "human", "scored": True, "coverage_lift": 0.1633,
             "auroc_SECONDARY_not_cross_rung_comparable": 0.5805,
             "n_essential": 681, "n_nonessential": 899}]
    reconcile_or_raise(fake)                                     # control: matching values pass
    for bad in ({"coverage_lift": 0.3126}, {"auroc_SECONDARY_not_cross_rung_comparable": 0.70},
                {"n_essential": 350}, {"scored": False}):
        perturbed = [dict(fake[0], **bad), fake[1]]
        try:
            reconcile_or_raise(perturbed)
        except ReconcileMismatch:
            pass
        else:
            raise AssertionError("reconcile gate did NOT raise on %s -- it is decorative" % bad)

    w = walled_record(tl.RUNGS_BY_KEY["saureus"], "no W0 record", "WALL_SCHEMA_UNVERIFIED")
    assert "coverage_lift" not in w and "auroc_SECONDARY_not_cross_rung_comparable" not in w
    n = permutation_null([0] * 90 + [1] * 10, 50, draws=40)
    assert abs(n["null_mean"]) < 0.08, n
    print("self-check OK: reconcile gate raises on all 4 perturbations, a walled rung carries no "
          "coverage key, the permutation null centres on ~0")
    return 0


# ---------------------------------------------------------------------------------------------------
# Step 7 — verdict rendering. The five honest limits live IN the template file, not here, so a render
# cannot omit them; a test asserts all five headings are present in the template and in every output.
# ---------------------------------------------------------------------------------------------------
MEMO_TEMPLATE = WIKI / "essentiality_transfer_ladder_memo_template.md"

HONEST_LIMIT_HEADINGS = (
    "### Label technology is confounded with distance across the full ladder",
    "### Cross-rung AUROC is a DIFFERENT SAMPLING FRAME and is refused here",
    "### Roughly 10% of the human miss is a catalogue PHRASING gap, not phylogeny",
    "### A rung is ONE organism",
    "### The bar is an asserted, derived consistency lock — not an endpoint test",
)


def bacterial_subladder_verdict(records, **thresholds):
    """The PRIMARY test: only the technology-matched bacterial rungs."""
    bact = [r for r in records if r["key"] in tl.BACTERIAL_RUNGS]
    return tl.classify_ladder(bact, **thresholds)


def _rung_table(records) -> str:
    head = ("| rung | organism | depth | tech | cov(ess) | cov(non) | **coverage_lift** | adjusted | "
            "null p95 / MAX | state |\n|---|---|---|---|---|---|---|---|---|---|")
    lines = [head]
    for r in sorted(records, key=lambda x: -x["depth"]):
        if r.get("scored"):
            lines.append("| `%s` | *%s* | %d | %s | %.4f | %.4f | **%.4f** | %.4f | %.4f / %.4f | scored |"
                         % (r["key"], r["organism"], r["depth"], r["technology"],
                            r["coverage_essential"], r["coverage_nonessential"], r["coverage_lift"],
                            r["coverage_lift_adjusted"], r["null_p95"], r["null_max"]))
        else:
            lines.append("| `%s` | *%s* | %d | %s | — | — | — | — | — | **%s** |"
                         % (r["key"], r["organism"], r["depth"], r["technology"], r["wall"]))
    return "\n".join(lines)


def render_memo(art: dict) -> str:
    """Fill the template. Post-hoc readings never replace the mechanical verdict (separate key)."""
    tpl = MEMO_TEMPLATE.read_text(encoding="utf-8")
    recs = art["rungs"]
    prim = bacterial_subladder_verdict(recs, **art["frozen_thresholds"])
    classes = "\n".join("- depth **%d**: %s%s" % (
        c["depth"], ", ".join("`%s`" % k for k in c["rungs"]),
        "  ← TIED: one distance class, so their spread is a consistency check, not a rung ordering"
        if len(c["rungs"]) > 1 else "") for c in art["distance_classes"])
    rec = "\n".join("- `%s`: coverage_lift **%.4f**, AUROC %.4f (both reproduced)" %
                    (k, v["coverage_lift"], v["auroc"]) for k, v in art["reconciled"].items())
    deriv = "\n".join("- **%s:** %s" % (k, v) for k, v in art["threshold_derivation"].items()
                      if isinstance(v, (str, int, float)))
    return tpl.format(
        date=art["date"], primary_test=art["primary_test"],
        primary_verdict=prim.verdict, primary_reason=prim.reason,
        secondary_verdict=art["verdict"]["verdict"],
        n_scored=art["n_scored"], n_walled=art["n_walled"],
        frozen=json.dumps(art["frozen_thresholds"]),
        rung_table=_rung_table(recs), distance_classes=classes,
        reconciled=rec, threshold_derivation=deriv)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--self-check", action="store_true")
    ap.add_argument("--emit", action="store_true")
    ap.add_argument("--emit-memo", action="store_true", help="also render the .md memo (Step 7)")
    a = ap.parse_args(argv)
    if a.self_check:
        return _self_check()

    # SMOKE FIRST: one rung before the full loop, so a broken path trips in seconds not after five fetches
    if _have("CEGv2.txt", "NEGv1.txt", "Homo_sapiens.gene_info.gz"):
        smoke = score_rung(tl.RUNGS_BY_KEY["human"], *load_human_rows())
        if not smoke.get("scored"):
            print("SMOKE FAILED on the human rung: %s" % smoke.get("wall"))
            return 2
        print("smoke OK (human rung scored: coverage_lift %.4f)" % smoke["coverage_lift"])

    records = build_records()
    try:
        reconciled = reconcile_or_raise(records)
    except ReconcileMismatch as e:
        print("RECONCILE FAILED -- nothing downstream is trustworthy: %s" % e)
        return 2
    print("reconcile OK: %s" % json.dumps(reconciled))

    v = tl.classify_ladder(records, **FROZEN_THRESHOLDS)
    art = {
        "record": "essentiality-transfer-ladder-v1",
        "date": str(date.today()),
        "question": tl.PREREGISTERED["question"],
        "primary_metric": tl.PREREGISTERED["primary_metric"],
        "primary_test": tl.PRIMARY_TEST,
        "preregistered": tl.PREREGISTERED,
        "frozen_thresholds": FROZEN_THRESHOLDS,
        "threshold_derivation": THRESHOLD_DERIVATION,
        "reconciled": reconciled,
        "distance_classes": [{"depth": d, "rungs": k} for d, k in tl.distance_classes()],
        "rungs": records,
        "verdict": v.as_dict(),
        "n_scored": v.n_scored,
        "n_walled": sum(1 for r in records if not r.get("scored")),
    }
    print("\n%-14s %5s %9s %9s %9s %9s  %s" % ("rung", "depth", "cov(ess)", "cov(non)", "LIFT",
                                               "adj", "state"))
    for r in sorted(records, key=lambda x: -x["depth"]):
        if r.get("scored"):
            print("%-14s %5d %9.4f %9.4f %9.4f %9.4f  scored (null_max %.4f)"
                  % (r["key"], r["depth"], r["coverage_essential"], r["coverage_nonessential"],
                     r["coverage_lift"], r["coverage_lift_adjusted"], r["null_max"]))
        else:
            print("%-14s %5d %9s %9s %9s %9s  %s" % (r["key"], r["depth"], "-", "-", "-", "-",
                                                     r["wall"]))
    print("\nVERDICT: %s\n  %s" % (v.verdict, v.reason))

    if a.emit or a.emit_memo:
        p = WIKI / ("essentiality_transfer_ladder_%s.json" % date.today())
        p.write_text(json.dumps(art, indent=2), encoding="utf-8")
        print("[-> %s]" % p.relative_to(ROOT))
    if a.emit_memo:
        m = WIKI / ("essentiality_transfer_ladder_%s.md" % date.today())
        m.write_text(render_memo(art), encoding="utf-8")
        print("[-> %s]" % m.relative_to(ROOT))

    return 3 if v.verdict.startswith("INDETERMINATE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
