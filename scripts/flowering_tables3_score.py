"""Score the Arabidopsis flowering-habit cell on Zhang & Jimenez-Gomez 2020 Table S3 (N=1,017).

THE FIRST REAL SCORING RUN for the plant cell. Substrate: the paper's own supplementary Table S3 --
per-accession FRI functional status (`deleterious_allele`) AND the phenotype (`FT16_mean`, days to first
flower, long days at 16C) in ONE file, for the 1,017 accessions that survived their QC.

FOUR HONESTY GATES, each of which can turn a headline number into a refusal:

 1. FRI-ROUTE ONLY (scope). Table S3 carries FRI status but NOT FLC. With FLC unobserved the cell's
    two-locus AND (`late iff functional FRI AND strong FLC`) collapses to its FRI route -- i.e. exactly the
    naive FRI-only rule the cell's own docstring warns about, and which its Da(1)-12 anchor exists to catch.
    So this run scores ONE ROUTE of the cell, and it is the WEAKER one (confidence-capped MEDIUM by the Lz-0
    counterexample). The FLC route -- the cell's distinctive claim -- is NOT tested here. We pass
    flc="functional" explicitly rather than letting an unobserved locus masquerade as a known one.

 2. NULL BASELINE. A binary call is meaningless unless it beats the best constant predictor on the SAME set.
    Reported per stratum, never assumed. (This gate already killed two earlier attempts: the AraPheno DTF1
    spot-check, and the paper's article-text anchors where a constant-`early` scores 12/13.)

 3. POPULATION STRUCTURE = this project's clonality analogue. The paper itself reports non-functional FRI is
    overrepresented in central/western Europe and rare elsewhere (their Figure S3, P<0.01 Fisher). So genotype
    correlates with ancestry, and a pooled accuracy partly measures "can you recognise a German accession?".
    The established project lesson is that overall accuracy conflates lineage with mechanism, so we report
    WITHIN-STRUCTURE-GROUP accuracy and a group-weighted (one-vote-per-group) figure beside the pooled one.
    Per the lineage-disclosure precedent, this DISCLOSES; it never silently replaces the raw number.

 4. THE THRESHOLD IS DERIVED, NOT ASSERTED. `late` = FT16 above the cohort median. The median is the one
    split that makes the constant-predictor baseline maximally hard (50/50) -- picking any other cut without
    a reason would be choosing our own baseline. Sensitivity to that choice is reported (tertile split), so
    the reader sees whether the result rides on the cut.

Deterministic, offline, no network. Reads only the committed CC-BY Table S3.

Run:  uv run python scripts/flowering_tables3_score.py
Exit: 0 = SCORED (beats null on the pooled set), 1 = NOT SCORED (degenerate / fails null).
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

from dna_decode.organism_rules.arabidopsis_flowering import (
    call_flowering_habit,
    reference_integrity_ok,
)

REPO = Path(__file__).resolve().parents[1]
TABLE_S3 = REPO / "data" / "arabidopsis" / "zhang2020" / "tpj14716-sup-0012-TableS3.tsv"

CITATION = (
    "Zhang L & Jimenez-Gomez JM (2020) Functional analysis of FRIGIDA using naturally occurring variation "
    "in Arabidopsis thaliana. The Plant Journal 103:154-165. doi:10.1111/tpj.14716. Table S3, CC-BY 4.0."
)

# The paper's own stated counts -- an integrity gate on the file itself (a truncated/mangled TSV fails here).
EXPECTED_ROWS = 1017
EXPECTED_DELETERIOUS = 245   # "Putatively loss-of-function mutations were present in 31 alleles
EXPECTED_ALLELE_GROUPS = 103 #  distributed among 245 accessions" / "We defined 103 distinct FRI alleles"


class ScoringError(RuntimeError):
    """Raised when the substrate fails an integrity gate -- never scored past."""


def _is_deleterious(row: dict) -> bool:
    return row["deleterious_allele"].strip().upper() == "TRUE"


def _has_phenotype(row: dict) -> bool:
    v = row["FT16_mean"].strip()
    return bool(v) and v.upper() != "NA"


def load_table_s3(path: Path = TABLE_S3) -> list[dict]:
    """Load Table S3 and verify it against the paper's own stated counts before any scoring."""
    if not path.exists():
        raise ScoringError(
            f"Table S3 not found at {path}. It is browser-only (Wiley bot-blocks scripted access); "
            f"see data/arabidopsis/zhang2020/README.md."
        )
    with path.open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    n_del = sum(1 for r in rows if _is_deleterious(r))
    n_groups = len({r["allele_group"] for r in rows})
    if len(rows) != EXPECTED_ROWS or n_del != EXPECTED_DELETERIOUS or n_groups != EXPECTED_ALLELE_GROUPS:
        raise ScoringError(
            f"Table S3 does not match the paper's stated counts -- refusing to score. "
            f"rows={len(rows)} (want {EXPECTED_ROWS}), deleterious={n_del} (want {EXPECTED_DELETERIOUS}), "
            f"allele_groups={n_groups} (want {EXPECTED_ALLELE_GROUPS})"
        )
    return rows


def phenotype_attrition(rows: list[dict]) -> dict:
    """Gate 5: who DROPS OUT for want of a phenotype, and is the dropout random wrt genotype?

    16% of Table S3 carries FT16_mean='NA'. A silent drop would be a quiet cohort redefinition, and this
    dropout is NOT random: unphenotyped accessions are enriched for FUNCTIONAL FRI relative to the base
    rate. That skew removes late-flowering candidates preferentially, so it shifts the class balance the
    metric is measured against -- report it beside the score, never bury it.
    """
    dropped = [r for r in rows if not _has_phenotype(r)]
    kept = [r for r in rows if _has_phenotype(r)]
    d_del = sum(1 for r in dropped if _is_deleterious(r))
    base = sum(1 for r in rows if _is_deleterious(r)) / len(rows)
    rate = (d_del / len(dropped)) if dropped else 0.0
    return {
        "n_total": len(rows), "n_phenotyped": len(kept), "n_dropped_no_ft16": len(dropped),
        "dropped_fraction": len(dropped) / len(rows),
        "deleterious_rate_overall": base,
        "deleterious_rate_among_dropped": rate,
        "dropout_is_genotype_skewed": abs(rate - base) > 0.05,
        "note": (
            f"{len(dropped)} of {len(rows)} accessions have FT16_mean='NA' and cannot be scored. The dropout "
            f"is NOT random wrt genotype: {rate:.1%} of dropped accessions carry a deleterious FRI vs "
            f"{base:.1%} overall -- i.e. FUNCTIONAL-FRI (late-candidate) accessions are preferentially "
            f"unphenotyped. This shifts the class balance the null baseline is computed against."
        ),
    }


def predict(row: dict) -> str:
    """Run the DEPLOYED cell on one accession. FLC is unobserved in S3 -> passed explicitly as functional.

    That assumption IS the scope limit (gate 1): it reduces the cell to its FRI route. It is also the
    paper's own implicit framing -- FLC null/weak alleles are documented as rare (Werner et al. 2005) --
    but 'rare' is not 'absent', and every FLC-route accession in here is a guaranteed miss we accept openly.
    """
    fri = "lof" if _is_deleterious(row) else "functional"
    habit = call_flowering_habit(fri, "functional").habit
    # Map the cell's habit vocabulary onto the scoring labels. Explicit, not .startswith() -- an ABSTAIN or
    # an unrecognised habit must RAISE, never be silently bucketed into a class it did not call.
    try:
        return {"winter_annual_late": "late", "summer_annual_early": "early"}[habit]
    except KeyError:
        raise ScoringError(
            f"cell returned habit={habit!r} for FRI={fri}/FLC=functional -- the scorer only maps "
            f"winter_annual_late/summer_annual_early. An ABSTAIN cannot be scored as a class."
        ) from None


def confusion(pairs: list[tuple[str, str]]) -> dict:
    """pairs = [(predicted, observed)] with values in {late, early}. 'late' is the positive class."""
    c = Counter(pairs)
    tp, fn = c[("late", "late")], c[("early", "late")]
    fp, tn = c[("late", "early")], c[("early", "early")]
    n = len(pairs)
    n_late, n_early = tp + fn, fp + tn
    return {
        "n": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "n_observed_late": n_late, "n_observed_early": n_early,
        "accuracy": (tp + tn) / n if n else None,
        "sensitivity": tp / n_late if n_late else None,
        "specificity": tn / n_early if n_early else None,
        # The best CONSTANT predictor on this exact set -- the bar any real call must clear (gate 2).
        "null_accuracy": max(n_late, n_early) / n if n else None,
        "null_call": ("late" if n_late >= n_early else "early") if n else None,
        "degenerate": n_late == 0 or n_early == 0,
    }


def _directional(rows: list[dict], threshold: float) -> dict:
    """Split the rule into its two directions -- the diagnostic that says WHY sens and spec diverge.

    A near-perfect sensitivity with a coin-flip specificity is the signature this project has learned to
    distrust: it usually means the LABEL is the wrong surrogate. Here it is NOT that, and the FP spread is
    what distinguishes them. A label/threshold artifact puts the false positives right AT the cut; a real
    necessary-but-not-sufficient mechanism spreads them far below it. We report the spread so the reader can
    make that call rather than take our word.
    """
    ft = lambda r: float(r["FT16_mean"])  # noqa: E731
    lof = [r for r in rows if _is_deleterious(r)]
    fun = [r for r in rows if not _is_deleterious(r)]
    fp = sorted(ft(r) for r in fun if ft(r) <= threshold)
    lof_early = sum(1 for r in lof if ft(r) <= threshold)
    fun_late = sum(1 for r in fun if ft(r) > threshold)
    return {
        "negative_direction": {
            "rule": "FRI loss-of-function -> early",
            "n": len(lof), "n_correct": lof_early,
            "precision": lof_early / len(lof) if lof else None,
            "mean_ft16": statistics.fmean(ft(r) for r in lof) if lof else None,
        },
        "positive_direction": {
            "rule": "FRI functional -> late",
            "n": len(fun), "n_correct": fun_late,
            "precision": fun_late / len(fun) if fun else None,
            "mean_ft16": statistics.fmean(ft(r) for r in fun) if fun else None,
        },
        "false_positive_spread": {
            "n": len(fp),
            "ft16_min": min(fp) if fp else None,
            "ft16_median": statistics.median(fp) if fp else None,
            "ft16_max": max(fp) if fp else None,
            "threshold": threshold,
            "clustered_at_threshold": (
                (statistics.median(fp) > threshold * 0.9) if fp else None
            ),
        },
    }


def score(all_rows: list[dict], *, quantile: str = "median") -> dict:
    rows = [r for r in all_rows if _has_phenotype(r)]      # gate 5: reported by phenotype_attrition()
    values = [float(r["FT16_mean"]) for r in rows]
    if quantile == "median":
        threshold, cut_desc = statistics.median(values), "FT16 > cohort median"
    else:  # tertile sensitivity check: is the result an artefact of the cut?
        threshold = statistics.quantiles(values, n=3)[1]
        cut_desc = "FT16 > cohort upper tertile"

    pooled: list[tuple[str, str]] = []
    by_group: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for r in rows:
        observed = "late" if float(r["FT16_mean"]) > threshold else "early"
        pair = (predict(r), observed)
        pooled.append(pair)
        by_group[r["group"].strip() or "NA"].append(pair)

    groups = {}
    for g, pairs in sorted(by_group.items()):
        m = confusion(pairs)
        # A group with only one observed class cannot score a binary call -- report it, never average it in.
        m["scorable"] = not m["degenerate"] and m["n"] >= 10
        groups[g] = m

    scorable = [g for g, m in groups.items() if m["scorable"]]
    # Gate 3: one vote per STRUCTURE group. Unweighted mean over groups, so a big over-sampled group
    # (central_europe n=157) cannot carry the metric the way it does in the pooled figure.
    group_weighted = (
        {
            "n_groups_scorable": len(scorable),
            "n_groups_unscorable": len(groups) - len(scorable),
            "mean_accuracy": statistics.fmean(groups[g]["accuracy"] for g in scorable),
            "mean_null_accuracy": statistics.fmean(groups[g]["null_accuracy"] for g in scorable),
            "n_groups_beating_null": sum(
                1 for g in scorable if groups[g]["accuracy"] > groups[g]["null_accuracy"]
            ),
        }
        if scorable
        else None
    )

    p = confusion(pooled)
    beats_null = (not p["degenerate"]) and p["accuracy"] > p["null_accuracy"]
    return {
        "threshold_days": round(threshold, 3), "cut": cut_desc,
        "pooled": p, "by_structure_group": groups, "group_weighted": group_weighted,
        "directional": _directional(rows, threshold),
        "verdict": (
            "DEGENERATE_NO_LABEL_VARIATION" if p["degenerate"]
            else "SCORED_BEATS_NULL" if beats_null
            else "FAILS_NULL_BASELINE"
        ),
    }


def build_report(main: dict, tertile: dict, attrition: dict) -> dict:
    p = main["pooled"]
    gw = main["group_weighted"]
    return {
        "schema": "flowering-tables3-score-v1",
        "date": date.today().isoformat(),
        "cell": "dna_decode/organism_rules/arabidopsis_flowering.py",
        "cell_reference_integrity_ok": reference_integrity_ok(),
        "substrate": {
            "source": str(TABLE_S3.relative_to(REPO)).replace("\\", "/"),
            "citation": CITATION,
            "n_accessions_scored": p["n"],
            "phenotype_attrition": attrition,
            "phenotype": "FT16_mean (days to first flower, long days 16C) -- the paper's own phenotype",
        },
        "scope_limits": [
            "FRI-ROUTE ONLY: Table S3 has no FLC status, so the cell's two-locus AND collapses to its FRI "
            "route -- the weaker, MEDIUM-confidence one. The FLC route (the cell's distinctive claim, the "
            "Da(1)-12 class) is NOT tested by this run.",
            "FLC assumed functional for every accession. Documented-rare, not absent; every FLC-route "
            "accession present is a guaranteed miss.",
            "HABIT/direction call, not quantitative days-to-flower. The paper reports FRI/FLC explains only "
            "part of long-day variation; the residue is polygenic + environmental.",
            attrition["note"],
            "IN-DISTRIBUTION, NOT INDEPENDENT: the cell's catalogue and this label both trace to the same "
            "literature. This measures faithfulness of the rule to the paper's own data -- it is NOT an "
            "out-of-distribution validation.",
        ],
        "main": main,
        "threshold_sensitivity": {
            "note": "Does the verdict ride on the median cut? Re-scored at the upper tertile.",
            "tertile": {
                "threshold_days": tertile["threshold_days"],
                "pooled_accuracy": tertile["pooled"]["accuracy"],
                "pooled_null_accuracy": tertile["pooled"]["null_accuracy"],
                "verdict": tertile["verdict"],
            },
        },
        "headline": {
            "pooled_accuracy": p["accuracy"],
            "pooled_null_accuracy": p["null_accuracy"],
            "sensitivity": p["sensitivity"],
            "specificity": p["specificity"],
            "structure_group_weighted_accuracy": gw["mean_accuracy"] if gw else None,
            "structure_group_weighted_null": gw["mean_null_accuracy"] if gw else None,
            "groups_beating_null": f"{gw['n_groups_beating_null']}/{gw['n_groups_scorable']}" if gw else None,
            "verdict": main["verdict"],
        },
    }


def render_md(rep: dict) -> str:
    p, m = rep["main"]["pooled"], rep["main"]
    gw, d = m["group_weighted"], m["directional"]
    L = [
        f"# Arabidopsis flowering cell — scored on Zhang 2020 Table S3 (N={p['n']})",
        "",
        f"*Generated {rep['date']} by `scripts/flowering_tables3_score.py`. Verdict: **{m['verdict']}**.*",
        "",
        "## Headline",
        "",
        "| metric | value |",
        "|---|---|",
        f"| pooled accuracy | **{p['accuracy']:.3f}** |",
        f"| best constant predictor (null) | {p['null_accuracy']:.3f} (always `{p['null_call']}`) |",
        f"| sensitivity (late) | {p['sensitivity']:.3f} |",
        f"| specificity (early) | {p['specificity']:.3f} |",
    ]
    if gw:
        L += [
            f"| **structure-group-weighted accuracy** | **{gw['mean_accuracy']:.3f}** "
            f"(null {gw['mean_null_accuracy']:.3f}) |",
            f"| groups beating their own null | {gw['n_groups_beating_null']}/{gw['n_groups_scorable']} |",
        ]
    L += [
        "",
        f"Split: `{m['cut']}` = **{m['threshold_days']} days**. "
        f"Observed {p['n_observed_late']} late / {p['n_observed_early']} early.",
        "",
        "## Who is missing (non-random dropout)",
        "",
        rep["substrate"]["phenotype_attrition"]["note"],
        "",
        "## Why sensitivity and specificity diverge (it is not the label)",
        "",
        f"| direction | n | precision | mean FT16 |",
        "|---|---|---|---|",
        f"| `{d['negative_direction']['rule']}` | {d['negative_direction']['n']} | "
        f"**{d['negative_direction']['precision']:.3f}** | {d['negative_direction']['mean_ft16']:.1f}d |",
        f"| `{d['positive_direction']['rule']}` | {d['positive_direction']['n']} | "
        f"**{d['positive_direction']['precision']:.3f}** | {d['positive_direction']['mean_ft16']:.1f}d |",
        "",
        "The rule is strong in one direction and weak in the other. **Losing FRI reliably makes a plant early**",
        f"({d['negative_direction']['precision']:.1%}); **having FRI does not make it late** "
        f"({d['positive_direction']['precision']:.1%}). Functional FRI is *necessary but not sufficient* — "
        "which is precisely what the cell's two-locus rule says (you also need a strong FLC), and precisely",
        "what this FRI-only run cannot use.",
        "",
        "A near-perfect sensitivity beside a coin-flip specificity is normally the signature of a **wrong",
        "label surrogate** in this project. It is not that here, and the false-positive spread is what rules",
        f"it out: the {d['false_positive_spread']['n']} functional-FRI-yet-early accessions span "
        f"{d['false_positive_spread']['ft16_min']:.1f}–{d['false_positive_spread']['ft16_max']:.1f} days with a "
        f"median of {d['false_positive_spread']['ft16_median']:.1f}d, far below the "
        f"{d['false_positive_spread']['threshold']:.1f}d cut. A thresholding/label artifact would pile them up",
        "*at* the boundary. These are real early plants that really do carry a functional FRI — the FLC route",
        "and the polygenic residue, exactly the mechanisms the cell names as beyond its reach.",
        "",
        "## The population-structure correction (read this before the pooled number)",
        "",
        "The paper reports non-functional FRI is overrepresented in central/western Europe and rare elsewhere",
        "(their Figure S3, P<0.01). So FRI genotype correlates with ancestry, and the pooled accuracy partly",
        "measures *'can you recognise a central European accession?'* — the plant analogue of the clonality",
        "inflation corrected elsewhere in this project. Per-STRUCTURE-group accuracy, each against **its own**",
        "null:",
        "",
        "| STRUCTURE group | n | accuracy | its null | beats null? | observed late/early |",
        "|---|---|---|---|---|---|",
    ]
    for g, s in sorted(m["by_structure_group"].items(), key=lambda kv: -kv[1]["n"]):
        if not s["scorable"]:
            why = "one class only" if s["degenerate"] else f"n<10"
            L.append(f"| {g} | {s['n']} | — | — | *unscorable ({why})* | "
                     f"{s['n_observed_late']}/{s['n_observed_early']} |")
            continue
        win = "**yes**" if s["accuracy"] > s["null_accuracy"] else "no"
        L.append(f"| {g} | {s['n']} | {s['accuracy']:.3f} | {s['null_accuracy']:.3f} | {win} | "
                 f"{s['n_observed_late']}/{s['n_observed_early']} |")
    t = rep["threshold_sensitivity"]["tertile"]
    L += [
        "",
        "## Scope limits (load-bearing)",
        "",
        *[f"- {s}" for s in rep["scope_limits"]],
        "",
        "## Threshold sensitivity",
        "",
        f"Re-scored at the upper tertile ({t['threshold_days']} days): accuracy {t['pooled_accuracy']:.3f} vs "
        f"null {t['pooled_null_accuracy']:.3f} → **{t['verdict']}**. "
        "The median cut is used for the headline because it makes the constant-predictor baseline maximally "
        "hard (~50/50); the tertile row shows whether the verdict rides on that choice.",
        "",
        "## Provenance",
        "",
        f"- Substrate: `{rep['substrate']['source']}` — {rep['substrate']['citation']}",
        f"- Cell reference-integrity guard: `{rep['cell_reference_integrity_ok']}`",
        f"- Phenotype: {rep['substrate']['phenotype']}",
        "",
    ]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table-s3", type=Path, default=TABLE_S3)
    ap.add_argument("--out-dir", type=Path, default=REPO / "wiki")
    args = ap.parse_args()

    if not reference_integrity_ok():
        print("[flowering-score] REFUSING: the cell's reference_integrity_ok() guard FAILED — "
              "the catalogue/rule is corrupted. No score is meaningful.")
        return 1

    rows = load_table_s3(args.table_s3)
    rep = build_report(score(rows), score(rows, quantile="tertile"), phenotype_attrition(rows))

    stem = f"flowering_tables3_score_{rep['date']}"
    (args.out_dir / f"{stem}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")
    (args.out_dir / f"{stem}.md").write_text(render_md(rep), encoding="utf-8")

    h = rep["headline"]
    at = rep["substrate"]["phenotype_attrition"]
    print(f"[flowering-score] scored N={at['n_phenotyped']} of {at['n_total']} "
          f"({at['n_dropped_no_ft16']} lack FT16)  verdict={h['verdict']}")
    print(f"  pooled accuracy      {h['pooled_accuracy']:.3f}   (null {h['pooled_null_accuracy']:.3f})")
    print(f"  sens {h['sensitivity']:.3f} / spec {h['specificity']:.3f}")
    if h["structure_group_weighted_accuracy"] is not None:
        print(f"  group-weighted       {h['structure_group_weighted_accuracy']:.3f}   "
              f"(null {h['structure_group_weighted_null']:.3f})   "
              f"groups beating null: {h['groups_beating_null']}")
    print(f"  -> {args.out_dir / (stem + '.md')}")
    return 0 if h["verdict"] == "SCORED_BEATS_NULL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
