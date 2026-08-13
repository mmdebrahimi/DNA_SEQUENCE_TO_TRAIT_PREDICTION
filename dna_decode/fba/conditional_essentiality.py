"""CONDITIONAL gene essentiality — the same gene, four media, and the model has to get the switch right.

Every essentiality number in this repo so far has been **single-condition**: one medium, one E/N call per
gene. That measures whether the model knows a gene is needed *at all*. It cannot measure the thing a strain
designer actually relies on — that a gene is dispensable on glucose and required on succinate, so deleting
it costs nothing in one process and kills the strain in another.

The Orth 2011 iJO1366 screen supplies exactly that: **1,075 E. coli K-12 genes × 4 minimal media**
(glucose aerobic / glucose anaerobic / lactate aerobic / succinate aerobic), each an experimental E or N,
with **68 genes conditionally essential** — essential in at least one medium and dispensable in at least
one other. Those 68 are the signal; the rest is background either way.

Two properties make this a better substrate than the single-condition sets:

1. **Two-sided by construction.** A conditionally-essential gene is its own control: the same gene, the
   same model, the same GPR — only the medium changes. A model that scores well by calling everything
   dispensable cannot fake this.
2. **It carries its own reproduction gate.** The supplement also ships the paper's OWN iJO1366 FBA
   predictions per condition, so the pipeline can be checked against a published result BEFORE any new
   number is trusted — the discipline that has caught two unit/wiring errors in this repo already.

The paper's FBA columns are a GATE, never a label. Labels are the experimental columns only.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# The four media, as (carbon-source exchange, aerobic?) against the model's own M9 mineral background.
# iML1515 / iJO1366 ship glucose-M9-aerobic as their default medium, so only the carbon source and the
# oxygen bound move. Everything else (the mineral panel) is left exactly as the reconstruction defines it.
CONDITIONS: dict[str, tuple[str, bool]] = {
    "glucose_aerobic": ("EX_glc__D_e", True),
    "glucose_anaerobic": ("EX_glc__D_e", False),
    "lactate_aerobic": ("EX_lac__L_e", True),
    "succinate_aerobic": ("EX_succ_e", True),
}

# Every carbon-source exchange any condition may use -- all are closed before one is opened, so a residual
# glucose uptake can never leak into the lactate or succinate condition (the failure mode that would make
# every medium silently score as glucose).
_ALL_CARBON = ("EX_glc__D_e", "EX_lac__L_e", "EX_lac__D_e", "EX_succ_e")

DEFAULT_LABELS = (Path(__file__).resolve().parents[2] / "data" / "raw" /
                  "ecoli_conditional_essentiality" /
                  "orth2011_table_s1_conditional_essentiality.tsv")


@dataclass(frozen=True)
class GeneRecord:
    """One gene's experimental calls + the paper's own FBA calls, across the four media."""
    gene_id: str
    gene: str
    experimental: dict[str, bool]      # condition -> essential?
    paper_fba: dict[str, bool]
    conditionally_essential: bool

    @property
    def n_conditions_essential(self) -> int:
        return sum(self.experimental.values())


def load_labels(path: str | Path | None = None) -> list[GeneRecord]:
    """Parse the committed gold standard. Comment lines start with '#'."""
    p = Path(path) if path else DEFAULT_LABELS
    rows: list[GeneRecord] = []
    header: list[str] | None = None
    for line in p.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if header is None:
            header = parts
            continue
        d = dict(zip(header, parts, strict=False))
        rows.append(GeneRecord(
            gene_id=d["gene_id"], gene=d["gene"],
            experimental={c: d[f"exp_{c}"] == "E" for c in CONDITIONS},
            paper_fba={c: d[f"fba_{c}"] == "E" for c in CONDITIONS},
            conditionally_essential=d["conditionally_essential"] == "YES",
        ))
    return rows


def conditionally_essential_genes(records: list[GeneRecord]) -> list[GeneRecord]:
    """Genes essential in >=1 medium AND dispensable in >=1 medium -- the two-sided subset.

    Recomputed from the experimental columns rather than trusting the supplement's own YES/NO flag, so a
    disagreement between the two is visible instead of inherited.
    """
    out = []
    for r in records:
        vals = list(r.experimental.values())
        if any(vals) and not all(vals):
            out.append(r)
    return out


def apply_condition(model, condition: str, uptake: float = 10.0) -> None:
    """Set the model's medium to one of the four conditions, IN PLACE.

    Closes every candidate carbon exchange first. Anaerobic sets the oxygen uptake to 0 rather than
    removing the exchange, so the reaction stays present and the only difference between the two glucose
    conditions is that single bound.
    """
    if condition not in CONDITIONS:
        raise KeyError(f"unknown condition {condition!r}; known: {sorted(CONDITIONS)}")
    carbon, aerobic = CONDITIONS[condition]
    have = {r.id for r in model.exchanges}
    if carbon not in have:
        raise KeyError(f"model {model.id} has no exchange {carbon!r} for condition {condition!r}")

    medium = dict(model.medium)
    for ex in _ALL_CARBON:
        medium.pop(ex, None)
    medium[carbon] = uptake
    if "EX_o2_e" in have:
        medium["EX_o2_e"] = 1000.0 if aerobic else 0.0
    model.medium = medium


def confusion_from_calls(experimental: dict[str, bool], predicted: dict[str, bool]) -> dict[str, int]:
    """TP/FP/FN/TN over a shared key set (positive = essential)."""
    keys = sorted(set(experimental) & set(predicted))
    tp = sum(1 for k in keys if experimental[k] and predicted[k])
    fp = sum(1 for k in keys if not experimental[k] and predicted[k])
    fn = sum(1 for k in keys if experimental[k] and not predicted[k])
    tn = sum(1 for k in keys if not experimental[k] and not predicted[k])
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "n": len(keys)}


def mcc(cm: dict[str, int]) -> float:
    """Matthews correlation. 0.0 on a degenerate margin (a single-class prediction), never a ZeroDivision."""
    tp, fp, fn, tn = cm["tp"], cm["fp"], cm["fn"], cm["tn"]
    den = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    return 0.0 if den == 0 else (tp * tn - fp * fn) / den


def continuous_readout(records: list[GeneRecord], ratios: dict[str, dict[str, float]]) -> dict:
    """Is the conditional signal ABSENT, or is the binary CUTOFF discarding it?

    FBA computes a continuous knockout growth ratio (mutant / wild type) per condition, and the deployed
    call thresholds it at 1% of wild type. That throws away everything in between, so this scores the raw
    ratio as a RANKING over every gene x condition cell on the two-sided subset (lower growth = more
    essential) and reports the best threshold the data could support.

    Measured 2026-08-12 on iML1515 / 268 cells (119 essential):
      * AUROC **0.5907** -- weak but above chance, so the sub-threshold variation is not pure noise.
      * deployed cutoff (<=0.01): MCC 0.0918, TP 10 / FN 109.
      * ORACLE best cutoff (<=0.3249): MCC 0.2544, TP 24, with the SAME 6 false positives.

    So retuning would roughly TRIPLE the conditional MCC at no precision cost -- but 64% of these genes have
    a perfectly FLAT ratio across all four media, so even the oracle leaves 95 of 119 essential cells
    missed. The readout costs real signal; most of the deficit is still the model.

    **`oracle_*` is fitted ON the evaluation set and is an UPPER BOUND, never a deployable number** -- the
    same rail the Track B deltaG arm carries. `deployable_threshold` is the honest counterpart.

    **AUROC carries run-to-run variation of about +-0.01** (0.598 / 0.6110 / 0.6099 over three identical
    invocations) because degenerate LP optima shift mid-range growth ratios between processes. The
    THRESHOLDED numbers are byte-stable across the same runs, since the shifts never cross a cutoff. Quote
    the AUROC as ~0.60, not to four decimals.
    """
    keys = sorted(CONDITIONS)
    y: list[int] = []
    score: list[float] = []
    for r in records:
        for c in keys:
            if r.gene_id in ratios.get(c, {}):
                y.append(1 if r.experimental[c] else 0)
                score.append(ratios[c][r.gene_id])
    n_pos, n_neg = sum(y), len(y) - sum(y)
    if not y or n_pos == 0 or n_neg == 0:
        return {"n_cells": len(y), "auroc": None, "note": "degenerate: one class only"}

    # rank-based AUROC with average ranks for ties (ties are common -- many ratios are exactly 1.0)
    paired = sorted(range(len(score)), key=lambda i: -score[i])   # ascending essentiality score
    ranks = [0.0] * len(score)
    for pos, i in enumerate(paired, start=1):
        ranks[i] = float(pos)
    by_val: dict[float, list[int]] = {}
    for i, v in enumerate(score):
        by_val.setdefault(v, []).append(i)
    for idxs in by_val.values():
        if len(idxs) > 1:
            avg = sum(ranks[i] for i in idxs) / len(idxs)
            for i in idxs:
                ranks[i] = avg
    pos_rank_sum = sum(ranks[i] for i in range(len(y)) if y[i] == 1)
    auroc = (pos_rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

    def at(thr: float) -> dict[str, int]:
        tp = sum(1 for i in range(len(y)) if score[i] <= thr and y[i] == 1)
        fp = sum(1 for i in range(len(y)) if score[i] <= thr and y[i] == 0)
        fn = sum(1 for i in range(len(y)) if score[i] > thr and y[i] == 1)
        tn = sum(1 for i in range(len(y)) if score[i] > thr and y[i] == 0)
        return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "n": len(y)}

    best_thr, best_cm, best_mcc = None, None, -2.0
    for thr in sorted(set(score)):
        cm = at(thr)
        v = mcc(cm)
        if v > best_mcc:
            best_thr, best_cm, best_mcc = thr, cm, v
    deployed = at(0.01)
    return {
        "n_cells": len(y), "n_essential_cells": n_pos,
        "auroc": round(auroc, 4),
        "deployed_threshold": 0.01,
        "deployed_confusion": deployed,
        "deployed_mcc": round(mcc(deployed), 4),
        "oracle_threshold": round(best_thr, 4),
        "oracle_confusion": best_cm,
        "oracle_mcc": round(best_mcc, 4),
        "oracle_note": ("the oracle threshold is fitted ON the evaluation set -- an UPPER BOUND, not a "
                        "deployable number. A deployable one needs a disjoint tuning split."),
    }


def deployable_threshold(records: list[GeneRecord], ratios: dict[str, dict[str, float]],
                         n_folds: int = 5) -> dict:
    """The honest version of `continuous_readout`'s oracle: fit the cutoff on a DISJOINT gene split.

    The oracle threshold is chosen on the same cells it is scored on, so it cannot be deployed. Here the
    conditionally-essential genes are split into `n_folds` groups **by gene** (never by cell -- splitting
    cells would leak, since the four cells of one gene share its ratio profile). For each fold the cutoff
    that maximises MCC on the OTHER folds is applied to the held-out one, and the held-out predictions are
    pooled into a single confusion matrix.

    If the held-out MCC lands near the oracle, the retune is real and deployable. If it collapses toward
    the deployed cutoff's score, the oracle was fitting noise -- which is exactly what a 268-cell set with
    AUROC 0.598 might do, and the reason this function exists rather than shipping the oracle.

    Folds are assigned by sorted gene id (deterministic, no RNG) so the number is reproducible.
    """
    keys = sorted(CONDITIONS)
    subset = sorted(conditionally_essential_genes(records), key=lambda r: r.gene_id)
    cells: list[tuple[str, int, float]] = []          # (gene_id, truth, ratio)
    for r in subset:
        for c in keys:
            if r.gene_id in ratios.get(c, {}):
                cells.append((r.gene_id, 1 if r.experimental[c] else 0, ratios[c][r.gene_id]))
    if not cells:
        return {"n_cells": 0, "held_out_mcc": None, "note": "no cells"}

    fold_of = {r.gene_id: i % n_folds for i, r in enumerate(subset)}

    def best_threshold(train: list[tuple[str, int, float]]) -> float:
        best, best_v = 0.01, -2.0
        for thr in sorted({c[2] for c in train}):
            tp = sum(1 for _, y, s in train if s <= thr and y == 1)
            fp = sum(1 for _, y, s in train if s <= thr and y == 0)
            fn = sum(1 for _, y, s in train if s > thr and y == 1)
            tn = sum(1 for _, y, s in train if s > thr and y == 0)
            v = mcc({"tp": tp, "fp": fp, "fn": fn, "tn": tn})
            if v > best_v:
                best, best_v = thr, v
        return best

    tp = fp = fn = tn = 0
    chosen = []
    for f in range(n_folds):
        train = [c for c in cells if fold_of[c[0]] != f]
        test = [c for c in cells if fold_of[c[0]] == f]
        if not train or not test:
            continue
        thr = best_threshold(train)
        chosen.append(round(thr, 4))
        for _, y, s in test:
            pred = s <= thr
            if pred and y == 1:
                tp += 1
            elif pred and y == 0:
                fp += 1
            elif not pred and y == 1:
                fn += 1
            else:
                tn += 1
    cm = {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "n": tp + fp + fn + tn}
    return {
        "n_cells": len(cells), "n_folds": n_folds,
        "thresholds_per_fold": chosen,
        "held_out_confusion": cm,
        "held_out_mcc": round(mcc(cm), 4),
        "note": ("cutoff fitted on disjoint gene folds and scored on held-out genes -- a DEPLOYABLE "
                 "estimate, unlike the oracle. Folds split BY GENE because the four cells of one gene "
                 "share a ratio profile and splitting cells would leak."),
    }


def pattern_distribution(records: list[GeneRecord],
                         predicted: dict[str, dict[str, bool]] | None = None,
                         conditions: tuple[str, ...] | None = None) -> dict:
    """WHY the switch score is low: what shape are the predictions, on the two-sided subset?

    Each gene becomes a 4-character pattern over the sorted conditions ('E' essential, '.' not), so the
    true patterns can be compared to the predicted ones directly.

    Measured 2026-08-12 on the 68 conditionally-essential genes: the TRUE patterns span 12 distinct shapes,
    while the paper's own iJO1366 predicts a CONSTANT pattern ('....' or 'EEEE') for **62 of 68 (91%)**.
    The model is not getting the switch wrong -- it is not making a switch at all, which is exactly why it
    lands within a point of the constant-predictor null.

    `predicted=None` describes the experimental labels themselves.
    """
    subset = conditionally_essential_genes(records)
    # MUST accept the caller's condition set. Defaulting to the 4-media CONDITIONS against a
    # 25-carbon-source prediction made every lookup miss -> every gene read as '....' -> a FALSE
    # '100% constant across 1 shape', which contradicted a positive exact-set count in the same run.
    keys = sorted(conditions) if conditions is not None else sorted(CONDITIONS)

    def pat(d: dict[str, bool]) -> str:
        return "".join("E" if d.get(c) else "." for c in keys)

    counts: dict[str, int] = {}
    n_constant = 0
    for r in subset:
        if predicted is None:
            p = pat(r.experimental)
        else:
            p = pat({c: predicted.get(c, {}).get(r.gene_id, False) for c in keys})
        counts[p] = counts.get(p, 0) + 1
        # A pattern is CONSTANT iff every position agrees -- length-agnostic. The literal ("....", "EEEE")
        # test that lived here was hardcoded to FOUR conditions, so on the 25-carbon-source panel a
        # 25-character all-dispensable or all-essential pattern matched NEITHER and the run reported
        # "0.0% constant" when the true figure was 184/217 = 84.8%. That is the SECOND hardcoded-4
        # assumption found in this one function; generalising `keys` alone was not enough.
        if len(set(p)) == 1:
            n_constant += 1
    return {
        "conditions_order": keys,
        "n_genes": len(subset),
        "patterns": dict(sorted(counts.items(), key=lambda kv: -kv[1])),
        "n_distinct_patterns": len(counts),
        "n_constant_pattern": n_constant,
        "constant_pattern_fraction": round(n_constant / len(subset), 4) if subset else None,
    }


def constant_baselines(records: list[GeneRecord],
                       conditions: tuple[str, ...] | None = None) -> dict[str, dict]:
    """The NULL controls the switch metric must always be reported against.

    A per-cell agreement of 0.57 sounds like signal until you notice that predicting **dispensable for
    everything** scores 0.5588 on the same subset -- because most conditionally-essential genes are
    essential in only one or two of the four media. Measured 2026-08-12: iJO1366 0.5735 and iML1515 0.5709
    against that 0.5588 null, i.e. roughly one point of real conditional resolution between them.

    Never quote `per_condition_agreement` without this.
    """
    keys = tuple(conditions) if conditions is not None else tuple(CONDITIONS)
    subset = conditionally_essential_genes(records)
    out = {}
    for name, const in (("always_essential", True), ("always_dispensable", False)):
        pred = {c: {r.gene_id: const for r in subset} for c in keys}
        out[name] = switch_accuracy(records, pred, conditions=keys)
    return out


def switch_accuracy(records: list[GeneRecord], predicted: dict[str, dict[str, bool]],
                    conditions: tuple[str, ...] | None = None) -> dict:
    """THE conditional metric: on the two-sided subset, does the model reproduce the SWITCH?

    For each conditionally-essential gene, compare the exact SET of media in which it is essential.
    `exact_set_match` is the strict version (the whole pattern right); `per_condition_agreement` is the
    lenient one (fraction of gene x condition cells right on that subset).

    A model that calls a gene essential in all four media, or in none, scores 0 on the strict metric even
    though a single-condition metric might score it well -- which is the entire point of measuring here.
    """
    keys = tuple(conditions) if conditions is not None else tuple(CONDITIONS)
    subset = conditionally_essential_genes(records)
    exact = 0
    cells_right = 0
    cells_total = 0
    for r in subset:
        pred_set = {c for c in keys if predicted.get(c, {}).get(r.gene_id, False)}
        true_set = {c for c in keys if r.experimental.get(c, False)}
        if pred_set == true_set:
            exact += 1
        for c in keys:
            cells_total += 1
            if predicted.get(c, {}).get(r.gene_id, False) == r.experimental.get(c, False):
                cells_right += 1
    return {
        "n_conditionally_essential": len(subset),
        "exact_set_match": exact,
        "exact_set_match_rate": round(exact / len(subset), 4) if subset else None,
        "per_condition_agreement": round(cells_right / cells_total, 4) if cells_total else None,
    }
