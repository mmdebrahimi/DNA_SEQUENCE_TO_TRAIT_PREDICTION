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


def constant_baselines(records: list[GeneRecord]) -> dict[str, dict]:
    """The NULL controls the switch metric must always be reported against.

    A per-cell agreement of 0.57 sounds like signal until you notice that predicting **dispensable for
    everything** scores 0.5588 on the same subset -- because most conditionally-essential genes are
    essential in only one or two of the four media. Measured 2026-08-12: iJO1366 0.5735 and iML1515 0.5709
    against that 0.5588 null, i.e. roughly one point of real conditional resolution between them.

    Never quote `per_condition_agreement` without this.
    """
    subset = conditionally_essential_genes(records)
    out = {}
    for name, const in (("always_essential", True), ("always_dispensable", False)):
        pred = {c: {r.gene_id: const for r in subset} for c in CONDITIONS}
        out[name] = switch_accuracy(records, pred)
    return out


def switch_accuracy(records: list[GeneRecord], predicted: dict[str, dict[str, bool]]) -> dict:
    """THE conditional metric: on the two-sided subset, does the model reproduce the SWITCH?

    For each conditionally-essential gene, compare the exact SET of media in which it is essential.
    `exact_set_match` is the strict version (the whole pattern right); `per_condition_agreement` is the
    lenient one (fraction of gene x condition cells right on that subset).

    A model that calls a gene essential in all four media, or in none, scores 0 on the strict metric even
    though a single-condition metric might score it well -- which is the entire point of measuring here.
    """
    subset = conditionally_essential_genes(records)
    exact = 0
    cells_right = 0
    cells_total = 0
    for r in subset:
        pred_set = {c for c in CONDITIONS if predicted.get(c, {}).get(r.gene_id, False)}
        true_set = {c for c in CONDITIONS if r.experimental[c]}
        if pred_set == true_set:
            exact += 1
        for c in CONDITIONS:
            cells_total += 1
            if predicted.get(c, {}).get(r.gene_id, False) == r.experimental[c]:
                cells_right += 1
    return {
        "n_conditionally_essential": len(subset),
        "exact_set_match": exact,
        "exact_set_match_rate": round(exact / len(subset), 4) if subset else None,
        "per_condition_agreement": round(cells_right / cells_total, 4) if cells_total else None,
    }
