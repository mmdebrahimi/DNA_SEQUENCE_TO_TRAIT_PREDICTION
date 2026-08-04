"""Keio-collection experimental essentiality labels (the free gold-standard validation label).

Source: RB-TnSeq mutant fitness of the Keio-parent strain BW25113 (Wetmore/Price fitness browser),
as vendored + used by Bernstein et al. 2023 (Mol Syst Biol) to validate iML1515 gene essentiality.
`fit_organism_Keio.tsv` gives per-gene fitness across carbon sources, keyed by b-number (`sysName`).

Experimental essentiality on a carbon source (Bernstein method): a gene with mean fitness on that
source **below `threshold` (default -2.0)** is a growth/no-growth "essential" call. Genes with NO
viable mutant are ABSENT from the table entirely — an independent essentiality signal (no mutant
could be made), reported separately by the validation harness.

All functions here are PURE (parse text -> labels); the network fetch lives in the validation script.
"""
from __future__ import annotations

import statistics

KEIO_FITNESS_URL = (
    "https://raw.githubusercontent.com/dbernste/E_coli_GEM_validation/HEAD/"
    "Fitness_Data/E_coli_BW25113/fit_organism_Keio.tsv"
)
DEFAULT_THRESHOLD = -2.0  # Bernstein 2023 growth/no-growth fitness cutoff


def parse_keio_fitness(
    text: str, carbon: str = "D-Glucose", threshold: float = DEFAULT_THRESHOLD
) -> dict[str, dict]:
    """Parse the fitness TSV -> {b_number: {"fitness": mean, "essential": bool}} for one carbon source.

    Only genes with >=1 usable fitness value on `carbon` are included (the assayable set). A gene is
    experimentally essential-on-<carbon> iff its mean fitness < `threshold`.
    """
    lines = text.splitlines()
    if not lines:
        return {}
    hdr = lines[0].split("\t")
    try:
        syscol = hdr.index("sysName")
    except ValueError:
        return {}
    cols = [i for i, h in enumerate(hdr) if carbon in h]
    if not cols:
        return {}
    out: dict[str, dict] = {}
    for ln in lines[1:]:
        c = ln.split("\t")
        need = max([syscol, *cols])
        if len(c) <= need:
            continue
        b = c[syscol].strip()
        if not b.startswith("b"):
            continue
        vals = []
        for i in cols:
            v = c[i].strip()
            if v in ("", "NA"):
                continue
            try:
                vals.append(float(v))
            except ValueError:
                continue
        if not vals:
            continue
        mean = statistics.fmean(vals)
        out[b] = {"fitness": mean, "essential": mean < threshold}
    return out


def confusion(
    labels: dict[str, bool], predictions: dict[str, bool]
) -> dict[str, int]:
    """PURE: confusion counts over the intersection of labelled + predicted genes.

    positive class = essential. Returns {tp, fp, tn, fn, n}.
    """
    keys = set(labels) & set(predictions)
    tp = fp = tn = fn = 0
    for k in keys:
        exp, pred = labels[k], predictions[k]
        if exp and pred:
            tp += 1
        elif not exp and pred:
            fp += 1
        elif not exp and not pred:
            tn += 1
        else:
            fn += 1
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn, "n": len(keys)}


def metrics_from_confusion(cm: dict[str, int]) -> dict[str, float]:
    """PURE: accuracy, precision, recall, F1, MCC from a confusion dict."""
    tp, fp, tn, fn = cm["tp"], cm["fp"], cm["tn"], cm["fn"]
    n = tp + fp + tn + fn
    acc = (tp + tn) / n if n else 0.0
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    mcc = ((tp * tn) - (fp * fn)) / denom if denom else 0.0
    return {
        "accuracy": acc,
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "mcc": mcc,
        "essential_prevalence": (tp + fn) / n if n else 0.0,
    }
