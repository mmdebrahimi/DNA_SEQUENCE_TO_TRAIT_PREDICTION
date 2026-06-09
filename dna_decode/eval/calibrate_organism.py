"""calibrate_organism — auto-select the per-organism AMR rule config from a labeled cohort.

Motivation (wider-AMR validation thread, 2026-06-08): the deterministic AMR rule's correct configuration
is ORGANISM-SPECIFIC, and not just the threshold. The validation samples showed three boundary flavors:

- TUNING   — right genes, wrong integer. Campylobacter cipro needs threshold 1 (single gyrA T86I),
             E. coli/Klebsiella need threshold 2 (gyrA+parC double mutant).
- CONTENT  — wrong genes counted. Acinetobacter meropenem over-calls because the intrinsic blaOXA-51
             family is in every isolate; the fix is to EXCLUDE intrinsic gene families. Salmonella cipro
             UNDER-calls because the QRDR-point-only counter excludes the legitimate plasmid qnr genes;
             the fix is a BROADER counter. (The same QRDR-point counter that makes Klebsiella perfect —
             excludes intrinsic OqxAB — makes Salmonella fail — excludes real qnr. So the COUNTER itself
             is organism-specific.)
- EXPRESSION — right genes, can't see their regulation (derepressed AmpC, efflux overexpression). This is
             the hard floor: NO presence-based config can fix it. calibrate_organism cannot cross it and
             reports the residual ceiling honestly.

This module does TUNING + CONTENT auto-calibration from a >=~15R/15S labeled cohort:
  1. choose the COUNTER (qrdr_point vs broad drug-class) and THRESHOLD by leave-one-out balanced accuracy,
  2. auto-flag INTRINSIC gene families (present in >=90% of BOTH R and S) and exclude them from counts.
It returns a CalibratedRule + the honest LOO estimate; it never claims to cross the EXPRESSION floor.

Validated 2026-06-08 falsifier (wiki/self_calibration_falsifier_2026-06-08.md): auto-threshold recovers
Campylobacter->1 / Klebsiella->2; family-level intrinsic-flag recovers blaOXA-51-family + blaADC.

Pure core (operates on per-strain determinant SYMBOL lists, no file I/O) so it is testable on synthetic
fixtures; `features_from_main_tsv` is the thin AMRFinder-main.tsv adapter.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from dna_decode.eval.amr_rules import cipro_determinants_from_main, qrdr_point_determinants

COUNTERS = ("qrdr_point", "broad")
THRESHOLDS = (1, 2, 3)

# qrdr_point counts fluoroquinolone QRDR target-site mutations — meaningful ONLY for fluoroquinolones.
# For other drugs it is all-zero (a degenerate predict-all-S), so restrict the counter grid by drug.
_FLUOROQUINOLONES = ("ciprofloxacin", "levofloxacin", "moxifloxacin", "norfloxacin", "ofloxacin")


def counters_for(drug: str) -> tuple[str, ...]:
    """Drug-appropriate counter grid. qrdr_point only for fluoroquinolones; broad for everything."""
    return ("qrdr_point", "broad") if drug.lower() in _FLUOROQUINOLONES else ("broad",)

# Beta-lactamase intrinsic FAMILY groupings: alleles of one intrinsic family appear under many distinct
# AMRFinder symbols (blaOXA-68/69/66/100/312 are all OXA-51-family), so per-SYMBOL prevalence never reaches
# the intrinsic threshold even when the FAMILY is universal. Collapse to family tokens before prevalence.
_OXA51_FAMILY = (
    "OXA-51", "OXA-64", "OXA-65", "OXA-66", "OXA-68", "OXA-69", "OXA-71", "OXA-90", "OXA-100",
    "OXA-113", "OXA-200", "OXA-202", "OXA-203", "OXA-208", "OXA-219", "OXA-312",
)


def is_point_mutation(symbol: str) -> bool:
    """True if the AMRFinder symbol is a target-site POINT mutation (e.g. gyrA_S83L, nalC_G71E).

    Point mutations are NEVER 'intrinsic genes' in the over-call sense — they are discriminative by COUNT
    (a wild-type strain simply has no such row), so they must be excluded from intrinsic-family flagging.
    AMRFinder POINT symbols carry an underscore-separated substitution suffix (gyrA_S83L, nalC_G71E); this
    proxy holds for every AMRFinder symbol observed (acquired genes like qnrB19/blaCTX-M-15 carry no '_').
    DEFERRED hardening (M2): prefer the authoritative AMRFinder `Method`=POINT column over this proxy —
    requires carrying Method through `features_from_main_tsv`, a feature-shape change tracked for follow-up."""
    return "_" in symbol


def family_token(symbol: str) -> str:
    """Collapse an AMRFinder element symbol to its gene-family token for prevalence grouping.

    blaOXA-51-family alleles -> 'blaOXA-51-family'; otherwise strip a trailing allele number
    (blaADC-76 -> blaADC, blaTEM-1 -> blaTEM, qnrB19 -> qnrB) and a POINT-mutation suffix
    (gyrA_S83L -> gyrA). Idempotent."""
    s = symbol.strip()
    if any(t in s for t in _OXA51_FAMILY):
        return "blaOXA-51-family"
    s = s.split("_", 1)[0]            # drop POINT-mutation suffix (gyrA_S83L -> gyrA)
    s = re.sub(r"-?\d+$", "", s)      # drop trailing allele number (blaADC-76 -> blaADC, qnrB19 -> qnrB)
    return s


def balanced_accuracy(counts_labels: list[tuple[int, str]], threshold: int) -> float:
    """Balanced accuracy of the rule 'R iff count >= threshold' over (count, label) pairs."""
    tp = fp = tn = fn = 0
    for n, lab in counts_labels:
        pred = "R" if n >= threshold else "S"
        if (lab, pred) == ("R", "R"): tp += 1
        elif (lab, pred) == ("S", "S"): tn += 1
        elif (lab, pred) == ("S", "R"): fp += 1
        else: fn += 1
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    return 0.5 * (sens + spec)


def confusion(counts_labels: list[tuple[int, str]], threshold: int) -> dict:
    tp = fp = tn = fn = 0
    for n, lab in counts_labels:
        pred = "R" if n >= threshold else "S"
        if (lab, pred) == ("R", "R"): tp += 1
        elif (lab, pred) == ("S", "S"): tn += 1
        elif (lab, pred) == ("S", "R"): fp += 1
        else: fn += 1
    N = len(counts_labels)
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    return {"acc": (tp + tn) / N if N else 0.0, "sens": sens, "spec": spec,
            "bal_acc": 0.5 * (sens + spec), "tp": tp, "tn": tn, "fp": fp, "fn": fn}


def intrinsic_families(per_strain_families: list[tuple[set[str], str]], *, min_prev: float = 0.90) -> list[str]:
    """Gene families present in >= min_prev of BOTH R and S strains => intrinsic / non-discriminating.

    `per_strain_families` is a list of (set-of-family-tokens, label). Returns sorted family tokens to
    exclude. A family equally common in both classes carries ~0 bits about the label (it is intrinsic or a
    ubiquitous background gene), so counting it only adds over-call noise."""
    R = [fams for fams, lab in per_strain_families if lab == "R"]
    S = [fams for fams, lab in per_strain_families if lab == "S"]
    if not R or not S:
        return []
    out = []
    allfams = set().union(*[f for f, _ in per_strain_families]) if per_strain_families else set()
    for fam in allfams:
        pr = sum(1 for f in R if fam in f) / len(R)
        ps = sum(1 for f in S if fam in f) / len(S)
        if pr >= min_prev and ps >= min_prev:
            out.append(fam)
    return sorted(out)


def _counts_for(strains: list[dict], counter: str, exclude: set[str]) -> list[int]:
    """Per-strain determinant count under `counter`, excluding intrinsic families."""
    return [sum(1 for sym in s[counter] if family_token(sym) not in exclude) for s in strains]


# LOO balanced-accuracy floor below which NO presence-based config is trustworthy — the organism×drug is
# EXPRESSION-driven (derepression/efflux/porin) and calibrate_organism honestly ABSTAINS rather than ship a
# bad rule. Implements the synthesis's "abstain, don't predict" for out-of-trust-zone organism×drugs.
# NOTE: the floor applies ONLY once the cohort is balanced + powered (see MIN_CLASS_COUNT); weak cohorts
# route to INSUFFICIENT_EVIDENCE instead, so EXPRESSION_FLOOR stays a mechanistic claim (not a dumping
# ground for under-powered data). 0.70 is a documented post-hoc constant (chosen from the gap between the
# 1.0 calibrated cohorts and the floor cohorts); revisit if a per-drug floor is ever warranted.
CALIBRATION_FLOOR = 0.70

# Minimum number of strains-with-runs per class for a calibration verdict to be meaningful. Below this
# (or one class entirely absent) -> INSUFFICIENT_EVIDENCE. A one-class cohort cannot distinguish
# "expression floor" from "no resistant examples loaded" (the Pseudomonas degenerate-cohort lesson).
MIN_CLASS_COUNT = 5


@dataclass
class CalibratedRule:
    drug: str
    counter: str
    threshold: int
    intrinsic_families_excluded: list[str]
    loo_balanced_accuracy: float
    full_balanced_accuracy: float
    full_confusion: dict
    per_config: dict                      # {"counter@thr": full_bal_acc}
    fold_selections: dict                 # {"counter@thr": times_selected_in_LOO}
    n: int
    n_R: int
    n_S: int
    verdict: str = "CALIBRATED"           # CALIBRATED | EXPRESSION_FLOOR (abstain) | INSUFFICIENT_EVIDENCE
    note: str = ""

    def predict(self, strain_features: dict) -> str:
        """Apply the calibrated rule to one strain's {counter: [symbols]} features -> 'R'/'S'."""
        excl = set(self.intrinsic_families_excluded)
        n = sum(1 for sym in strain_features.get(self.counter, []) if family_token(sym) not in excl)
        return "R" if n >= self.threshold else "S"


def _config_grid(counters, thresholds):
    return [(c, t) for c in counters for t in thresholds]


def _select_best_config(strains, labels, exclude, counters, thresholds):
    """Return (counter, threshold) maximizing full-cohort balanced accuracy.
    Tie-break (deterministic, anti-overfit): higher bal_acc, then LOWER threshold (more sensitive),
    then prefer the MORE-SPECIFIC counter 'qrdr_point' over 'broad' (qrdr_point excludes efflux noise)."""
    counter_rank = {"qrdr_point": 0, "broad": 1}
    best = None
    for c, t in _config_grid(counters, thresholds):
        counts = list(zip(_counts_for(strains, c, exclude), labels))
        ba = balanced_accuracy(counts, t)
        score = (ba, -t, -counter_rank.get(c, 9))   # maximize ba, then lower t, then lower counter_rank
        if best is None or score > best[0]:
            best = (score, c, t, ba)
    return best[1], best[2], best[3]


def calibrate(strains: list[dict], labels: list[str], drug: str, *,
              counters: tuple[str, ...] | None = None, thresholds: tuple[int, ...] = THRESHOLDS,
              intrinsic_min_prev: float = 0.90) -> CalibratedRule:
    """Auto-select (counter, threshold) + intrinsic exclusions for `drug` from a labeled cohort.

    `strains`: list of per-strain feature dicts {counter_name: [determinant_symbols]} (see
    `features_from_main_tsv`). `labels`: aligned list of 'R'/'S'. The DEPLOYED config is the deterministic
    full-cohort `_select_best_config` pick; `loo_balanced_accuracy` separately estimates the selection
    PROCEDURE's generalization (per-fold re-selection). Verdict: INSUFFICIENT_EVIDENCE (a class <
    MIN_CLASS_COUNT) / CALIBRATED (LOO balanced-acc >= floor) / EXPRESSION_FLOOR (powered but sub-floor).
    Intrinsic families computed once on the full cohort for the deployed rule (documented in-sample)."""
    if len(strains) != len(labels):
        raise ValueError(f"strains ({len(strains)}) and labels ({len(labels)}) length mismatch")
    if not strains:
        raise ValueError("empty cohort")
    if counters is None:
        counters = counters_for(drug)
    n_R = labels.count("R"); n_S = labels.count("S")

    # intrinsic families (full-cohort) — gene-PRESENCE symbols only (point mutations excluded: they are
    # discriminative-by-count, never intrinsic genes) -> family tokens
    per_strain_fams = [({family_token(sym) for c in s for sym in s[c] if not is_point_mutation(sym)}, lab)
                       for s, lab in zip(strains, labels)]
    excl = set(intrinsic_families(per_strain_fams, min_prev=intrinsic_min_prev))

    # full-cohort bal-acc per config (transparency)
    per_config = {}
    for c, t in _config_grid(counters, thresholds):
        counts = list(zip(_counts_for(strains, c, excl), labels))
        per_config[f"{c}@{t}"] = round(balanced_accuracy(counts, t), 4)

    # DEPLOYED config = deterministic full-cohort selection (NOT a modal LOO pick — that introduced a
    # selection-procedure-vs-deployed-rule mismatch + a Counter.most_common tie-break ambiguity). The LOO
    # number below estimates the SELECTION PROCEDURE's generalization; the deployed rule is this.
    mc, mt, _ = _select_best_config(strains, labels, excl, counters, thresholds)
    full_conf = confusion(list(zip(_counts_for(strains, mc, excl), labels)), mt)

    # Min-class guard: a one-class / under-powered cohort cannot yield a meaningful verdict. Short-circuit
    # to INSUFFICIENT_EVIDENCE BEFORE the floor check so EXPRESSION_FLOOR stays a mechanistic claim.
    if min(n_R, n_S) < MIN_CLASS_COUNT:
        return CalibratedRule(
            drug=drug, counter=mc, threshold=mt, intrinsic_families_excluded=sorted(excl),
            loo_balanced_accuracy=None, full_balanced_accuracy=round(full_conf["bal_acc"], 4),
            full_confusion=full_conf, per_config=per_config, fold_selections={},
            n=len(strains), n_R=n_R, n_S=n_S, verdict="INSUFFICIENT_EVIDENCE",
            note=(f"INSUFFICIENT_EVIDENCE: {n_R}R/{n_S}S — fewer than MIN_CLASS_COUNT={MIN_CLASS_COUNT} in "
                  f"a class (or one class absent). A one-class/under-powered cohort cannot distinguish an "
                  f"expression floor from missing examples. No verdict; load more strains-with-runs."))

    # LOO: per held-out strain select best config on the REST, predict held-out; collect (pred, true) so we
    # can compute a TRUE balanced accuracy (the field used to be plain accuracy — misnamed).
    fold_sel: Counter = Counter()
    loo_preds: list[tuple[str, str]] = []
    idx = list(range(len(strains)))
    for i in idx:
        train_s = [strains[j] for j in idx if j != i]
        train_l = [labels[j] for j in idx if j != i]
        tf = [({family_token(sym) for c in s for sym in s[c] if not is_point_mutation(sym)}, lab)
              for s, lab in zip(train_s, train_l)]
        fold_excl = set(intrinsic_families(tf, min_prev=intrinsic_min_prev))
        c, t, _ = _select_best_config(train_s, train_l, fold_excl, counters, thresholds)
        fold_sel[f"{c}@{t}"] += 1
        n_held = sum(1 for sym in strains[i].get(c, []) if family_token(sym) not in fold_excl)
        loo_preds.append(("R" if n_held >= t else "S", labels[i]))
    tp = sum(1 for p, y in loo_preds if p == "R" and y == "R")
    fn = sum(1 for p, y in loo_preds if p == "S" and y == "R")
    tn = sum(1 for p, y in loo_preds if p == "S" and y == "S")
    fp = sum(1 for p, y in loo_preds if p == "R" and y == "S")
    loo_sens = tp / (tp + fn) if (tp + fn) else None
    loo_spec = tn / (tn + fp) if (tn + fp) else None
    loo_ba = round(0.5 * (loo_sens + loo_spec), 4) if (loo_sens is not None and loo_spec is not None) else None

    deployed_cfg = f"{mc}@{mt}"
    if loo_ba is None:
        verdict = "INSUFFICIENT_EVIDENCE"
        note = "INSUFFICIENT_EVIDENCE: a class was absent in LOO predictions; no balanced accuracy."
    elif loo_ba >= CALIBRATION_FLOOR:
        verdict = "CALIBRATED"
        cmp = "qrdr_point@2" if drug.lower() == "ciprofloxacin" else None
        note = (f"deployed {deployed_cfg}; LOO balanced-acc {loo_ba:.3f}."
                + (f" (differs from cipro default {cmp} — organism-specific.)" if cmp and deployed_cfg != cmp else ""))
    else:
        verdict = "EXPRESSION_FLOOR"
        note = (f"NO presence-based config clears the floor (LOO balanced-acc {loo_ba:.3f} < "
                f"{CALIBRATION_FLOOR}) on a balanced+powered cohort ({n_R}R/{n_S}S) — EXPRESSION-driven "
                f"(derepression/efflux overexpression/porin loss); gene-presence cannot decode it. ABSTAIN. "
                f"Intrinsic families excluded ({sorted(excl)}) but the residual is the uncrossable floor.")
    return CalibratedRule(
        drug=drug, counter=mc, threshold=mt, intrinsic_families_excluded=sorted(excl),
        loo_balanced_accuracy=loo_ba, full_balanced_accuracy=round(full_conf["bal_acc"], 4),
        full_confusion=full_conf, per_config=per_config, fold_selections=dict(fold_sel),
        n=len(strains), n_R=n_R, n_S=n_S, verdict=verdict, note=note)


def features_from_main_tsv(main_tsv: Path, drug: str) -> dict:
    """Build one strain's {counter: [symbols]} features from an AMRFinder main.tsv.

    'qrdr_point' = QRDR target-alteration POINT mutations; 'broad' = all drug-class determinants. Absent
    main.tsv -> empty lists (a determinant-free strain, predicted S)."""
    qp = [d["symbol"] for d in qrdr_point_determinants(main_tsv)]
    broad = [d["symbol"] for d in cipro_determinants_from_main(main_tsv, drug)]
    return {"qrdr_point": qp, "broad": broad}
