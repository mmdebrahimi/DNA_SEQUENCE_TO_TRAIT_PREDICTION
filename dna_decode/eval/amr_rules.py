"""Deterministic, interpretable AMR caller from AMRFinderPlus curated determinants.

The Phase-2 decision (2026-06-05, `plans/AMR_embedding_niche_decision_2026-06-05.md`): for E. coli AMR,
the honest decoder is mechanism features, NOT frozen-pool embeddings (which lost to the QRDR-POINT
knowledge baseline on cipro). This module is the deterministic, fully-transparent caller that decision
points to — it reads AMRFinderPlus's CURATED resistance determinants (main.tsv, not the --mutation_all
noise) and calls R/S per drug, reporting exactly which determinants drove the call.

Rule (v0): predict RESISTANT iff the genome carries >=1 curated determinant whose AMRFinder Class/Subclass
matches the drug (via mic_tiers.amrfinder_classes_for). Pure, interpretable, no trained model. Documented
operating characteristics on the cipro N=147 cohort: accuracy 0.85 / sensitivity 0.96 / specificity 0.75
(the simple "any determinant" rule is sensitivity-favoring — it over-calls strains carrying a single
low-level determinant; multi-determinant strains are high-confidence). See tests/test_amr_rules.py.
"""
from __future__ import annotations

import csv
from pathlib import Path

from dna_decode.data.mic_tiers import amrfinder_classes_for


def cipro_determinants_from_main(main_tsv: Path, drug: str) -> list[dict]:
    """Return the curated AMRFinder determinants (main.tsv rows) relevant to `drug`.

    A row is drug-relevant iff any token of its Class or Subclass matches the drug's AMRFinder class set
    (mic_tiers.amrfinder_classes_for) — case-insensitive substring on the '/'-joined multi-class strings
    AMRFinder emits (e.g. 'AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN'). Each returned dict carries the
    element symbol + name + class for transparent reporting."""
    classes = {c.upper() for c in amrfinder_classes_for(drug)}
    out: list[dict] = []
    p = Path(main_tsv)
    if not p.exists():
        return out
    for r in csv.DictReader(p.open(encoding="utf-8"), delimiter="\t"):
        cls = (r.get("Class") or "").upper()
        sub = (r.get("Subclass") or "").upper()
        if any(c in cls or c in sub for c in classes):
            out.append({
                "symbol": (r.get("Element symbol") or "").strip(),
                "name": (r.get("Element name") or "").strip(),
                "class": r.get("Class") or "", "subclass": r.get("Subclass") or "",
                "pct_identity": r.get("% Identity to reference"),
            })
    return out


_DEFAULT_THRESHOLD = 2  # cipro-validated: clinical fluoroquinolone-R needs multiple QRDR hits (see below)


def call_resistance(main_tsv: Path, drug: str, resistance_threshold: int = _DEFAULT_THRESHOLD) -> dict:
    """Deterministic, interpretable R/S call for one genome's AMRFinder main.tsv.

    Rule (tiered, count-based): R iff #curated drug-relevant determinants >= `resistance_threshold`.
    Default threshold 2 — cipro-validated (N=147: acc 0.939 / sens 0.931 / spec 0.947); clinical
    fluoroquinolone resistance typically needs >=2 QRDR hits (e.g. gyrA + parC), so a single determinant
    is usually low-level/susceptible (cohort: 15/17 single-determinant strains were S). The earlier
    >=1 rule over-called (spec 0.75). NOTE: acquired-gene-dominant drugs (e.g. cef plasmid beta-lactamases,
    where ONE blaCTX-M = resistance) should pass `resistance_threshold=1`; 2 is the QRDR/point-mutation
    default, not universal. Confidence: HIGH if n>=threshold+1 OR n==0; MODERATE near the boundary
    (threshold-1 <= n <= threshold). Absent main.tsv → 'INDETERMINATE'."""
    p = Path(main_tsv)
    if not p.exists():
        return {"prediction": "INDETERMINATE", "confidence": None, "n_determinants": 0,
                "determinants": [], "rule": "amrfinder_curated_determinant_v1", "caveat":
                "AMRFinder main.tsv not found for this genome"}
    dets = cipro_determinants_from_main(p, drug)
    n = len(dets)
    pred = "R" if n >= resistance_threshold else "S"
    if n == 0 or n >= resistance_threshold + 1:
        conf = "HIGH"
    else:                                    # boundary zone (threshold-1 .. threshold)
        conf = "MODERATE"
    boundary_note = ""
    if pred == "S" and n == resistance_threshold - 1 and n >= 1:
        boundary_note = (f" {n} determinant present (below the >={resistance_threshold} threshold; "
                         f"single-determinant strains are usually low-level/susceptible but ~12% were R).")
    return {
        "prediction": pred,
        "confidence": conf,
        "n_determinants": n,
        "determinants": dets,
        "resistance_threshold": resistance_threshold,
        "rule": f"amrfinder_curated_determinant_v1 (R iff >={resistance_threshold} curated drug-class determinants)",
        "caveat": ("interpretable deterministic call from AMRFinder's curated database. cipro N=147 "
                   "op-chars at threshold=2: acc 0.939 / sens 0.931 / spec 0.947 (≈ the POINT-XGB 0.943 "
                   "classifier, but a transparent rule). ~7% of R strains carry <2 detected determinants "
                   "(under-call). Acquired-gene drugs may need threshold=1." + boundary_note),
    }


def evaluate_cohort(runs_root: Path, accession_label_pairs: list[tuple[str, int]], drug: str,
                    resistance_threshold: int = _DEFAULT_THRESHOLD) -> dict:
    """Operating characteristics of the deterministic rule over a labelled cohort.
    `accession_label_pairs` = [(assembly_accession, 0|1), ...]. Returns counts + acc/sens/spec."""
    tp = fp = tn = fn = na = 0
    for acc, y in accession_label_pairs:
        c = call_resistance(Path(runs_root) / acc / "main.tsv", drug, resistance_threshold)
        if c["prediction"] == "INDETERMINATE":
            na += 1
            continue
        r = c["prediction"] == "R"
        if r and y == 1:
            tp += 1
        elif r and y == 0:
            fp += 1
        elif (not r) and y == 0:
            tn += 1
        else:
            fn += 1
    n = tp + fp + tn + fn
    return {
        "n": n, "na": na, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "accuracy": round((tp + tn) / n, 3) if n else None,
        "sensitivity": round(tp / (tp + fn), 3) if (tp + fn) else None,
        "specificity": round(tn / (tn + fp), 3) if (tn + fp) else None,
    }
