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


def cipro_determinants_from_main(main_tsv: Path, drug: str,
                                 subclass_any: frozenset | None = None) -> list[dict]:
    """Return the curated AMRFinder determinants (main.tsv rows) relevant to `drug`.

    A row is drug-relevant iff any token of its Class or Subclass matches the drug's AMRFinder class set
    (mic_tiers.amrfinder_classes_for) — case-insensitive substring on the '/'-joined multi-class strings
    AMRFinder emits (e.g. 'AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN'). Each returned dict carries the
    element symbol + name + class for transparent reporting.

    `subclass_any` (drug-specific REFINEMENT, e.g. ceftriaxone): when provided, a class-matched row is
    kept ONLY if its Subclass contains one of these tokens. This is how ceftriaxone separates EXTENDED-
    spectrum determinants (Subclass CEPHALOSPORIN/CARBAPENEM — real 3rd-gen-cephalosporin resistance:
    blaCTX-M, blaCMY, blaNDM…) from intrinsic/narrow-spectrum BETA-LACTAM (blaTEM-1, blaEC) that confer
    ampicillin-R but NOT ceftriaxone-R. Validated: the broad class match gave cef spec 0.41; the
    extended-spectrum refinement gives acc 0.933 / sens 0.962 / spec 0.912 on N=60 (see DRUG_RULE)."""
    classes = {c.upper() for c in amrfinder_classes_for(drug)}
    refine = {t.upper() for t in subclass_any} if subclass_any else None
    out: list[dict] = []
    p = Path(main_tsv)
    if not p.exists():
        return out
    for r in csv.DictReader(p.open(encoding="utf-8"), delimiter="\t"):
        cls = (r.get("Class") or "").upper()
        sub = (r.get("Subclass") or "").upper()
        if not any(c in cls or c in sub for c in classes):
            continue
        if refine is not None and not any(t in sub for t in refine):
            continue
        out.append({
            "symbol": (r.get("Element symbol") or "").strip(),
            "name": (r.get("Element name") or "").strip(),
            "class": r.get("Class") or "", "subclass": r.get("Subclass") or "",
            "pct_identity": r.get("% Identity to reference"),
        })
    return out


_DEFAULT_THRESHOLD = 2  # cipro-validated: clinical fluoroquinolone-R needs multiple QRDR hits (see below)

# Per-drug deterministic-rule config (validated on cached BV-BRC cohorts, 2026-06-06). Each entry:
#   threshold          : min #curated determinants for an R call
#   subclass_any       : Subclass-refinement tokens (None = match on the broad drug class set)
#   validated          : op-char provenance string
# Drugs absent here fall back to (_DEFAULT_THRESHOLD, no refinement).
DRUG_RULE: dict[str, dict] = {
    "ciprofloxacin": {"threshold": 2, "subclass_any": None,
                      "validated": "N=147 acc 0.939/sens 0.931/spec 0.947 (QRDR point-mutations: >=2 hits)"},
    "ceftriaxone":   {"threshold": 1, "subclass_any": frozenset({"CEPHALOSPORIN", "CARBAPENEM"}),
                      "validated": "N=60 acc 0.933/sens 0.962/spec 0.912 (extended-spectrum bla only)"},
    "tetracycline":  {"threshold": 1, "subclass_any": None,
                      "validated": "N=12 acc 0.833/sens 1.0/spec 0.667 (acquired tet genes; small N)"},
    "gentamicin":    {"threshold": 1, "subclass_any": None,
                      "validated": "NOT VALIDATED — no cohort yet (acquired aminoglycoside-modifying genes; threshold=1 by mechanism analogy)"},
}


def rule_for(drug: str) -> dict:
    """Return the per-drug rule config (threshold + subclass refinement + provenance)."""
    return DRUG_RULE.get(drug.lower(), {"threshold": _DEFAULT_THRESHOLD, "subclass_any": None,
                                        "validated": "unconfigured drug — default threshold=2, no refinement"})


def call_resistance(main_tsv: Path, drug: str, resistance_threshold: int | None = None) -> dict:
    """Deterministic, interpretable R/S call for one genome's AMRFinder main.tsv.

    Rule (tiered, count-based): R iff #curated drug-relevant determinants >= threshold. When
    `resistance_threshold` is None the PER-DRUG validated config (DRUG_RULE) supplies both the threshold
    and any Subclass refinement; pass an int to override the threshold explicitly. Per-drug op-chars:
      - ciprofloxacin: threshold 2 (QRDR point-mutations need >=2 hits) — N=147 acc 0.939.
      - ceftriaxone:   threshold 1 + EXTENDED-SPECTRUM subclass refinement (CEPHALOSPORIN/CARBAPENEM only;
                       plain BETA-LACTAM like blaTEM-1 is ampicillin-R not ceftriaxone-R) — N=60 acc 0.933.
      - tetracycline:  threshold 1 (acquired tet genes) — N=12 acc 0.833 (small N).
      - gentamicin:    threshold 1 by mechanism analogy — NOT yet cohort-validated.
    Confidence: HIGH if n>=threshold+1 OR n==0; MODERATE near the boundary. Absent main.tsv → INDETERMINATE."""
    p = Path(main_tsv)
    if not p.exists():
        return {"prediction": "INDETERMINATE", "confidence": None, "n_determinants": 0,
                "determinants": [], "rule": "amrfinder_curated_determinant_v1", "caveat":
                "AMRFinder main.tsv not found for this genome"}
    cfg = rule_for(drug)
    if resistance_threshold is None:
        resistance_threshold = cfg["threshold"]
    dets = cipro_determinants_from_main(p, drug, subclass_any=cfg["subclass_any"])
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
    refine = cfg["subclass_any"]
    refine_note = (f" ceftriaxone uses an extended-spectrum subclass refinement "
                   f"({'/'.join(sorted(refine))}); plain BETA-LACTAM (blaTEM-1) is excluded as "
                   f"ampicillin-R-not-ceftriaxone-R." if refine else "")
    return {
        "prediction": pred,
        "confidence": conf,
        "n_determinants": n,
        "determinants": dets,
        "resistance_threshold": resistance_threshold,
        "rule": (f"amrfinder_curated_determinant_v1 (R iff >={resistance_threshold} curated "
                 f"{'extended-spectrum ' if refine else ''}drug-class determinants)"),
        "caveat": ("interpretable deterministic call from AMRFinder's curated database. Per-drug validated "
                   f"op-chars: {cfg['validated']}." + refine_note + boundary_note),
    }


def evaluate_cohort(runs_root: Path, accession_label_pairs: list[tuple[str, int]], drug: str,
                    resistance_threshold: int | None = None) -> dict:
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
