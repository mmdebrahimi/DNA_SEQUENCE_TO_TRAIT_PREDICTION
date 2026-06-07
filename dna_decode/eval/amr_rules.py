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

from dna_decode.data.mic_tiers import CO_RESISTANCE_MECHANISMS, amrfinder_classes_for

# The pinned AMRFinderPlus image (tag encodes the DB date 2026-03-24.1) that produced the curated
# determinants this caller reads. Stamped into every call's provenance so an R/S verdict is auditable +
# reproducible against a known determinant source. Mirror of scripts/drug_mechanism_audit.AMRFINDER_IMAGE.
AMRFINDER_IMAGE_PINNED = "ncbi/amr:4.2.7-2026-03-24.1"

# Resistance mechanisms the curated-determinant rule CANNOT see: efflux overexpression, regulatory
# changes, and porin loss are expression/regulatory phenotypes, not acquired genes/point mutations in
# AMRFinder's curated set. A SUSCEPTIBLE (S) call therefore cannot RULE OUT resistance via these —
# they are the rule's structural blind spots (and the dominant cause of its false negatives).
UNDETECTABLE_MECHANISMS = sorted(CO_RESISTANCE_MECHANISMS)  # ['efflux', 'porin_loss', 'regulatory']


# Fluoroquinolone QRDR target-alteration genes (point mutations in these = the canonical cipro mechanism).
# Counting ONLY these POINT mutations is the CROSS-ORGANISM-robust cipro rule: it excludes intrinsic
# chromosomal efflux genes (e.g. K. pneumoniae oqxAB) that AMRFinder tags QUINOLONE-class but that are
# present in susceptible isolates too — those saturate the broad-determinant count and break transfer.
# Validated 2026-06-07: QRDR-POINT>=2 → Klebsiella cipro acc 1.000 (vs 0.577 broad) AND E. coli 0.925.
QRDR_GENES = ("gyrA", "gyrB", "parC", "parE")


def qrdr_point_determinants(main_tsv: Path) -> list[dict]:
    """Return the fluoroquinolone QRDR target-alteration POINT mutations (gyrA/gyrB/parC/parE_<mut>,
    AMRFinder Method POINT*) as determinant dicts for transparent reporting. The cross-organism-robust
    cipro signal — excludes intrinsic efflux/acquired genes that inflate the broad QUINOLONE-class count."""
    p = Path(main_tsv)
    out: list[dict] = []
    if not p.exists():
        return out
    for r in csv.DictReader(p.open(encoding="utf-8"), delimiter="\t"):
        sym = (r.get("Element symbol") or "")
        meth = (r.get("Method") or "").upper()
        if "POINT" in meth and any(sym == g or sym.startswith(g + "_") for g in QRDR_GENES):
            out.append({
                "symbol": sym.strip(),
                "name": (r.get("Element name") or "").strip(),
                "class": r.get("Class") or "", "subclass": r.get("Subclass") or "",
                "pct_identity": r.get("% Identity to reference"),
            })
    return out


def qrdr_point_count(main_tsv: Path) -> int:
    """Count fluoroquinolone QRDR target-alteration POINT mutations. Returns 0 if main.tsv absent."""
    return len(qrdr_point_determinants(main_tsv))


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
    "ciprofloxacin": {"threshold": 2, "subclass_any": None, "counter": "qrdr_point",
                      "validated": "E.coli N=147 acc 0.925/sens 0.875/spec 0.973; Klebsiella N=30 acc 1.0 "
                                   "(QRDR target POINT mutations gyrA/parC/parE >=2 — cross-organism-robust)"},
    "ceftriaxone":   {"threshold": 1, "subclass_any": frozenset({"CEPHALOSPORIN", "CARBAPENEM"}),
                      "validated": "N=60 acc 0.933/sens 0.962/spec 0.912 (extended-spectrum bla only)"},
    "tetracycline":  {"threshold": 1, "subclass_any": None, "gene_prefixes": ("tet",),
                      "validated": "E.coli N=12 acc 0.917/sens 1.0/spec 0.833; Klebsiella N=30 acc 0.8/spec 1.0/sens 0.6 (acquired tet* genes only — excludes intrinsic oqxAB efflux; sens-limited by efflux-mediated tet-R, a curated-determinant blind spot)"},
    "gentamicin":    {"threshold": 1, "subclass_any": frozenset({"GENTAMICIN"}),
                      "validated": "N=128 acc 0.945/sens 0.893/spec 0.96 (GENTAMICIN-subclass only; aph/aadA streptomycin-kanamycin genes excluded)"},
    "meropenem":     {"threshold": 1, "subclass_any": frozenset({"CARBAPENEM"}),
                      "validated": "Klebsiella N=30 acc 0.867/sens 1.0/spec 0.733 (acquired carbapenemase, CARBAPENEM-subclass: blaKPC/NDM/OXA-48; vs naive 0.533). Excludes ESBL/AmpC; blind to porin-loss-mediated R (expected FN mode)."},
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
      - ciprofloxacin: threshold 2 over QRDR TARGET POINT mutations (gyrA/parC/parE; counter='qrdr_point').
                       Cross-organism-robust — excludes intrinsic chromosomal efflux (e.g. Klebsiella OqxAB)
                       that broad QUINOLONE-class counting sweeps in. E.coli N=147 acc 0.925; Klebsiella 1.0.
      - ceftriaxone:   threshold 1 + EXTENDED-SPECTRUM subclass refinement (CEPHALOSPORIN/CARBAPENEM only;
                       plain BETA-LACTAM like blaTEM-1 is ampicillin-R not ceftriaxone-R) — N=60 acc 0.933.
      - tetracycline:  threshold 1 (acquired tet genes) — N=12 acc 0.833 (small N).
      - gentamicin:    threshold 1 + GENTAMICIN-subclass refinement — N=128 acc 0.945.
    Confidence: HIGH if n>=threshold+1 OR n==0; MODERATE near the boundary. Absent main.tsv → INDETERMINATE."""
    p = Path(main_tsv)
    if not p.exists():
        return {"prediction": "INDETERMINATE", "confidence": None, "n_determinants": 0,
                "determinants": [], "rule": "amrfinder_curated_determinant_v1", "caveat":
                "AMRFinder main.tsv not found for this genome"}
    cfg = rule_for(drug)
    if resistance_threshold is None:
        resistance_threshold = cfg["threshold"]
    # Counter selection: cipro uses the cross-organism-robust QRDR target-POINT-mutation count
    # (excludes intrinsic efflux/acquired genes); other drugs use the broad drug-class determinant
    # count with optional Subclass refinement.
    if cfg.get("counter") == "qrdr_point":
        dets = qrdr_point_determinants(p)
    else:
        dets = cipro_determinants_from_main(p, drug, subclass_any=cfg.get("subclass_any"))
        # gene_prefixes: keep only determinants whose SYMBOL starts with one of these (acquired-gene
        # refinement) — excludes intrinsic chromosomal genes that share the drug's broad Subclass but
        # don't confer acquired resistance (e.g. tet excludes K. pneumoniae oqxAB efflux). Cross-organism.
        prefixes = cfg.get("gene_prefixes")
        if prefixes:
            dets = [d for d in dets if d["symbol"].lower().startswith(tuple(p_.lower() for p_ in prefixes))]
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
    refine_note = (f" {drug} counts only determinants whose AMRFinder Subclass indicates "
                   f"{'/'.join(sorted(refine))} activity (broader-class genes that don't confer "
                   f"{drug} resistance are excluded)." if refine else "")
    # Blind-spot honesty: a SUSCEPTIBLE call cannot rule out resistance via mechanisms the curated-
    # determinant rule can't see (efflux overexpression / porin loss / regulatory). Surface them so a
    # negative is read as "no curated determinant found", NOT "definitely susceptible".
    undetectable = UNDETECTABLE_MECHANISMS if pred == "S" else []
    return {
        "prediction": pred,
        "confidence": conf,
        "n_determinants": n,
        "determinants": dets,
        "resistance_threshold": resistance_threshold,
        "undetectable_mechanisms": undetectable,
        "rule": (f"amrfinder_curated_determinant_v1 (R iff >={resistance_threshold} curated "
                 f"{'extended-spectrum ' if refine else ''}drug-class determinants)"),
        "caveat": ("interpretable deterministic call from AMRFinder's curated database. Per-drug validated "
                   f"op-chars: {cfg['validated']}." + refine_note + boundary_note
                   + (f" An S call cannot rule out resistance via {', '.join(UNDETECTABLE_MECHANISMS)} "
                      "(expression/regulatory mechanisms absent from AMRFinder's curated determinants)."
                      if pred == "S" else "")),
    }


# Genotype-vs-phenotype discordance taxonomy — the honest failure-mode buckets. Used by evaluate_cohort
# (and any caller with a ground-truth label) to characterize WHERE the deterministic rule fails, not just
# how often. This is the "failure-tolerant tool" deliverable: name the failure mode, don't hide it.
FN_UNDETECTED_MECHANISM = "FN_undetected_mechanism"          # true R, called S — efflux/porin/regulatory/low-level
FP_DETERMINANT_NO_PHENOTYPE = "FP_determinant_without_phenotype"  # true S, called R — label-noise/silent-or-low-expr/borderline-MIC


def discordance_bucket(prediction: str, true_label: int) -> str | None:
    """Classify a genotype-vs-phenotype mismatch into an honest failure-mode bucket.

    Returns None when concordant or INDETERMINATE. FN (R phenotype missed) → the rule's structural blind
    spots (UNDETECTABLE_MECHANISMS); FP (R called on a susceptible isolate) → determinant present but not
    phenotypically expressed (label noise / silent or low-expression gene / MIC near the breakpoint)."""
    if prediction == "INDETERMINATE":
        return None
    pred_r = prediction == "R"
    if pred_r == bool(true_label):
        return None
    return FN_UNDETECTED_MECHANISM if (not pred_r and true_label == 1) else FP_DETERMINANT_NO_PHENOTYPE


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
        # honest failure-mode breakdown (discordance taxonomy)
        "discordance": {
            FN_UNDETECTED_MECHANISM: fn,        # R missed — likely efflux/porin-loss/regulatory/low-level
            FP_DETERMINANT_NO_PHENOTYPE: fp,    # called R but S — likely label-noise/low-expression/borderline-MIC
        },
    }
