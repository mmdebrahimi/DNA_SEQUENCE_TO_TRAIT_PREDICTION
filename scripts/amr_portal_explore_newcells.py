"""EXPLORATORY new-cell curation on the EBI AMR Portal — Tier 3/4 of the unscored-cell triage.

Per-cell curation, NOT a sweep (`wiki/amr_portal_unscored_triage_2026-06-28.md`). The intrinsic-gene
guardrail: do NOT spray the generic E. coli `DRUG_RULE` across organisms. The ONE legitimate cross-organism
transfer is **ciprofloxacin via `qrdr_point`** — its rule counts gyrA/parC/parE TARGET-SITE point mutations,
a CONSERVED mechanism across Gram-negatives (the cipro DRUG_RULE's own note says "cross-organism-robust";
Campylobacter cipro already scored 1.0 this way). Acquired-gene drugs (tet/cef/gent/mero) are NOT transferred
here — their determinants vary with each organism's intrinsic flora (a separate per-cell curation, or ABSTAIN).

This scores cipro-via-qrdr_point on a CURATED candidate set + honestly STRESS-TESTS where it should NOT
transfer (Pseudomonas efflux-driven R; Staph Gram-positive grlA gene-naming). Output is a NAMESPACE-SEPARATE
EXPLORATORY artifact (`wiki/amr_portal_newcell_exploratory_<date>.json`) — these organisms have NO deployed
claim, so they do NOT enter the deployed `amr_portal_independent` card (scope contract + shared-key trap).
A cell that scores cleanly is a CANDIDATE for a frozen-registry promotion (a user-ratified freeze amendment),
NOT auto-promoted. FROZEN AMR surface untouched (read-only `call_resistance`).
"""
from __future__ import annotations

import json
import sys
import tempfile
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.eval.amr_rules import call_resistance  # noqa: E402
from scripts.amr_portal_score_independent import (  # noqa: E402
    DEFAULT_CRYPTIC, DEFAULT_PARQUET, GENO_PARQUET,
    _load_genotype, _load_leak_set, _load_phenotype_cells, confusion, genotype_to_main_tsv, wilson_ci,
)

# Curated candidate cells: (AMR-Portal organism, drug, call_resistance organism=, mechanism rationale).
# organism=None -> the DRUG_RULE default (cipro = qrdr_point). All cipro here (the conserved-QRDR transfer).
CANDIDATES = [
    ("Neisseria gonorrhoeae", "ciprofloxacin", None, "QRDR gyrA S91F/D95 + parC conserved; huge N"),
    ("Enterobacter cloacae", "ciprofloxacin", None, "Enterobacterale; QRDR gyrA/parC conserved (E.coli-like)"),
    ("Enterobacter hormaechei", "ciprofloxacin", None, "Enterobacterale; QRDR conserved"),
    ("Serratia marcescens", "ciprofloxacin", None, "Enterobacterale; QRDR conserved (small N)"),
    # --- honest guardrail stress-tests (qrdr_point expected to UNDER-call) ---
    ("Acinetobacter baumannii", "ciprofloxacin", None, "QRDR gyrA S83L conserved BUT S-poor cohort (221S)"),
    ("Pseudomonas aeruginosa", "ciprofloxacin", None, "STRESS: efflux(MexAB)-driven R w/o 2 QRDR muts -> under-call expected"),
    ("Staphylococcus aureus", "ciprofloxacin", None, "STRESS: Gram-positive grlA/gyrA gene-naming differs -> may under-call"),
]


def _verdict(sens, spec) -> str:
    if sens is None or spec is None:
        return "UNSCORABLE"
    balacc = (sens + spec) / 2
    if balacc >= 0.85 and sens >= 0.70 and spec >= 0.70:
        return "TRANSFERS_CLEAN"        # the conserved rule transfers -> promotion candidate
    if balacc >= 0.70:
        return "TRANSFERS_PARTIAL"       # imperfect (under-call) -> organism-specific curation could lift
    return "NEEDS_CURATION"              # conserved rule does NOT transfer -> curate or ABSTAIN


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pheno", type=Path, default=DEFAULT_PARQUET)
    ap.add_argument("--geno", type=Path, default=GENO_PARQUET)
    ap.add_argument("--cryptic", type=Path, default=DEFAULT_CRYPTIC)
    a = ap.parse_args(argv)
    for p in (a.pheno, a.geno):
        if not p.exists():
            print(f"ERROR: AMR Portal parquet not found: {p}", file=sys.stderr)
            return 2
    cells = [(o, d) for o, d, _, _ in CANDIDATES]
    routing = {(o, d): reg for o, d, reg, _ in CANDIDATES}
    rationale = {(o, d): why for o, d, _, why in CANDIDATES}
    leak = _load_leak_set(a.cryptic)
    geno = _load_genotype(a.geno, {o for o, _ in cells})
    pheno = _load_phenotype_cells(a.pheno, cells, leak)
    tmp = Path(tempfile.gettempdir()) / "_amrportal_explore.tsv"
    out_cells = {}
    for (org, drug), isolates in pheno.items():
        conf = {"TP": 0, "FP": 0, "TN": 0, "FN": 0}
        n_indet = 0
        for bs, leaked, sir in isolates:
            if leaked:
                continue
            tmp.write_text(genotype_to_main_tsv(geno.get(bs, [])), encoding="utf-8")
            pred = call_resistance(tmp, drug, organism=routing[(org, drug)])["prediction"]
            if pred not in ("R", "S"):
                n_indet += 1
                continue
            c = confusion(pred, sir)
            if c:
                conf[c] += 1
        tp, fp, tn, fn = conf["TP"], conf["FP"], conf["TN"], conf["FN"]
        nR, nS = tp + fn, tn + fp
        sens = round(tp / nR, 3) if nR else None
        spec = round(tn / nS, 3) if nS else None
        acc = round((tp + tn) / (tp + tn + fp + fn), 3) if (tp + tn + fp + fn) else None
        out_cells[f"{org}|{drug}"] = {
            "organism": org, "drug": drug, "rule": "qrdr_point (DRUG_RULE default)",
            "rationale": rationale[(org, drug)], "n_R": nR, "n_S": nS, "n_indeterminate": n_indet,
            "sens": sens, "spec": spec, "accuracy": acc,
            "sens_wilson95": wilson_ci(tp, nR), "spec_wilson95": wilson_ci(tn, nS),
            "powered": nR >= 10 and nS >= 10, "verdict": _verdict(sens, spec),
        }
    clean = [k for k, c in out_cells.items() if c["verdict"] == "TRANSFERS_CLEAN" and c["powered"]]
    out = {
        "_schema": "amr-portal-newcell-exploratory-v1", "date": _date.today().isoformat(),
        "tier": "exploratory-new-cell-curation (NOT a deployed claim; NOT in the deployed independent card)",
        "transfer_basis": ("ciprofloxacin qrdr_point counts CONSERVED gyrA/parC/parE target-site point "
                           "mutations -> the one drug whose rule legitimately transfers across Gram-negatives "
                           "(the intrinsic-gene guardrail blocks acquired-gene drugs from such transfer)"),
        "promotion_note": ("a TRANSFERS_CLEAN powered cell is a CANDIDATE for a calibrated-registry promotion "
                           "= a USER-ratified frozen-surface amendment (sha-pinned file); NOT auto-promoted"),
        "n_candidates": len(out_cells), "n_transfers_clean_powered": len(clean),
        "cells": out_cells, "frozen_surface_changed": False,
    }
    outp = REPO / "wiki" / f"amr_portal_newcell_exploratory_{_date.today().isoformat()}.json"
    outp.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"{'cell':<40} {'nR':>5} {'nS':>5} {'sens':>6} {'spec':>6} {'acc':>6} verdict")
    for k, c in out_cells.items():
        s = lambda x: f"{x:.3f}" if isinstance(x, float) else "  -  "
        print(f"{k:<40} {c['n_R']:>5} {c['n_S']:>5} {s(c['sens']):>6} {s(c['spec']):>6} {s(c['accuracy']):>6} "
              f"{'POW ' if c['powered'] else 'underp '}{c['verdict']}")
    print(f"\n{len(clean)}/{len(out_cells)} TRANSFERS_CLEAN+powered -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
