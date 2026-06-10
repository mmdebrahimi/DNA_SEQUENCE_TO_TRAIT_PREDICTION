"""Step 6 of the expression_context plan — independent-cohort eval overlay + promotion gate.

Runs the PRIMARY expression_context detector (the exact frozen falsifier rule) on the INDEPENDENT
Acinetobacter meropenem cohort built by build_acinetobacter_indep_cohort.py, resolving each assembly FASTA
from the refseq cache BY VERSIONED ACCESSION (provenance-safe, not run-dir adjacency). Tabulates the
right-shaped endpoint for an ABSTAIN-override: rescue-rate + S-upgrades + Wilson 95% upper bound on the
false-upgrade rate. For Acinetobacter|meropenem every baseline call is ABSTAIN (EXPRESSION_FLOOR), so the
whole cohort IS the prior-ABSTAIN set.

PROMOTE iff s_upgrades == 0 AND r_rescues >= 1 AND n_S >= 15 (the r_rescues clause guards against an inert
never-firing detector trivially passing the zero-S-upgrade test). Eval-only: does NOT flip the registry.

Usage: .venv/Scripts/python.exe scripts/expression_context_validate.py
Exit 0 = PROMOTE, 1 = HOLD, 2 = cohort/refs missing.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data import refseq
from dna_decode.eval.expression_context import detect_is_upstream_junction

IS_REF = Path("data/isaba1_ref/ISAba1_ref.fna")
TARGET_REF = Path("data/isaba1_ref/OXA51fam_ref.fna")
COHORT = Path("data/raw/acinetobacter_meropenem_indep/selected.tsv")
GCACHE = Path("data/raw/acinetobacter_meropenem_indep/refseq")
MIN_S = 15


def wilson_upper(successes: int, n: int, z: float = 1.96) -> float | None:
    """Wilson score interval UPPER bound for a binomial proportion. None if n==0."""
    if n == 0:
        return None
    phat = successes / n
    denom = 1 + z * z / n
    center = phat + z * z / (2 * n)
    margin = z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n))
    return round(min(1.0, (center + margin) / denom), 4)


def compute_rescue_endpoint(signals: list[bool], labels: list[str]) -> dict:
    """Pure endpoint from per-strain (signal, label). labels are 'R'/'S' (the true AST phenotype).

    Baseline for every strain is ABSTAIN (EXPRESSION_FLOOR), so a signal=True is an ABSTAIN->R upgrade.
    r_rescues = true-R strains the detector fires on (correct upgrades);
    s_upgrades = true-S strains the detector fires on (false upgrades — must be 0).
    PROMOTE iff s_upgrades==0 AND r_rescues>=1 AND n_S>=MIN_S, else HOLD with the binding reason.
    """
    assert len(signals) == len(labels), "signals/labels length mismatch"
    n_R = sum(1 for y in labels if y == "R")
    n_S = sum(1 for y in labels if y == "S")
    r_rescues = sum(1 for s, y in zip(signals, labels) if s and y == "R")
    s_upgrades = sum(1 for s, y in zip(signals, labels) if s and y == "S")
    reasons = []
    if s_upgrades != 0:
        reasons.append(f"s_upgrades={s_upgrades} (must be 0)")
    if r_rescues < 1:
        reasons.append("r_rescues=0 (detector inert — no true-R upgrade)")
    if n_S < MIN_S:
        reasons.append(f"n_S={n_S} (<{MIN_S})")
    verdict = "PROMOTE" if not reasons else "HOLD"
    return {
        "n_R": n_R, "n_S": n_S, "r_rescues": r_rescues, "s_upgrades": s_upgrades,
        "abstain_rescue_rate": round(r_rescues / n_R, 3) if n_R else None,
        "false_upgrade_rate": round(s_upgrades / n_S, 3) if n_S else None,
        "false_upgrade_wilson95_upper": wilson_upper(s_upgrades, n_S),
        "verdict": verdict, "hold_reasons": reasons,
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def main() -> int:
    if not COHORT.exists():
        print(f"cohort not found at {COHORT}; run build_acinetobacter_indep_cohort.py first")
        return 2
    if not IS_REF.exists() or not TARGET_REF.exists():
        print("ISAba1 / OXA refs missing"); return 2
    labels_by_acc = {}
    for ln in COHORT.read_text().splitlines():
        if "\t" in ln:
            a, rs = ln.split("\t"); labels_by_acc[a] = rs.strip()

    rows, signals, labels = [], [], []
    for acc, lab in labels_by_acc.items():
        try:
            fasta = refseq.fasta_path(acc, GCACHE)
        except Exception:
            fasta = None
        if not fasta or not Path(fasta).exists():
            rows.append({"acc": acc, "label": lab, "skipped": "no assembly FASTA"}); continue
        res = detect_is_upstream_junction(fasta, is_ref=IS_REF, target_ref=TARGET_REF)
        if res["status"] != "ok":
            rows.append({"acc": acc, "label": lab, "skipped": res.get("reason")}); continue
        ev = res["evidence"]
        signals.append(bool(res["signal"])); labels.append(lab)
        rows.append({"acc": acc, "label": lab, "signal": res["signal"],
                     "n_is_hits": ev["n_is_hits"], "n_target_hits": ev["n_target_hits"],
                     "junction": ev["junction"], "fasta_sha256": _sha256(fasta)})

    endpoint = compute_rescue_endpoint(signals, labels)
    artifact = {
        "_schema": "expression-context-validation-v1", "date": _date.today().isoformat(),
        "organism": "Acinetobacter_baumannii", "drug": "meropenem",
        "detector": "expression_context primary (frozen falsifier rule; same-contig + OXA-upstream 400bp; no IS-orientation)",
        "refs": {"is_ref": str(IS_REF), "is_ref_sha256": _sha256(IS_REF),
                 "target_ref": str(TARGET_REF), "target_ref_sha256": _sha256(TARGET_REF)},
        "n_evaluated": len(labels), "endpoint": endpoint, "per_strain": rows,
    }
    out_json = Path(f"wiki/expression_context_acinetobacter_validation_{_date.today().isoformat()}.json")
    out_md = out_json.with_suffix(".md")
    out_json.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    ep = endpoint
    md = (f"# expression_context independent validation — Acinetobacter x meropenem — {_date.today().isoformat()}\n\n"
          f"PRIMARY detector (frozen falsifier rule) on the INDEPENDENT cohort (disjoint from the cached 30).\n\n"
          f"## Verdict: {ep['verdict']}\n\n"
          f"| metric | value |\n|---|---|\n"
          f"| n evaluated | {len(labels)} ({ep['n_R']}R/{ep['n_S']}S) |\n"
          f"| r_rescues (true-R upgraded ABSTAIN->R) | {ep['r_rescues']} |\n"
          f"| s_upgrades (false-R; must be 0) | {ep['s_upgrades']} |\n"
          f"| abstain_rescue_rate | {ep['abstain_rescue_rate']} |\n"
          f"| false-upgrade rate | {ep['false_upgrade_rate']} |\n"
          f"| false-upgrade Wilson95 upper | {ep['false_upgrade_wilson95_upper']} |\n\n"
          f"Gate: PROMOTE iff s_upgrades==0 AND r_rescues>=1 AND n_S>=15. "
          f"{'PROMOTE -> eligible for experimental opt-in (NOT default-on).' if ep['verdict']=='PROMOTE' else 'HOLD: ' + '; '.join(ep['hold_reasons'])}\n\n"
          f"Promotion is opt-in only; default-on deferred until n_S materially larger (the Wilson upper bound "
          f"is still wide at small n_S). Raw BLAST hits + FASTA/ref SHA256 in the JSON sidecar for reproducibility.\n")
    out_md.write_text(md, encoding="utf-8")
    print(f"VERDICT: {ep['verdict']}  (rescues={ep['r_rescues']}, s_upgrades={ep['s_upgrades']}, "
          f"n_S={ep['n_S']}, Wilson95-upper={ep['false_upgrade_wilson95_upper']})")
    print(f"Artifacts: {out_md} + {out_json}")
    return 0 if ep["verdict"] == "PROMOTE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
