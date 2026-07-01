"""Validate the IrisPlex 6-SNP v0.1 (deployed Walsh/HIrisPlex rule) on OpenSNP — and quantify its
value-OVER-v0 (rs12913832 alone).

The wrapper-vs-underlying-tool rail in reverse: v0 (single-SNP) is the cheap baseline; v0.1 (6-SNP
deployed model) must EARN its complexity by recovering the heterozygotes v0 abstains on WITHOUT losing
accuracy. Headline = (a) accuracy/coverage on the full blue/brown set, (b) the rescue: how v0.1 does on
exactly the rs12913832-heterozygote subset v0 calls 'intermediate'.

Reuses the zip-native OpenSNP helpers (no 21GB extraction, no network) + the sourced IrisPlex coefficients.
HONESTY: same self-reported label as v0; complete-case (all 6 SNPs callable) — coverage reported, not hidden;
ancestry-confounded (European-calibrated model). Emits wiki/eye_colour_irisplex_v01_validation_<date>.json.
"""
from __future__ import annotations

import json
import sys
import zipfile
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.eye_colour import call_eye_colour  # noqa: E402
from dna_decode.data.eye_colour_irisplex import _SNP_ORDER, predict_irisplex  # noqa: E402
from scripts.eye_colour_opensnp_ingest import (  # noqa: E402
    DEFAULT_ZIP, _eye_colour_by_user, _find_phenotype_member, _genotype_members_by_user,
    _pick_member, rsids_from_member,
)
from scripts.eye_colour_opensnp_validate import bin_eye_colour  # noqa: E402


def _confusion(pairs: list[tuple[str, str]]) -> dict:
    """pairs = (label, pred) over blue/brown only; brown=positive."""
    tp = sum(1 for l, p in pairs if l == "brown" and p == "brown")
    fn = sum(1 for l, p in pairs if l == "brown" and p == "blue")
    tn = sum(1 for l, p in pairs if l == "blue" and p == "blue")
    fp = sum(1 for l, p in pairs if l == "blue" and p == "brown")
    n = tp + fn + tn + fp
    return {"TP": tp, "FP": fp, "TN": tn, "FN": fn, "n": n,
            "accuracy": round((tp + tn) / n, 3) if n else None,
            "brown_sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "blue_spec": round(tn / (tn + fp), 3) if (tn + fp) else None}


def run(zip_path: Path, limit: int | None = None) -> dict:
    if not zip_path.exists():
        return {"status": "ZIP_NOT_PRESENT", "zip": str(zip_path)}
    zf = zipfile.ZipFile(str(zip_path))
    pheno = _find_phenotype_member(zf)
    if not pheno:
        return {"status": "NO_PHENOTYPE_CSV"}
    eye = _eye_colour_by_user(zf, pheno)
    geno_by_uid = _genotype_members_by_user(zf)
    want = set(_SNP_ORDER)

    n_other = n_nogeno = n_incomplete = 0
    # full blue/brown set under v0.1 argmax (blue/brown only; intermediate-pred excluded from binary)
    v01_pairs: list[tuple[str, str]] = []
    v01_intermediate_pred = 0
    # paired comparison on COMPLETE-CASE users: v0 vs v0.1
    paired_v0: list[tuple[str, str]] = []
    paired_v01: list[tuple[str, str]] = []
    # the rescue subset: users v0 abstains on (rs12913832 heterozygote)
    rescue_total = 0
    rescue_v01_correct = 0
    rescue_breakdown = {"brown": {"brown": 0, "blue": 0, "intermediate": 0},
                        "blue": {"brown": 0, "blue": 0, "intermediate": 0}}
    # the DEPLOYED forensic operating mode: 0.7-threshold category (blue/brown/intermediate/undefined).
    # argmax always emits a colour; the 0.7 threshold is the rule the real IrisPlex workflow uses.
    v01_thresh_pairs: list[tuple[str, str]] = []
    n_thresh_undefined = 0
    scored = 0

    for uid, raw_colour in eye.items():
        label = bin_eye_colour(raw_colour)
        if label is None or label == "other":
            n_other += 1
            continue
        member = _pick_member(geno_by_uid.get(uid, []))
        if not member:
            n_nogeno += 1
            continue
        gts = rsids_from_member(zf, member, want)
        pred = predict_irisplex(gts)
        if pred["status"] != "PREDICTED":
            n_incomplete += 1
            continue
        scored += 1
        v0_call = call_eye_colour(gts["rs12913832"])["prediction"]   # single-SNP baseline
        v01_call = pred["prediction"]

        # full-set v0.1 binary (argmax)
        if v01_call in ("blue", "brown"):
            v01_pairs.append((label, v01_call))
        else:
            v01_intermediate_pred += 1
        # DEPLOYED 0.7-threshold mode: count only confident blue/brown; intermediate/undefined -> abstain
        cat = pred.get("category_at_0.7")
        if cat in ("blue", "brown"):
            v01_thresh_pairs.append((label, cat))
        else:
            n_thresh_undefined += 1
        # paired (complete-case) — both models on the same users; binary only when each gives blue/brown
        if v0_call in ("blue", "brown"):
            paired_v0.append((label, v0_call))
        if v01_call in ("blue", "brown"):
            paired_v01.append((label, v01_call))
        # the rescue: v0 abstained (intermediate)
        if v0_call == "intermediate":
            rescue_total += 1
            rescue_breakdown[label][v01_call if v01_call in ("blue", "brown", "intermediate") else "intermediate"] += 1
            if v01_call == label:
                rescue_v01_correct += 1
        if limit and scored >= limit:
            break

    return {
        "status": "SCORED" if scored else "NO_USERS_SCORED",
        "schema": "eye-colour-irisplex-v01-v1", "date": _date.today().isoformat(),
        "source": "OpenSNP archive dump (archive.org/opensnp_data_dumps, 2017-12-08)",
        "rule": "IrisPlex 6-SNP model form (Walsh/HIrisPlex); coefficients published-model-consistent, "
                "sourced via open impl (primary Walsh table NOT directly cross-checked)",
        "label_tier": "self-reported (near-independent, non-circular, noisy)",
        "n_complete_case_scored": scored,
        "n_other_excluded": n_other, "n_no_genotype_file": n_nogeno,
        "n_incomplete_6snp": n_incomplete,
        "v01_full_set_binary_argmax": _confusion(v01_pairs),
        "v01_predicted_intermediate": v01_intermediate_pred,
        "v01_deployed_threshold_0.7": {
            **_confusion(v01_thresh_pairs),
            "n_undefined_or_intermediate_abstained": n_thresh_undefined,
            "note": "the DEPLOYED forensic decision rule (0.7 threshold -> blue/brown/undefined); "
                    "abstains below 0.7. This is the like-for-like validation of the deployed rule; "
                    "argmax above is the more-permissive always-call variant.",
        },
        "paired_complete_case": {
            "v0_single_snp": _confusion(paired_v0),
            "v01_irisplex": _confusion(paired_v01),
            "note": "same complete-case users; coverage = n binary-called. v0.1 calls more (resolves hets).",
        },
        "heterozygote_rescue": {
            "definition": "users v0 abstains on (rs12913832 AG -> 'intermediate')",
            "n_v0_abstained": rescue_total,
            "n_v01_correct": rescue_v01_correct,
            "v01_recall_on_rescue": round(rescue_v01_correct / rescue_total, 3) if rescue_total else None,
            "breakdown_label_to_v01pred": rescue_breakdown,
        },
        "caveats": [
            "self-reported label (not a lab assay)",
            "complete-case: requires all 6 SNPs callable; coverage in n_complete_case_scored",
            "ancestry-confounded: IrisPlex European-calibrated (within-ancestry split = a further step)",
            "rs16891982 is C/G palindromic -> literal forward-strand assumption (bounded; named)",
            "coefficients SOURCED from brianbhsu/eye-color (published Walsh/HIrisPlex model), not fabricated",
        ],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args(argv)
    res = run(a.zip, limit=a.limit)
    if res.get("status") == "SCORED":
        out = REPO / "wiki" / f"eye_colour_irisplex_v01_validation_{_date.today().isoformat()}.json"
        out.write_text(json.dumps(res, indent=2), encoding="utf-8")
        print(f"[wrote {out}]")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
