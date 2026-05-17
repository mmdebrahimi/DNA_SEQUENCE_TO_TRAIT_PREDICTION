"""Merge cipro mechanism audit (Experiment 1) x MIC audit (Experiment 2) into
one cipro_mechanism_phenotype_audit table.

Codex's cross-question integration plan step 3: "Merge into one
cipro_mechanism_phenotype_audit table" — the load-bearing artifact for
deciding whether the N=38 labels are clean enough for any model conclusion.

Produces a per-strain row with:
- cohort_binary_label (R/S)
- MIC tier (HIGH_R / BORDERLINE / NO_MIC / ...)
- median MIC
- primary AMRFinder mechanism class
- list of mechanisms found
- decisive flag (HIGH_R + QRDR or HIGH_S + no_silent_mechanism)
- noise flag (BORDERLINE label OR R-label without any known mechanism OR S-label with silent mechanism)
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from datetime import date as _date
from pathlib import Path


def _classify_noise(row: dict) -> str:
    """Single-row noise classification."""
    label = row.get("cohort_binary_label")
    tier = row.get("mic_tier", "?")
    mechs = row.get("mechanisms_present", [])

    if label == 1:
        if tier == "HIGH_R" and "QRDR_target_alteration" in mechs:
            return "CLEAN_R_QRDR"
        if tier == "HIGH_R" and mechs:
            return "CLEAN_R_nonQRDR"
        if tier == "HIGH_R" and not mechs:
            return "SUSPECT_R_no_mechanism"  # high MIC but AMRFinder finds nothing
        if tier in {"BORDERLINE", "AMBIGUOUS", "CONFLICT"}:
            return "NOISY_R_borderline"
        if tier == "NO_MIC":
            return "NOISY_R_no_mic"
        return "OTHER_R"
    else:
        if tier == "HIGH_S" and not mechs:
            return "CLEAN_S_no_mechanism"
        if tier == "HIGH_S" and mechs:
            return "SUSPECT_S_silent_mechanism"
        if mechs and tier != "HIGH_S":
            return "SUSPECT_S_borderline_mechanism"  # S with mechanism + borderline MIC = likely mislabeled
        if tier in {"BORDERLINE", "AMBIGUOUS", "CONFLICT"}:
            return "NOISY_S_borderline"
        if tier == "NO_MIC":
            return "NOISY_S_no_mic"
        return "OTHER_S"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mech-audit", type=Path, default=Path("wiki/cipro_mechanism_audit_2026-05-17.json"))
    parser.add_argument("--mic-audit", type=Path, default=Path("wiki/cipro_mic_audit_2026-05-17.json"))
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = Path(f"wiki/cipro_mechanism_phenotype_audit_{_date.today().isoformat()}.md")

    if not args.mech_audit.exists():
        print(f"[merge] FAIL mech audit not found: {args.mech_audit}")
        return 2
    if not args.mic_audit.exists():
        print(f"[merge] FAIL mic audit not found: {args.mic_audit}")
        return 2

    mech = json.loads(args.mech_audit.read_text(encoding="utf-8"))
    mic = json.loads(args.mic_audit.read_text(encoding="utf-8"))

    # Build per-strain merged rows
    mic_by_strain = {r["strain_id"]: r for r in mic["per_strain"]}
    mech_by_strain = {r["strain_id"]: r for r in mech["per_strain"]}

    all_ids = sorted(set(mic_by_strain) | set(mech_by_strain))
    merged: list[dict] = []
    for sid in all_ids:
        mic_r = mic_by_strain.get(sid, {})
        mech_r = mech_by_strain.get(sid, {})
        row = {
            "strain_id": sid,
            "accession": mech_r.get("accession") or mic_r.get("accession", ""),
            "mlst": mech_r.get("mlst") or mic_r.get("mlst", ""),
            "cohort_binary_label": mic_r.get("cohort_binary_label", mech_r.get("cohort_label")),
            "cohort_binary_RS": "R" if (mic_r.get("cohort_binary_label", mech_r.get("cohort_label")) == 1) else "S",
            "mic_tier": mic_r.get("tier", "?"),
            "median_mic": (mic_r.get("detail", {}) or {}).get("median_mic"),
            "n_mic_rows": mic_r.get("n_mic_rows", 0),
            "n_ast_rows": mic_r.get("n_ast_rows", 0),
            "clsi_call": (mic_r.get("detail", {}) or {}).get("clsi_call"),
            "eucast_call": (mic_r.get("detail", {}) or {}).get("eucast_call"),
            "mechanism_status": mech_r.get("status", "MISSING"),
            "primary_mechanism": mech_r.get("primary_mechanism_class", "MISSING"),
            "mechanisms_present": mech_r.get("mechanisms_present", []),
            "mech_hits": mech_r.get("mech_hits", {}),
            "n_mech_hits": mech_r.get("n_hits", 0),
        }
        row["noise_class"] = _classify_noise(row)
        merged.append(row)

    # Aggregates
    noise_counts: dict[str, int] = defaultdict(int)
    for r in merged:
        noise_counts[r["noise_class"]] += 1
    r_strains = [r for r in merged if r["cohort_binary_label"] == 1]
    s_strains = [r for r in merged if r["cohort_binary_label"] == 0]

    clean_r = sum(1 for r in r_strains if r["noise_class"].startswith("CLEAN_R"))
    clean_s = sum(1 for r in s_strains if r["noise_class"].startswith("CLEAN_S"))
    suspect_r = sum(1 for r in r_strains if r["noise_class"].startswith("SUSPECT_R"))
    suspect_s = sum(1 for r in s_strains if r["noise_class"].startswith("SUSPECT_S"))
    noisy_r = sum(1 for r in r_strains if r["noise_class"].startswith("NOISY_R"))
    noisy_s = sum(1 for r in s_strains if r["noise_class"].startswith("NOISY_S"))

    cohort_signal_quality = (clean_r + clean_s) / max(1, len(merged))

    verdict = (
        "SIGNAL_DOMINATES" if cohort_signal_quality >= 0.7 else
        "MIXED" if cohort_signal_quality >= 0.4 else
        "NOISE_DOMINATES"
    )

    print(f"[merge] N={len(merged)} merged rows")
    print(f"\n=== Noise class distribution ===")
    for cls, n in sorted(noise_counts.items(), key=lambda x: -x[1]):
        print(f"  {cls:30s} {n}")
    print(f"\nClean (R+S): {clean_r}R + {clean_s}S = {clean_r + clean_s}")
    print(f"Suspect (R+S): {suspect_r}R + {suspect_s}S = {suspect_r + suspect_s}")
    print(f"Noisy (R+S):  {noisy_r}R + {noisy_s}S = {noisy_r + noisy_s}")
    print(f"Cohort signal quality: {cohort_signal_quality:.2f}")
    print(f"Verdict: {verdict}")

    # JSON sidecar
    json_path = args.output.with_suffix(".json")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "mech_audit": str(args.mech_audit),
        "mic_audit": str(args.mic_audit),
        "n_merged": len(merged),
        "noise_class_counts": dict(noise_counts),
        "clean_count": clean_r + clean_s,
        "suspect_count": suspect_r + suspect_s,
        "noisy_count": noisy_r + noisy_s,
        "cohort_signal_quality": cohort_signal_quality,
        "verdict": verdict,
        "per_strain": merged,
    }
    json_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(f"[merge] wrote JSON: {json_path}")

    # Markdown packet
    lines = [
        f"# Cipro mechanism x phenotype audit — N={len(merged)} cohort ({_date.today().isoformat()})",
        "",
        f"**Purpose:** the load-bearing combined audit table from Codex's cross-question integration plan. Joins AMRFinder mechanism classification (Experiment 1) x MIC noise tiering (Experiment 2) to decide whether the N=38 cipro labels are clean enough for any genomic predictor to fit.",
        "",
        f"## Verdict: **{verdict}**",
        f"- Cohort signal quality (CLEAN / total): **{cohort_signal_quality:.2f}**",
        f"- Clean strains: {clean_r}R + {clean_s}S = **{clean_r + clean_s}**",
        f"- Suspect strains: {suspect_r}R + {suspect_s}S = **{suspect_r + suspect_s}** (likely mislabeled or atypical biology)",
        f"- Noisy strains: {noisy_r}R + {noisy_s}S = **{noisy_r + noisy_s}** (label is structurally unreliable)",
        "",
        "## Noise class distribution",
        "",
        "| class | count | meaning |",
        "|---|---:|---|",
        f"| CLEAN_R_QRDR | {sum(1 for r in r_strains if r['noise_class']=='CLEAN_R_QRDR')} | HIGH_R MIC + textbook gyrA/parC/parE mutations |",
        f"| CLEAN_R_nonQRDR | {sum(1 for r in r_strains if r['noise_class']=='CLEAN_R_nonQRDR')} | HIGH_R MIC + plasmid/efflux/regulatory mechanism |",
        f"| SUSPECT_R_no_mechanism | {sum(1 for r in r_strains if r['noise_class']=='SUSPECT_R_no_mechanism')} | HIGH_R MIC but no AMRFinder hits — novel mechanism or AMRFinder gap |",
        f"| NOISY_R_borderline | {sum(1 for r in r_strains if r['noise_class']=='NOISY_R_borderline')} | borderline/ambiguous MIC — label may be wrong |",
        f"| NOISY_R_no_mic | {sum(1 for r in r_strains if r['noise_class']=='NOISY_R_no_mic')} | no MIC in BV-BRC — label is opaque |",
        f"| CLEAN_S_no_mechanism | {sum(1 for r in s_strains if r['noise_class']=='CLEAN_S_no_mechanism')} | HIGH_S MIC + no cipro mechanism — clean susceptible |",
        f"| SUSPECT_S_silent_mechanism | {sum(1 for r in s_strains if r['noise_class']=='SUSPECT_S_silent_mechanism')} | HIGH_S MIC but mechanism present — silent / non-functional |",
        f"| SUSPECT_S_borderline_mechanism | {sum(1 for r in s_strains if r['noise_class']=='SUSPECT_S_borderline_mechanism')} | borderline MIC + mechanism — likely mislabeled to S |",
        f"| NOISY_S_borderline | {sum(1 for r in s_strains if r['noise_class']=='NOISY_S_borderline')} | borderline MIC with no mechanism — label opaque |",
        f"| NOISY_S_no_mic | {sum(1 for r in s_strains if r['noise_class']=='NOISY_S_no_mic')} | no MIC in BV-BRC — label opaque |",
        "",
        "## Per-strain merged table",
        "",
        "| strain_id | accession | label | mic_tier | med MIC | primary mech | mechs | noise_class | mlst |",
        "|---|---|---|---|---:|---|---|---|---|",
    ]
    for r in sorted(merged, key=lambda x: (x["cohort_binary_label"], x["noise_class"], x["strain_id"])):
        med_str = f"{r['median_mic']:.3f}" if r['median_mic'] is not None else "-"
        mechs_str = ",".join(r["mechanisms_present"]) if r["mechanisms_present"] else "-"
        lines.append(
            f"| {r['strain_id']} | {r['accession']} | {r['cohort_binary_RS']} | {r['mic_tier']} | "
            f"{med_str} | {r['primary_mechanism']} | {mechs_str} | {r['noise_class']} | {r['mlst']} |"
        )
    lines.extend([
        "",
        "## How to use",
        "",
        "- **SIGNAL_DOMINATES (>=0.70 clean):** the cohort is clean enough that NT's FAIL on cipro is a genuine model issue, not label noise. Pivot to curated baseline / Bakta / per-gene NT.",
        "- **MIXED (0.40-0.70 clean):** label noise is a real confounder. Re-run Stage 1 + curated baseline on the CLEAN subset only; if AUROC improves substantially, label noise was a meaningful contributor.",
        "- **NOISE_DOMINATES (<0.40 clean):** the cohort is structurally too noisy for ANY N=38 classifier. The next leverage is either (a) curate labels manually, (b) expand cohort to N=150 with strict MIC filters, or (c) switch to a different ground-truth source than BV-BRC.",
        "",
        "## Recommended next experiments",
        "",
        "1. Run **curated AMR baseline** (`scripts/cipro_curated_baseline.py`) on full N=38 first to establish a ceiling.",
        "2. Re-run on `--restrict-to-decisive` subset to test whether label noise was the limiter.",
        "3. If both fail and PIVOT TRIGGER condition 4 is NOT met, the next experiment is **per-gene NT windows** on the CLEAN_R_QRDR + CLEAN_S_no_mechanism strains.",
        "",
        f"_JSON sidecar: `{json_path}`_",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"[merge] wrote packet: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
