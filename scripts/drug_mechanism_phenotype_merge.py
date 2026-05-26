"""Drug-agnostic mechanism x MIC merge for audit-aware prediction packets."""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from datetime import date as _date
from pathlib import Path

from dna_decode.data.mic_tiers import CO_RESISTANCE_MECHANISMS, primary_mechanisms_for, supported_drugs


def classify_noise(row: dict, drug: str) -> tuple[str, bool, list[str]]:
    label = row.get("cohort_binary_label")
    tier = row.get("mic_tier", "?")
    mechs = set(row.get("mechanisms_present", []))
    primary_set = set(primary_mechanisms_for(drug))
    primary_present = bool(mechs & primary_set)
    co_resistance = sorted(mechs & set(CO_RESISTANCE_MECHANISMS))
    opacity_flag = False

    if label == 1:
        if tier == "HIGH_R" and primary_present:
            noise = "CLEAN_R_primary_mechanism"
        elif tier == "HIGH_R" and not primary_present and mechs:
            noise = "OPAQUE_R_co_resistance_only"
            opacity_flag = True
        elif tier == "HIGH_R" and not mechs:
            noise = "OPAQUE_R_no_mechanism"
            opacity_flag = True
        elif tier in {"BORDERLINE", "AMBIGUOUS", "CONFLICT"}:
            noise = "NOISY_R_borderline"
        elif tier == "NO_MIC":
            noise = "NOISY_R_no_mic"
        else:
            noise = "OTHER_R"
    else:
        if tier == "HIGH_S" and not primary_present:
            noise = "CLEAN_S_no_primary_mechanism"
        elif tier == "HIGH_S" and primary_present:
            noise = "SUSPECT_S_silent_primary_mechanism"
        elif primary_present and tier != "HIGH_S":
            noise = "SUSPECT_S_borderline_primary_mechanism"
        elif tier in {"BORDERLINE", "AMBIGUOUS", "CONFLICT"}:
            noise = "NOISY_S_borderline"
        elif tier == "NO_MIC":
            noise = "NOISY_S_no_mic"
        else:
            noise = "OTHER_S"

    return noise, opacity_flag, co_resistance


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--drug", required=True, choices=supported_drugs())
    parser.add_argument("--mech-audit", type=Path, required=True)
    parser.add_argument("--mic-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--suspend-threshold",
        type=float,
        default=0.40,
        help=(
            "Below this signal_quality, fire SUSPEND_CONDITION_4. Default 0.40 matches "
            "scripts/cipro_mechanism_phenotype_merge.py NOISE_DOMINATES threshold (2026-05-17 "
            "calibration; see plans/Cef_Audit_Aware_Packet_Design.md edit #2 lock)."
        ),
    )
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = Path(f"wiki/{args.drug.lower()}_mechanism_phenotype_merge_{_date.today().isoformat()}.md")

    if not args.mech_audit.exists():
        print(f"[merge] mech audit missing: {args.mech_audit}", file=sys.stderr)
        return 2
    if not args.mic_audit.exists():
        print(f"[merge] mic audit missing: {args.mic_audit}", file=sys.stderr)
        return 2

    mech = json.loads(args.mech_audit.read_text(encoding="utf-8"))
    mic = json.loads(args.mic_audit.read_text(encoding="utf-8"))
    mic_by_strain = {r["strain_id"]: r for r in mic["per_strain"]}
    mech_by_strain = {r["strain_id"]: r for r in mech["per_strain"]}

    merged: list[dict] = []
    for sid in sorted(set(mic_by_strain) | set(mech_by_strain)):
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
        noise, opacity, co_res = classify_noise(row, args.drug)
        row["noise_class"] = noise
        row["mechanism_opacity_flag"] = opacity
        row["co_resistance_modifiers"] = co_res
        row["primary_mechanisms"] = sorted(
            set(row.get("mechanisms_present", [])) & set(primary_mechanisms_for(args.drug))
        )
        merged.append(row)

    noise_counts: dict[str, int] = defaultdict(int)
    for row in merged:
        noise_counts[row["noise_class"]] += 1
    clean_count = sum(1 for row in merged if row["noise_class"].startswith("CLEAN_"))
    suspect_count = sum(1 for row in merged if row["noise_class"].startswith("SUSPECT_"))
    opacity_count = sum(1 for row in merged if row["mechanism_opacity_flag"])
    # signal_quality formula LOCKED: clean / total. Matches scripts/cipro_mechanism_phenotype_merge.py
    # line 163 (`cohort_signal_quality = clean_count / max(1, len(merged))`). Earlier draft used
    # clean / (clean + suspect + opacity); this would bias the gate above the K/N null baseline
    # because the denominator excludes NOISY/OTHER buckets. See
    # plans/Cef_Audit_Aware_Packet_Design.md edit #1 lock + the K/N null sanity check feedback.
    signal_quality = clean_count / max(1, len(merged))
    gate_verdict = (
        "SUSPEND_CONDITION_4"
        if signal_quality < args.suspend_threshold
        else "RUN_FULL_AND_CLEAN"
        if signal_quality >= 0.70
        else "MIXED"
    )

    json_path = args.output.with_suffix(".json")
    payload = {
        "drug": args.drug,
        "mech_audit": str(args.mech_audit),
        "mic_audit": str(args.mic_audit),
        "gate_verdict": gate_verdict,
        "signal_quality": signal_quality,
        "noise_class_counts": dict(noise_counts),
        "clean_count": clean_count,
        "suspect_count": suspect_count,
        "opacity_count": opacity_count,
        "per_strain": merged,
    }
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")

    lines = [
        f"# {args.drug} mechanism x phenotype merge ({_date.today().isoformat()})",
        "",
        f"- Gate verdict: `{gate_verdict}`",
        f"- Signal quality: `{signal_quality:.3f}`",
        f"- Clean count: `{clean_count}`",
        f"- Suspect count: `{suspect_count}`",
        f"- Opacity count: `{opacity_count}`",
        "",
        "| class | count |",
        "|---|---:|",
    ]
    for cls, count in sorted(noise_counts.items(), key=lambda x: (-x[1], x[0])):
        lines.append(f"| {cls} | {count} |")
    lines.extend(
        [
            "",
            "| strain_id | accession | label | mic_tier | primary_mechanisms | co_resistance | noise_class | opacity |",
            "|---|---|---|---|---|---|---|---|",
        ]
    )
    for row in merged:
        lines.append(
            f"| {row['strain_id']} | {row['accession']} | {row['cohort_binary_RS']} | "
            f"{row['mic_tier']} | {','.join(row['primary_mechanisms']) or '-'} | "
            f"{','.join(row['co_resistance_modifiers']) or '-'} | {row['noise_class']} | "
            f"{'!' if row['mechanism_opacity_flag'] else ''} |"
        )
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[merge] wrote {args.output}")
    print(f"[merge] wrote {json_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
