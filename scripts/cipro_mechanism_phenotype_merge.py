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


PRIMARY_CIPRO_MECHANISMS = {"QRDR_target_alteration", "plasmid_protect_modify"}
CO_RESISTANCE_MECHANISMS = {"efflux", "regulatory", "porin_loss"}


def _classify_noise(row: dict) -> tuple[str, bool, list[str]]:
    """Per-row noise classification + opacity flag + co-resistance modifiers.

    Returns (noise_class, mechanism_opacity_flag, co_resistance_modifiers).

    - **Primary cipro mechanism** = QRDR target alteration OR plasmid quinolone
      protection/modification (qnr / aac6-Ib-cr). These are the textbook
      cipro-conferring mechanisms.
    - **Co-resistance modifiers** = efflux / regulatory / porin_loss. Real
      biology, but does NOT confer cipro-R on their own at clinically-meaningful
      MIC levels. Reported as a separate column; does NOT drive the noise class.
    - **Mechanism opacity flag** = True when AMRFinder mechanism status is
      MISSING or NO_MECHANISM despite a HIGH-tier MIC (i.e., the biology may be
      real but outside AMRFinder's catalog). Distinguishes "labels noisy" from
      "tool incomplete" — Codex round-2 critique 2026-05-17.
    """
    label = row.get("cohort_binary_label")
    tier = row.get("mic_tier", "?")
    mechs = set(row.get("mechanisms_present", []))
    primary_present = bool(mechs & PRIMARY_CIPRO_MECHANISMS)
    co_resistance = sorted(mechs & CO_RESISTANCE_MECHANISMS)
    mechanism_status = row.get("mechanism_status", "MISSING")
    opacity_flag = False

    if label == 1:
        if tier == "HIGH_R" and primary_present:
            noise = "CLEAN_R_primary_mechanism"
        elif tier == "HIGH_R" and not primary_present and mechs:
            # Co-resistance only — AMRFinder finds something but not a primary
            # cipro mechanism. Possibly a non-catalog primary OR a tool miss.
            noise = "OPAQUE_R_co_resistance_only"
            opacity_flag = True
        elif tier == "HIGH_R" and not mechs:
            # High MIC + nothing in AMRFinder catalog = tool/parser opacity
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
            # Primary mechanism present + S label that isn't decisively-S
            # is the strongest mislabeling signal.
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
        noise, opacity, co_res = _classify_noise(row)
        row["noise_class"] = noise
        row["mechanism_opacity_flag"] = opacity
        row["co_resistance_modifiers"] = co_res
        row["has_primary_cipro_mechanism"] = bool(set(row.get("mechanisms_present", [])) & PRIMARY_CIPRO_MECHANISMS)
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
    opaque_r = sum(1 for r in r_strains if r["noise_class"].startswith("OPAQUE_R"))
    noisy_r = sum(1 for r in r_strains if r["noise_class"].startswith("NOISY_R"))
    noisy_s = sum(1 for r in s_strains if r["noise_class"].startswith("NOISY_S"))
    opacity_count = sum(1 for r in merged if r.get("mechanism_opacity_flag"))
    clean_count = clean_r + clean_s

    cohort_signal_quality = clean_count / max(1, len(merged))

    verdict = (
        "SIGNAL_DOMINATES" if cohort_signal_quality >= 0.7 else
        "MIXED" if cohort_signal_quality >= 0.4 else
        "NOISE_DOMINATES"
    )

    # Pre-curated-baseline gate (Patch 4): clean_count drives the next-experiment
    # recommendation. Codex round-2 critique noted that clean_count=0 may be a
    # MECHANISM_OPACITY problem (AMRFinder gap) rather than a label problem —
    # so opacity_count is reported separately and acts as a tie-breaker.
    if clean_count >= 20:
        gate = "RUN_CURATED_BASELINE_FULL_AND_CLEAN"
        next_step = "Run cipro_curated_baseline.py at full N=38 + on clean-subset filter. Verdict is load-bearing."
    elif clean_count >= 10:
        gate = "RUN_CURATED_BASELINE_FULL_ONLY"
        next_step = "Run cipro_curated_baseline.py at full N=38 only. Verdict is descriptive — clearly note label noise upper bound."
    elif opacity_count >= 5:
        # High clean_count failure but high opacity = AMRFinder may be wrong,
        # not labels. Debug AMRFinder before declaring labels unusable.
        gate = "MECHANISM_DEBUG_BRANCH"
        next_step = (
            f"clean_count={clean_count} is low BUT opacity_count={opacity_count} is high — "
            "AMRFinder may be missing mechanisms, not labels being unusable. "
            "Run scripts/cipro_curated_baseline.py with --skip-gate flag for an INFORMATIONAL run; "
            "separately run mechanism debug (manual gyrA/parC inspection on HIGH_R opacity strains)."
        )
    else:
        gate = "SUSPEND_CONDITION_4"
        next_step = (
            f"clean_count={clean_count} too low AND opacity_count={opacity_count} not the bottleneck — "
            "the N=38 cohort is structurally unusable for PIVOT TRIGGER condition 4. "
            "Next experiments: (a) cohort expansion to N=150 with strict MIC filter, OR "
            "(b) per-gene NT windows on the small clean set as a diagnostic."
        )

    print(f"[merge] N={len(merged)} merged rows")
    print(f"\n=== Noise class distribution ===")
    for cls, n in sorted(noise_counts.items(), key=lambda x: -x[1]):
        print(f"  {cls:38s} {n}")
    print(f"\nClean (R+S): {clean_r}R + {clean_s}S = {clean_count}")
    print(f"Suspect (R+S): {suspect_r}R + {suspect_s}S = {suspect_r + suspect_s}")
    print(f"Opaque-R (mechanism missing): {opaque_r}")
    print(f"Mechanism opacity flag count: {opacity_count}")
    print(f"Noisy (R+S):  {noisy_r}R + {noisy_s}S = {noisy_r + noisy_s}")
    print(f"Cohort signal quality: {cohort_signal_quality:.2f}")
    print(f"Verdict: {verdict}")
    print(f"\nPre-curated-baseline gate: {gate}")
    print(f"Recommended next step: {next_step}")

    # JSON sidecar
    json_path = args.output.with_suffix(".json")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "mech_audit": str(args.mech_audit),
        "mic_audit": str(args.mic_audit),
        "n_merged": len(merged),
        "noise_class_counts": dict(noise_counts),
        "clean_count": clean_count,
        "clean_r": clean_r,
        "clean_s": clean_s,
        "suspect_count": suspect_r + suspect_s,
        "opaque_r": opaque_r,
        "opacity_count": opacity_count,
        "noisy_count": noisy_r + noisy_s,
        "cohort_signal_quality": cohort_signal_quality,
        "verdict": verdict,
        "pre_curated_gate": gate,
        "recommended_next_step": next_step,
        "primary_cipro_mechanisms_definition": sorted(PRIMARY_CIPRO_MECHANISMS),
        "co_resistance_mechanisms_definition": sorted(CO_RESISTANCE_MECHANISMS),
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
        f"- Clean strains: {clean_r}R + {clean_s}S = **{clean_count}**",
        f"- Suspect strains: {suspect_r}R + {suspect_s}S = **{suspect_r + suspect_s}** (likely mislabeled or atypical biology)",
        f"- Opaque-R strains (HIGH_R + no primary mechanism): **{opaque_r}** — tool/parser miss vs novel biology",
        f"- Total opacity flagged: **{opacity_count}** (informs whether to suspect labels OR AMRFinder)",
        f"- Noisy strains: {noisy_r}R + {noisy_s}S = **{noisy_r + noisy_s}** (label is structurally unreliable)",
        "",
        f"## Pre-curated-baseline gate: **{gate}**",
        f"- {next_step}",
        "",
        "Definitions:",
        f"- **Primary cipro mechanism** = QRDR target alteration (gyrA/parC/parE) OR plasmid quinolone protection/modification (qnr* / aac6-Ib-cr). Textbook cipro-conferring.",
        f"- **Co-resistance modifiers** = efflux + regulatory + porin_loss. Real biology but do not confer cipro-R alone; reported separately, do not drive noise classification.",
        f"- **Mechanism opacity flag** = HIGH_R + no primary mechanism. AMRFinder may have a catalog gap; distinct from label noise.",
        "",
        "## Noise class distribution",
        "",
        "| class | count | meaning |",
        "|---|---:|---|",
        f"| CLEAN_R_primary_mechanism | {sum(1 for r in r_strains if r['noise_class']=='CLEAN_R_primary_mechanism')} | HIGH_R MIC + QRDR or qnr/aac present |",
        f"| OPAQUE_R_co_resistance_only | {sum(1 for r in r_strains if r['noise_class']=='OPAQUE_R_co_resistance_only')} | HIGH_R + only efflux/regulatory/porin found — tool gap likely |",
        f"| OPAQUE_R_no_mechanism | {sum(1 for r in r_strains if r['noise_class']=='OPAQUE_R_no_mechanism')} | HIGH_R but no AMRFinder hits at all |",
        f"| NOISY_R_borderline | {sum(1 for r in r_strains if r['noise_class']=='NOISY_R_borderline')} | borderline/ambiguous MIC — label may be wrong |",
        f"| NOISY_R_no_mic | {sum(1 for r in r_strains if r['noise_class']=='NOISY_R_no_mic')} | no MIC in BV-BRC — label is opaque |",
        f"| CLEAN_S_no_primary_mechanism | {sum(1 for r in s_strains if r['noise_class']=='CLEAN_S_no_primary_mechanism')} | HIGH_S MIC + no QRDR / no qnr — clean susceptible |",
        f"| SUSPECT_S_silent_primary_mechanism | {sum(1 for r in s_strains if r['noise_class']=='SUSPECT_S_silent_primary_mechanism')} | HIGH_S MIC but primary mechanism present — silent / non-functional |",
        f"| SUSPECT_S_borderline_primary_mechanism | {sum(1 for r in s_strains if r['noise_class']=='SUSPECT_S_borderline_primary_mechanism')} | borderline MIC + primary mech — likely mislabeled to S |",
        f"| NOISY_S_borderline | {sum(1 for r in s_strains if r['noise_class']=='NOISY_S_borderline')} | borderline MIC with no primary mech — label opaque |",
        f"| NOISY_S_no_mic | {sum(1 for r in s_strains if r['noise_class']=='NOISY_S_no_mic')} | no MIC in BV-BRC — label opaque |",
        "",
        "## Per-strain merged table",
        "",
        "| strain_id | accession | label | mic_tier | med MIC | primary? | mechs | co-res | noise_class | opacity | mlst |",
        "|---|---|---|---|---:|---|---|---|---|---|---|",
    ]
    for r in sorted(merged, key=lambda x: (x["cohort_binary_label"], x["noise_class"], x["strain_id"])):
        med_str = f"{r['median_mic']:.3f}" if r['median_mic'] is not None else "-"
        mechs_str = ",".join(r["mechanisms_present"]) if r["mechanisms_present"] else "-"
        co_res_str = ",".join(r.get("co_resistance_modifiers", []) or []) or "-"
        primary_str = "yes" if r.get("has_primary_cipro_mechanism") else "no"
        opacity_str = "!" if r.get("mechanism_opacity_flag") else " "
        lines.append(
            f"| {r['strain_id']} | {r['accession']} | {r['cohort_binary_RS']} | {r['mic_tier']} | "
            f"{med_str} | {primary_str} | {mechs_str} | {co_res_str} | {r['noise_class']} | {opacity_str} | {r['mlst']} |"
        )
    lines.extend([
        "",
        "## How to use (gated)",
        "",
        f"Pre-curated-baseline gate fired: **{gate}**.",
        f"- {next_step}",
        "",
        "Verdict-level interpretation:",
        "- **SIGNAL_DOMINATES (>=0.70 clean):** the cohort is clean enough that NT's FAIL on cipro is a genuine model issue, not label noise. Curated baseline verdict is load-bearing.",
        "- **MIXED (0.40-0.70 clean):** label noise is a real confounder. Curated baseline at full N is descriptive; consider a clean-subset rerun if `clean_count >= 20`.",
        "- **NOISE_DOMINATES (<0.40 clean):** the N=38 cohort is structurally noisy. Two branches: if `opacity_count >= 5` -> AMRFinder may be missing mechanisms, debug before declaring labels unusable; otherwise -> N=150 expansion or per-gene NT diagnostic.",
        "",
        f"_JSON sidecar: `{json_path}`_",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"[merge] wrote packet: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
