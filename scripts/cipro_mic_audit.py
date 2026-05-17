"""Cipro AST/MIC rejoin audit on the N=38 cohort.

Experiment 2 from the Open Questions brainstorm. Loaded cohort persists binary
ast_ciprofloxacin labels only; raw BV-BRC AST CSV carries MIC values. Rejoin
per-strain MIC + method + source + duplicate/conflict status. Flag rows
near CLSI/EUCAST breakpoints — those are the inherently-noisy phenotype labels
a genomic predictor can't recover.

Cipro breakpoints (E. coli):
- CLSI 2024: S ≤0.5 / I 1.0 / R ≥2.0 ug/mL
- EUCAST 14.0: S ≤0.25 / I 0.5 / R ≥1.0 ug/mL

Confidence tiers (per Codex synthesis):
- HIGH_R: MIC >= 4x CLSI-R breakpoint (>= 8.0 ug/mL)
- HIGH_S: MIC <= 1/4x CLSI-S breakpoint (<= 0.125 ug/mL)
- BORDERLINE: within 2x of either breakpoint
- AMBIGUOUS: CLSI and EUCAST disagree on R/S call
- NO_MIC: no MIC value in any AST row for this strain
- CONFLICT: multiple AST rows for the strain disagree on R/S

Output: wiki/cipro_mic_audit_<date>.md + .json sidecar.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from datetime import date as _date
from pathlib import Path
from statistics import median

import pandas as pd

from dna_decode.data.cohort import load_cohort

CLSI_R = 2.0
CLSI_S = 0.5
EUCAST_R = 1.0
EUCAST_S = 0.25


def _parse_mic(raw: str, sign: str) -> float | None:
    s = (raw or "").strip()
    if not s or s.upper() in {"NA", "N/A", "NULL", "NONE", "-"}:
        return None
    s = s.lstrip("<>=").strip()
    try:
        return float(s)
    except ValueError:
        return None


def _confidence_tier(mics: list[float], r_calls: list[str]) -> tuple[str, dict]:
    """Classify a strain's confidence tier from rejoined MIC + R/S calls."""
    distinct_calls = {c.upper() for c in r_calls if c}
    conflict = len(distinct_calls & {"R", "RESISTANT"}) > 0 and \
               len(distinct_calls & {"S", "SUSCEPTIBLE"}) > 0
    if conflict:
        return "CONFLICT", {"distinct_calls": sorted(distinct_calls)}
    if not mics:
        return "NO_MIC", {}
    med = median(mics)
    clsi_call = "R" if med >= CLSI_R else ("S" if med <= CLSI_S else "I")
    eucast_call = "R" if med >= EUCAST_R else ("S" if med <= EUCAST_S else "I")
    detail = {
        "median_mic": med,
        "min_mic": min(mics),
        "max_mic": max(mics),
        "n_mic_rows": len(mics),
        "clsi_call": clsi_call,
        "eucast_call": eucast_call,
        "distance_to_clsi_r": med / CLSI_R,
        "distance_to_clsi_s": med / CLSI_S,
    }
    if clsi_call != eucast_call:
        return "AMBIGUOUS", detail
    if med >= 4 * CLSI_R:  # >= 8.0
        return "HIGH_R", detail
    if med <= CLSI_S / 4:  # <= 0.125
        return "HIGH_S", detail
    if CLSI_S / 2 <= med <= 2 * CLSI_R:  # 0.25 - 4.0
        return "BORDERLINE", detail
    return "DECISIVE", detail


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=Path("data/processed/gate_b_n40_cipro_cohort.parquet"))
    parser.add_argument("--ast-csv", type=Path, default=Path("C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv"))
    parser.add_argument("--drug", default="ciprofloxacin")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = Path(f"wiki/cipro_mic_audit_{_date.today().isoformat()}.md")

    cohort = load_cohort(args.cohort)
    drug_lower = args.drug.lower()
    strains_with_label = [s for s in cohort.strains if drug_lower in s.ast_labels]
    print(f"[mic_audit] cohort: {len(cohort.strains)} total; {len(strains_with_label)} with {args.drug} label")

    cohort_ids = {s.strain_id for s in strains_with_label}
    bin_label = {s.strain_id: int(s.ast_labels[drug_lower]) for s in strains_with_label}
    accession = {s.strain_id: s.assembly_accession for s in strains_with_label}
    mlst = {s.strain_id: s.mlst for s in strains_with_label}

    print(f"[mic_audit] loading raw AST: {args.ast_csv}")
    df = pd.read_csv(args.ast_csv, sep=None, engine="python", dtype=str, keep_default_na=False)
    df = df[df["Antibiotic"].str.contains(args.drug, case=False, na=False)]
    df = df[df["Genome ID"].isin(cohort_ids)]
    print(f"[mic_audit] raw AST rows for {args.drug} on cohort: {len(df)}")

    per_strain_rows: dict[str, list[dict]] = defaultdict(list)
    for _, row in df.iterrows():
        gid = row["Genome ID"]
        mic = _parse_mic(row.get("Measurement Value", "") or row.get("Measurement", ""), row.get("Measurement Sign", ""))
        per_strain_rows[gid].append({
            "mic_value": mic,
            "mic_sign": row.get("Measurement Sign", ""),
            "mic_unit": row.get("Measurement Unit", ""),
            "method": row.get("Laboratory Typing Method", ""),
            "testing_standard": row.get("Testing Standard", ""),
            "phenotype": row.get("Resistant Phenotype", ""),
            "source": row.get("Source", ""),
        })

    # Tier each strain
    audit_rows: list[dict] = []
    tier_counts: dict[str, int] = defaultdict(int)
    for sid in sorted(cohort_ids):
        rows = per_strain_rows.get(sid, [])
        mics = [r["mic_value"] for r in rows if r["mic_value"] is not None]
        calls = [r["phenotype"] for r in rows if r["phenotype"]]
        tier, detail = _confidence_tier(mics, calls)
        tier_counts[tier] += 1
        # Per-class summary
        ast_label = bin_label[sid]
        audit_rows.append({
            "strain_id": sid,
            "accession": accession[sid],
            "mlst": mlst[sid],
            "cohort_binary_label": ast_label,
            "cohort_binary_R_or_S": "R" if ast_label == 1 else "S",
            "n_ast_rows": len(rows),
            "n_mic_rows": len(mics),
            "tier": tier,
            "detail": detail,
            "rows": rows,
        })

    # R-set vs S-set tier distribution
    r_tier_counts: dict[str, int] = defaultdict(int)
    s_tier_counts: dict[str, int] = defaultdict(int)
    for r in audit_rows:
        bucket = r_tier_counts if r["cohort_binary_label"] == 1 else s_tier_counts
        bucket[r["tier"]] += 1

    # Decisively-R subset (HIGH_R only) — the "clean" cipro-R cohort
    decisive_R = [r["strain_id"] for r in audit_rows if r["cohort_binary_label"] == 1 and r["tier"] == "HIGH_R"]
    decisive_S = [r["strain_id"] for r in audit_rows if r["cohort_binary_label"] == 0 and r["tier"] == "HIGH_S"]
    decisive_total = len(decisive_R) + len(decisive_S)

    # Verdict
    n_borderline_R = r_tier_counts.get("BORDERLINE", 0) + r_tier_counts.get("AMBIGUOUS", 0)
    n_no_mic_R = r_tier_counts.get("NO_MIC", 0)
    n_total_R = sum(r_tier_counts.values())
    pct_clean_R = (r_tier_counts.get("HIGH_R", 0) / n_total_R * 100) if n_total_R else 0.0
    verdict = (
        "CLEAN" if pct_clean_R >= 70 and n_no_mic_R <= 2 else
        "MIXED" if pct_clean_R >= 40 else
        "NOISY"
    )

    print(f"\n=== Tier distribution (R set, n={n_total_R}) ===")
    for t, n in sorted(r_tier_counts.items(), key=lambda x: -x[1]):
        print(f"  {t:12s} {n}")
    print(f"\n=== Tier distribution (S set, n={sum(s_tier_counts.values())}) ===")
    for t, n in sorted(s_tier_counts.items(), key=lambda x: -x[1]):
        print(f"  {t:12s} {n}")
    print(f"\nDecisive-R: {len(decisive_R)} / Decisive-S: {len(decisive_S)} / Total decisive: {decisive_total}")
    print(f"Verdict: {verdict}")

    # JSON sidecar
    json_path = args.output.with_suffix(".json")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "cohort_path": str(args.cohort),
        "ast_csv": str(args.ast_csv),
        "drug": args.drug,
        "breakpoints": {
            "clsi_r": CLSI_R, "clsi_s": CLSI_S,
            "eucast_r": EUCAST_R, "eucast_s": EUCAST_S,
        },
        "n_cohort": len(strains_with_label),
        "r_tier_counts": dict(r_tier_counts),
        "s_tier_counts": dict(s_tier_counts),
        "decisive_R_ids": decisive_R,
        "decisive_S_ids": decisive_S,
        "verdict": verdict,
        "per_strain": audit_rows,
    }
    json_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(f"[mic_audit] wrote JSON: {json_path}")

    # Markdown packet
    lines = [
        f"# Cipro AST/MIC audit — N=38 cohort ({_date.today().isoformat()})",
        "",
        f"**Purpose:** rejoin raw BV-BRC AST MIC values onto the cohort to surface borderline/no-MIC strains whose binary R/S label is structurally noisy.",
        f"**Source:** `{args.ast_csv}`",
        f"**Breakpoints:** CLSI 2024 (S ≤{CLSI_S} / R ≥{CLSI_R}), EUCAST 14.0 (S ≤{EUCAST_S} / R ≥{EUCAST_R}) ug/mL.",
        "",
        f"## Verdict: **{verdict}**",
        f"- Decisive-R (HIGH_R only): **{len(decisive_R)}** / {n_total_R}",
        f"- Decisive-S (HIGH_S only): **{len(decisive_S)}** / {sum(s_tier_counts.values())}",
        f"- Total decisive subset: **{decisive_total}** strains",
        f"- Clean-R fraction: **{pct_clean_R:.1f}%**",
        "",
        "## Tier distribution",
        "",
        "| tier | R strains | S strains | description |",
        "|---|---:|---:|---|",
        f"| HIGH_R | {r_tier_counts.get('HIGH_R', 0)} | {s_tier_counts.get('HIGH_R', 0)} | MIC ≥8.0 (4× CLSI-R) — decisive R |",
        f"| HIGH_S | {r_tier_counts.get('HIGH_S', 0)} | {s_tier_counts.get('HIGH_S', 0)} | MIC ≤0.125 (1/4× CLSI-S) — decisive S |",
        f"| DECISIVE | {r_tier_counts.get('DECISIVE', 0)} | {s_tier_counts.get('DECISIVE', 0)} | CLSI/EUCAST agree, not borderline |",
        f"| BORDERLINE | {r_tier_counts.get('BORDERLINE', 0)} | {s_tier_counts.get('BORDERLINE', 0)} | MIC in [0.25, 4.0] — within 2× breakpoint |",
        f"| AMBIGUOUS | {r_tier_counts.get('AMBIGUOUS', 0)} | {s_tier_counts.get('AMBIGUOUS', 0)} | CLSI and EUCAST disagree on call |",
        f"| CONFLICT | {r_tier_counts.get('CONFLICT', 0)} | {s_tier_counts.get('CONFLICT', 0)} | multiple AST rows disagree R vs S |",
        f"| NO_MIC | {r_tier_counts.get('NO_MIC', 0)} | {s_tier_counts.get('NO_MIC', 0)} | no MIC value parsed |",
        "",
        "## Per-strain audit",
        "",
        "| strain_id | accession | label | tier | median MIC | n_rows | n_mic | mlst |",
        "|---|---|---|---|---:|---:|---:|---|",
    ]
    for r in sorted(audit_rows, key=lambda x: (x["cohort_binary_label"], x["tier"], x["strain_id"])):
        med = r["detail"].get("median_mic", "")
        med_str = f"{med:.3f}" if isinstance(med, (int, float)) else "-"
        lines.append(
            f"| {r['strain_id']} | {r['accession']} | {r['cohort_binary_R_or_S']} | {r['tier']} | "
            f"{med_str} | {r['n_ast_rows']} | {r['n_mic_rows']} | {r['mlst']} |"
        )
    lines.extend([
        "",
        "## Decisive subset (recommended for clean re-run)",
        "",
        f"- **Decisive-R (HIGH_R):** {sorted(decisive_R)}",
        f"- **Decisive-S (HIGH_S):** {sorted(decisive_S)}",
        f"- **N = {decisive_total}** (vs cohort N=38)",
        "",
        "## How to use",
        "",
        "1. If verdict is CLEAN: cohort labels are trustworthy; downstream NT/k-mer FAIL is a model issue, not label noise.",
        "2. If verdict is MIXED or NOISY: re-run Stage 1 on the decisive subset; if AUROC ≥0.10 better, label noise was a real confounder.",
        "3. Use the AMBIGUOUS / CONFLICT / NO_MIC strains as a flag list to exclude or re-label.",
        "4. Feed `decisive_R_ids` + `decisive_S_ids` to a curated-baseline experiment for a cleaner ground truth.",
        "",
        f"_JSON sidecar: `{json_path}`_",
    ])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"[mic_audit] wrote packet: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
