"""Decoder-suite provenance-disjoint VALIDATION REPORT CARD (v0 roll-up of the standing capability).

Anchor-4 (idea-anchor 2026-06-10): turn the proven one-off provenance-disjoint NCBI-PD validation into a
STANDING, suite-wide trust surface. This v0 ROLLS UP what already exists on disk into one maintained report
card; it does NOT score or census (those are the Stage-1 census + Stage-2 validator scripts). It is pure,
read-only over its inputs, and writes exactly ONE pair of output files nothing else produces.

Inputs (all already on disk; no network, no Docker):
  - wiki/provenance_disjoint_validation_*.json   -> SCORED cells (provenance-disjoint-validation-v1 schema)
  - wiki/provdisjoint_census_results.json        -> powering verdicts (provdisjoint-census-results-v1)
  - dna_decode/data/calibrated_amr_rules.json    -> ABSTAINS_BY_DESIGN cells (EXPRESSION_FLOOR verdict)

Cell-state machine (the probe's honest-tiering requirement — "validated" is per-cell, never a suite headline):
  SCORED                    a Stage-2 provdisjoint JSON exists -> show acc/sens/spec/n + its honest tier
  POWERED_UNSCORED          censused >= MIN/class both classes, not yet scored
  UNDERPOWERED              censused < MIN/class (e.g. surveillance-dominated organisms)
  ABSTAINS_BY_DESIGN        registry verdict EXPRESSION_FLOOR (rule refuses expression-driven R it can't decode)
  NOT_CENSUSED              bacterial + census-able, but no census yet
  NO_FREE_PHENOTYPE_SOURCE  fungal/antiviral/antimalarial -> no free isolate-level AST source (structural non-cell)

HONEST TIER (do NOT inflate): every SCORED cell is PROVENANCE-disjoint (different submitter/lab/country),
NOT methodology-independent (most NCBI submitters use CLSI broth microdilution) and NOT external clinical
validation. There is deliberately NO aggregate "X% validated" headline number.

Usage: .venv/Scripts/python.exe scripts/build_validation_report_card.py
Exit 0 always (it's a report, not a gate).
"""
from __future__ import annotations

import glob
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WIKI = ROOT / "wiki"

# Other-kingdom decoders: shipped, but NOT census-able on NCBI-PD AST_phenotypes (structural non-cells).
# (organism, drug, engine). Ledger row 58: no free isolate-level phenotype source exists for these.
OTHER_KINGDOM = [
    ("Candida_auris", "fluconazole", "fungal ERG11 BLAST target-site"),
    ("Candida_auris", "voriconazole", "fungal ERG11 BLAST target-site"),
    ("Candida_auris", "caspofungin", "fungal FKS1 BLAST target-site"),
    ("Candida_auris", "micafungin", "fungal FKS1 BLAST target-site"),
    ("Influenza_A", "oseltamivir", "influenza NA-inhibitor target-site"),
    ("Plasmodium_falciparum", "artemisinin", "P. falciparum kelch13 target-site"),
]

PROV_TIER = ("provenance-disjoint (different submitter/lab/country); NOT methodology-independent "
             "(most submitters use CLSI broth microdilution); NOT external clinical validation")


def _key(org: str, drug: str) -> tuple[str, str]:
    return (org.lower(), drug.lower())


def load_scored() -> dict:
    cells = {}
    for f in sorted(glob.glob(str(WIKI / "provenance_disjoint_validation_*.json"))):
        try:
            d = json.loads(Path(f).read_text(encoding="utf-8"))
        except Exception:
            continue
        org, drug = d.get("organism"), d.get("drug")
        if org and drug:
            cells[_key(org, drug)] = {**d, "_file": Path(f).name}
    return cells


def load_census() -> dict:
    p = WIKI / "provdisjoint_census_results.json"
    out = {}
    if p.exists():
        for r in json.loads(p.read_text(encoding="utf-8")).get("results", []):
            out[_key(r["organism"], r["drug"])] = r
    return out


def load_registry() -> dict:
    p = ROOT / "dna_decode" / "data" / "calibrated_amr_rules.json"
    out = {}
    if p.exists():
        for k, v in json.loads(p.read_text(encoding="utf-8")).get("rules", {}).items():
            org, _, drug = k.partition("|")
            if drug:
                out[_key(org, drug)] = v
    return out


def classify(key, scored, census, registry) -> dict:
    org, drug = key
    if key in scored:
        m = scored[key].get("metrics", {})
        return {"state": "SCORED",
                "acc": m.get("acc"), "sens": m.get("sens"), "spec": m.get("spec"),
                "n": m.get("n_scored"),
                "tp": m.get("tp"), "fp": m.get("fp"), "tn": m.get("tn"), "fn": m.get("fn"),
                "abstain": m.get("abstain"),
                "tier": scored[key].get("independence_tier", PROV_TIER),
                "file": scored[key].get("_file")}
    if key in registry and str(registry[key].get("verdict", "")).upper() == "EXPRESSION_FLOOR":
        return {"state": "ABSTAINS_BY_DESIGN",
                "note": f"registry verdict EXPRESSION_FLOOR ({registry[key].get('counter')}@{registry[key].get('threshold')}) "
                        "— rule refuses expression-driven R it cannot decode"}
    if key in census:
        c = census[key]
        if c.get("powered"):
            return {"state": "POWERED_UNSCORED",
                    "note": f"censused {c.get('other_R')}R/{c.get('other_S')}S provenance-disjoint (>=MIN/class); not yet scored"}
        return {"state": "UNDERPOWERED",
                "note": f"censused {c.get('other_R')}R/{c.get('other_S')}S provenance-disjoint (< MIN/class) — surveillance-dominated"}
    return {"state": "NOT_CENSUSED",
            "note": "bacterial + census-able; no provenance census yet"}


def main() -> int:
    scored, census, registry = load_scored(), load_census(), load_registry()

    rows = []
    seen = set()
    # 1) bacterial cells: union of scored + censused + registry keys
    for key in sorted(set(scored) | set(census) | set(registry)):
        rows.append((key, classify(key, scored, census, registry)))
        seen.add(key)
    # 2) other-kingdom structural non-cells
    for org, drug, engine in OTHER_KINGDOM:
        key = _key(org, drug)
        if key in seen:
            continue
        rows.append((key, {"state": "NO_FREE_PHENOTYPE_SOURCE",
                           "note": f"{engine}; no free isolate-level AST source on NCBI-PD (structural non-cell)"}))

    counts = {}
    for _, c in rows:
        counts[c["state"]] = counts.get(c["state"], 0) + 1

    today = _date.today().isoformat()
    artifact = {"_schema": "decoder-validation-report-card-v0", "date": today,
                "honest_tier": PROV_TIER, "no_aggregate_headline": True,
                "state_counts": counts,
                "cells": [{"organism": k[0], "drug": k[1], **c} for k, c in rows]}
    (WIKI / "decoder_validation_report_card.json").write_text(json.dumps(artifact, indent=2), encoding="utf-8")

    # markdown
    L = []
    L.append(f"# Decoder-suite provenance-disjoint validation report card — {today}\n")
    L.append("Standing trust surface for the shipped deterministic AMR decoders (Anchor-4). Each cell is the "
             "DEPLOYED `call_resistance(organism, drug)` rule scored on a FRESH, leakage-checked, "
             "**provenance-disjoint** NCBI-PD cohort (submitters OUTSIDE NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA).\n")
    L.append("> **Honest tier (do NOT inflate):** every SCORED cell is provenance-disjoint (different "
             "submitter/lab/country), **NOT** methodology-independent (most submitters use CLSI broth "
             "microdilution) and **NOT** external clinical validation. There is deliberately **no aggregate "
             "“X% validated” number** — read the grid cell by cell.\n")
    L.append("## State legend\n")
    L.append("| state | meaning |\n|---|---|")
    L.append("| `SCORED` | Stage-2 provdisjoint run exists — acc/sens/spec shown |")
    L.append("| `POWERED_UNSCORED` | censused ≥ 20/class both classes; not yet scored |")
    L.append("| `UNDERPOWERED` | censused < 20/class (surveillance-dominated organism) |")
    L.append("| `ABSTAINS_BY_DESIGN` | registry EXPRESSION_FLOOR — rule refuses what it can't decode |")
    L.append("| `NOT_CENSUSED` | bacterial + census-able; no census yet |")
    L.append("| `NO_FREE_PHENOTYPE_SOURCE` | fungal/antiviral/antimalarial — no free isolate-level AST (structural non-cell) |\n")
    L.append("## State counts\n")
    L.append("| state | cells |\n|---|---|")
    for s in ("SCORED", "POWERED_UNSCORED", "UNDERPOWERED", "ABSTAINS_BY_DESIGN", "NOT_CENSUSED", "NO_FREE_PHENOTYPE_SOURCE"):
        if s in counts:
            L.append(f"| `{s}` | {counts[s]} |")
    L.append("\n## Cells\n")
    L.append("| organism | drug | state | acc | sens | spec | n | detail |\n|---|---|---|---|---|---|---|---|")
    for k, c in rows:
        org, drug = k
        if c["state"] == "SCORED":
            detail = c.get("tier", "")[:48]
            L.append(f"| {org} | {drug} | `SCORED` | {c.get('acc')} | {c.get('sens')} | {c.get('spec')} | "
                     f"{c.get('n')} | TP{c.get('tp')} FP{c.get('fp')} TN{c.get('tn')} FN{c.get('fn')} |")
        else:
            L.append(f"| {org} | {drug} | `{c['state']}` | — | — | — | — | {c.get('note','')} |")
    L.append("\n## Provenance\n")
    L.append("- SCORED cells: `wiki/provenance_disjoint_validation_*.json` (Stage-2 `provenance_disjoint_validate.py`).")
    L.append("- Powering: `wiki/provdisjoint_census_results.json` (Stage-1 `ncbi_pd_provenance_census.py`).")
    L.append("- ABSTAINS: `dna_decode/data/calibrated_amr_rules.json` (EXPRESSION_FLOOR verdicts).")
    L.append("- Rebuild: `.venv/Scripts/python.exe scripts/build_validation_report_card.py` (read-only roll-up; re-run as cells land).")
    (WIKI / "decoder_validation_report_card.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    print(f"report card written: wiki/decoder_validation_report_card.md")
    print(f"state counts: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
