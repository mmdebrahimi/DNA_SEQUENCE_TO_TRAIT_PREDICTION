"""Decoder-suite provenance-disjoint VALIDATION REPORT CARD (standing roll-up).

Anchor-4: a standing, suite-wide trust surface. Pure read-only roll-up of what already exists on disk into
one maintained report card; it does NOT score or census (those are the Stage-1 census + Stage-2 validator).

Inputs (all on disk; no network, no Docker):
  - dna_decode/data/shipped_decoder_surface.py    -> the authoritative DEPLOYED-CLAIM row set (the grid)
  - wiki/provenance_disjoint_validation_*.json     -> SCORED cells (provenance-disjoint-validation-v1)
  - wiki/provdisjoint_census_results.json          -> powering verdicts (provdisjoint-census-results-v1)
  - dna_decode/data/calibrated_amr_rules.json      -> ABSTAINS_BY_DESIGN cells (EXPRESSION_FLOOR verdict)

Rows = the shipped-decoder surface UNION the observed scored/census/registry keys, so un-censused shipped
decoders still render (NOT_CENSUSED) and a new decoder cannot ship invisibly (the surface is coverage-tested
against the CLI drug catalogs).

Cell-state machine (the probe's + brainstorm's honest-tiering requirement — per-cell, never a suite headline):
  SCORED                    a Stage-2 provdisjoint JSON exists -> acc/sens/spec/n + its honest tier
  POWERED_UNSCORED          censused >= MIN/class both classes, not yet scored
  UNDERPOWERED              censused < MIN/class (surveillance-dominated organisms)
  ABSTAINS_BY_DESIGN        registry verdict EXPRESSION_FLOOR (rule refuses what it can't decode)
  NOT_CENSUSED              bacterial + census-able, no census yet
  LABEL_CONFOUNDED          phenotype LABEL is an unreliable surrogate (oxacillin AST vs mecA / cefoxitin)
  NO_FREE_PHENOTYPE_SOURCE  fungal/antiviral/antimalarial -> no free isolate-level AST source (non-cell)

Surface `phenotype_source_status` (no_free_source / label_confounded) is a STRUCTURAL label property and
takes precedence over observations (we never present a misleading clean SCORED on a confounded label, nor a
NCBI-PD cell where no free phenotype exists).

HONEST TIER (do NOT inflate): every SCORED cell is PROVENANCE-disjoint (different submitter/lab/country),
NOT methodology-independent (most submitters use CLSI broth microdilution) and NOT external clinical
validation. There is deliberately NO aggregate "X% validated" headline.

Usage: .venv/Scripts/python.exe scripts/build_validation_report_card.py   (exit 0 always — a report, not a gate)
"""
from __future__ import annotations

import glob
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WIKI = ROOT / "wiki"
sys.path.insert(0, str(ROOT))

from dna_decode.data.shipped_decoder_surface import surface_index  # noqa: E402
from dna_decode.data.cell_key import canonical_cell_key  # noqa: E402

# Reframed tier (lineage-disclosure layer): the headline must say BOTH that this is provenance-disjoint
# AND that the R classes are clonally dominated, with the lineage-effective N + cluster-weighted metrics
# (with CI) disclosed in the lineage table. It is NOT lineage-independent external validation.
PROV_TIER = ("isolate-level provenance-disjoint stress test (different submitter/lab/country); R classes "
             "clonally dominated — lineage-effective N + cluster-weighted metrics (with Wilson CI) disclosed "
             "in the lineage table; NOT methodology-independent (most submitters use CLSI broth "
             "microdilution) and NOT lineage-independent external clinical validation")

LINEAGE_SIDECAR = "provdisjoint_lineage_metrics.json"

# M2: the canonical (organism, drug) join key — shared with the lineage sidecar + scored JSONs.
_key = canonical_cell_key


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


def load_lineage_metrics() -> dict:
    """Read the lineage-disclosure sidecar -> {canonical_key: cell}. Empty if absent."""
    p = WIKI / LINEAGE_SIDECAR
    out = {}
    if p.exists():
        try:
            for c in json.loads(p.read_text(encoding="utf-8")).get("cells", []):
                out[_key(c["organism"], c["drug"])] = c
        except Exception:  # noqa: BLE001 — a malformed sidecar must not break the read-only roll-up
            pass
    return out


def _assert_weighted_renderable(w: dict) -> None:
    """C3 emitter guard: a cluster-weighted point estimate may NEVER be rendered without its Wilson CI
    + effective-N. A bare weighted sens/spec is a honesty inversion (tiny-N point with no uncertainty)."""
    for metric in ("sens", "spec"):
        if w.get(metric) is not None:
            ci = w.get(f"{metric}_ci")
            assert isinstance(ci, (list, tuple)) and len(ci) == 2, \
                f"cluster-weighted {metric}={w.get(metric)} rendered without a Wilson CI (C3 violation)"
            assert w.get(f"{metric}_eff_n") is not None, \
                f"cluster-weighted {metric}={w.get(metric)} rendered without effective-N (C3 violation)"


def build_lineage_block(scell: dict | None) -> dict:
    """Project a lineage-sidecar cell into the report-card lineage block.

    A SCORED cell with no sidecar row -> status 'not_computed'; a partial cohort (genomes missing) ->
    status 'incomplete' with k/N; a complete cohort -> status 'scored' with per-threshold effective-N +
    CI-bearing cluster-weighted metrics. Never silently blank (M1/M2)."""
    if scell is None:
        return {"status": "not_computed"}
    if not scell.get("lineage_tier_emitted"):
        return {"status": "incomplete" if scell.get("partial") else "not_computed",
                "raw_N": scell.get("raw_N"),
                "n_genomes_missing": scell.get("n_genomes_missing", 0)}
    eff, cw = {}, {}
    for t, blk in scell.get("thresholds", {}).items():
        eff[t] = {"R": blk["effective_lineage_N_R"], "S": blk["effective_lineage_N_S"]}
        w = blk["cluster_weighted"]
        _assert_weighted_renderable(w)
        cw[t] = {"sens": w["sens"], "sens_ci": w["sens_ci"], "sens_eff_n": w["sens_eff_n"],
                 "spec": w["spec"], "spec_ci": w["spec_ci"], "spec_eff_n": w["spec_eff_n"],
                 "n_discordant": w["n_discordant"]}
    return {"status": "scored", "raw_N": scell.get("raw_N"),
            "effective_lineage_N": eff, "cluster_weighted": cw, "grade": scell.get("lineage_grade")}


def classify(key, scored, census, registry, surface=None) -> dict:
    """Resolve one cell's state. Surface structural-label properties (no_free_source / label_confounded)
    take precedence over observations; otherwise SCORED > ABSTAINS > census > NOT_CENSUSED."""
    status = (surface or {}).get("phenotype_source_status")
    if status == "no_free_source":
        return {"state": "NO_FREE_PHENOTYPE_SOURCE",
                "note": f"{(surface or {}).get('engine','')}; no free isolate-level AST source (structural non-cell)"}
    if status == "label_confounded":
        return {"state": "LABEL_CONFOUNDED",
                "note": "phenotype LABEL is an unreliable surrogate (oxacillin AST vs mecA; cefoxitin is the CLSI surrogate)"}

    if key in scored:
        m = scored[key].get("metrics", {})
        return {"state": "SCORED",
                "acc": m.get("acc"), "sens": m.get("sens"), "spec": m.get("spec"), "n": m.get("n_scored"),
                "tp": m.get("tp"), "fp": m.get("fp"), "tn": m.get("tn"), "fn": m.get("fn"),
                "tier": scored[key].get("independence_tier", PROV_TIER), "file": scored[key].get("_file")}
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
    return {"state": "NOT_CENSUSED", "note": "bacterial + census-able; no provenance census yet"}


def main() -> int:
    scored, census, registry = load_scored(), load_census(), load_registry()
    lineage = load_lineage_metrics()
    surface = surface_index()

    rows = []
    for key in sorted(set(surface) | set(scored) | set(census) | set(registry)):
        c = classify(key, scored, census, registry, surface.get(key))
        # Lineage augments (never demotes) the state machine: only SCORED cells carry a lineage block.
        if c["state"] == "SCORED":
            c["lineage"] = build_lineage_block(lineage.get(key))
        rows.append((key, c))

    counts = {}
    for _, c in rows:
        counts[c["state"]] = counts.get(c["state"], 0) + 1

    today = _date.today().isoformat()
    artifact = {"_schema": "decoder-validation-report-card-v0", "date": today,
                "honest_tier": PROV_TIER, "no_aggregate_headline": True,
                "state_counts": counts,
                "cells": [{"organism": k[0], "drug": k[1], **c} for k, c in rows]}
    (WIKI / "decoder_validation_report_card.json").write_text(json.dumps(artifact, indent=2), encoding="utf-8")

    STATES = ("SCORED", "POWERED_UNSCORED", "UNDERPOWERED", "ABSTAINS_BY_DESIGN",
              "NOT_CENSUSED", "LABEL_CONFOUNDED", "NO_FREE_PHENOTYPE_SOURCE")
    L = []
    L.append(f"# Decoder-suite provenance-disjoint validation report card — {today}\n")
    L.append("Standing trust surface for the shipped deterministic AMR decoders (Anchor-4). Rows are the "
             "DEPLOYED-CLAIM surface (`dna_decode/data/shipped_decoder_surface.py`) unioned with observed "
             "scored/census cells. Each cell is the DEPLOYED `call_resistance(organism, drug)` rule scored on "
             "a FRESH, leakage-checked, **provenance-disjoint** NCBI-PD cohort (submitters OUTSIDE "
             "NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA).\n")
    L.append("> **Honest tier (do NOT inflate):** every SCORED cell is an isolate-level provenance-disjoint "
             "stress test (different submitter/lab/country). The R classes are **clonally dominated** — the "
             "raw-isolate sens/spec is inflated by over-sampled clones, so the **Lineage disclosure** table "
             "below reports lineage-effective N + cluster-weighted sens/spec (one vote per lineage) with a "
             "Wilson CI. It is **NOT** methodology-independent (most submitters use CLSI broth microdilution) "
             "and **NOT** lineage-independent external clinical validation. There is deliberately **no "
             "aggregate “X% validated” number** — read the grid cell by cell.\n")
    L.append("## State legend\n")
    L.append("| state | meaning |\n|---|---|")
    L.append("| `SCORED` | Stage-2 provdisjoint run exists — acc/sens/spec shown |")
    L.append("| `POWERED_UNSCORED` | censused ≥ 20/class both classes; not yet scored |")
    L.append("| `UNDERPOWERED` | censused < 20/class (surveillance-dominated organism) |")
    L.append("| `ABSTAINS_BY_DESIGN` | registry EXPRESSION_FLOOR — rule refuses what it can't decode |")
    L.append("| `NOT_CENSUSED` | bacterial + census-able; no census yet |")
    L.append("| `LABEL_CONFOUNDED` | phenotype label is an unreliable surrogate (oxacillin AST vs mecA) |")
    L.append("| `NO_FREE_PHENOTYPE_SOURCE` | fungal/antiviral/antimalarial — no free isolate-level AST (structural non-cell) |\n")
    L.append("## State counts\n")
    L.append("| state | cells |\n|---|---|")
    for s in STATES:
        if s in counts:
            L.append(f"| `{s}` | {counts[s]} |")
    L.append("\n## Cells\n")
    L.append("| organism | drug | state | acc | sens | spec | n | detail |\n|---|---|---|---|---|---|---|---|")
    for k, c in rows:
        org, drug = k
        if c["state"] == "SCORED":
            L.append(f"| {org} | {drug} | `SCORED` | {c.get('acc')} | {c.get('sens')} | {c.get('spec')} | "
                     f"{c.get('n')} | TP{c.get('tp')} FP{c.get('fp')} TN{c.get('tn')} FN{c.get('fn')} |")
        else:
            L.append(f"| {org} | {drug} | `{c['state']}` | — | — | — | — | {c.get('note','')} |")
    # ---- Lineage disclosure (clonality-corrected) ----
    scored_rows = [(k, c) for k, c in rows if c["state"] == "SCORED"]
    L.append("\n## Lineage disclosure (clonality-corrected)\n")
    L.append("Raw sens/spec counts one vote per ISOLATE; clones inflate it. Below: lineage-effective N "
             "(greedy-representative Mash clustering — chaining-resistant, NOT single-linkage) + "
             "cluster-weighted sens/spec (one vote per same-label lineage; mixed-label clones are "
             "DISCORDANT, never majority-voted) with a 95% Wilson CI. Weighted N is tiny — the CI is the "
             "point. Weighted metrics shown at Mash 0.005 (conservative); the JSON carries 0.001 too.\n")
    L.append("| organism | drug | raw N | eff lineages R/S @.001 | eff lineages R/S @.005 | "
             "wtd sens [95% CI] (n) | wtd spec [95% CI] (n) | discordant | grade |\n"
             "|---|---|---|---|---|---|---|---|---|")
    for k, c in scored_rows:
        org, drug = k
        lin = c.get("lineage", {"status": "not_computed"})
        if lin.get("status") != "scored":
            note = ("lineage: incomplete "
                    f"({lin.get('n_genomes_missing', '?')} genomes missing)"
                    if lin.get("status") == "incomplete" else "lineage: not computed")
            L.append(f"| {org} | {drug} | {lin.get('raw_N', '—')} | — | — | — | — | — | {note} |")
            continue
        eff = lin["effective_lineage_N"]
        cw = lin["cluster_weighted"]

        def _eff(t):
            e = eff.get(t)
            return f"{e['R']}/{e['S']}" if e else "—"

        def _wtd(metric):
            w = cw.get("0.005")
            if not w or w.get(metric) is None:
                return "—"
            lo, hi = w[f"{metric}_ci"]
            return f"{w[metric]} [{lo}–{hi}] (n={w[f'{metric}_eff_n']})"

        disc = cw.get("0.005", {}).get("n_discordant", 0)
        L.append(f"| {org} | {drug} | {lin.get('raw_N')} | {_eff('0.001')} | {_eff('0.005')} | "
                 f"{_wtd('sens')} | {_wtd('spec')} | {disc} | {lin.get('grade', '—')} |")

    L.append("\n## Provenance\n")
    L.append("- Row set: `dna_decode/data/shipped_decoder_surface.py` (deployed-claim surface) ∪ observed cells.")
    L.append("- SCORED cells: `wiki/provenance_disjoint_validation_*.json` (Stage-2 `provenance_disjoint_validate.py`).")
    L.append("- Powering: `wiki/provdisjoint_census_results.json` (Stage-1 `ncbi_pd_provenance_census.py`).")
    L.append("- ABSTAINS: `dna_decode/data/calibrated_amr_rules.json` (EXPRESSION_FLOOR verdicts).")
    L.append("- Lineage disclosure: `wiki/provdisjoint_lineage_metrics.json` (`scripts/compute_lineage_metrics.py`).")
    L.append("- Rebuild: `.venv/Scripts/python.exe scripts/build_validation_report_card.py` (read-only roll-up; re-run as cells land).")
    (WIKI / "decoder_validation_report_card.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    print("report card written: wiki/decoder_validation_report_card.md")
    print(f"state counts: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
