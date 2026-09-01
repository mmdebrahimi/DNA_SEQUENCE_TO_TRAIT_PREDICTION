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

# The report card reads its AMR grid from the Evidence-Contract Registry (v0.1). registry.surface_index()
# re-exports the frozen-surface-shaped dict FROM the registry's AMR cells (== shipped_decoder_surface by
# construction; pinned by tests/test_cell_registry.py) so the registry is the single source the card reads.
from dna_decode.data.cell_registry import surface_index  # noqa: E402
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


def invisible_fraction_from_metrics(metrics: dict) -> float | None:
    """Determinant-invisible fraction for a SCORED cell = of the measured-R isolates the cell SCORED (not
    abstained), the fraction it calls non-R (fn / (tp + fn) = 1 − sens). The honest 'how much resistance
    this cell structurally misses' number, DESCRIPTIVE (not an endorsement input). None when no measured-R
    was scored. NOTE: the truly-invisible vs rule-limited split (per-isolate determinants) is available only
    in the NCBI-PD atlas (`wiki/determinant_blindness_atlas.*`); the frozen provdisjoint JSONs store only the
    confusion counts, so this is the aggregate fraction."""
    tp, fn = metrics.get("tp"), metrics.get("fn")
    if tp is None or fn is None:
        return None
    denom = tp + fn
    return round(fn / denom, 3) if denom else None


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


def load_source_concentration() -> dict:
    """Read the source-concentration measurement -> {canonical_key: cell}. Empty if absent.

    NAMESPACE-SEPARATE for the same reason as `load_prospective`: this shares its (organism, drug) key
    with a provenance-disjoint cell, so merging it would silently overwrite one number with the other --
    the shared-key trap. It AUGMENTS a cell; it never replaces a metric and never changes a state.

    WHAT IT ADDS, and why it is not redundant with the lineage layer. Lineage discloses clonal domination
    WITHIN a cohort. This discloses how many independent SOURCES the cohort draws on at all, which bounds
    what the cohort could ever have detected. Measured case: the e.coli x gentamicin cell is 95% one
    BioProject and contains ZERO carriers of the `rmt` determinant family, so it reported sens 0.893 while
    two source-diverse measurements of the same cell reported 0.429 and 0.523. The number was never wrong
    about its cohort; it was a statement about one hospital's isolates.
    """
    f = WIKI / "provdisjoint_source_concentration.json"
    if not f.exists():
        return {}
    try:
        d = json.loads(f.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001 — a malformed sidecar must not break the read-only roll-up
        return {}
    if not d.get("complete"):
        # A partial provenance sweep cannot support a concentration claim; render nothing rather than a
        # floor that reads like a measurement.
        return {}
    return {_key(c["organism"], c["drug"]): c for c in d.get("cells", []) if c.get("organism")}


def load_doubt_layer() -> dict:
    """Read the L2 doubt-layer measurement -> {drug: cell}. Empty if absent.

    NAMESPACE-SEPARATE, and keyed by DRUG rather than by (organism, drug) — which is the honest shape,
    because the completeness screen runs across the whole cached index and is NOT a per-cohort
    measurement. `build_doubt_block` stamps that scope so a drug-level result can never be read as a
    statement about one cell's cohort.

    WHAT IT ADDS, and why it is not redundant with the three layers above. Lineage discloses clonal
    domination inside a cohort; source concentration discloses how few sources the cohort drew on;
    prospective discloses a temporally-clean re-score. All three ask *how good is the evidence for this
    cell's number*. This asks a different question: **does the deployed RULE fail to represent a
    determinant family present in the data at all?** — the completeness failure that produced the
    gentamicin `rmt` blind spot, which no amount of better cohort evidence would have surfaced.
    """
    hits = sorted(WIKI.glob("doubt_layer_per_cell_*.json"))
    if not hits:
        return {}
    try:
        d = json.loads(hits[-1].read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001 — a malformed sidecar must not break the read-only roll-up
        return {}
    return {str(c["drug"]).strip().lower(): c
            for c in d.get("determinant_completeness_arm", []) if c.get("drug")}


def build_doubt_block(dcell: dict | None) -> dict:
    """Project the doubt measurement into a report-card block. Never silently blank."""
    if dcell is None:
        return {"status": "not_measured",
                "note": "the completeness screen has not run for this drug — unassessable, not clean"}
    if dcell.get("status") != "scored":
        return {"status": dcell.get("status", "unknown"), "note": dcell.get("note", "")}
    strong = dcell.get("strong") or []
    return {
        "status": "measured",
        "scope": "drug-level across the whole cached determinant index — NOT this cell's cohort",
        "n_families_uncounted": dcell.get("n_families_uncounted"),
        "n_strong": dcell.get("n_strong", 0),
        "n_raw_signature": dcell.get("n_raw_signature"),
        "strong_families": [s["evidence"]["symbol"] for s in strong],
        "known_gap_recovered": dcell.get("known_gap_recovered", False),
        "note": ("families the deployed rule cannot represent whose labelled carriers are uniformly "
                 "resistant, family-wise corrected. Changes no metric and no cell state."),
    }


SINGLE_SOURCE_SHARE = 0.80


def build_source_block(scell: dict | None) -> dict:
    """Project the concentration measurement into a report-card block. Never silently blank."""
    if scell is None:
        return {"status": "not_measured"}
    bp = scell.get("bioproject") or {}
    share = bp.get("largest_share")
    blk = {"status": "measured",
           "n_cohort": scell.get("n_cohort"),
           "distinct_bioprojects": bp.get("distinct"),
           "largest_share": share,
           "dominant": (bp.get("largest") or [None])[0],
           "n_unknown_provenance": bp.get("n_unknown"),
           "distinct_centers": (scell.get("sra_center") or {}).get("distinct")}
    # A DISCLOSURE flag, not a demotion: the cell's state and metrics are untouched. It says the estimate
    # rests on one source, which is a different fact from the estimate being wrong.
    blk["single_source"] = bool(share is not None and share >= SINGLE_SOURCE_SHARE)
    return blk


def load_prospective() -> dict:
    """Read the prospective-lock scored artifacts -> {canonical_key: cell}. Empty if none.

    NAMESPACE-SEPARATE from `load_scored()` on purpose. A prospective cell shares its (organism, drug)
    key with a provenance-disjoint cell, so merging them would silently overwrite one with the other --
    the documented shared-key trap that Fix C in the external-cohort arm exists to avoid. Prospective
    AUGMENTS a cell (its own block + its own table); it never replaces the provdisjoint number, and never
    changes a cell's state.

    Keeps the NEWEST artifact per cell (they accrue over time and are re-scored as N grows).
    """
    out: dict = {}
    for f in sorted(glob.glob(str(WIKI / "prospective_lock_validation_*.json"))):
        try:
            d = json.loads(Path(f).read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001 — a malformed artifact must not break the read-only roll-up
            continue
        org, drug = d.get("organism"), d.get("drug")
        if not (org and drug):
            continue
        k = _key(org, drug)
        prev = out.get(k)
        if prev is None or str(d.get("generated", "")) >= str(prev.get("generated", "")):
            out[k] = d
    return out


def build_prospective_block(pcell: dict | None) -> dict:
    """Project a prospective-lock artifact into the report-card block. Never silently blank."""
    if pcell is None:
        return {"status": "not_accrued"}
    if not pcell.get("prospective_lock_verified"):
        # the scorer hard-fails on drift, so this should be unreachable -- but a stale artifact from a
        # drifted decoder must never be rendered as if it validated the CURRENT one.
        return {"status": "lock_unverified", "generated": pcell.get("generated")}

    # SUPERSEDED-BY-SURFACE-CHANGE (added 2026-08-31 with the gentamicin v2 lock).
    # `prospective_lock_verified` records that the lock held WHEN THE ARTIFACT WAS WRITTEN. It is not
    # re-checked afterwards, so once the decoder is revised every prior prospective score keeps saying
    # True while silently describing a RETIRED rule. That is the whole failure the lock exists to
    # prevent, arriving through the back door. Re-verify against the LIVE surface and withhold the
    # numbers on mismatch.
    #
    # This resets the prospective clock rather than merely relabelling: the v2 cutoff is later, so the
    # isolates that were prospective for v1 are PRE-lock for v2 and cannot be re-scored into evidence.
    # That cost is real and is the honest price of the revision -- it applies to the ciprofloxacin cell
    # too, whose RULE did not change but whose pinned surface did; behavioural sameness is an argument,
    # and hash-pinning exists so evidence never rests on one.
    stamped = ((pcell.get("lock_manifest") or {}).get("surface_sha256")) or {}
    if stamped:
        try:
            from dna_decode.eval.prospective_lock import surface_hashes
            live = surface_hashes()
            drifted = sorted(f for f, h in stamped.items() if live.get(f) != h)
        except Exception:  # noqa: BLE001 — a read-only roll-up must not break on an import failure
            drifted = []
        if drifted:
            return {"status": "superseded_by_surface_change",
                    "generated": pcell.get("generated"),
                    "lock_date": (pcell.get("lock_manifest") or {}).get("lock_date"),
                    "drifted_files": drifted,
                    "note": ("scored against a decoder that has since been revised; its numbers describe "
                             "the RETIRED rule and are withheld. The prospective clock restarts at the "
                             "new lock date -- these isolates are pre-lock for the current decoder.")}
    conf = pcell.get("confusion") or {}
    pw = pcell.get("powering") or {}
    blk = {"status": "scored" if conf.get("n_scored") else "not_accrued",
           "generated": pcell.get("generated"),
           "lock_date": (pcell.get("lock_manifest") or {}).get("lock_date"),
           "n_scored": conf.get("n_scored"), "R": pw.get("scored_R"), "S": pw.get("scored_S"),
           "acc": conf.get("acc"), "sens": conf.get("sens"), "spec": conf.get("spec"),
           "abstain": conf.get("abstain"), "powering": pw.get("status")}
    blk["regression"] = _prospective_regression(blk)
    return blk


# CONVENTION, not a derived truth: a POWERED prospective cell whose sensitivity is below this detects
# fewer than half of resistant isolates. It is chosen because it needs NO comparator and no distributional
# assumption -- a threshold tuned against one cohort's provdisjoint delta would be exactly the
# single-cohort inference this project's own gates reject. Stated here so it is arguable and changeable,
# not buried.
PROSPECTIVE_SENS_FLOOR = 0.5


def _prospective_regression(blk: dict) -> bool:
    """Does this prospective result contradict the cell's SCORED standing? Read-only, no state change."""
    if blk.get("status") != "scored" or blk.get("powering") != "POWERED":
        return False        # an accruing / underpowered cell makes no claim either way
    sens = blk.get("sens")
    return sens is not None and sens < PROSPECTIVE_SENS_FLOOR


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
                "invisible_fraction": invisible_fraction_from_metrics(m),
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


def load_naive_value_add() -> list[dict]:
    """Curated-vs-naive value-add rows from the committed naive-comparator JSONs (2026-06-27).

    The wrapper-vs-underlying-tool rail: the deployed call_resistance rule must BEAT naive AMRFinder use
    ('any drug-class determinant -> R') on balanced accuracy. Two surfaces:
      - Oxford independent measured-MIC  (external_validation_oxford_naive_comparator_*.json)
      - the 10-cell provenance-disjoint  (provdisjoint_naive_comparator_*.json)
    Read-only; returns [] if the JSONs are absent.
    """
    rows: list[dict] = []
    ox = sorted(glob.glob(str(WIKI / "external_validation_oxford_naive_comparator_*.json")))
    if ox:
        try:
            d = json.loads(Path(ox[-1]).read_text(encoding="utf-8"))
            for drug, m in (d.get("drugs") or {}).items():
                rows.append({"surface": "Oxford ext. MIC", "organism": "Escherichia_coli_Shigella",
                             "drug": drug, "frozen_balacc": (m.get("frozen") or {}).get("balacc"),
                             "naive_balacc": (m.get("naive") or {}).get("balacc"),
                             "delta": (m.get("delta_frozen_minus_naive") or {}).get("balacc"),
                             "verdict": m.get("value_add_verdict")})
        except (OSError, json.JSONDecodeError):
            pass
    pd = sorted(glob.glob(str(WIKI / "provdisjoint_naive_comparator_*.json")))
    if pd:
        try:
            d = json.loads(Path(pd[-1]).read_text(encoding="utf-8"))
            for key, m in (d.get("cells") or {}).items():
                org, _, drug = key.partition(":")
                rows.append({"surface": "provdisjoint", "organism": org, "drug": drug,
                             "frozen_balacc": m.get("frozen_balacc"), "naive_balacc": m.get("naive_balacc"),
                             "delta": m.get("delta_balacc"), "verdict": m.get("value_add_verdict")})
        except (OSError, json.JSONDecodeError):
            pass
    return rows


def main() -> int:
    scored, census, registry = load_scored(), load_census(), load_registry()
    lineage = load_lineage_metrics()
    prospective = load_prospective()
    source_conc = load_source_concentration()
    doubt = load_doubt_layer()
    surface = surface_index()

    rows = []
    for key in sorted(set(surface) | set(scored) | set(census) | set(registry) | set(prospective)):
        c = classify(key, scored, census, registry, surface.get(key))
        # Lineage augments (never demotes) the state machine: only SCORED cells carry a lineage block.
        if c["state"] == "SCORED":
            c["lineage"] = build_lineage_block(lineage.get(key))
        # Prospective likewise AUGMENTS: attached wherever an artifact exists, regardless of state, so a
        # prospective result can never be hidden by the cell's provdisjoint state.
        # Source concentration AUGMENTS too, and only where it was measured. It never demotes: a
        # single-source cell keeps its state and its published metrics.
        if key in source_conc:
            c["source_concentration"] = build_source_block(source_conc.get(key))
        # L2 DOUBT augments too, under its OWN key. Keyed by DRUG (the screen is index-wide, not
        # per-cohort), and attached only where the drug was actually screened -- an unmeasured drug
        # gets no block rather than one reading as a clean bill.
        if key[1] in doubt:
            c["doubt_layer"] = build_doubt_block(doubt.get(key[1]))
        if key in prospective:
            c["prospective"] = build_prospective_block(prospective.get(key))
            # TOP-LEVEL so a consumer filtering on `state == SCORED` cannot miss it. The state itself is
            # deliberately NOT demoted -- the provenance-disjoint result is still what it was; this says
            # the cell has a CONTRADICTING prospective result, which is a different fact.
            if c["prospective"].get("regression"):
                c["prospective_regression"] = True
                c["deployment_caveat"] = (
                    f"prospective sens {c['prospective'].get('sens')} on {c['prospective'].get('n_scored')} "
                    f"post-lock isolates -- the frozen rule under-calls this cell; see the "
                    f"Prospective-lock disclosure table")
        rows.append((key, c))

    counts = {}
    for _, c in rows:
        counts[c["state"]] = counts.get(c["state"], 0) + 1

    today = _date.today().isoformat()
    artifact = {"_schema": "decoder-validation-report-card-v0", "date": today,
                "honest_tier": PROV_TIER, "no_aggregate_headline": True,
                "state_counts": counts,
                "curated_vs_naive_value_add": load_naive_value_add(),
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
    L.append("`blind.` = determinant-invisible fraction (of the scored measured-R, the fraction the cell "
             "calls non-R = FN/(TP+FN) = 1−sens) — the honest 'how much resistance this cell structurally "
             "misses'; DESCRIPTIVE, not an endorsement input. The truly-invisible vs rule-limited split is "
             "in `wiki/determinant_blindness_atlas.md` (NCBI-PD cells).\n")
    L.append("| organism | drug | state | acc | sens | spec | n | blind. | detail |"
             "\n|---|---|---|---|---|---|---|---|---|")
    for k, c in rows:
        org, drug = k
        if c["state"] == "SCORED":
            inv = c.get("invisible_fraction")
            L.append(f"| {org} | {drug} | `SCORED` | {c.get('acc')} | {c.get('sens')} | {c.get('spec')} | "
                     f"{c.get('n')} | {'—' if inv is None else inv} | "
                     f"TP{c.get('tp')} FP{c.get('fp')} TN{c.get('tn')} FN{c.get('fn')} |")
        else:
            L.append(f"| {org} | {drug} | `{c['state']}` | — | — | — | — | — | {c.get('note','')} |")
    # ---- Source-concentration disclosure (how many independent sources back each number) ----
    src_rows = [(k, c) for k, c in rows if c.get("source_concentration", {}).get("status") == "measured"]
    if src_rows:
        L.append("\n## Source-concentration disclosure (how many sources back each SCORED number)\n")
        L.append("The lineage table above corrects for clonal domination WITHIN a cohort. This asks the "
                 "question one level up: how many independent SOURCES does the cohort draw on at all? "
                 "Every cell here is provenance-DISJOINT from the tuning data — that was the design goal "
                 "and it was met. Provenance-DIVERSE is a different property, was never claimed, and is "
                 "what this measures.\n")
        L.append("WHY IT MATTERS, measured: `escherichia_coli_shigella x gentamicin` is 95% one "
                 "BioProject and contains ZERO carriers of the `rmt` determinant family. It reports "
                 "sens 0.893; two source-diverse measurements of the same cell with the same frozen rule "
                 "report 0.429 and 0.523. A cohort with no carriers of a determinant family cannot detect "
                 "a rule blind to that family.\n")
        L.append("The error is NOT directional: a 97%-single-source cipro cell reads PESSIMISTIC "
                 "(spec 0.700 vs 0.988 on an 8-BioProject set). A single-site estimate is an estimate of "
                 "that site. **These rows change no metric and no cell state.**\n")
        L.append("| organism | drug | N | BioProjects | largest share | dominant | unknown provenance |\n"
                 "|---|---|---|---|---|---|---|")
        for k, c in src_rows:
            org, drug = k
            s = c["source_concentration"]
            share = s.get("largest_share")
            share_s = "—" if share is None else f"{share:.0%}"
            flag = "  **SINGLE-SOURCE**" if s.get("single_source") else ""
            L.append(f"| {org} | {drug} | {s.get('n_cohort', '—')} | {s.get('distinct_bioprojects', '—')} "
                     f"| {share_s}{flag} | {s.get('dominant') or '—'} | {s.get('n_unknown_provenance', 0)} |")
        n_single = sum(1 for _, c in src_rows if c["source_concentration"].get("single_source"))
        L.append(f"\n**{n_single} of {len(src_rows)}** cells rest on ONE BioProject holding "
                 f"≥{int(SINGLE_SOURCE_SHARE * 100)}% of the cohort.\n")

    # ---- L2 doubt disclosure (does the RULE fail to represent a determinant family?) ----
    doubt_rows = [(k, c) for k, c in rows if c.get("doubt_layer", {}).get("status") == "measured"]
    if doubt_rows:
        L.append("\n## Catalog-completeness disclosure (L2 doubt — can the RULE even represent it?)\n")
        L.append("The three tables above all ask *how good is the evidence for this cell's number*. "
                 "This asks a different question: **does the deployed rule fail to represent a "
                 "determinant family that is present in the data at all?** That is the failure which "
                 "produced the gentamicin `rmt` blind spot, and no amount of better cohort evidence "
                 "would have surfaced it.\n")
        L.append("Rows are **drug-level** — the screen runs across the whole cached determinant index, "
                 "NOT this cell's cohort. `STRONG` means a family the rule cannot represent whose "
                 "labelled carriers are uniformly resistant, **after a family-wise correction** over "
                 "the families screened for that drug. The correction is load-bearing: the raw purity "
                 "signature fires on 5 families and 4 are coincidences (cipro `qnrA1` at 4R/0S is "
                 "p=0.030 against ~125 families screened).\n")
        L.append("**These rows change no metric and no cell state.**\n")
        L.append("| organism | drug | families rule can't represent | raw signature | **STRONG** | "
                 "families |\n|---|---|---|---|---|---|")

        for k, c in doubt_rows:
            org, drug = k
            d = c["doubt_layer"]
            fams = ", ".join(f"`{s}`" for s in d.get("strong_families") or []) or "—"
            mark = "  **KNOWN GAP**" if d.get("known_gap_recovered") else ""
            L.append(f"| {org} | {drug} | {d.get('n_families_uncounted', '—')} | "
                     f"{d.get('n_raw_signature', '—')} | {d.get('n_strong', 0)}{mark} | {fams} |")

        # Per DRUG, not per row: a drug spans several organism cells and summing rows would multiply
        # one screen result by however many cells happen to share that drug.
        per_drug = {k[1]: c["doubt_layer"].get("n_strong", 0) for k, c in doubt_rows}
        n_strong = sum(per_drug.values())
        noun = "family survives" if n_strong == 1 else "families survive"
        L.append(f"\nAcross {len(per_drug)} screened drugs, **{n_strong or 'no'}** determinant {noun} "
                 "the correction. **Honest limit:** this project has exactly ONE independently "
                 "confirmed completeness gap, so recovering it is a single case and **not a rate** — "
                 "it bounds nothing about gaps never confirmed.\n")

    # ---- Prospective-lock disclosure (temporal) ----
    prosp_rows = [(k, c) for k, c in rows if c.get("prospective")]
    if prosp_rows:
        L.append("\n## Prospective-lock disclosure (temporal — leakage-free BY CONSTRUCTION)\n")
        L.append("A SEPARATE arm from the provenance-disjoint numbers above, not a replacement for them. "
                 "Every isolate here became public STRICTLY AFTER the decoder was frozen and sha256-pinned "
                 "(`wiki/prospective_lock_manifest_2026-06-22.json`), so the decoder cannot have been tuned "
                 "to it — the leakage argument is temporal, not statistical. `verify_lock` re-hashes the "
                 "live decoder on every scoring run and hard-fails on drift.\n")
        L.append("HONEST SCOPE: N is small and ACCRUES over time; this is a temporal stress test, NOT "
                 "lineage-independent clinical validation, and these rows are NOT clonality-corrected "
                 "(the lineage table above applies to the provdisjoint cohorts only).\n")
        L.append("| organism | drug | lock date | N (R/S) | acc | sens | spec | abstain | powering | as of |\n"
                 "|---|---|---|---|---|---|---|---|---|---|")
        for k, c in prosp_rows:
            org, drug = k
            p = c["prospective"]
            if p.get("status") != "scored":
                L.append(f"| {org} | {drug} | — | — | — | — | — | — | — | {p.get('status')} |")
                continue

            def _f(v):
                return "—" if v is None else f"{v:.3f}"

            flag = "  **REGRESSION**" if p.get("regression") else ""
            L.append(f"| {org} | {drug} | {p.get('lock_date', '—')} | "
                     f"{p.get('n_scored')} ({p.get('R')}R/{p.get('S')}S) | {_f(p.get('acc'))} | "
                     f"{_f(p.get('sens'))}{flag} | {_f(p.get('spec'))} | {p.get('abstain')} | "
                     f"{p.get('powering')} | {p.get('generated')} |")
        L.append("\nA LOW prospective sens with HIGH spec means the rule under-calls — it is missing "
                 "determinants, not mislabelling. Diagnose the false negatives' features before reading it "
                 "as decay; see `wiki/prospective_lock_first_accrual_2026-08-24.md`, where exactly that "
                 "diagnosis located a real catalog gap rather than drift.\n")

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

    # ---- Curated layer value-over-naive-baseline (wrapper-vs-underlying-tool rail) ----
    nva = load_naive_value_add()
    if nva:
        L.append("\n## Curated layer value-over-naive-baseline\n")
        L.append("The deployed `call_resistance` rule vs NAIVE AMRFinder use ('any drug-class determinant → R', "
                 "no subclass/point/threshold refinement) on the SAME labels, balanced accuracy. The curated "
                 "layer must BEAT naive tool use on INDEPENDENT data, else the number only proves the tool "
                 "works (the validate-wrapper-vs-underlying-tool rail). Reconciled cells only.\n")
        L.append("| surface | organism | drug | frozen balacc | naive balacc | Δ | verdict |\n"
                 "|---|---|---|---|---|---|---|")
        for r in nva:
            if r.get("verdict") == "RECONCILE_MISMATCH":
                continue
            L.append(f"| {r['surface']} | {r['organism']} | {r['drug']} | {r.get('frozen_balacc')} | "
                     f"{r.get('naive_balacc')} | {r.get('delta')} | {r.get('verdict')} |")
        adds = sum(1 for r in nva if r.get("verdict") == "CURATED_LAYER_ADDS_VALUE")
        ties = sum(1 for r in nva if r.get("verdict") == "NAIVE_TIES_CURATED")
        L.append(f"\n_{adds} cells the curated layer adds value, {ties} ties, "
                 f"{sum(1 for r in nva if r.get('verdict') == 'NAIVE_BEATS_CURATED')} naive-beats. "
                 "Sources: `wiki/external_validation_oxford_naive_comparator_*.json` + "
                 "`wiki/provdisjoint_naive_comparator_*.json`; full synthesis "
                 "`wiki/curated_vs_naive_value_add_synthesis_2026-06-27.md`._")
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
