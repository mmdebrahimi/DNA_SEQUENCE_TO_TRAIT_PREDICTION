"""CABBAGE lineage metrics — a NAMESPACE-SEPARATE sidecar feeding the existing G8 gate.

Why this is not written into `wiki/provdisjoint_lineage_metrics.json`
--------------------------------------------------------------------
The four scored CABBAGE cells key on their `amrfinder_organism`, and two of those keys
ALREADY EXIST in the provdisjoint lineage sidecar: `Campylobacter x ciprofloxacin` and
`Escherichia_coli_Shigella x gentamicin`. The collision is MEASURED, not hypothetical.
`scripts/build_validation_report_card.py` reads that sidecar by an exact filename constant
(`LINEAGE_SIDECAR = "provdisjoint_lineage_metrics.json"`), so the real guard is PATH
INEQUALITY — there is no glob to avoid. Writing CABBAGE cells under `canonical_cell_key`
into that file would silently overwrite two provdisjoint numbers; this module cannot
address that path at all.

The bar is IMPORTED, never restated
-----------------------------------
`MIN_EFFECTIVE_LINEAGES` comes from `rejection_gates`. A literal `3` here would let the
two drift apart, and the gate — not this module — owns the verdict.

Merge-or-refuse, defined before the first increment
---------------------------------------------------
Cells accrue one at a time (each needs its own assembly fetch + Mash run), so this sidecar
is a STAGED artifact. "Write what this run computed" is correct for a one-shot and silently
lossy for a staged one, and the lossy output is indistinguishable from a complete one. So a
write that would drop a cell already on disk RAISES rather than truncating.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from dna_decode.eval.rejection_gates import MIN_EFFECTIVE_LINEAGES

SCHEMA = "cabbage-lineage-metrics-v1"

#: This module's own sidecar. Relative to the project root.
SIDECAR_NAME = "cabbage_lineage_metrics.json"

#: The file this module must NEVER write. Named so the guard is explicit and testable.
FORBIDDEN_SIDECAR_NAME = "provdisjoint_lineage_metrics.json"

#: Inherited verbatim from `scripts/compute_lineage_metrics.py:50`. NOT searched or tuned:
#: `variance_ratio = intra/inter` tends to 0 as clusters tend to singletons, so minimising it
#: is monotonically biased toward over-splitting. A fixed grid is the documented fix.
THRESHOLDS: tuple[float, ...] = (0.001, 0.005)


class SidecarPathRefused(Exception):
    """Raised when a write would target the provdisjoint sidecar."""


class CellDropRefused(Exception):
    """Raised when a staged write would drop a cell already present on disk."""


class WeightedPointNotRenderable(Exception):
    """Raised when a cluster-weighted point is emitted without its interval."""


def cabbage_cell_key(species: str, drug: str) -> str:
    """Key a CABBAGE cell by its SPECIES (not `amrfinder_organism`) plus drug.

    Deliberately distinct in shape from `cell_key.canonical_cell_key`, which keys on the
    AMRFinder organism. Two of the four CABBAGE cells collide under that key; keying on the
    species string keeps `Salmonella enterica x ciprofloxacin` and
    `Salmonella enterica x ceftriaxone` distinguishable from each other and from the
    provdisjoint `Salmonella` rows, and makes an accidental cross-file merge fail loudly
    on an unrecognised key rather than quietly succeed.
    """
    return f"{species.strip()} x {drug.strip().lower()}"


def sidecar_path(root: str | Path = ".") -> Path:
    """Resolve this module's sidecar, refusing the forbidden filename by construction."""
    p = Path(root) / "wiki" / SIDECAR_NAME
    if p.name == FORBIDDEN_SIDECAR_NAME:
        raise SidecarPathRefused(f"refusing to target {FORBIDDEN_SIDECAR_NAME}")
    return p


def assert_not_forbidden(path: str | Path) -> None:
    """Guard any caller-supplied path against the provdisjoint sidecar."""
    if Path(path).name == FORBIDDEN_SIDECAR_NAME:
        raise SidecarPathRefused(
            f"{path} is the provdisjoint lineage sidecar; two CABBAGE cells share keys with "
            "it, so a write here would silently overwrite provdisjoint numbers"
        )


def load_cells(root: str | Path = ".") -> list[dict[str, Any]]:
    """Read the sidecar's cells, or [] when it does not exist yet."""
    p = sidecar_path(root)
    if not p.exists():
        return []
    doc = json.loads(p.read_text(encoding="utf-8"))
    return list(doc.get("cells", []))


def merge_cells(
    existing: Sequence[Mapping[str, Any]],
    incoming: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Merge by cell key — incoming wins per key, nothing is ever dropped.

    Idempotent on replay: merging the same incoming set twice yields the same result.
    """
    by_key: dict[str, dict[str, Any]] = {}
    for cell in list(existing) + list(incoming):
        by_key[cell["key"]] = dict(cell)
    return [by_key[k] for k in sorted(by_key)]


def assert_no_cell_dropped(
    prior: Sequence[Mapping[str, Any]],
    merged: Sequence[Mapping[str, Any]],
) -> None:
    """Refuse a staged write that loses a cell already on disk."""
    lost = {c["key"] for c in prior} - {c["key"] for c in merged}
    if lost:
        raise CellDropRefused(
            "refusing to write: would drop cell(s) already on disk: "
            + ", ".join(sorted(lost))
            + ". This sidecar accrues one cell per run, so a partial write is "
              "indistinguishable from a complete one."
        )


def write_cells(
    cells: Sequence[Mapping[str, Any]],
    root: str | Path = ".",
    *,
    generated: str | None = None,
) -> Path:
    """Merge the supplied cells into the sidecar and write it, refusing any cell drop."""
    p = sidecar_path(root)
    assert_not_forbidden(p)
    prior = load_cells(root)
    merged = merge_cells(prior, cells)
    assert_no_cell_dropped(prior, merged)
    doc = {
        "schema": SCHEMA,
        "generated": generated,
        "bar_source": "dna_decode.eval.rejection_gates.MIN_EFFECTIVE_LINEAGES",
        "min_effective_lineages_bar": MIN_EFFECTIVE_LINEAGES,
        "thresholds": list(THRESHOLDS),
        "cells": merged,
    }
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Cluster STRUCTURE, not just cluster counts.
#
# The TB arm reported 43 discordant clusters and a healthy-looking cluster count
# while those clusters hid 97% of the isolates, leaving sens/spec on a ~3%
# non-random residue. A count alone cannot show that; largest-cluster fraction and
# discordant-ISOLATE fraction can.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class WeightedPoint:
    """A cluster-weighted rate that cannot be rendered without its interval.

    Weighted N is small by construction (one vote per lineage), so the interval is
    the point. `value=None` means not computed — that is renderable; a value with
    no interval is not.
    """

    label: str
    value: float | None
    ci: tuple[float, float] | None = None
    k: int | None = None
    n: int | None = None


def assert_weighted_renderable(point: WeightedPoint) -> None:
    """Refuse a weighted point that carries a value but no interval."""
    if point.value is not None and point.ci is None:
        raise WeightedPointNotRenderable(
            f"weighted point {point.label!r} has value {point.value} but no Wilson "
            "interval; weighted N is small, so the interval is the point"
        )


def weighted_point(label: str, k: int, n: int) -> WeightedPoint:
    """Build a weighted point with its Wilson interval always attached."""
    from dna_decode.eval.clonality import wilson_ci

    if n <= 0:
        return WeightedPoint(label=label, value=None, ci=None, k=k, n=n)
    pt = WeightedPoint(label=label, value=round(k / n, 4), ci=wilson_ci(k, n), k=k, n=n)
    assert_weighted_renderable(pt)
    return pt


def minority_class(labels: Mapping[str, Any]) -> str:
    """Which class has fewer ISOLATES — the class G8's "smallest class" refers to.

    Ties resolve to "R", which is the minority in all four CABBAGE cells; the tie
    branch exists so the function is total, not because a tie is expected.
    """
    from dna_decode.eval.clonality import _normalize_label

    norm = [_normalize_label(v) for v in labels.values()]
    n_r = sum(1 for x in norm if x == "R")
    n_s = sum(1 for x in norm if x == "S")
    return "S" if n_s < n_r else "R"


def cluster_structure(
    clusters: Mapping[str, int],
    labels: Mapping[str, Any],
) -> dict[str, Any]:
    """Effective lineage counts PLUS the structure figures that make them readable.

    Discordant (mixed-label) clusters are excluded from both class counts and
    reported separately — `clonality.cluster_class` returns DISCORDANT and is never
    majority-voted.
    """
    from dna_decode.eval.clonality import cluster_class, cluster_members, effective_lineage_n

    if not clusters:
        return {
            "n_isolates": 0, "n_clusters": 0,
            "largest_cluster_fraction": None,
            "n_discordant_clusters": 0, "n_discordant_isolates": 0,
            "discordant_isolate_fraction": None,
            "effective_lineage_N_R": 0, "effective_lineage_N_S": 0,
        }

    members = cluster_members(dict(clusters))
    n_isolates = len(clusters)
    sizes = [len(v) for v in members.values()]

    discordant = {
        cid: sids for cid, sids in members.items()
        if cluster_class([labels[s] for s in sids]) == "DISCORDANT"
    }
    n_disc_iso = sum(len(v) for v in discordant.values())

    return {
        "n_isolates": n_isolates,
        "n_clusters": len(members),
        "largest_cluster_fraction": round(max(sizes) / n_isolates, 4),
        "n_discordant_clusters": len(discordant),
        "n_discordant_isolates": n_disc_iso,
        "discordant_isolate_fraction": round(n_disc_iso / n_isolates, 4),
        "effective_lineage_N_R": effective_lineage_n(dict(clusters), dict(labels), "R"),
        "effective_lineage_N_S": effective_lineage_n(dict(clusters), dict(labels), "S"),
    }


def min_effective_lineages(
    clusters: Mapping[str, int],
    labels: Mapping[str, Any],
) -> tuple[int, str]:
    """The evidence G8 consumes: the MINORITY class's effective lineage count.

    Returns (count, which_class). See `wiki/cabbage_g8_acceptance_bar.json` — this is
    a scoped reading of "smallest class", since the majority class is not clustered.
    """
    cls = minority_class(labels)
    struct = cluster_structure(clusters, labels)
    return struct[f"effective_lineage_N_{cls}"], cls
