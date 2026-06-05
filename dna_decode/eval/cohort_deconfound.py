"""Cohort de-confound gate — a PRECONDITION for any embedding-vs-classical falsifier.

Lesson (2026-06-04, EP-4/AMR): a sampling-independent label TYPE (e.g. AMR MIC) does NOT
guarantee a de-confounded COHORT. If the two classes are drawn from disjoint lineages /
geographies / studies, a classifier "predicting" the label actually predicts batch. Pathotype
died on this; cef gate_b was caught by it; cipro N=147 passed.

CONSTRUCT-VALID CRITERION (hardened 2026-06-04 after review): the right question is NOT "does a
provenance axis predict the label" — for AMR, resistance is *legitimately* clonal (resistant
lineages spread), so lineage always predicts resistance somewhat, and an "axis predicts label"
rule would over-block. The right question is **within-group label CONTRAST**: does the same
lineage / country / year contain BOTH classes, with enough matched strains to learn from? A
provenance axis with (near-)zero within-group contrast = the cohort is separable by that axis =
confounded. (A permutation / Cramér's-V association screen is a deferred refinement; it is the
wrong primary statistic here precisely because clonal resistance inflates association without
implying the verdict is batch.)

VERDICT is a 3-state PROMOTABILITY contract:
  - DE_CONFOUNDED  (= ADMIT)           → cohort screen PASSED; may back a promotable falsifier verdict.
                                          NECESSARY, NOT SUFFICIENT — it proves within-group contrast
                                          exists on every available axis, not that the model verdict is
                                          pure biology.
  - WARN           (= DIAGNOSTIC_ONLY) → run a falsifier for diagnosis but the result is NON-PROMOTABLE
                                          (primary contrast thin, OR a secondary provenance axis aliases).
  - CONFOUNDED     (= BLOCK)           → no within-lineage contrast → refuse; no verdict.

Pure functions over plain lists → exhaustively unit-testable; no I/O. `confound_report_for_cohort`
pulls (label, MLST, country, year) from a Cohort for a given drug.
"""
from __future__ import annotations

from collections import defaultdict

# verdict constants (values kept stable for back-compat; meaning is the 3-state contract above)
DE_CONFOUNDED = "DE_CONFOUNDED"   # ADMIT — promotable
WARN = "WARN"                     # DIAGNOSTIC_ONLY — non-promotable
CONFOUNDED = "CONFOUNDED"         # BLOCK — refuse

# thresholds (engineering defaults, NOT calibrated science at N~50-150 — see brainstorm).
MIN_SHARED_LINEAGES = 3           # primary axis must have >= this many both-class lineages
MIN_MATCHED_MINORITY = 5          # min of (minority-class strains inside shared groups); kills the
                                  # "6 lineages x 1R+1S + big tail" weak-pass (matched=6 ok; 1R+1S*3=3 not)
CLEAN_SHARED_LINEAGES = 5         # primary >= this (with matched support) → clean; 3-4 → WARN

_MISSING_TOKENS = {"", "none", "nan", "unknown", "na", "n/a", "null"}


def _norm(x):
    """Normalize a grouping value; missing/placeholder → None (excluded from contrast credit)."""
    if x is None:
        return None
    s = str(x).strip()
    return None if s.lower() in _MISSING_TOKENS else s


def axis_contrast(labels, groups) -> dict:
    """Within-group label contrast for one provenance axis. Missing groups excluded."""
    members = defaultdict(list)
    for y, g in zip(labels, groups):
        gn = _norm(g)
        if gn is not None:
            members[gn].append(int(y))
    n_groups = len(members)
    shared = {g: ys for g, ys in members.items() if len(set(ys)) > 1}   # groups with BOTH classes
    pos_in_shared = sum(1 for g in shared for y in members[g] if y == 1)
    neg_in_shared = sum(1 for g in shared for y in members[g] if y == 0)
    coverage = sum(len(v) for v in members.values())
    return {
        "n_groups": n_groups,
        "shared_groups": len(shared),
        "pos_in_shared": pos_in_shared,
        "neg_in_shared": neg_in_shared,
        "matched_minority": min(pos_in_shared, neg_in_shared),
        "coverage": coverage,            # non-missing strains on this axis
    }


def confound_report(labels, lineages, regions=None, years=None, *,
                    min_shared_lineages: int = MIN_SHARED_LINEAGES,
                    min_matched_minority: int = MIN_MATCHED_MINORITY,
                    clean_shared_lineages: int = CLEAN_SHARED_LINEAGES) -> dict:
    """labels: 0/1. lineages: primary axis (MLST). regions/years: optional secondary provenance axes.
    Returns a report dict with `verdict` (3-state) and `promotable` (True only for DE_CONFOUNDED)."""
    n = len(labels)
    for name, ax in (("lineages", lineages), ("regions", regions), ("years", years)):
        if ax is not None and len(ax) != n:
            raise ValueError(f"{name} length {len(ax)} != labels length {n}")
    if sum(1 for y in labels if int(y) == 1) == 0 or sum(1 for y in labels if int(y) == 0) == 0:
        return {"verdict": CONFOUNDED, "promotable": False, "reason": "degenerate: a class is empty",
                "primary": None, "secondary": {}}

    primary = axis_contrast(labels, lineages)

    # HARD BLOCK: primary lineage axis lacks within-lineage contrast (no signal beyond lineage to learn).
    if primary["shared_groups"] < min_shared_lineages or primary["matched_minority"] < min_matched_minority:
        return {"verdict": CONFOUNDED, "promotable": False,
                "reason": (f"primary lineage axis lacks within-group contrast: "
                           f"{primary['shared_groups']} shared lineage(s) (need >={min_shared_lineages}), "
                           f"matched minority {primary['matched_minority']} (need >={min_matched_minority})"),
                "primary": primary, "secondary": {}}

    # Secondary provenance axes: a low-cardinality axis (>=2 non-missing groups) with (near-)zero
    # within-group contrast ALIASES the label → downgrade to non-promotable WARN.
    secondary = {}
    aliasing = []
    for name, ax in (("country", regions), ("year", years)):
        if ax is None:
            continue
        c = axis_contrast(labels, ax)
        secondary[name] = c
        if c["n_groups"] >= 2 and c["matched_minority"] < min_matched_minority:
            aliasing.append(name)

    if aliasing:
        return {"verdict": WARN, "promotable": False,
                "reason": (f"primary lineage contrast OK ({primary['shared_groups']} shared, "
                           f"matched {primary['matched_minority']}), but secondary axis/axes "
                           f"{aliasing} alias the label (near-zero within-group contrast) "
                           f"→ non-promotable"),
                "primary": primary, "secondary": secondary}

    if primary["shared_groups"] >= clean_shared_lineages and primary["matched_minority"] >= min_matched_minority:
        return {"verdict": DE_CONFOUNDED, "promotable": True,
                "reason": (f"{primary['shared_groups']} shared lineages, matched minority "
                           f"{primary['matched_minority']}; no secondary axis aliases the label. "
                           f"(Screen passed — necessary, not sufficient.)"),
                "primary": primary, "secondary": secondary}

    return {"verdict": WARN, "promotable": False,
            "reason": (f"borderline primary contrast: {primary['shared_groups']} shared lineages "
                       f"(< {clean_shared_lineages} for a clean pass), matched {primary['matched_minority']}"),
            "primary": primary, "secondary": secondary}


def confound_report_for_cohort(cohort, drug: str, **kw) -> dict:
    labels, lin, reg, yr = [], [], [], []
    keys_lower = None
    for s in cohort.strains:
        keys_lower = {k.lower(): k for k in s.ast_labels}
        if drug.lower() not in keys_lower:
            continue
        labels.append(int(s.ast_labels[keys_lower[drug.lower()]]))
        lin.append(getattr(s, "mlst", None))
        reg.append(getattr(s, "country", None))
        yr.append(getattr(s, "year", None))
    return confound_report(labels, lin, reg, yr, **kw)


def render_report(rep: dict) -> str:
    p = rep.get("primary")
    base = f"[{rep['verdict']}{'·promotable' if rep.get('promotable') else '·NON-PROMOTABLE'}] "
    if p is None:
        return base + rep["reason"]
    sec = "; ".join(f"{k}: matched={v['matched_minority']}/{v['n_groups']}g"
                    for k, v in rep.get("secondary", {}).items())
    return (base + f"primary lineages: {p['shared_groups']} shared, "
            f"matched minority {p['matched_minority']} (R-in-shared {p['pos_in_shared']}, "
            f"S-in-shared {p['neg_in_shared']}); {sec}; {rep['reason']}")
