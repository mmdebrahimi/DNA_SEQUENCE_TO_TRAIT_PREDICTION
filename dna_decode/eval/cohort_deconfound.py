"""Cohort de-confounding gate — a PRECONDITION for any embedding-vs-classical falsifier.

Lesson (2026-06-04, EP-4/AMR): a sampling-independent label TYPE (e.g. AMR MIC) does NOT
guarantee a de-confounded COHORT. If the two classes are drawn from disjoint lineages /
geographies / studies, a classifier "predicting" the label actually predicts batch, not
biology. Pathotype died on this (label IS sampling context); cef gate_b_cohort was
compromised by it (R≈USA, S≈Africa/India, R/S sharing only 1 MLST). cipro N=147 passed
(6 shared R/S lineages). The difference is invisible without this check — so run it BEFORE
building a falsifier, not after.

Primary signal = within-LINEAGE class co-occurrence (always available via MLST): do R and S
strains actually share lineages, and is a meaningful minority fraction of each class inside
those shared lineages? Geography (country) is a secondary informational flag that can
downgrade an otherwise-clean verdict to WARN.

Pure functions over plain lists → exhaustively unit-testable; no I/O. A cohort adapter
(`confound_report_for_cohort`) pulls labels/MLST/country from a Cohort for a given drug.
"""
from __future__ import annotations

from collections import Counter

DE_CONFOUNDED = "DE_CONFOUNDED"   # R and S genuinely co-occur within lineages → verdict measures biology
WARN = "WARN"                     # borderline co-occurrence (or clean lineages but separable geography)
CONFOUNDED = "CONFOUNDED"         # classes near lineage/geography-separable → verdict would measure batch

# thresholds (grounded on the two known data points: cipro N=147 PASS vs cef gate_b FAIL)
MIN_SHARED_LINEAGES = 3           # cef had 1 → CONFOUNDED; cipro had 6 → ok
MIN_MINORITY_SHARED_FRAC = 0.10   # cef min(3%,4%)=3% → CONFOUNDED; cipro min(13%,49%)=13% → ok
CLEAN_SHARED_LINEAGES = 5         # >= this (with the frac) → DE_CONFOUNDED; 3-4 → WARN
GEO_DOMINANCE = 0.60              # a class is "region-dominant" if >= this share sits in one region


def confound_report(labels, lineages, regions=None, *,
                    min_shared_lineages: int = MIN_SHARED_LINEAGES,
                    min_minority_shared_frac: float = MIN_MINORITY_SHARED_FRAC,
                    clean_shared_lineages: int = CLEAN_SHARED_LINEAGES,
                    geo_dominance: float = GEO_DOMINANCE) -> dict:
    """labels: 0/1 ints. lineages: hashable lineage id per strain. regions: optional region/country
    per strain (blank/None ignored). Returns a report dict incl. `verdict`."""
    n = len(labels)
    if len(lineages) != n or (regions is not None and len(regions) != n):
        raise ValueError("labels, lineages, regions must be equal length")
    pos = [i for i, y in enumerate(labels) if int(y) == 1]
    neg = [i for i, y in enumerate(labels) if int(y) == 0]
    if not pos or not neg:
        return {"verdict": CONFOUNDED, "reason": "degenerate: a class is empty",
                "n_pos": len(pos), "n_neg": len(neg), "shared_lineages": 0,
                "pos_in_shared_frac": 0.0, "neg_in_shared_frac": 0.0, "geo": None}

    pos_lin = {lineages[i] for i in pos}
    neg_lin = {lineages[i] for i in neg}
    shared = pos_lin & neg_lin
    pos_in_shared = sum(1 for i in pos if lineages[i] in shared)
    neg_in_shared = sum(1 for i in neg if lineages[i] in shared)
    pos_frac = pos_in_shared / len(pos)
    neg_frac = neg_in_shared / len(neg)
    minority_frac = min(pos_frac, neg_frac)

    # geography (secondary, informational + WARN-downgrade)
    geo = None
    geo_separable = False
    if regions is not None:
        def top(idxs):
            c = Counter(str(regions[i]).strip() for i in idxs if str(regions[i]).strip())
            return c.most_common(1)[0] if c else (None, 0)
        pr, pc = top(pos); nr, nc = top(neg)
        pos_present = sum(1 for i in pos if str(regions[i]).strip())
        neg_present = sum(1 for i in neg if str(regions[i]).strip())
        pos_dom = (pc / pos_present) if pos_present else 0.0
        neg_dom = (nc / neg_present) if neg_present else 0.0
        geo_separable = bool(pr and nr and pr != nr and pos_dom >= geo_dominance and neg_dom >= geo_dominance)
        geo = {"pos_top_region": pr, "pos_top_share": round(pos_dom, 2),
               "neg_top_region": nr, "neg_top_share": round(neg_dom, 2),
               "separable": geo_separable}

    # verdict — lineage is the hard gate; geography can only downgrade a clean verdict to WARN.
    if len(shared) < min_shared_lineages or minority_frac < min_minority_shared_frac:
        verdict = CONFOUNDED
        reason = (f"only {len(shared)} shared lineage(s) (need >={min_shared_lineages}) / "
                  f"minority-class shared fraction {minority_frac:.2f} (need >={min_minority_shared_frac:.2f})")
    elif len(shared) >= clean_shared_lineages and minority_frac >= min_minority_shared_frac:
        verdict = DE_CONFOUNDED
        reason = f"{len(shared)} shared lineages; minority-class shared fraction {minority_frac:.2f}"
        if geo_separable:
            verdict = WARN
            reason += f" BUT geography near-separable ({geo['pos_top_region']} vs {geo['neg_top_region']})"
    else:
        verdict = WARN
        reason = (f"borderline: {len(shared)} shared lineages, minority shared fraction {minority_frac:.2f}")

    return {
        "verdict": verdict, "reason": reason,
        "n_pos": len(pos), "n_neg": len(neg),
        "pos_lineages": len(pos_lin), "neg_lineages": len(neg_lin),
        "shared_lineages": len(shared),
        "pos_in_shared": pos_in_shared, "neg_in_shared": neg_in_shared,
        "pos_in_shared_frac": round(pos_frac, 3), "neg_in_shared_frac": round(neg_frac, 3),
        "minority_shared_frac": round(minority_frac, 3),
        "geo": geo,
    }


def confound_report_for_cohort(cohort, drug: str, **kw) -> dict:
    """Adapter: pull (label, MLST, country) per strain that carries `drug` in ast_labels."""
    labels, lineages, regions = [], [], []
    for s in cohort.strains:
        if drug.lower() not in {k.lower() for k in s.ast_labels}:
            continue
        lab = s.ast_labels.get(drug) if drug in s.ast_labels else s.ast_labels.get(drug.lower())
        labels.append(int(lab))
        lineages.append(str(getattr(s, "mlst", None)))
        regions.append(getattr(s, "country", None))
    return confound_report(labels, lineages, regions, **kw)


def render_report(rep: dict) -> str:
    g = rep.get("geo")
    geo = (f"; geography {g['pos_top_region']}({g['pos_top_share']}) vs {g['neg_top_region']}"
           f"({g['neg_top_share']}), separable={g['separable']}" if g else "")
    return (f"[{rep['verdict']}] {rep['n_pos']}R/{rep['n_neg']}S; "
            f"shared lineages {rep['shared_lineages']} "
            f"(R-in-shared {rep['pos_in_shared_frac']}, S-in-shared {rep['neg_in_shared_frac']}); "
            f"{rep['reason']}{geo}")
