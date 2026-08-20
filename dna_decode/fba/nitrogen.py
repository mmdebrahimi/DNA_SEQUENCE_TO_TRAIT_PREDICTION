"""Nitrogen-source conditional essentiality — the first NON-carbon substrate axis.

Mirrors the carbon path (`fitness_browser.carbon_conditions` / `apply_carbon_condition`) but for
`expGroup='nitrogen source'`. Kept as a SEPARATE module rather than generalising the carbon functions in
place: the carbon path is load-bearing with many pinned tests, and this axis has a genuinely different
medium contract (a fixed carbon source is held while the nitrogen source varies).

Probed 2026-08-17: feba.db Keio has 32 nitrogen experiments over 16 distinct sources; 13 map to iML1515
exchanges. `Gly-DL-Asp` / `Gly-Glu` are dipeptides with no iML1515 exchange, and `casamino acids` is a
mixture — all three are unmappable by design, the same exclusion the carbon panel makes.
"""
from __future__ import annotations

import re
import sqlite3

from .conditional_essentiality import GeneRecord
from .fitness_browser import ESSENTIAL_FITNESS, ORG_ID

NITROGEN_EXP_GROUP = "nitrogen source"

#: Curated {Fitness-Browser label -> iML1515 exchange}. Explicit rather than fuzzy-matched: a nitrogen
#: label like "D-Serine" would fuzzy-match the carbon index too, and silently scoring the wrong compound
#: is the failure mode this whole module exists to avoid.
NITROGEN_EXCHANGES: dict[str, str] = {
    "Adenosine": "EX_adn_e",
    "Ammonium chloride": "EX_nh4_e",
    "Cytidine": "EX_cytd_e",
    "D-Alanine": "EX_ala__D_e",
    "D-Serine": "EX_ser__D_e",
    "Glycine": "EX_gly_e",
    "L-Alanine": "EX_ala__L_e",
    "L-Arginine": "EX_arg__L_e",
    "L-Asparagine": "EX_asn__L_e",
    "L-Aspartic Acid": "EX_asp__L_e",
    "L-Glutamine": "EX_gln__L_e",
    "L-Serine": "EX_ser__L_e",
    "Putrescine Dihydrochloride": "EX_ptrc_e",
}

#: Deliberately excluded, with the reason. Kept in code so the exclusion is auditable, not silent.
NITROGEN_UNMAPPABLE: dict[str, str] = {
    "Gly-DL-Asp": "dipeptide; no iML1515 exchange",
    "Gly-Glu": "dipeptide; no iML1515 exchange",
    "casamino acids": "mixture, not a single compound (same exclusion the carbon panel makes)",
}

DEFAULT_CARBON = "EX_glc__D_e"

#: Salt / hydrate forms stripped before an exact name comparison. Order matters: the longer forms must
#: come first so "dihydrochloride" is not left as "di" by an earlier "hydrochloride" replacement.
_SALT_FORMS = ("dihydrochloride", "hydrochloride", "monohydrate", "hexahydrate", "dihydrate",
               "hydrate", "disodium salt", "sodium salt", "potassium salt", "chloride", "salt")


def normalize_compound(name: str) -> str:
    """Lowercase, strip salt/hydrate decoration, drop non-alphanumerics.

    Used ONLY for EXACT equality against a model metabolite's own name. It is deliberately not a fuzzy
    matcher: a hand-guessed mapping in the stress probe paired "sodium fluoride" with EX_fe2_e (ferrous
    iron), which would have manufactured a data point out of an unrelated metabolite.
    """
    s = str(name).lower()
    for j in _SALT_FORMS:
        s = s.replace(j, " ")
    return re.sub(r"[^a-z0-9]", "", s)


def exchange_name_index(model) -> dict[str, str]:
    """{normalized metabolite name -> EX_ id}. First writer wins, so the map is deterministic."""
    idx: dict[str, str] = {}
    for r in model.reactions:
        if r.id.startswith("EX_") and r.id.endswith("_e"):
            mets = list(r.metabolites)
            if mets:
                idx.setdefault(normalize_compound(mets[0].name), r.id)
    return idx


def nitrogen_conditions_for_org(conn: sqlite3.Connection, model, org_id: str) -> dict[str, str]:
    """{nitrogen source -> exchange} for ANY organism: curated map first, then EXACT name match.

    The curated `NITROGEN_EXCHANGES` entries win where they apply (they were hand-verified against
    iML1515). Everything else is admitted only on exact post-normalization equality with the model's own
    metabolite name -- never a substring or fuzzy match.
    """
    have = {r.id for r in model.exchanges}
    idx = exchange_name_index(model)
    out: dict[str, str] = {}
    for (cond,) in conn.execute(
            "SELECT DISTINCT condition_1 FROM Experiment WHERE orgId=? AND expGroup=?",
            (org_id, NITROGEN_EXP_GROUP)):
        curated = NITROGEN_EXCHANGES.get(cond)
        if curated and curated in have:
            out[cond] = curated
            continue
        ex = idx.get(normalize_compound(cond))
        if ex and ex in have:
            out[cond] = ex
    return out


def nitrogen_conditions(conn: sqlite3.Connection, model) -> dict[str, str]:
    """{nitrogen-source label -> iML1515 exchange} for sources present in BOTH the assay and the model."""
    have = {r.id for r in model.exchanges}
    present = {c for (c,) in conn.execute(
        "SELECT DISTINCT condition_1 FROM Experiment WHERE orgId=? AND expGroup=?",
        (ORG_ID, NITROGEN_EXP_GROUP))}
    return {k: v for k, v in NITROGEN_EXCHANGES.items() if k in present and v in have}


def apply_nitrogen_condition(model, exchange: str, all_nitrogen: tuple[str, ...] = (),
                             uptake: float = 10.0, carbon: str = DEFAULT_CARBON,
                             carbon_uptake: float = 10.0) -> None:
    """Set `exchange` as the SOLE nitrogen source while holding `carbon` fixed. IN PLACE.

    Every other candidate nitrogen exchange -- crucially `EX_nh4_e`, which is in the default medium -- is
    closed FIRST. Without that, a residual ammonium uptake makes every condition silently score as
    ammonium; this is the exact failure mode `_ALL_CARBON` guards against on the carbon axis.

    NOTE (biology, not a bug): several of these sources (alanine, serine, aspartate, glutamine ...) also
    supply CARBON. That matches the real assay -- glucose minimal medium plus the test compound as sole N
    source -- so these are not nitrogen-only perturbations and must not be described as such.
    """
    have = {r.id for r in model.exchanges}
    if exchange not in have:
        raise KeyError(f"model {model.id} has no exchange {exchange!r}")
    if carbon not in have:
        raise KeyError(f"model {model.id} has no carbon exchange {carbon!r}")

    medium = dict(model.medium)
    for ex in set(all_nitrogen) | set(NITROGEN_EXCHANGES.values()) | {"EX_nh4_e"}:
        medium.pop(ex, None)
    medium[carbon] = carbon_uptake
    medium[exchange] = uptake
    model.medium = medium


#: Fields that are DERIVED FROM SOLVES and therefore may not be quoted if the panel is not reproducible.
#: `best_constant_null` is deliberately NOT here -- it comes from the labels, not from any solve.
REDACTED_ON_NONDETERMINISM = ("per_cell_agreement", "per_condition", "predictions")


def determinism_verdict(pass_a: dict[str, dict[str, float]], pass_b: dict[str, dict[str, float]],
                        keys, genes, frac: float = 0.01, min_safety_factor: float = 1000.0,
                        metric_a: float | None = None, metric_b: float | None = None) -> dict:
    """Are two passes reproducible AT THE LEVEL OF THE CLAIMS THEY SUPPORT?

    **Why this is not a bit-equality check.** The first version of this gate asked whether two passes
    agreed to within 1e-12 on every growth ratio, and reported 147 of 2,015 cells differing. Probing it
    (`scripts/fba_nitrogen_determinism_probe.py`) showed the largest disagreement anywhere was 3.2e-11 --
    ordinary float64 behaviour for an LP objective divided by another LP objective -- that ZERO cells
    crossed the `frac` call line, and that the headline metric was identical to six decimals across three
    passes including one on a freshly loaded model. The differing-cell COUNT was itself unstable across
    runs (147, then 58, then 66), which is the signature of a threshold sitting inside float noise.

    So bit-equality was measuring the wrong thing. What the pre-registration actually needs to rule out is
    the failure that forced today's retraction: **a conclusion that changes depending on which run you
    quote.** That is a claim-level property, and it is checked here as three conditions:

    1. **No call flips.** No cell may land on opposite sides of `frac` in the two passes. This is the only
       drift that can change an essentiality call, and therefore the only one that can change a claim.
    2. **The headline metric is identical** at reported precision (when supplied).
    3. **A DERIVED margin, not an asserted tolerance.** The gate compares the largest observed drift
       against `min_margin_to_threshold` -- how close the nearest cell in the whole panel actually gets to
       the decision line. `safety_factor = min_margin / max_abs_delta` is the number of orders of magnitude
       of headroom between noise and the line. This adapts to the data instead of being picked: if a
       future panel's ratios crowd the threshold, the margin shrinks and the gate FAILS on its own. A
       fixed tolerance could not do that.

    `min_safety_factor=1000` demands three orders of magnitude of headroom. It is not tuned to let this
    panel through -- the measured factor here is ~1e8, five orders of magnitude clear of the bar -- and
    `test_determinism_gate_fails_when_a_cell_sits_near_the_line` pins that it still fails when it should.

    This gate is STRICTER than the one it replaces in two ways the original never checked: the metric
    must match, and a failure now REDACTS the numbers (`REDACTED_ON_NONDETERMINISM`) instead of printing
    them beside a `deterministic: false` flag.
    """
    deltas: list[float] = []
    flips: list[tuple[str, str, float, float]] = []
    margins: list[float] = []
    for c in keys:
        for g in genes:
            a = pass_a[c].get(g, 0.0)
            b = pass_b[c].get(g, 0.0)
            margins.append(abs(a - frac))
            d = abs(a - b)
            if d > 0.0:
                deltas.append(d)
            if (a <= frac) != (b <= frac):
                flips.append((c, g, a, b))

    max_delta = max(deltas) if deltas else 0.0
    min_margin = min(margins) if margins else 0.0
    safety = float("inf") if max_delta == 0.0 else min_margin / max_delta
    metric_ok = (metric_a is None and metric_b is None) or (metric_a == metric_b)
    ok = (not flips) and metric_ok and safety >= min_safety_factor

    return {
        "deterministic_at_claim_level": ok,
        "n_cells": len(keys) * len(genes),
        "n_call_flips": len(flips),
        "call_flip_examples": [(c, g, a, b) for c, g, a, b in flips[:5]],
        "max_abs_delta": max_delta,
        "min_margin_to_threshold": min_margin,
        "safety_factor": safety,
        "min_safety_factor": min_safety_factor,
        "headline_metric_matches": metric_ok,
        "headline_metric": [metric_a, metric_b],
        "basis": ("claim-level reproducibility: no cell crosses the call line, the headline metric is "
                  "identical, and the nearest cell to the threshold is >= min_safety_factor times "
                  "further from it than the largest observed numerical drift"),
    }


def redact_unverified(payload: dict, deterministic: bool) -> dict:
    """Blank every solve-derived claim when the panel did not clear the determinism gate.

    **This is the control that failed.** The first nitrogen run wrote `"deterministic": false` into its
    artifact and then reported `per_cell_agreement`, `per_condition` and all three P1/P2/P3 verdicts
    beside it -- numbers a reader could quote in good faith, despite the pre-registration saying "no
    verdict is reported" on a determinism failure. A flag that sits next to the numbers it invalidates is
    not a control. Removing the numbers is.
    """
    if deterministic:
        return payload
    out = dict(payload)
    for f in REDACTED_ON_NONDETERMINISM:
        out[f] = None
    out["verdict"] = "NON_DETERMINISTIC_NO_VERDICT"
    out["redaction_note"] = (
        "Withheld " + ", ".join(REDACTED_ON_NONDETERMINISM) + " because the panel failed its "
        "pre-registered determinism gate. Per wiki/fba_nitrogen_prereg_2026-08-17.md these numbers may "
        "not be reported, so they are removed rather than printed beside a false flag.")
    return out


def load_nitrogen_records(conn: sqlite3.Connection, conditions: dict[str, str],
                          gene_filter: set[str] | None = None,
                          threshold: float = ESSENTIAL_FITNESS,
                          org_id: str = ORG_ID) -> list[GeneRecord]:
    """Experimental conditional essentiality per (gene, nitrogen source), as GeneRecords.

    Mirrors `fitness_browser.load_records` but scoped to `expGroup='nitrogen source'`. Replicates for the
    same source are averaged. A gene is kept only if it has a value in EVERY condition -- a partial row
    would make a missing measurement look like a dispensability call.
    """
    keys = tuple(sorted(conditions))
    exp_to_cond = {}
    for name, cond in conn.execute(
            "SELECT expName, condition_1 FROM Experiment WHERE orgId=? AND expGroup=?",
            (org_id, NITROGEN_EXP_GROUP)):
        if cond in conditions:
            exp_to_cond[name] = cond

    agg: dict[tuple[str, str], list[float]] = {}
    q = ("SELECT g.sysName, f.expName, f.fit FROM GeneFitness f "
         "JOIN Gene g ON g.orgId=f.orgId AND g.locusId=f.locusId WHERE f.orgId=?")
    for sysname, expname, fit in conn.execute(q, (org_id,)):
        cond = exp_to_cond.get(expname)
        if cond is None or not sysname:
            continue
        if gene_filter is not None and sysname not in gene_filter:
            continue
        agg.setdefault((sysname, cond), []).append(float(fit))

    by_gene: dict[str, dict[str, float]] = {}
    for (gene, cond), vals in agg.items():
        by_gene.setdefault(gene, {})[cond] = sum(vals) / len(vals)

    out: list[GeneRecord] = []
    for gene, per_cond in by_gene.items():
        if len(per_cond) != len(keys):
            continue
        exp = {c: (per_cond[c] < threshold) for c in keys}
        out.append(GeneRecord(gene_id=gene, gene=gene, experimental=exp,
                              paper_fba={c: False for c in keys},
                              conditionally_essential=any(exp.values()) and not all(exp.values())))
    return out
