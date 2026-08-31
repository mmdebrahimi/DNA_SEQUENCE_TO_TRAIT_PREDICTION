"""Inline trust-surface — a decoder call's own "how well is this validated" badge.

The project's load-bearing honesty discipline lives in the standing report cards
(`wiki/*_report_card.json`): each (organism, drug) cell carries a DIFFERENT, explicitly-labelled
independence tier (provenance-disjoint measured AST / free wet-lab fold-change / in-distribution
knowledge baseline / abstains-by-design / no-free-source / not-censused). Until now that tier was
visible only by reading the wiki — NOT in the tool's own output. This module surfaces it INLINE:
given a call's `(drug, organism)`, return the cell's honest tier + headline metric + the source card,
so every prediction carries its trust badge instead of leaving it buried.

PURE-PYTHON, NO NETWORK: reads the committed standing report-card JSONs at lookup time (cached). It
NEVER fabricates a number — a cell with no validation evidence returns NOT_CENSUSED / NO_FREE_SOURCE /
UNKNOWN, never a borrowed metric. It also NEVER averages across tiers (the project's no-aggregate-
headline discipline) — it reports the single best-evidence cell that matches.

Resolution order (strongest evidence first): HIV free wet-lab → TB independent measured → bacterial
EBI AMR Portal independent measured → bacterial NCBI-PD provenance-disjoint → SARS-CoV-2 in-distribution
→ shipped-surface structural fallback (no-free-source / label-confounded / not-censused) → UNKNOWN.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

# Report cards load from the PACKAGED copy first (a built wheel force-includes them at
# dna_decode/report_cards/ -- see pyproject), falling back to the repo-root wiki/ in an editable checkout.
# This is the packaging gate (2026-06-24): without the packaged copy, a wheel install silently degrades
# every trust badge because site-packages/wiki/ does not exist.
_PKG_CARDS = Path(__file__).resolve().parent.parent / "report_cards"   # installed wheel
_WIKI = Path(__file__).resolve().parent.parent.parent / "wiki"          # editable / source tree


def _card_path(name: str) -> Path:
    pkg = _PKG_CARDS / name
    return pkg if pkg.exists() else (_WIKI / name)

# --- tiers (ordered strongest -> weakest) ---
INDEPENDENT_WETLAB = "INDEPENDENT_WETLAB"            # HIV: free isolate-level wet-lab fold-change (non-circular)
INDEPENDENT_MEASURED = "INDEPENDENT_MEASURED"        # bacteria/TB: free measured AST, provenance/BioSample-disjoint
PROVENANCE_DISJOINT = "PROVENANCE_DISJOINT"          # NCBI-PD: submitter/lab/country-disjoint (not methodology-indep)
IN_DISTRIBUTION = "IN_DISTRIBUTION"                  # knowledge baseline (catalog + labels same source) — NOT independent
UNDERPOWERED = "UNDERPOWERED"
LABEL_CONFOUNDED = "LABEL_CONFOUNDED"                # phenotype label is an unreliable surrogate
ABSTAINS_BY_DESIGN = "ABSTAINS_BY_DESIGN"            # decoder refuses (expression-driven R it cannot decode)
NO_FREE_PHENOTYPE_SOURCE = "NO_FREE_PHENOTYPE_SOURCE"  # no free isolate-level phenotype source exists
NOT_CENSUSED = "NOT_CENSUSED"                        # shipped decoder, not yet scored on a validation cohort
UNKNOWN = "UNKNOWN"                                  # drug/organism not a recognised decoder cell

_CAVEAT = {
    INDEPENDENT_WETLAB: "free, independent, isolate-level wet-lab fold-change (Stanford HIVDB PhenoSense; non-circular)",
    INDEPENDENT_MEASURED: "independent measured-AST cohort (provenance/BioSample-disjoint); non-circular",
    PROVENANCE_DISJOINT: "provenance-disjoint (different submitter/lab/country); NOT methodology-independent",
    IN_DISTRIBUTION: "in-distribution knowledge baseline (catalog + labels share a source); NOT independent validation",
    UNDERPOWERED: "validation cohort is underpowered (too few isolates / one class); read the metric with caution",
    LABEL_CONFOUNDED: "the phenotype LABEL is an unreliable surrogate; the genotype call may be the more trustworthy output",
    ABSTAINS_BY_DESIGN: "the decoder ABSTAINS by design — expression/regulation-driven R it cannot decode from gene presence",
    NO_FREE_PHENOTYPE_SOURCE: "no free isolate-level measured-phenotype source exists; the catalog is curated but NOT validated here",
    NOT_CENSUSED: "shipped decoder not yet scored on a validation cohort — no independence claim",
    UNKNOWN: "not a recognised decoder cell on a standing report card",
}

# drugs whose validation lives on the SARS-CoV-2 in-distribution surface (separate card namespace)
_SARSCOV2_DRUGS = {"nirmatrelvir", "ensitrelvir", "lufotrelvir"}
_HCMV_DRUGS = {"ganciclovir", "valganciclovir", "cidofovir", "foscarnet", "letermovir"}


@lru_cache(maxsize=None)
def _load(name: str):
    p = _card_path(name)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


def _genus(organism: str | None) -> str:
    """Coarse genus token: first word, lowercased, underscores treated as spaces (the decoder is genus-routed,
    and the cards spell the same organism three ways — 'Escherichia' / 'Escherichia_coli_Shigella' /
    'Escherichia coli' all collapse to 'escherichia')."""
    if not organism:
        return ""
    return organism.strip().replace("_", " ").split()[0].lower() if organism.strip() else ""


def _norm_org(s: str | None) -> str:
    """Full normalized organism string (lowercased, underscores->spaces) for EXACT cell matching."""
    return (s or "").strip().lower().replace("_", " ")


def _pick_bacterial_cell(cells: list[dict], d: str, organism: str | None):
    """Pick a bacterial card cell for (drug, organism). EXACT normalized-organism match wins; genus fallback
    fires ONLY when the genus resolves to a SINGLE distinct organism for this drug. A bare/under-specified
    genus that spans >=2 distinct species cells (e.g. 'Shigella' -> flexneri vs sonnei, different metrics) is
    AMBIGUOUS -> never silently borrows one species' number. Returns (cell, status) with status in
    {'exact', 'genus', 'ambiguous', 'none'}."""
    no, g = _norm_org(organism), _genus(organism)
    drug_cells = [c for c in cells if str(c.get("drug", "")).strip().lower() == d]
    for c in drug_cells:                                   # exact normalized-organism match first
        if no and _norm_org(c.get("organism")) == no:
            return c, "exact"
    gmatch = [c for c in drug_cells if g and _genus(c.get("organism")) == g]
    distinct = {_norm_org(c.get("organism")) for c in gmatch}
    if len(distinct) == 1:
        return gmatch[0], "genus"
    if len(distinct) >= 2:
        return None, "ambiguous"
    return None, "none"


# namespace genus tokens (the guard): an HIV/TB/SARS card may lend evidence ONLY when the requested organism
# is ABSENT (drug-only lookup) or normalizes into that namespace -- never to a contradictory organism.
_HIV_GENUS = {"hiv", "hiv-1", "hiv1"}
_TB_GENUS = {"mycobacterium", "m.tuberculosis", "mtb", "tuberculosis"}
_SARS_GENUS = {"sars-cov-2", "sarscov2", "sars", "betacoronavirus"}
_HCMV_GENUS = {"hcmv", "cmv", "cytomegalovirus", "human betaherpesvirus 5", "betaherpesvirus", "human herpesvirus 5"}

_MISMATCH_CAVEAT = ("drug recognised in another organism's namespace but the requested organism does NOT "
                    "match it -- refusing to lend that cell's evidence (no fabricated tier; see evidence_cell)")
_AMBIGUOUS_CAVEAT = ("the requested genus spans >=2 distinct species cells with DIFFERENT metrics -- refusing "
                     "to borrow one species' number; pass the exact species (e.g. 'Shigella sonnei')")
_REASON_CAVEAT = {"namespace_mismatch": _MISMATCH_CAVEAT, "ambiguous_genus": _AMBIGUOUS_CAVEAT}


def _rec(tier: str, source_card: str, headline: str = "", metric: float | None = None,
         n: int | None = None, cell: str = "", *, reason: str | None = None,
         requested_cell: str = "", evidence_cell: str | None = None) -> dict:
    ev = evidence_cell if evidence_cell is not None else (cell or None)
    return {
        "tier": tier,
        "independent": tier in (INDEPENDENT_WETLAB, INDEPENDENT_MEASURED),
        "headline": headline,
        "metric": metric,
        "n": n,
        "cell": cell,                          # backward-compat alias of evidence_cell
        "requested_cell": requested_cell or cell or None,
        "evidence_cell": ev,
        "reason": reason,
        "source_card": source_card,
        "caveat": _REASON_CAVEAT.get(reason, _CAVEAT[tier]),
    }


def lookup_trust(drug: str, organism: str | None = None) -> dict:
    """Best-evidence honest validation badge for a (drug, organism) decoder cell. Always returns a dict
    (tier UNKNOWN if nothing matches); never fabricates a metric, never averages across tiers.

    NAMESPACE GUARD: the HIV/TB/SARS cards are matched by drug, but only LEND their evidence when the
    requested organism is absent (drug-only) OR normalizes into that namespace. A contradictory organism
    (e.g. rifampicin + Escherichia) is REFUSED -> tier UNKNOWN + reason='namespace_mismatch', with the
    rejected candidate exposed in evidence_cell so the borrowing is auditable (never silently applied)."""
    d = (drug or "").strip().lower()
    g = _genus(organism)
    req = f"{organism or '?'}|{d}"
    rejected: tuple[str, str] | None = None   # (evidence_cell, source_card) of a drug-match the organism rejected

    def _compatible(ns: set[str]) -> bool:
        return g == "" or g in ns

    # 1. HIV — free wet-lab fold-change (no organism axis in the card; guarded by namespace)
    hiv = _load("hiv_decoder_report_card.json")
    if hiv:
        for c in hiv.get("cells", []):
            if str(c.get("drug", "")).strip().lower() == d:
                if _compatible(_HIV_GENUS):
                    auc = c.get("auc_call_separates_fold")
                    return _rec(INDEPENDENT_WETLAB, "wiki/hiv_decoder_report_card.md",
                                headline=f"AUC {auc} (N={c.get('n')}, {c.get('drug_class')})" if auc else "scored",
                                metric=auc, n=c.get("n"), cell=f"HIV-1|{d}", requested_cell=req)
                rejected = rejected or (f"HIV-1|{d}", "wiki/hiv_decoder_report_card.md")
                break

    # 2. TB — independent measured AST (guarded by namespace)
    tb = _load("tb_report_card.json")
    if tb:
        for c in tb.get("independent", []):
            if str(c.get("drug", "")).strip().lower() == d:
                if _compatible(_TB_GENUS):
                    acc = c.get("raw_acc")
                    n = (c.get("n_R") or 0) + (c.get("n_S") or 0)
                    return _rec(INDEPENDENT_MEASURED, "wiki/tb_report_card.md",
                                headline=f"acc {round(acc, 3)} (N={n})" if acc is not None else "scored",
                                metric=acc, n=n, cell=f"M.tuberculosis|{d}", requested_cell=req)
                rejected = rejected or (f"M.tuberculosis|{d}", "wiki/tb_report_card.md")
                break

    # 3. bacteria — EBI AMR Portal INDEPENDENT measured AST (exact organism > unambiguous genus)
    portal = _load("amr_portal_independent_report_card.json")
    if portal:
        c, status = _pick_bacterial_cell(portal.get("cells", []), d, organism)
        if status == "ambiguous":
            return _rec(UNKNOWN, "wiki/amr_portal_independent_report_card.md",
                        reason="ambiguous_genus", requested_cell=req, evidence_cell=None)
        if c:
            acc = c.get("accuracy")
            n = (c.get("n_R") or 0) + (c.get("n_S") or 0)
            tier = UNDERPOWERED if str(c.get("tier", "")).upper().startswith("UNDERPOWER") else INDEPENDENT_MEASURED
            return _rec(tier, "wiki/amr_portal_independent_report_card.md",
                        headline=f"acc {round(acc, 3)} (N={n})" if acc is not None else "scored",
                        metric=acc, n=n, cell=f"{c.get('organism')}|{d}", requested_cell=req)

    # 4. bacteria — NCBI-PD provenance-disjoint card (exact > unambiguous genus); also the structural non-cells
    deck = _load("decoder_validation_report_card.json")
    if deck:
        c, status = _pick_bacterial_cell(deck.get("cells", []), d, organism)
        if status == "ambiguous":
            return _rec(UNKNOWN, "wiki/decoder_validation_report_card.md",
                        reason="ambiguous_genus", requested_cell=req, evidence_cell=None)
        if c:
            st = str(c.get("state", "")).upper()
            if st == "SCORED":
                acc = c.get("acc")
                return _rec(PROVENANCE_DISJOINT, "wiki/decoder_validation_report_card.md",
                            headline=f"acc {acc} (N={c.get('n')})" if acc is not None else "scored",
                            metric=acc, n=c.get("n"), cell=f"{c.get('organism')}|{d}", requested_cell=req)
            _STATE_TIER = {"ABSTAINS_BY_DESIGN": ABSTAINS_BY_DESIGN, "LABEL_CONFOUNDED": LABEL_CONFOUNDED,
                           "UNDERPOWERED": UNDERPOWERED, "NO_FREE_PHENOTYPE_SOURCE": NO_FREE_PHENOTYPE_SOURCE,
                           "NOT_CENSUSED": NOT_CENSUSED}
            tier = _STATE_TIER.get(st, NOT_CENSUSED)
            return _rec(tier, "wiki/decoder_validation_report_card.md",
                        cell=f"{c.get('organism')}|{d}", requested_cell=req)

    # 5. SARS-CoV-2 — in-distribution knowledge baseline (separate namespace, underpowered; guarded)
    if d in _SARSCOV2_DRUGS:
        if _compatible(_SARS_GENUS):
            return _rec(IN_DISTRIBUTION, "wiki/sarscov2_mpro_validation_result_2026-06-23.md",
                        headline="in-distribution (CoV-RDB), underpowered", cell=f"SARS-CoV-2|{d}",
                        requested_cell=req)
        rejected = rejected or (f"SARS-CoV-2|{d}", "wiki/sarscov2_mpro_validation_result_2026-06-23.md")

    # 5b. HCMV — in-distribution knowledge baseline (herpesvirus; catalog curated from Chou recombinant
    # fold-change; NO free held-out per-isolate phenotype exists -> independent is a CLOSED negative). Guarded.
    if d in _HCMV_DRUGS:
        if _compatible(_HCMV_GENUS):
            return _rec(IN_DISTRIBUTION, "wiki/hcmv_decoder_report_card.md",
                        headline="in-distribution (Chou recombinant fold-change); independent = closed (no free held-out)",
                        cell=f"HCMV|{d}", requested_cell=req)
        rejected = rejected or (f"HCMV|{d}", "wiki/hcmv_decoder_report_card.md")

    # 6. shipped-surface structural fallback (no card cell yet)
    try:
        from dna_decode.data.shipped_decoder_surface import SHIPPED_DECODER_SURFACE
        for (org, drg, _eng, _scope, status, _grp) in SHIPPED_DECODER_SURFACE:
            if drg.strip().lower() == d and (not g or _genus(org) == g):
                tier = {"no_free_source": NO_FREE_PHENOTYPE_SOURCE, "label_confounded": LABEL_CONFOUNDED,
                        "ncbi_pd": NOT_CENSUSED}.get(status, NOT_CENSUSED)
                return _rec(tier, "dna_decode/data/shipped_decoder_surface.py",
                            cell=f"{org}|{d}", requested_cell=req)
    except Exception:
        pass

    # nothing matched. If a namespace card recognised the drug but the organism rejected it, say so.
    if rejected:
        return _rec(UNKNOWN, rejected[1], reason="namespace_mismatch",
                    requested_cell=req, evidence_cell=rejected[0])
    return _rec(UNKNOWN, "", requested_cell=req, evidence_cell=None)


def prospective_regressions() -> list[dict]:
    """Every cell whose PROSPECTIVE-lock result contradicts its standing validation. Read-only.

    Sourced from the report card's top-level `prospective_regression` flag (a POWERED post-lock cohort
    whose sensitivity is below the card's floor). Empty when the card is absent or nothing regressed.
    """
    deck = _load("decoder_validation_report_card.json")
    out = []
    for c in (deck or {}).get("cells", []):
        if not c.get("prospective_regression"):
            continue
        p = c.get("prospective") or {}
        out.append({"organism": c.get("organism"), "drug": c.get("drug"),
                    "sens": p.get("sens"), "n": p.get("n_scored"), "lock_date": p.get("lock_date"),
                    "caveat": c.get("deployment_caveat", "")})
    return out


def prospective_regression_for(drug: str, organism: str | None = None) -> dict | None:
    """The prospective regression on THIS cell, if any. Reuses the card's own organism matcher."""
    deck = _load("decoder_validation_report_card.json")
    if not deck:
        return None
    c, status = _pick_bacterial_cell(deck.get("cells", []), str(drug).strip().lower(), organism)
    if not c or status == "ambiguous" or not c.get("prospective_regression"):
        return None
    p = c.get("prospective") or {}
    return {"organism": c.get("organism"), "drug": c.get("drug"),
            "sens": p.get("sens"), "n": p.get("n_scored"), "lock_date": p.get("lock_date"),
            "caveat": c.get("deployment_caveat", "")}


@lru_cache(maxsize=1)
def _doubt_card() -> dict | None:
    """Newest committed per-cell doubt-layer measurement, or None. Read-only, never fabricated."""
    for base in (_PKG_CARDS, _WIKI):
        hits = sorted(base.glob("doubt_layer_per_cell_*.json")) if base.exists() else []
        if hits:
            try:
                return json.loads(hits[-1].read_text(encoding="utf-8"))
            except (OSError, ValueError):
                return None
    return None


def doubt_layer_for(drug: str) -> dict | None:
    """This cell's standing L2 doubt-layer status, or None when it has never been measured.

    AUGMENT-ONLY, and that is the whole discipline: this reports whether the cell has a known
    CATALOG-COMPLETENESS gap, which is a different question from how well the cell is VALIDATED. It
    must never move a tier or a metric — the same rule the lineage, prospective and
    source-concentration layers hold. Absence returns None (never-measured), never a clean bill.
    """
    card = _doubt_card()
    if not card:
        return None
    d = str(drug).strip().lower()
    for cell in card.get("determinant_completeness_arm", []):
        if str(cell.get("drug", "")).lower() != d:
            continue
        if cell.get("status") != "scored":
            return {"arm": "determinant_completeness", "status": cell.get("status"),
                    "note": cell.get("note", ""), "source": "doubt_layer_per_cell"}
        strong = cell.get("strong") or []
        return {"arm": "determinant_completeness", "status": "scored",
                "n_strong_completeness_signals": cell.get("n_strong", 0),
                "n_families_uncounted": cell.get("n_families_uncounted"),
                "known_gap_recovered": cell.get("known_gap_recovered", False),
                "strong_families": [s["evidence"]["symbol"] for s in strong],
                "source": "doubt_layer_per_cell"}
    for cell in card.get("position_novelty_arm", []):
        if str(cell.get("drug", "")).lower() == d:
            return {"arm": "position_novelty", "status": "scored",
                    "sens_on_blindspot": cell.get("sens_on_blindspot"),
                    "fp_on_catalog_negative_S": cell.get("fp_on_catalog_negative_S"),
                    "lift": cell.get("lift"), "source": "doubt_layer_per_cell"}
    return None


def _cell_layer_for(drug: str, organism: str | None, key: str) -> dict | None:
    """A named disclosure block from this cell's report-card row, or None. Read-only.

    Reuses the card's own organism matcher so a layer can never be attached to the wrong cell, and
    returns None on an ambiguous match rather than guessing.
    """
    deck = _load("decoder_validation_report_card.json")
    if not deck:
        return None
    c, status = _pick_bacterial_cell(deck.get("cells", []), str(drug).strip().lower(), organism)
    if not c or status == "ambiguous":
        return None
    block = c.get(key)
    return block if isinstance(block, dict) else None


# Every per-cell disclosure block the report card can carry. DERIVED consumers key off this rather
# than hand-listing, so a new layer cannot silently stay CLI-unreachable -- which is exactly what
# happened to `lineage` and `source_concentration`: both rendered on the card for cells whose calls
# never mentioned them, and the gap was invisible until the layers were enumerated against the CLI.
DISCLOSURE_LAYERS = ("lineage", "source_concentration", "prospective", "doubt_layer")


def trust_block(drug: str, organism: str | None = None) -> dict:
    """Public always-safe accessor for embedding in an amr-mechanism-call-v1 record's `validation` field.

    CROSS-CUTTING prospective annotation: a cell's badge tier may resolve from ANY card (AMR-Portal, TB,
    HIV, ...), so the prospective-regression note is attached HERE rather than inside one card's branch --
    E. coli x gentamicin resolves at the AMR-Portal card and would otherwise never learn that its own
    post-lock cohort contradicts it. The tier and metric are UNCHANGED; this only adds the contradiction.
    """
    badge = lookup_trust(drug, organism)
    pr = prospective_regression_for(drug, organism)
    if pr:
        badge["prospective_regression"] = pr
        badge["caveat"] = (
            f"{badge.get('caveat', '')} || PROSPECTIVE REGRESSION: post-lock sens {pr['sens']} on "
            f"N={pr['n']} isolates public after {pr['lock_date']} -- the frozen rule UNDER-CALLS this "
            f"cell; the accuracy above is in-distribution and does not describe it")
    # L2 DOUBT, augment-only. Attached under its OWN key and NEVER folded into `tier`/`headline`/the
    # metrics: "does this cell have a known catalog-completeness gap" is a different question from
    # "how well is this cell validated", and merging them would silently overwrite one with the
    # other -- the shared-key trap this project has hit before. Absent -> the key is simply absent.
    dl = doubt_layer_for(drug)
    if dl:
        badge["doubt_layer"] = dl
    # LINEAGE + SOURCE CONCENTRATION + the full PROSPECTIVE block, each under its own key. These
    # rendered on the standing report card while every CLI call stayed silent about them, which
    # matters most where it is worst: `escherichia_coli_shigella x gentamicin` is 95% ONE BioProject
    # and holds zero carriers of the `rmt` family, so its 0.893 is a statement about one hospital's
    # isolates (source-diverse measurements of the same cell report 0.523). A caller deciding on that
    # number could not see the caveat that explains it.
    #
    # `prospective` is attached WHENEVER it exists, not only when it regresses: the pre-existing
    # `prospective_regression` key surfaces bad news only, so a caller could not distinguish "no
    # post-lock data" from "post-lock data that AGREED". Both keys ship; neither replaces the other.
    for _layer in ("lineage", "source_concentration", "prospective"):
        _b = _cell_layer_for(drug, organism, _layer)
        if _b:
            badge[_layer] = _b
    return badge


def one_line(badge: dict) -> str:
    """Compact human-readable badge for CLI output."""
    head = f" -- {badge['headline']}" if badge.get("headline") else ""
    return f"validation: {badge['tier']}{head}  ({badge['caveat']}; see {badge['source_card'] or 'n/a'})"


def lineage_one_line(badge: dict) -> str | None:
    """The clonality correction, compactly. None when there is nothing renderable.

    A cell's raw sens/spec counts one vote per ISOLATE, so an over-sampled clone carries the metric;
    the lineage layer collapses each same-label lineage to one vote. Every SCORED R class here is
    clonally dominated, which is why this belongs beside the headline rather than only in the wiki.

    NEVER prints a weighted point without its interval — the effective N is tiny, so the CI IS the
    result. That mirrors `build_validation_report_card._assert_weighted_renderable`, which refuses to
    render such a point at all; a compact renderer that quietly dropped the CI would undo the guard.
    """
    lin = badge.get("lineage")
    if not isinstance(lin, dict) or lin.get("status") != "scored":
        return None
    weighted = lin.get("cluster_weighted") or {}
    if not weighted:
        return None
    thr = min(weighted, key=lambda t: float(t))       # finest available: least collapsed
    w = weighted[thr] or {}
    eff = (lin.get("effective_lineage_N") or {}).get(thr) or {}
    head = (f"clonality: raw N={lin.get('raw_N')} collapses to {eff.get('R', '?')}R/"
            f"{eff.get('S', '?')}S effective lineages at Mash {thr}")
    sens, s_ci = w.get("sens"), w.get("sens_ci")
    spec, p_ci = w.get("spec"), w.get("spec_ci")
    if sens is not None and s_ci and spec is not None and p_ci:
        head += (f"; cluster-weighted sens {sens} [{s_ci[0]}-{s_ci[1]}], "
                 f"spec {spec} [{p_ci[0]}-{p_ci[1]}]")
    grade = lin.get("grade")
    return f"{head}{f' -- {grade}' if grade else ''}"


def concentration_one_line(badge: dict) -> str | None:
    """The source-concentration caveat, when it is decision-relevant. None otherwise.

    Printed only for a SINGLE-SOURCE cell, because that is where the headline metric stops describing
    the population a caller is applying it to. The measured case: `escherichia_coli_shigella x
    gentamicin` reports sens 0.893 from a cohort that is 95% one BioProject and contains no carriers
    of the `rmt` family; source-diverse measurements of the same cell with the same frozen rule report
    0.523. The number was never wrong about its cohort -- it was a statement about one hospital.

    The direction is NOT always optimistic: a 97%-single-source ciprofloxacin cell reads PESSIMISTIC
    (spec 0.700 vs 0.988 on an 8-BioProject set). So this says the estimate is narrow, never that it
    is inflated.
    """
    s = badge.get("source_concentration")
    if not isinstance(s, dict) or s.get("status") != "measured" or not s.get("single_source"):
        return None
    share = s.get("largest_share")
    share_s = "" if share is None else f" ({share:.0%} of the cohort)"
    n_bp = s.get("distinct_bioprojects")
    return (f"source concentration: SINGLE-SOURCE -- {n_bp} BioProject(s), one dominant{share_s}. "
            "The metric above describes that source's isolates; it is a narrow estimate, in either "
            "direction, not necessarily an inflated one")
