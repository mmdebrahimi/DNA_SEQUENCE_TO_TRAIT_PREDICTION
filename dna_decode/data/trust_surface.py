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


def trust_block(drug: str, organism: str | None = None) -> dict:
    """Public always-safe accessor for embedding in an amr-mechanism-call-v1 record's `validation` field."""
    return lookup_trust(drug, organism)


def one_line(badge: dict) -> str:
    """Compact human-readable badge for CLI output."""
    head = f" -- {badge['headline']}" if badge.get("headline") else ""
    return f"validation: {badge['tier']}{head}  ({badge['caveat']}; see {badge['source_card'] or 'n/a'})"
