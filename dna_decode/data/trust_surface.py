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

_WIKI = Path(__file__).resolve().parent.parent.parent / "wiki"

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
    p = _WIKI / name
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


def _rec(tier: str, source_card: str, headline: str = "", metric: float | None = None,
         n: int | None = None, cell: str = "") -> dict:
    return {
        "tier": tier,
        "independent": tier in (INDEPENDENT_WETLAB, INDEPENDENT_MEASURED),
        "headline": headline,
        "metric": metric,
        "n": n,
        "cell": cell,
        "source_card": source_card,
        "caveat": _CAVEAT[tier],
    }


def lookup_trust(drug: str, organism: str | None = None) -> dict:
    """Best-evidence honest validation badge for a (drug, organism) decoder cell. Always returns a dict
    (tier UNKNOWN if nothing matches); never fabricates a metric, never averages across tiers."""
    d = (drug or "").strip().lower()
    g = _genus(organism)

    # 1. HIV — free wet-lab fold-change (matched by drug; HIV has no organism axis in the surface)
    hiv = _load("hiv_decoder_report_card.json")
    if hiv:
        for c in hiv.get("cells", []):
            if str(c.get("drug", "")).strip().lower() == d:
                auc = c.get("auc_call_separates_fold")
                return _rec(INDEPENDENT_WETLAB, "wiki/hiv_decoder_report_card.md",
                            headline=f"AUC {auc} (N={c.get('n')}, {c.get('drug_class')})" if auc else "scored",
                            metric=auc, n=c.get("n"), cell=f"HIV-1|{d}")

    # 2. TB — independent measured AST (matched by drug)
    tb = _load("tb_report_card.json")
    if tb:
        for c in tb.get("independent", []):
            if str(c.get("drug", "")).strip().lower() == d:
                acc = c.get("raw_acc")
                n = (c.get("n_R") or 0) + (c.get("n_S") or 0)
                return _rec(INDEPENDENT_MEASURED, "wiki/tb_report_card.md",
                            headline=f"acc {round(acc, 3)} (N={n})" if acc is not None else "scored",
                            metric=acc, n=n, cell=f"M.tuberculosis|{d}")

    # 3. bacteria — EBI AMR Portal INDEPENDENT measured AST (genus + drug)
    portal = _load("amr_portal_independent_report_card.json")
    if portal:
        for c in portal.get("cells", []):
            if str(c.get("drug", "")).strip().lower() == d and _genus(c.get("organism")) == g and g:
                acc = c.get("accuracy")
                n = (c.get("n_R") or 0) + (c.get("n_S") or 0)
                tier = UNDERPOWERED if str(c.get("tier", "")).upper().startswith("UNDERPOWER") else INDEPENDENT_MEASURED
                return _rec(tier, "wiki/amr_portal_independent_report_card.md",
                            headline=f"acc {round(acc, 3)} (N={n})" if acc is not None else "scored",
                            metric=acc, n=n, cell=f"{c.get('organism')}|{d}")

    # 4. bacteria — NCBI-PD provenance-disjoint card (genus + drug); also carries the structural non-cells
    deck = _load("decoder_validation_report_card.json")
    if deck:
        for c in deck.get("cells", []):
            if str(c.get("drug", "")).strip().lower() == d and _genus(c.get("organism")) == g and g:
                st = str(c.get("state", "")).upper()
                if st == "SCORED":
                    acc = c.get("acc")
                    return _rec(PROVENANCE_DISJOINT, "wiki/decoder_validation_report_card.md",
                                headline=f"acc {acc} (N={c.get('n')})" if acc is not None else "scored",
                                metric=acc, n=c.get("n"), cell=f"{c.get('organism')}|{d}")
                _STATE_TIER = {"ABSTAINS_BY_DESIGN": ABSTAINS_BY_DESIGN, "LABEL_CONFOUNDED": LABEL_CONFOUNDED,
                               "UNDERPOWERED": UNDERPOWERED, "NO_FREE_PHENOTYPE_SOURCE": NO_FREE_PHENOTYPE_SOURCE,
                               "NOT_CENSUSED": NOT_CENSUSED}
                tier = _STATE_TIER.get(st, NOT_CENSUSED)
                return _rec(tier, "wiki/decoder_validation_report_card.md", cell=f"{c.get('organism')}|{d}")

    # 5. SARS-CoV-2 — in-distribution knowledge baseline (separate namespace, underpowered)
    if d in _SARSCOV2_DRUGS:
        return _rec(IN_DISTRIBUTION, "wiki/sarscov2_mpro_validation_result_2026-06-23.md",
                    headline="in-distribution (CoV-RDB), underpowered", cell=f"SARS-CoV-2|{d}")

    # 6. shipped-surface structural fallback (no card cell yet)
    try:
        from dna_decode.data.shipped_decoder_surface import SHIPPED_DECODER_SURFACE
        for (org, drg, _eng, _scope, status, _grp) in SHIPPED_DECODER_SURFACE:
            if drg.strip().lower() == d and (not g or _genus(org) == g):
                tier = {"no_free_source": NO_FREE_PHENOTYPE_SOURCE, "label_confounded": LABEL_CONFOUNDED,
                        "ncbi_pd": NOT_CENSUSED}.get(status, NOT_CENSUSED)
                return _rec(tier, "dna_decode/data/shipped_decoder_surface.py", cell=f"{org}|{d}")
    except Exception:
        pass

    return _rec(UNKNOWN, "", cell=f"{organism or '?'}|{d}")


def trust_block(drug: str, organism: str | None = None) -> dict:
    """Public always-safe accessor for embedding in an amr-mechanism-call-v1 record's `validation` field."""
    return lookup_trust(drug, organism)


def one_line(badge: dict) -> str:
    """Compact human-readable badge for CLI output."""
    head = f" -- {badge['headline']}" if badge.get("headline") else ""
    return f"validation: {badge['tier']}{head}  ({badge['caveat']}; see {badge['source_card'] or 'n/a'})"
