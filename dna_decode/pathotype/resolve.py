"""v0 deterministic decision table: cluster profile -> derived pathotype call.

Implements the ledger-locked 11-class decision table + abstention rules. Pure
function of a boolean cluster profile (+ optional partial-hit / QC context); no I/O,
no detection — so it is exhaustively unit-testable on synthetic profiles.

Compatibility-resolver framing (ledger v5): every call is `*_COMPATIBLE` — a
genotype compatibility statement with abstention, NOT a clinical prediction. Calls
in SUPPORTED_CLASSES have external-validity support (ExPEC/EPEC/ETEC); the rest are
marked scope-limited.
"""
from __future__ import annotations

from dna_decode.pathotype.markers import (
    EXPEC_STRONG, EXPEC_SUPPORT, SUPPORTED_CLASSES, RULES_VERSION,
)

AAF_CLUSTERS = ["AAF_I", "AAF_II", "AAF_III", "AAF_IV", "AAF_V"]


def resolve_call(profile: dict[str, bool], *,
                 partial_clusters: frozenset[str] = frozenset(),
                 qc_pass: bool = True) -> dict:
    """Resolve a cluster presence profile to a derived pathotype call.

    Args:
        profile: cluster -> bool (confidently present). Missing key == False.
        partial_clusters: clusters with sub-threshold (PARTIAL) hits -> trigger
            AMBIGUOUS when they are the only evidence for a primary call.
        qc_pass: assembly QC verdict. False -> AMBIGUOUS_LOW_QC (never commensal).

    Returns:
        derived_call dict (primary, secondary, confidence_tier, rule_id, reason,
        external_validity).
    """
    g = lambda c: bool(profile.get(c, False))

    if not qc_pass:
        return _call("AMBIGUOUS_LOW_QC", [], "AMBIGUOUS", "RULE_QC_001",
                     "assembly QC failed; marker absence is not interpretable", partial_clusters)

    stx = g("STX1") or g("STX2")
    lee = g("LEE")
    bfp = g("BFP_EAF")
    etec = g("LT") or g("ST")
    aaf = any(g(c) for c in AAF_CLUSTERS)
    aggR, aat, aai = g("EAEC_REG"), g("EAEC_TRANSPORT"), g("EAEC_T6SS")
    eaec_confident = (aggR and (aaf or aat or aai)) or (aaf and (aat or aai))
    eaec_partial = (aggR or aaf or aat or aai) and not eaec_confident

    expec_strong = sum(g(c) for c in EXPEC_STRONG)
    expec_support = sum(g(c) for c in EXPEC_SUPPORT)
    upec = expec_strong >= 2

    eiec = g("EIEC_FLAG")
    daec = g("AFA_DRA") and expec_strong == g("AFA_DRA") and not (stx or lee or etec or eaec_confident) \
        and not (g("P_FIMBRIAE") or g("S_FIMBRIAE") or g("HEMOLYSIN") or g("CNF1"))
    daec = daec or g("DAEC_FLAG")

    # --- primary DEC modules (intestinal pathotypes) confidently present ---
    module_calls: list[str] = []
    if stx and lee:
        module_calls.append("EHEC_COMPATIBLE")
    elif stx:
        module_calls.append("STEC_NON_LEE")
    if lee and not stx:
        module_calls.append("tEPEC_COMPATIBLE" if bfp else "aEPEC_COMPATIBLE")
    if etec:
        module_calls.append("ETEC_COMPATIBLE")
    if eaec_confident:
        module_calls.append("EAEC_COMPATIBLE")

    # distinct DEC pathotype modules (EHEC/STEC count as one stx module for hybrid count)
    dec_modules = set()
    if stx:
        dec_modules.add("STX")
    if lee and not stx:
        dec_modules.add("EPEC")
    if etec:
        dec_modules.add("ETEC")
    if eaec_confident:
        dec_modules.add("EAEC")

    # --- out-of-scope evidence with no in-scope DEC module ---
    if (eiec or daec) and not module_calls and not upec:
        which = "ipaH/EIEC" if eiec else "afa/dra-only DAEC"
        return _call("UNCLASSIFIED", [], "UNCLASSIFIED", "RULE_OOS_001",
                     f"out-of-v0-scope evidence ({which}); not a six-class DEC pathotype", partial_clusters)

    # --- HYBRID: >=2 primary DEC modules confidently present ---
    if len(dec_modules) >= 2:
        return _call("HYBRID", sorted(module_calls), "HYBRID", "RULE_HYBRID_001",
                     f"{len(dec_modules)} primary DEC modules present: {sorted(dec_modules)}; "
                     f"reported as multilabel, no forced primary", partial_clusters)

    # --- single confident DEC module ---
    if len(module_calls) == 1:
        primary = module_calls[0]
        secondary = ["UPEC_COMPATIBLE"] if upec else []
        return _call(primary, secondary, "CONFIDENT", _rule_for(primary),
                     _reason_for(primary, stx, lee, bfp, etec), partial_clusters)

    # --- no DEC module: ExPEC / ambiguous / commensal ---
    if upec:
        return _call("UPEC_COMPATIBLE", [], "CONFIDENT", "RULE_UPEC_001",
                     f"{expec_strong} strong ExPEC markers (+{expec_support} support); "
                     f"extraintestinal-compatible (definitive UPEC needs clinical metadata)",
                     partial_clusters)

    # ambiguous: partial EAEC, single strong ExPEC, or partial primary hits
    if eaec_partial or expec_strong == 1 or (partial_clusters & _primary_clusters()):
        bits = []
        if eaec_partial:
            bits.append("incomplete EAEC module (regulator/transport/AAF not co-present)")
        if expec_strong == 1:
            bits.append("single strong ExPEC marker (UPEC needs >=2)")
        if partial_clusters & _primary_clusters():
            bits.append(f"sub-threshold primary hits: {sorted(partial_clusters & _primary_clusters())}")
        return _call("AMBIGUOUS", [], "AMBIGUOUS", "RULE_AMB_001",
                     "; ".join(bits), partial_clusters)

    # nothing fired, QC ok
    return _call("COMMENSAL_LOW_MARKER_BURDEN", [], "CONFIDENT", "RULE_COMM_001",
                 "no DEC module and sub-threshold ExPEC markers (absence of v0 markers, "
                 "not biological commensal truth)", partial_clusters)


def _primary_clusters() -> frozenset[str]:
    return frozenset(["STX1", "STX2", "LEE", "LT", "ST",
                      "EAEC_REG", "EAEC_TRANSPORT", "EAEC_T6SS"] + AAF_CLUSTERS)


def _rule_for(call: str) -> str:
    return {
        "EHEC_COMPATIBLE": "RULE_EHEC_001", "STEC_NON_LEE": "RULE_STEC_001",
        "tEPEC_COMPATIBLE": "RULE_tEPEC_001", "aEPEC_COMPATIBLE": "RULE_aEPEC_001",
        "ETEC_COMPATIBLE": "RULE_ETEC_001", "EAEC_COMPATIBLE": "RULE_EAEC_001",
    }.get(call, "RULE_GEN_001")


def _reason_for(call, stx, lee, bfp, etec) -> str:
    return {
        "EHEC_COMPATIBLE": "STX + LEE detected",
        "STEC_NON_LEE": "STX without LEE (LEE-negative STEC)",
        "tEPEC_COMPATIBLE": "LEE + BFP/EAF, no STX (typical EPEC)",
        "aEPEC_COMPATIBLE": "LEE without BFP/EAF, no STX (atypical EPEC)",
        "ETEC_COMPATIBLE": "LT and/or ST enterotoxin detected",
        "EAEC_COMPATIBLE": "aggR/AAF + transport/T6SS support",
    }.get(call, "single DEC module")


def _call(primary, secondary, tier, rule_id, reason, partial_clusters) -> dict:
    return {
        "primary": primary,
        "secondary": secondary,
        "confidence_tier": tier,
        "rule_id": rule_id,
        "rule_version": RULES_VERSION,
        "reason": reason,
        "external_validity": "supported" if primary in SUPPORTED_CLASSES else "scope_limited",
        "partial_clusters": sorted(partial_clusters),
    }
