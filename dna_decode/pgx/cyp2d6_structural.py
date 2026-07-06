"""CYP2D6 STRUCTURAL surface (v0.2) — read-depth copy-number caller (the "structural caller" the SNP
surface cannot be).

The SNP cell (`cyp2d6_catalog` + `cyp2d6_caller`) is honestly blind to structural alleles (gene deletion
*5, duplications *xN, CYP2D6-CYP2D7 hybrids *13/*36/*68) — those need read-depth on a BAM/CRAM, not a phased
VCF (`cnv_hybrid_unassessed=true`). This module is that missing surface: from a CYP2D6-region vs
single-copy-control DEPTH RATIO (computed off a real CRAM/BAM), it calls integer COPY NUMBER and a
structural status — DELETION (*5-consistent) / NORMAL / DUPLICATION (*xN-consistent).

EMPIRICAL CALIBRATION (real 1000G 30x CRAMs, verified 2026-07-06 — see wiki/cyp2d6_structural_*):
  ratio = mean_depth(chr22:42126000-42130000, CYP2D6 body) / mean_depth(chr22:41000000-41050000, control).
  Confirmed truth classes separate cleanly:
    DELETION  (1 copy: *X/*5): ratio ~0.58-0.66
    NORMAL    (2 copy):        ratio ~1.26   <- the NORMAL_BASELINE
    DUPLICATION (>=3 copy):    ratio ~1.73-2.29
  So per-copy depth-fraction = NORMAL_BASELINE / 2, and CN = round(ratio / per_copy) = round(2*ratio/BASE).
  (BASE > 1.0 for a 2-copy region because CYP2D7 paralog reads leak into the CYP2D6 body + region GC/
  mappability — an empirical baseline, NOT a physical 1.0.)

HONEST SCOPE (load-bearing — this is COPY NUMBER, not full structural typing):
  * RESOLVES: *5 gene deletion (CN drop) + *xN duplication (CN gain). A het *5 -> CN 1; a duplication -> CN>=3.
  * DOES NOT RESOLVE: hybrid ALLELE IDENTITY (*13 vs *36 vs *68) — that needs PSV (paralogous-sequence-
    variant) analysis of CYP2D6-vs-CYP2D7 reads (the full Cyrius/Aldy/StellarPGx approach). A hybrid that
    changes copy number (e.g. *36x2) shows as a CN change; its allele label is NOT called here.
  * The raw single-region depth ratio is paralog-confounded at the edges; this v0.2 is a COARSE integer-CN
    caller, not a base-pair-precise breakpoint caller. Reference tool: Cyrius (Illumina).
"""
from __future__ import annotations

from dataclasses import dataclass

# CYP2D6 gene body + a nearby single-copy control (GRCh38 chr22). The gene body overlaps the CYP2D6/CYP2D7
# homology, so the baseline is empirical (below), not a physical 1.0.
CYP2D6_REGION = ("chr22", 42126000, 42130000)
CONTROL_REGION = ("chr22", 41000000, 41050000)

# Empirically calibrated on confirmed 2-copy 1000G samples (median of the NORMAL-truth ratios). A 2-copy
# CYP2D6 body reads ~1.26x the single-copy control (paralog leak + GC/mappability). Refined by
# scripts/cyp2d6_structural_probe.py from the real depth run; pinned here as the deployed constant.
NORMAL_BASELINE = 1.26
MAX_COPY_NUMBER = 6   # clamp (defensive; CYP2D6 copy number rarely exceeds ~4-5)

# --- HYBRID DETECTION (v0.3) — the CYP2D7 paralog read-depth signal ------------------------------------
# A CYP2D6-CYP2D7 HYBRID allele (*13/*36/*68) carries extra CYP2D7-derived sequence, so it ELEVATES the
# CYP2D7 read depth. Measuring depth over the CYP2D7 paralog (just telomeric of CYP2D6) vs the control gives
# a hybrid signal that DISTINGUISHES a hybrid from a pure CYP2D6 duplication (a pure *xN dup elevates CYP2D6
# but NOT CYP2D7). This is the read-depth-tractable half of hybrid resolution; the exact IDENTITY (*13 vs
# *36 vs *68) still needs PSV analysis of CYP2D6-vs-CYP2D7 reads (Cyrius-class, future).
CYP2D7_REGION = ("chr22", 42139500, 42144300)
# Threshold calibrated on the 38-sample structural set (scripts/cyp2d6_hybrid_validate): d7/ctl >= 1.25 gives
# spec 1.00 / sens 0.62 — a HIGH-SPECIFICITY detector (a positive is trustworthy; it never fired on a pure
# dup or normal). Sensitivity is partial: the *68 family (common, non-functional) detects cleanly (>1.5);
# some subtle *36 + the opposite-signature *13 (LOW d7) are missed.
HYBRID_D7_THRESHOLD = 1.25


def hybrid_suspected(cyp2d7_depth_ratio: float, threshold: float = HYBRID_D7_THRESHOLD) -> bool:
    """True if the CYP2D7-paralog/control depth ratio is elevated -> a CYP2D6-CYP2D7 hybrid is present.
    HIGH-SPECIFICITY: a True is trustworthy (spec 1.0 in validation); a False does NOT rule out a subtle
    *36/*13 hybrid (sens ~0.62)."""
    return cyp2d7_depth_ratio >= threshold


def hybrid_detection(cyp2d7_depth_ratio: float, threshold: float = HYBRID_D7_THRESHOLD) -> dict:
    """CYP2D7 depth ratio -> hybrid-presence call (NOT identity). Detects that a CYP2D6-CYP2D7 hybrid
    (*13/*36/*68) is present; never resolves WHICH (that needs PSV analysis, Cyrius-class)."""
    susp = hybrid_suspected(cyp2d7_depth_ratio, threshold)
    note = ("elevated CYP2D7 paralog depth -> a CYP2D6-CYP2D7 HYBRID allele (*13/*36/*68) is PRESENT "
            "(high-specificity detector; a positive is trustworthy). Detects PRESENCE only — the exact "
            "identity (*13 vs *36 vs *68) needs CYP2D6-vs-CYP2D7 PSV analysis (Cyrius-class)."
            if susp else
            "CYP2D7 depth not elevated — no hybrid detected (a subtle *36 or the opposite-signature *13 "
            "hybrid can still be missed; sens ~0.62 in validation).")
    return {"cyp2d7_depth_ratio": round(cyp2d7_depth_ratio, 3), "hybrid_suspected": susp,
            "hybrid_identity_resolved": False, "note": note}


@dataclass(frozen=True)
class StructuralCall:
    copy_number: int
    status: str            # "deletion" | "normal_copy_number" | "duplication" | "homozygous_deletion"
    star_consistent: str   # the star-allele class this CN is consistent with (NOT a specific allele)
    depth_ratio: float
    hybrid_identity_unresolved: bool   # ALWAYS True — this surface never resolves *13/*36/*68 identity
    note: str


def copy_number_from_ratio(ratio: float, normal_baseline: float = NORMAL_BASELINE) -> int:
    """Integer CYP2D6 copy number from the depth ratio. CN = round(2*ratio/baseline), clamped to [0, MAX].
    A 2-copy sample (ratio==baseline) -> CN 2; a het deletion (ratio~baseline/2) -> CN 1."""
    if ratio is None or ratio < 0 or normal_baseline <= 0:
        raise ValueError(f"invalid ratio/baseline: {ratio!r}/{normal_baseline!r}")
    cn = round(2.0 * ratio / normal_baseline)
    return max(0, min(MAX_COPY_NUMBER, cn))


def structural_call(ratio: float, normal_baseline: float = NORMAL_BASELINE) -> StructuralCall:
    """Depth ratio -> copy-number structural call (deletion / normal / duplication). NEVER resolves hybrid
    allele identity (that needs PSV analysis) -> hybrid_identity_unresolved is always True."""
    cn = copy_number_from_ratio(ratio, normal_baseline)
    if cn == 0:
        status, star = "homozygous_deletion", "*5/*5 (both copies deleted)"
    elif cn == 1:
        status, star = "deletion", "*X/*5 (one gene deletion) — SNP call is the single remaining allele"
    elif cn == 2:
        status, star = "normal_copy_number", "two CYP2D6 copies (no deletion/duplication)"
    else:
        status, star = "duplication", f"{cn} copies — *xN duplication (or a copy-adding hybrid, e.g. *36xN)"
    note = ("Copy number from CYP2D6/control read-depth ratio (real CRAM/BAM). Resolves *5 deletion + *xN "
            "duplication ONLY; hybrid IDENTITY (*13/*36/*68) is NOT resolved (needs CYP2D6-vs-CYP2D7 PSV "
            "analysis — Cyrius-class). Coarse integer-CN, not a base-pair breakpoint caller.")
    return StructuralCall(copy_number=cn, status=status, star_consistent=star, depth_ratio=round(ratio, 3),
                          hybrid_identity_unresolved=True, note=note)


_HYBRID_STARS = ("*13", "*36", "*61", "*63", "*68")


def truth_copy_number_class(getrm_truth: str) -> str | None:
    """Map a GeT-RM CYP2D6 truth string to its COPY-NUMBER class (deletion/normal/duplication) for
    validation. Returns None (EXCLUDED from CN scoring) when the net copy number is NOT determinable from
    the star label alone — i.e. any truth containing a CYP2D6-CYP2D7 HYBRID (*13/*36/*61/*63/*68; its copy
    contribution depends on tandem arrangement) or a del+dup mix. This is the *structural copy-number* axis,
    distinct from SNP-allele identity; hybrid samples are still MEASURED (depth), just not CN-scored."""
    import re
    s = getrm_truth.strip()
    if any(h in s for h in _HYBRID_STARS):
        return None          # hybrid copy contribution unresolved by the star label -> not CN-scored
    has_del = "*5" in re.split(r"[/|+()]", s)
    has_dup = bool(re.search(r"[xX]\d", s))
    if has_del and has_dup:
        return None          # del + dup -> ambiguous net CN
    if has_dup:
        return "duplication"
    if has_del:
        # *5/*5 -> homozygous_deletion; *X/*5 -> deletion (one copy)
        alleles = [a for a in re.split(r"[/|]", s) if a.strip()]
        return "homozygous_deletion" if all("*5" in a for a in alleles) else "deletion"
    return "normal_copy_number"
