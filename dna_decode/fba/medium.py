"""Growth media for FBA validation — match the medium to how the LABELS were measured.

A gene-essentiality label is only meaningful relative to a growth condition: an amino-acid biosynthesis
gene is essential on minimal medium and dispensable on rich, and that is biology, not model error. So a
validation that scores a MINIMAL-medium model against RICH-medium labels reports a mismatch as a defect.

Measured on yeast/iMM904 vs the SGD gold standard (which comes from the deletion collection on **YPD**,
a rich medium) while iMM904 ships a glucose-MINIMAL default:

| medium | TP | FP | FN | TN | MCC | precision |
|---|---|---|---|---|---|---|
| minimal (model default) | 43 | **67** | 92 | 703 | **0.2524** | 0.391 |
| YPD-like rich (label-matched) | 34 | **13** | 101 | 757 | **0.3773** | **0.723** |

False positives collapse by 81% and MCC rises ~50% relative, from a config change with no new biochemistry
— exactly the direction the mismatch predicts, since a minimal-medium model must call biosynthesis genes
essential that are dispensable when the nutrient is supplied.
"""
from __future__ import annotations

# Exchange reactions supplied by a rich (YPD-like) medium: the 20 proteinogenic amino acids plus the
# nucleobases/nucleosides a rich yeast-extract/peptone medium provides. BiGG ids; any absent from a given
# model are skipped rather than raising, so this is portable across reconstructions.
RICH_MEDIUM_EXCHANGES: tuple[str, ...] = (
    "EX_ala__L_e", "EX_arg__L_e", "EX_asn__L_e", "EX_asp__L_e", "EX_cys__L_e", "EX_gln__L_e",
    "EX_glu__L_e", "EX_gly_e", "EX_his__L_e", "EX_ile__L_e", "EX_leu__L_e", "EX_lys__L_e",
    "EX_met__L_e", "EX_phe__L_e", "EX_pro__L_e", "EX_ser__L_e", "EX_thr__L_e", "EX_trp__L_e",
    "EX_tyr__L_e", "EX_val__L_e",
    "EX_ade_e", "EX_gua_e", "EX_ura_e", "EX_csn_e", "EX_thym_e", "EX_ins_e",
)

DEFAULT_SUPPLEMENT_UPTAKE = 10.0


def rich_medium(model, uptake: float = DEFAULT_SUPPLEMENT_UPTAKE) -> dict[str, float]:
    """The model's own default medium PLUS the rich supplements it actually has exchanges for.

    Additive by design: the carbon source, oxygen and mineral bounds a reconstruction ships are kept, and
    only the supplements are opened. Returns the medium dict; the caller assigns it (so a caller can diff
    or log it before applying).
    """
    have = {r.id for r in model.exchanges}
    med = dict(model.medium)
    for ex in RICH_MEDIUM_EXCHANGES:
        if ex in have:
            med[ex] = uptake
    return med


def apply_rich_medium(model, uptake: float = DEFAULT_SUPPLEMENT_UPTAKE) -> list[str]:
    """Set the model to a rich medium in place. Returns the exchange ids actually opened."""
    before = set(model.medium)
    model.medium = rich_medium(model, uptake=uptake)
    return sorted(set(model.medium) - before)
