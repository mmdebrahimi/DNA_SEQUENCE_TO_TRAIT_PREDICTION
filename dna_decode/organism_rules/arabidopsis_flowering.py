"""Arabidopsis thaliana flowering-HABIT decoder — the plant analog of the AMR/pigmentation curated cells.

The DETERMINISTIC counterpart to a closed negative: the Arabidopsis flowering-time EMBEDDING test failed under
de-confounding (2026-06-12, PlantCaduceus within-group r2 -0.13 — it learned population structure, not the
causal signal; the 3rd de-confounded embedding failure across the kingdom boundary). But the causal signal is
CATALOGUED: natural flowering-habit variation is driven predominantly by discrete loss-of-function variants at
two interacting loci. This cell reads that mechanism directly.

THE BIOLOGY (the rule this encodes):
    FLC (FLOWERING LOCUS C) is a dosage-dependent floral REPRESSOR; FRI (FRIGIDA) is required for HIGH FLC.
    So the WINTER-ANNUAL habit (late, vernalization-requiring) needs BOTH: functional FRI *AND* a strong FLC.
    Losing EITHER gives the SUMMER-ANNUAL / rapid-cycling habit (early, no vernalization requirement).

    => late  iff  FRI functional AND FLC strong          (a multi-locus AND)
    => early iff  FRI LoF  OR  FLC weak/null             (either route)

    FLC is DOWNSTREAM of FRI, so an FLC null/weak allele calls EARLY regardless of FRI (higher confidence);
    the FRI-LoF route is confidence-capped by the Lz-0 counterexample (FRI deletion yet LATE, via FRI-
    INDEPENDENT FLC upregulation — FRL1/FES1-class activators the two-locus rule cannot see).

WHY THIS IS A NON-FROZEN organism_rules CELL: the rule is a multi-locus AND with an epistasis caveat — a
shape the frozen count/OR `amr_rules.DRUG_RULE` engine cannot express (same reason the TMP-SMX `sul AND dfr`
overlay and the TB cell live outside it). The frozen decoder surface is untouched.

VALIDATED (both routes, 2026-07-16/17 — the rule's shape is measured, not asserted):
    FRI route  — Zhang 2020 Table S3, N=854 (`scripts/flowering_tables3_score.py`). The honest figure is
                 the population-structure-weighted 0.710 vs its own 0.676 null (+3.4pp), NOT the pooled
                 0.733-vs-0.502: FRI genotype tracks ancestry, so pooled accuracy partly measures "can you
                 recognise a central European accession?". Directional: FRI-LoF->early 93.9% (strong) but
                 FRI-functional->late only 65.8% (weak) — functional FRI is NECESSARY, NOT SUFFICIENT,
                 which is exactly what the AND below says and why the FLC route matters.
    FLC route  — the distinctive claim, validated on n=106 (`scripts/flowering_flc_route_test.py`) by
                 joining measured FLC EXPRESSION (AraPheno phenotype 29, Atwell 2010) to S3. ALL FOUR
                 cells of the AND call their majority correctly:
                     functional + strong -> 85% late   (LATE  ✓)
                     functional + weak   -> 39% late   (EARLY ✓)  <- the Da(1)-12 class, 46pp separation
                     lof + strong        -> 17% late   (EARLY ✓)  <- the Lz-0 class: REAL but 1/6 = RARE,
                                                                     which JUSTIFIES the MEDIUM cap below
                     lof + weak          -> 10% late   (EARLY ✓)
                 FLC earns its place: net +5 calls fixed (14 rescued / 9 broken) on the 70 functional-FRI
                 accessions a FRI-only rule calls ALL late; **within-ancestry 0.803 vs FRI-only 0.767 vs
                 null 0.751 — the FLC route roughly TRIPLES the within-ancestry advantage.**
                 CAVEAT that travels with that number: the gain RIDES ON THE THRESHOLD (q30 +0.066 / q50
                 +0.047 / q60 +0.000 / q70 -0.085), holding only in the biologically plausible low-quantile
                 range (weak FLC alleles are RARE — Werner 2005 — which a median split cannot represent).
                 And FLC expression is a PROXY for allele status, not the same measurement.
    The premise (FRI drives FLC) is confirmed too: median FLC 1.265 with functional FRI vs 0.164 with
    FRI-LoF — a 7.7x ratio.

HONEST SCOPE (load-bearing):
  - PARTIAL. FRI/FLC explains ~40-70% of long-day flowering-time variation; the rest is polygenic +
    environment (photoperiod, vernalization, temperature). This decodes the HABIT / DIRECTION
    (early / summer-annual vs late / winter-annual), NOT quantitative days-to-flower. Anything outside the
    catalogued two-locus mechanism ABSTAINS rather than guessing (mirrors AMR R/S vs MIC-continuous).
  - v0 input = ALLELE CALLS (the wheel-only `--observed` pattern the HIV/fungal cells use). Genome-mode
    (detect the FRI-Col premature stop / FRI-Ler start-codon deletion from a sequence) needs VERIFIED variant
    coordinates + a reference — a v0.1 follow-on, deliberately NOT fabricated here.
  - Faithful-to-literature: this applies published functional-allele assignments; it is not a new model.

Pure-python, wheel-only, offline, deterministic.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# ---- curated causal-allele catalog (functional status per named natural allele) -------------------------
# status: "functional" (FRI: activates FLC / FLC: strong repressor) | "lof" (null) | "weak" (FLC only:
# reduced steady-state mRNA, not a null) | "unknown".
# Sources are named per entry; see the module docstring + wiki/deterministic_flowering_scoping_2026-07-16.md.

FRI_ALLELES: dict[str, dict] = {
    "Col": {"status": "lof", "lesion": "premature stop codon",
            "source": "Johanson et al. 2000 (FRI molecular analysis); the Col-0 reference is a rapid cycler"},
    "Ler": {"status": "lof", "lesion": "start-codon deletion",
            "source": "Johanson et al. 2000; NB the FRI-Ler allele is not a true null (it retains some FLC "
                      "induction) but the accession is early — see Ler's additional weak FLC"},
    "Cvi": {"status": "lof", "lesion": "first-intron substitution + in-frame stop",
            "source": "Gazzani et al. 2003 / Shindo et al. 2005 (substitution alleles, not the Col/Ler deletions)"},
    "Wil-2": {"status": "lof", "lesion": "first-intron substitution",
              "source": "Gazzani et al. 2003 / Shindo et al. 2005"},
    "Sf-2": {"status": "functional", "lesion": None,
             "source": "the functional FRI introgressed into Col to build the classic late-flowering Col-FRI "
                       "winter-annual line"},
    "H51": {"status": "functional", "lesion": None, "source": "functional-FRI winter-annual reference"},
}

FLC_ALLELES: dict[str, dict] = {
    "Col": {"status": "functional", "lesion": None,
            "source": "strong FLC; Col is early only because its FRI is LoF"},
    "Ler": {"status": "weak", "lesion": "naturally-occurring weak allele (reduced steady-state mRNA)",
            "source": "Michaels et al. 2003 PNAS (attenuation of FLC activity)"},
    "Van-0": {"status": "lof", "lesion": "nonsense mutation",
              "source": "Werner et al. 2005 (FRI-independent variation) — early despite functional FRI"},
    "Bur-0": {"status": "lof", "lesion": "aberrant splicing -> behaves as a null",
              "source": "Werner et al. 2005 (XAM mapping)"},
    "Da(1)-12": {"status": "weak", "lesion": "weak allele (reduced steady-state mRNA)",
                 "source": "Michaels et al. 2003 PNAS — functional FRI yet summer-annual"},
    "Shakhdara": {"status": "weak", "lesion": "weak allele",
                  "source": "Michaels et al. 2003 PNAS — functional FRI yet summer-annual"},
}

_FRI_OK = {"functional"}
_FLC_STRONG = {"functional"}
_FLC_REDUCED = {"weak", "lof"}
_VALID = {"functional", "lof", "weak", "unknown"}

# The mechanism this two-locus rule CANNOT see (named, not hidden) — the Lz-0 class.
UNSEEN_MECHANISMS = (
    "FRI-INDEPENDENT FLC upregulation (e.g. the Lz-0 accession: an FRI deletion yet LATE, via FRL1/FES1-class "
    "activators) — an FRI-LoF genome can still be late",
    "the polygenic + environment-dependent residue (~30-60% of long-day variation): photoperiod, temperature, "
    "vernalization history, FT/CO-pathway variation",
)


class FloweringInputError(ValueError):
    """Raised on an unknown allele name / invalid status (never a silent wrong call)."""


@dataclass
class FloweringCall:
    habit: str                # "summer_annual_early" | "winter_annual_late" | "ABSTAIN"
    vernalization_required: bool | None
    confidence: str           # "high" | "medium" | "low"
    fri_status: str
    flc_status: str
    rule: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "organism": "Arabidopsis_thaliana", "trait": "flowering_habit",
            "regime": "A_curated_catalog", "rule": self.rule,
            "habit": self.habit, "vernalization_required": self.vernalization_required,
            "confidence": self.confidence,
            "fri_status": self.fri_status, "flc_status": self.flc_status,
            "notes": self.notes,
            "scope_limit": ("PARTIAL: FRI/FLC explains ~40-70% of long-day flowering-time variation; this is a "
                            "HABIT/direction call, NOT quantitative days-to-flower"),
            "undetectable_mechanisms": list(UNSEEN_MECHANISMS),
        }


def status_for(locus: str, allele: str) -> str:
    """Resolve a named natural allele -> functional status. Unknown name -> FloweringInputError."""
    table = {"FRI": FRI_ALLELES, "FLC": FLC_ALLELES}.get(locus.upper())
    if table is None:
        raise FloweringInputError(f"unknown locus {locus!r}; expected FRI or FLC")
    if allele in table:
        return table[allele]["status"]
    if allele.lower() in _VALID:                  # allow passing a status directly
        return allele.lower()
    raise FloweringInputError(
        f"unknown {locus.upper()} allele {allele!r}; known: {sorted(table)} (or a status: {sorted(_VALID)})")


def call_flowering_habit(fri: str, flc: str) -> FloweringCall:
    """Deterministic flowering-habit call from FRI + FLC allele names (or statuses).

    late (winter-annual) iff FRI functional AND FLC strong; early (summer-annual) if EITHER is lost/weak;
    ABSTAIN if either locus is unknown (the polygenic residue is not guessed).
    """
    fri_s, flc_s = status_for("FRI", fri), status_for("FLC", flc)
    notes: list[str] = []

    if fri_s == "unknown" or flc_s == "unknown":
        return FloweringCall("ABSTAIN", None, "low", fri_s, flc_s, "arabidopsis_flowering_habit_v0",
                             ["an FRI/FLC locus is uncalled — the two-locus mechanism cannot be evaluated; "
                              "abstaining rather than guessing the polygenic residue"])

    if flc_s in _FLC_REDUCED:
        # FLC is DOWNSTREAM of FRI: a weak/null repressor gives the early habit regardless of FRI.
        if fri_s in _FRI_OK:
            notes.append("functional FRI but reduced FLC -> early via the FLC route (the Da(1)-12 / Shakhdara "
                         "class); a naive FRI-only rule would mis-call this LATE")
        habit, conf = "summer_annual_early", "high"
        notes.append(f"FLC {flc_s} -> the dosage-dependent repressor is lost/attenuated (downstream of FRI)")
    elif fri_s not in _FRI_OK:
        # FRI-LoF with a strong FLC: the classic summer-annual route, capped by the Lz-0 counterexample.
        habit, conf = "summer_annual_early", "medium"
        notes.append("FRI LoF with a strong FLC -> early via the FRI route (the Col/Ler class)")
        notes.append("confidence capped MEDIUM: FRI-INDEPENDENT FLC upregulation can still give a LATE habit "
                     "(the Lz-0 counterexample) — this two-locus rule cannot see it")
    else:
        habit, conf = "winter_annual_late", "high"
        notes.append("functional FRI + strong FLC -> high FLC -> the winter-annual habit (the Col-FRI[Sf-2] "
                     "class); vernalization required for rapid spring flowering")

    return FloweringCall(habit, habit == "winter_annual_late", conf, fri_s, flc_s,
                         "arabidopsis_flowering_habit_v0", notes)


def reference_integrity_ok() -> bool:
    """Biology contract guard — a corrupted catalog/rule fails this. Pins the four literature anchors,
    including the one a naive FRI-only rule gets WRONG (Da(1)-12: functional FRI yet summer-annual)."""
    col = call_flowering_habit("Col", "Col")             # FRI-LoF, strong FLC -> early
    col_fri = call_flowering_habit("Sf-2", "Col")        # functional FRI + strong FLC -> LATE (Col-FRI)
    ler = call_flowering_habit("Ler", "Ler")             # FRI-LoF + weak FLC -> early (double hit)
    da12 = call_flowering_habit("Sf-2", "Da(1)-12")      # functional FRI, weak FLC -> EARLY (the FLC route)
    van0 = call_flowering_habit("Sf-2", "Van-0")         # functional FRI, FLC null -> EARLY
    return (col.habit == "summer_annual_early"
            and col_fri.habit == "winter_annual_late" and col_fri.vernalization_required is True
            and ler.habit == "summer_annual_early"
            and da12.habit == "summer_annual_early"      # a naive FRI-only rule would say LATE -> WRONG
            and van0.habit == "summer_annual_early")
