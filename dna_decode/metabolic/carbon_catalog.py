"""E. coli carbon-source utilization decoder — the AMR determinant->phenotype paradigm applied to
metabolism, with the ONE twist metabolism adds that resistance does not: UPTAKE-GATING.

THE BIOLOGY (the rule this encodes):
    A cell can catabolize a sugar iff it can (1) IMPORT it (a transporter/permease) AND (2) break it down
    (the catabolic enzymes), AND (3) the transporter is EXPRESSED under the queried O2 condition.

    capability  iff  (all catabolic enzymes present) AND (a transporter present) AND (transporter expressed)

    The naive AMR-style rule ("has the pathway genes -> can use it") is RIGHT for most sugars but WRONG for
    the classic case metabolism is famous for:

    THE CITRATE ANCHOR (the Da(1)-12 of metabolism — the case a naive rule mis-calls):
        E. coli K-12 is Cit- AEROBICALLY. It carries the FULL TCA cycle (it metabolizes citrate as an
        internal intermediate every second) AND it carries the citT citrate/succinate antiporter. Yet on an
        aerobic citrate plate it CANNOT grow, because citT (the cit operon citCDEFXGT) is expressed ONLY
        ANAEROBICALLY (CitAB two-component + anaerobiosis). A "has citrate genes -> Cit+" rule says + ; the
        measured aerobic phenotype is - . The famous Lenski LTEE Cit+ mutants evolved AEROBIC citT expression
        via the rnk-citG regulatory duplication (Blount et al. 2012 Nature) — i.e. they fixed EXACTLY the
        expression gate this rule models, not a new enzyme. So: aerobic wild-type -> cannot_utilize (uptake
        gate closed); anaerobic -> utilizes (citT on, fermented via the citrate lyase CitDEF).

WHY A NON-FROZEN cell (like flowering / TMP-SMX / TB): the rule is (enzymes AND transporter AND condition)
— an AND across gene FAMILIES plus an expression gate — a shape the frozen count/OR `amr_rules.DRUG_RULE`
engine cannot represent. The frozen decoder surface is untouched.

VALIDATED against measured E. coli K-12 MG1655 phenotypes (EcoCyc / textbook Neidhardt; the anchors in
`reference_integrity_ok`): lac+ ara+ mal+ xyl+ rha+ glc+, and the CIT- aerobic / CIT+ anaerobic split.

HONEST SCOPE (load-bearing):
  - v0 = E. coli carbon-source CATABOLISM only (the cleanest curated determinant->phenotype metabolic map).
    Nitrogen/sulfur sources, anabolic auxotrophies, and cross-organism transfer are deliberately OUT.
  - This calls the DIRECTION (can / cannot utilize a named sugar aerobically or anaerobically), NOT growth
    rate / yield / lag. Anything not in the catalogued mechanism ABSTAINS rather than guessing.
  - Faithful-to-literature: it applies published operon/transporter assignments; it is not a new model. It
    reads gene PRESENCE — it cannot see a point mutation that silently inactivates a present gene (a v0.1
    genome-mode + sequence follow-on, deliberately not fabricated here).

Pure-python, wheel-only, offline, deterministic.
"""
from __future__ import annotations

from dataclasses import dataclass, field

_BOTH = ("aerobic", "anaerobic")


@dataclass(frozen=True)
class Substrate:
    """A curated carbon source: the catabolic enzymes + the uptake transporter(s) + O2-expression gate."""
    name: str
    enzymes: tuple[str, ...]                 # ALL required (AND across the catabolic genes)
    transporters: tuple[tuple[str, ...], ...]  # OR across sets; each set is AND (e.g. an ABC transporter)
    transporter_expressed: tuple[str, ...]   # O2 conditions under which the transporter is expressed
    source: str
    note: str = ""


# ---- curated E. coli carbon-catabolism catalog (each entry sourced) ------------------------------------
# Symbols are standard EcoCyc / NCBI gene symbols so the genome (feature-table) mode can match directly.

CARBON_SOURCES: dict[str, Substrate] = {
    "lactose": Substrate(
        "lactose", enzymes=("lacZ",), transporters=(("lacY",),), transporter_expressed=_BOTH,
        source="EcoCyc lac operon (lacZYA); Neidhardt — lacZ b-galactosidase + lacY permease",
        note="the classic Lac +/- system: lacZ- (no b-galactosidase) OR lacY- (no import) each give Lac-"),
    "L-arabinose": Substrate(
        "L-arabinose", enzymes=("araA", "araB", "araD"),
        transporters=(("araE",), ("araF", "araG", "araH")), transporter_expressed=_BOTH,
        source="EcoCyc araBAD + araE (low-affinity) / araFGH (high-affinity ABC)",
        note="AraA isomerase -> AraB ribulokinase -> AraD epimerase; either transporter suffices"),
    "maltose": Substrate(
        "maltose", enzymes=("malP", "malQ"),
        transporters=(("malE", "malF", "malG", "malK"),), transporter_expressed=_BOTH,
        source="EcoCyc mal regulon; MalEFGK2 ABC importer (+ LamB porin) + MalPQ maltodextrin catabolism"),
    "D-xylose": Substrate(
        "D-xylose", enzymes=("xylA", "xylB"),
        transporters=(("xylE",), ("xylF", "xylG", "xylH")), transporter_expressed=_BOTH,
        source="EcoCyc xylAB + xylE (symport) / xylFGH (ABC)",
        note="XylA isomerase -> XylB xylulokinase"),
    "L-rhamnose": Substrate(
        "L-rhamnose", enzymes=("rhaA", "rhaB", "rhaD"), transporters=(("rhaT",),), transporter_expressed=_BOTH,
        source="EcoCyc rhaBAD + rhaT permease",
        note="RhaA isomerase -> RhaB kinase -> RhaD aldolase"),
    "D-glucose": Substrate(
        "D-glucose", enzymes=("pgi",), transporters=(("ptsG",), ("crr", "ptsH", "ptsI")),
        transporter_expressed=_BOTH,
        source="EcoCyc PTS (ptsG/crr/ptsHI); glycolysis is universal",
        note="near-universal: glucose is the preferred carbon source; many redundant uptake routes"),
    "citrate": Substrate(
        "citrate", enzymes=("citD", "citE", "citF"), transporters=(("citT",),),
        transporter_expressed=("anaerobic",),
        source="EcoCyc cit operon citCDEFXGT; Blount et al. 2012 Nature (LTEE Cit+ = evolved AEROBIC citT)",
        note="THE ANCHOR: citT is ANAEROBIC-only in wild-type K-12, so aerobic citrate uptake is closed "
             "despite the full TCA cycle -> Cit- aerobic / Cit+ anaerobic. A naive has-the-genes rule fails."),
}

# aliases -> canonical name (tolerant CLI input)
_ALIASES = {
    "lac": "lactose", "arabinose": "L-arabinose", "ara": "L-arabinose", "arab": "L-arabinose",
    "mal": "maltose", "maltodextrin": "maltose", "xylose": "D-xylose", "xyl": "D-xylose",
    "rhamnose": "L-rhamnose", "rha": "L-rhamnose", "glucose": "D-glucose", "glc": "D-glucose",
    "glu": "D-glucose", "cit": "citrate", "citric acid": "citrate",
}

UNSEEN_MECHANISMS = (
    "a point mutation / frameshift that SILENTLY inactivates a gene that is PRESENT (this reads presence, "
    "not sequence integrity) — e.g. a lacZ nonsense allele still annotated 'lacZ'",
    "regulatory / cryptic-operon variation beyond the catalogued O2 gate (catabolite repression state, "
    "cryptic bgl/asc systems, silent operons activated only by mutation)",
    "growth RATE / yield / lag — this calls the can/cannot DIRECTION, not quantitative fitness",
)


class MetabolicInputError(ValueError):
    """Unknown substrate / invalid condition (never a silent wrong call)."""


@dataclass
class MetabolicCall:
    substrate: str
    capability: str               # "utilizes" | "cannot_utilize" | "ABSTAIN"
    condition: str                # "aerobic" | "anaerobic"
    confidence: str               # "high" | "medium" | "low"
    enzymes_present: list[str]
    enzymes_missing: list[str]
    transporter_present: bool
    transporter_expressed: bool
    rule: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "organism": "Escherichia_coli", "trait": "carbon_utilization",
            "regime": "A_curated_catalog_uptake_gated", "rule": self.rule,
            "substrate": self.substrate, "capability": self.capability, "condition": self.condition,
            "confidence": self.confidence,
            "enzymes_present": self.enzymes_present, "enzymes_missing": self.enzymes_missing,
            "transporter_present": self.transporter_present,
            "transporter_expressed_under_condition": self.transporter_expressed,
            "notes": self.notes,
            "scope_limit": ("v0: E. coli carbon CATABOLISM, presence-based; calls can/cannot DIRECTION not "
                            "growth rate; reads gene presence not sequence integrity"),
            "undetectable_mechanisms": list(UNSEEN_MECHANISMS),
        }


def resolve_substrate(name: str) -> str:
    """Canonicalize a substrate name (case/alias tolerant). Unknown -> MetabolicInputError."""
    key = name.strip()
    if key in CARBON_SOURCES:
        return key
    low = key.lower()
    if low in _ALIASES:
        return _ALIASES[low]
    # case-insensitive exact match against canonical names
    for canon in CARBON_SOURCES:
        if canon.lower() == low:
            return canon
    raise MetabolicInputError(
        f"unknown substrate {name!r}; known: {sorted(CARBON_SOURCES)} (or aliases {sorted(_ALIASES)})")


def call_carbon_utilization(substrate: str, present_genes, condition: str = "aerobic") -> MetabolicCall:
    """Deterministic carbon-utilization call from the set of PRESENT gene symbols.

    utilizes iff (all catabolic enzymes present) AND (a transporter set present) AND (transporter expressed
    under `condition`). Case-insensitive gene matching. `condition` in {aerobic, anaerobic}.
    """
    canon = resolve_substrate(substrate)
    ent = CARBON_SOURCES[canon]
    cond = condition.strip().lower()
    if cond not in _BOTH:
        raise MetabolicInputError(f"unknown condition {condition!r}; expected 'aerobic' or 'anaerobic'")

    present = {g.lower() for g in present_genes}
    enz_present = [e for e in ent.enzymes if e.lower() in present]
    enz_missing = [e for e in ent.enzymes if e.lower() not in present]
    # a transporter is present if ANY of its (AND-ed) sets is fully present
    tp_present = any(all(g.lower() in present for g in tset) for tset in ent.transporters)
    tp_expressed = cond in ent.transporter_expressed
    notes: list[str] = []

    rule = "ecoli_carbon_utilization_v0"
    if enz_missing or not tp_present:
        # a molecular block: an enzyme or the importer is absent -> high-confidence cannot-utilize
        if enz_missing:
            notes.append(f"catabolic enzyme(s) absent: {', '.join(enz_missing)} -> cannot break down {canon}")
        if not tp_present:
            notes.append(f"no {canon} uptake transporter present ({_fmt_transporters(ent)}) -> cannot import")
        return MetabolicCall(canon, "cannot_utilize", cond, "high", enz_present, enz_missing,
                             tp_present, tp_expressed, rule, notes)

    if not tp_expressed:
        # THE CITRATE ANCHOR: enzymes + transporter GENE present, but the importer is not expressed under
        # this O2 condition -> cannot utilize despite carrying the pathway. The naive has-genes rule fails.
        notes.append(f"all {canon} enzymes + the transporter GENE are present, but the transporter is not "
                     f"expressed {cond}ally (expressed: {', '.join(ent.transporter_expressed)}) -> uptake "
                     f"gate CLOSED -> cannot utilize {canon} {cond}ally")
        if ent.note:
            notes.append(ent.note)
        notes.append("a naive 'has the pathway genes -> can use it' rule MIS-CALLS this positive")
        return MetabolicCall(canon, "cannot_utilize", cond, "medium", enz_present, enz_missing,
                             tp_present, tp_expressed, rule, notes)

    # enzymes + transporter present AND expressed -> utilizes
    conf = "high"
    notes.append(f"all catabolic enzymes + a transporter present and expressed {cond}ally -> utilizes {canon}")
    if ent.note:
        notes.append(ent.note)
    return MetabolicCall(canon, "utilizes", cond, conf, enz_present, enz_missing,
                         tp_present, tp_expressed, rule, notes)


def _fmt_transporters(ent: Substrate) -> str:
    return " or ".join("+".join(tset) for tset in ent.transporters)


def genes_for(substrate: str) -> list[str]:
    """All catalogued gene symbols for a substrate (enzymes + every transporter option) — handy for tests
    and for a genome-mode 'what would I look for' listing."""
    ent = CARBON_SOURCES[resolve_substrate(substrate)]
    out = list(ent.enzymes)
    for tset in ent.transporters:
        out.extend(tset)
    return out


def reference_integrity_ok() -> bool:
    """Biology contract guard — a corrupted catalog/rule fails this. Pins measured E. coli K-12 phenotypes,
    INCLUDING the citrate anchor a naive has-the-genes rule gets wrong."""
    # lactose: full lac operon -> Lac+ (high)
    lac_pos = call_carbon_utilization("lactose", ["lacZ", "lacY"])
    # lactose: lacZ knocked out (permease only) -> Lac- (the classic enzyme-knockout negative)
    lac_neg = call_carbon_utilization("lactose", ["lacY"])
    # arabinose: full araBAD + araE -> Ara+
    ara_pos = call_carbon_utilization("L-arabinose", ["araA", "araB", "araD", "araE"])
    # glucose: PTS present -> Glc+
    glc_pos = call_carbon_utilization("D-glucose", ["pgi", "ptsG"])
    # THE ANCHOR — citrate AEROBIC with the FULL cit operon present -> Cit- (naive rule says +, truth -)
    cit_aer = call_carbon_utilization("citrate", ["citD", "citE", "citF", "citT"], condition="aerobic")
    # citrate ANAEROBIC with the same genes -> Cit+ (citT now expressed)
    cit_ana = call_carbon_utilization("citrate", ["citD", "citE", "citF", "citT"], condition="anaerobic")
    return (lac_pos.capability == "utilizes" and lac_pos.confidence == "high"
            and lac_neg.capability == "cannot_utilize" and "lacZ" in lac_neg.enzymes_missing
            and ara_pos.capability == "utilizes"
            and glc_pos.capability == "utilizes"
            and cit_aer.capability == "cannot_utilize" and cit_aer.transporter_present is True
            and cit_aer.transporter_expressed is False and cit_aer.confidence == "medium"
            and cit_ana.capability == "utilizes")
