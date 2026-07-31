"""Horse (Equus caballus) coat-colour decoder — the curated-catalog + EPISTASIS paradigm applied to the
best-characterised animal coat-colour system (the sibling of dog_coat.py, a different organism's visible trait).

THE BIOLOGY (the rule this encodes — OMIA-curated causal variants; the standard equine colour model):

  Coat colour is set by a small number of well-characterised loci acting in a FIXED order. Two pigments:
  EUMELANIN (black) and PHAEOMELANIN (red). The loci:

    E  (MC1R)     pigment-type SWITCH. e/e -> NO black pigment -> CHESTNUT (red), RECESSIVE-EPISTATIC over A
                  (a chestnut horse shows red even if genetically A/A bay). Causal e: MC1R p.Ser83Phe
                  (Marklund 1996, C901T). OMIA 001199-9796.
    A  (ASIP)     black DISTRIBUTION, expressed ONLY if E- (black-capable). A_ = BAY (black to points),
                  a/a = uniform BLACK. Causal a: ASIP 11-bp deletion (Rieder 2001). OMIA 001249... .
    CR (SLC45A2)  CREAM dilution, INCOMPLETE-DOMINANT (dose matters). 1 copy dilutes red -> palomino
                  (chestnut) / buckskin (bay) / smoky-black (black); 2 copies dilute both -> cremello /
                  perlino / smoky-cream (blue eyes). Causal: SLC45A2 c.457G>A p.Asp153Asn (Mariat 2003).
                  OMIA 001344-9796.
    D  (TBX3)     DUN dilution + primitive markings (dorsal stripe). 3-allele series D > nd1 > nd2: D =
                  dun-diluted; nd1/nd2 = non-dun (nd2 = 1617-bp deletion, no primitive markings). Dun on
                  chestnut = red dun, on bay = (bay) dun, on black = grullo/grulla. Imsland 2016 Nat Genet.
    G  (STX17)    GREY, autosomal DOMINANT + progressive: a G_ horse is BORN its base colour and greys to
                  white with age -> EPISTATIC for the ADULT colour (a bay G/n horse greys out). Causal:
                  4.6-kb intron-6 duplication (Rosengren Pielberg 2008). OMIA 001356-9796.

  Resolution order: BASE (E x A) -> dilutions (CREAM, DUN) -> GREY (epistatic override for the adult coat).

  THE EPISTASIS ANCHORS (the cases a naive has-the-allele rule mis-calls):
    (1) e/e + A/A -> CHESTNUT, not bay (E recessive-epistatic; no black for A to distribute).
    (2) G/n + any base -> GREY as an adult (born the base colour); a base-only rule mis-reports the base.
  `reference_integrity_ok` pins exactly these.

WHY A NON-FROZEN cell (like dog_coat / flowering / TB): a fixed-order epistasis across loci with a
dose-dependent dilution + a dominant progressive override — not the frozen count/OR amr_rules shape. Imports
nothing from the frozen decoder surface.

HONEST SCOPE (load-bearing):
  - v0 = the FIVE named loci (E/A/CR/D/G). Other dilutions/patterns (Champagne SLC36A1, Silver PMEL, Pearl,
    Roan, Tobiano/Overo/Sabino white-spotting, Appaloosa LP, Flaxen, sooty) are DELIBERATELY OUT; a horse
    carrying them ABSTAINS on the affected axis (via --present) rather than a confident wrong call.
  - Input is PER-LOCUS allele CALLS (the VGL/Etalon report shape), e.g. E=E/e,A=A/a,CR=Cr/N,D=nd1/nd1,G=n/n.
    A genome/VCF caller for the causal variants is a v0.1 follow-on.
  - Calls the COLOUR (base + cream/dun dilution + grey), NOT shade (sooty/flaxen), spotting extent, or the
    born-vs-current stage of a grey (it flags "greying with age"). Anything uncatalogued ABSTAINS.
  - Faithful-to-literature: applies the published OMIA locus/allele assignments; it is not a new model.
  - KNOWLEDGE_BASELINE: no free per-individual validation substrate (unlike the dog cell's Darwin's Ark).

Pure-python, wheel-only, offline, deterministic. Regime-A curated catalog. Companion/livestock visible-trait
genetics — NOT any human/forensic application.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from dna_decode.data.horse_coat import INDETERMINATE, call_horse_base_colour

RULES_VERSION = "horse-coat-colour-v0.1.0"


@dataclass(frozen=True)
class Locus:
    name: str
    gene: str
    alleles: tuple[str, ...]        # dominance order, most-dominant first
    source: str
    note: str = ""


LOCI: dict[str, Locus] = {
    "E": Locus("E", "MC1R", ("E", "e"),
               source="OMIA 001199-9796 (MC1R p.Ser83Phe = recessive red `e`, Marklund 1996 C901T)",
               note="pigment-type switch: e/e -> chestnut (red), RECESSIVE-EPISTATIC over A"),
    "A": Locus("A", "ASIP", ("A", "a"),
               source="OMIA (Agouti/ASIP; a = 11-bp deletion, Rieder 2001)",
               note="black distribution, expressed only if E-: A_ = bay, a/a = uniform black"),
    "CR": Locus("CR", "SLC45A2", ("Cr", "N"),
                source="OMIA 001344-9796 (SLC45A2/MATP c.457G>A p.Asp153Asn, Mariat 2003)",
                note="cream dilution, INCOMPLETE-DOMINANT: 1 copy dilutes red, 2 copies dilute red+black"),
    "D": Locus("D", "TBX3", ("D", "nd1", "nd2"),
               source="OMIA (TBX3 regulatory; Imsland 2016 Nat Genet; nd2 = 1617-bp deletion)",
               note="dun dilution + primitive markings; D > nd1 > nd2 (nd1/nd2 non-dun)"),
    "G": Locus("G", "STX17", ("G", "n"),
               source="OMIA 001356-9796 (STX17 4.6-kb intron-6 duplication, Rosengren Pielberg 2008)",
               note="grey, DOMINANT + progressive: born base colour, greys to white with age (EPISTATIC adult)"),
}

# dilution / pattern loci NOT modelled in v0 -> the affected appearance axis ABSTAINS when declared present.
UNMODELLED_LOCI = {
    "CH": "champagne (SLC36A1) — dilution not modelled in v0",
    "Z": "silver (PMEL) — dilutes black to chocolate/flaxen mane, not modelled",
    "PRL": "pearl (SLC45A2 prl) — recessive dilution, interacts with cream; not modelled in v0",
    "RN": "roan — white hairs intermixed, not a base dilution",
    "TO": "tobiano / overo / sabino white-spotting — changes how much coat is coloured",
    "LP": "appaloosa (LP/PATN) leopard-complex spotting",
    "STY": "sooty / flaxen — shade modifiers, v0 calls the base colour not the shade",
}

_ALLELE_ALIASES = {
    "cr": "Cr", "n": "N", "nd1": "nd1", "nd2": "nd2", "g": "G",
}
_CREAM = "Cr"
_GREY = "G"


class HorseInputError(ValueError):
    """Unknown locus/allele or malformed genotype (never a silent wrong call)."""


@dataclass
class HorseCoatCall:
    coat_color: str                 # human label, e.g. "bay", "palomino", "grullo (black dun)", "grey (born bay)"
    base_color: str                 # chestnut | bay | black
    pigment_type: str               # "eumelanin+phaeomelanin" | "phaeomelanin"
    dilutions: list[str]            # e.g. ["cream x1", "dun"]
    greying: bool                   # True if G_ -> greys with age
    confidence: str                 # high | medium | low
    abstains_on: list[str]
    per_locus: dict
    rule: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "organism": "Equus_caballus", "trait": "coat_colour",
            "regime": "A_curated_catalog_epistatic", "rule": self.rule,
            "coat_color": self.coat_color, "base_color": self.base_color,
            "pigment_type": self.pigment_type, "dilutions": self.dilutions, "greying_with_age": self.greying,
            "confidence": self.confidence, "abstains_on": self.abstains_on, "per_locus": self.per_locus,
            "notes": self.notes,
            "scope_limit": ("v0: five loci E/A/CR/D/G; calls base + cream/dun dilution + grey; NOT sooty/flaxen "
                            "shade, spotting extent, champagne/silver/pearl, or the current stage of greying"),
            "evidence_tier": "knowledge_baseline (curated OMIA catalog; no free per-individual validation substrate)",
        }


def _canon_allele(locus: str, tok: str) -> str:
    t = tok.strip()
    if t in LOCI[locus].alleles:
        return t
    low = t.lower()
    if low in _ALLELE_ALIASES and _ALLELE_ALIASES[low] in LOCI[locus].alleles:
        return _ALLELE_ALIASES[low]
    raise HorseInputError(
        f"unknown {locus}-locus allele {tok!r}; recognised: {list(LOCI[locus].alleles)}")


def parse_genotype(locus: str, spec: str) -> tuple[str, str]:
    if locus not in LOCI:
        raise HorseInputError(f"unknown locus {locus!r}; v0 loci: {list(LOCI)}")
    parts = [p for p in spec.replace("|", "/").split("/") if p.strip()]
    if len(parts) != 2:
        raise HorseInputError(f"{locus} genotype {spec!r} is not diploid (expected a1/a2, e.g. E/e)")
    return _canon_allele(locus, parts[0]), _canon_allele(locus, parts[1])


def _has_dominant(geno: tuple[str, str], allele: str) -> bool:
    return allele in geno


def call_horse_coat(loci_genotypes: dict[str, str], present_loci: list[str] | None = None) -> HorseCoatCall:
    """Deterministic horse coat-colour call from per-locus allele genotypes.

    `loci_genotypes`: {locus -> 'a1/a2'} for any subset of E/A/CR/D/G (E required — the pigment switch;
    A needed to resolve bay-vs-black when E-). `present_loci`: v0-unmodelled loci present on the horse
    (e.g. ['Z','TO']) -> the affected appearance axis is added to abstains_on.
    """
    rule = f"horse_coat_colour_epistatic_v0 ({RULES_VERSION})"
    notes: list[str] = []
    per: dict = {}
    geno: dict[str, tuple[str, str]] = {}
    for loc, spec in loci_genotypes.items():
        L = loc.strip().upper()
        if L not in LOCI:
            if L in UNMODELLED_LOCI:
                raise HorseInputError(
                    f"locus {loc!r} is a v0-unmodelled locus ({UNMODELLED_LOCI[L]}); pass it via "
                    f"present_loci=[...] so the affected axis ABSTAINS instead of a wrong call")
            raise HorseInputError(f"unknown locus {loc!r}; v0 loci: {list(LOCI)}")
        geno[L] = parse_genotype(L, spec)
        per[L] = "/".join(geno[L])

    abstains: list[str] = []
    for pl in (present_loci or []):
        P = pl.strip().upper()
        if P in UNMODELLED_LOCI:
            abstains.append(f"{P}: {UNMODELLED_LOCI[P]}")
        elif P not in LOCI:
            raise HorseInputError(f"unknown present locus {pl!r}")

    if "E" not in geno:
        raise HorseInputError("E (MC1R) genotype is required — it is the pigment-type switch (chestnut vs "
                              "black-capable) at the top of the epistasis")

    # ---- 1. BASE colour (E x A) — REUSE the deployed base rule (dna_decode.data.horse_coat, Rieder 2001 /
    #         VGL), which this cell EXTENDS with cream/dun/grey. Don't duplicate the base epistasis. ----
    base_raw = call_horse_base_colour("".join(geno["E"]), "".join(geno["A"]) if "A" in geno else "")
    if base_raw == "chestnut":
        base = "chestnut"
        pigment = "phaeomelanin"
        if "A" in geno and _has_dominant(geno["A"], "A"):
            notes.append("A/_ present but coat is chestnut (e/e is RECESSIVE-EPISTATIC — no black pigment for "
                         "Agouti to distribute); a naive 'A -> bay' rule would MIS-CALL this horse")
    elif base_raw == INDETERMINATE:      # E- (E is a valid 2-allele call) but A absent/uncallable
        base = "bay_or_black_unknown"
        pigment = "eumelanin+phaeomelanin"
        notes.append("E- (black-capable) but A (ASIP) genotype absent -> bay-vs-black UNKNOWN; pass A")
        abstains.append("A: bay-vs-black not resolvable without the Agouti genotype")
    else:
        base = base_raw                  # "bay" | "black"
        pigment = "eumelanin+phaeomelanin"

    # ---- 2. dilutions: CREAM (dose) + DUN ----
    dilutions: list[str] = []
    cream = geno["CR"].count(_CREAM) if "CR" in geno else 0
    if cream:
        dilutions.append(f"cream x{cream}")
    dun = "D" in geno and _has_dominant(geno["D"], "D")
    if dun:
        dilutions.append("dun")

    color = _compose_base_with_dilution(base, cream, dun, notes)

    # ---- 3. GREY (epistatic for the ADULT coat) ----
    greying = "G" in geno and _has_dominant(geno["G"], _GREY)
    if greying:
        born = color
        color = f"grey (born {born}; greys to white with age)"
        notes.append("G/_ -> GREY: born the base colour, progressively whitens with age (EPISTATIC for the "
                     "adult coat); a base-only rule would report the born colour and miss the greying")

    # ---- confidence ----
    conf = "high"
    if base == "bay_or_black_unknown":
        conf = "medium"
    if abstains:
        conf = "medium" if conf == "high" else conf
    return HorseCoatCall(color, base if not base.endswith("unknown") else "undetermined", pigment,
                         dilutions, greying, conf, abstains, per, rule, notes)


def _compose_base_with_dilution(base: str, cream: int, dun: bool, notes: list) -> str:
    """Compose the human colour label from base + cream dose + dun (pre-grey)."""
    if base == "bay_or_black_unknown":
        core = "black-based (bay or black — A needed)"
    elif cream >= 2:
        core = {"chestnut": "cremello", "bay": "perlino", "black": "smoky cream"}[base] + " (double cream, blue-eyed)"
    elif cream == 1:
        core = {"chestnut": "palomino", "bay": "buckskin", "black": "smoky black"}[base]
    else:
        core = base
    if dun:
        if base == "chestnut" and cream == 0:
            core = "red dun"
        elif base == "black" and cream == 0:
            core = "grullo/grulla (black dun)"
        elif base == "bay" and cream == 0:
            core = "bay dun (classic dun)"
        else:
            core = f"{core} + dun (primitive markings)"
        notes.append("dun adds a dorsal stripe + primitive markings")
    return core


def reference_integrity_ok() -> bool:
    """Biology contract guard — pins known horse genotypes -> colours, INCLUDING the two epistasis anchors a
    naive has-the-allele rule gets wrong (e/e hides Agouti; Grey is epistatic for the adult coat)."""
    chestnut = call_horse_coat({"E": "e/e", "A": "A/A"})              # ANCHOR 1: e/e hides A/A -> chestnut
    bay = call_horse_coat({"E": "E/e", "A": "A/a"})
    black = call_horse_coat({"E": "E/E", "A": "a/a"})
    palomino = call_horse_coat({"E": "e/e", "A": "a/a", "CR": "Cr/N"})
    buckskin = call_horse_coat({"E": "E/e", "A": "A/a", "CR": "Cr/N"})
    cremello = call_horse_coat({"E": "e/e", "A": "a/a", "CR": "Cr/Cr"})
    grullo = call_horse_coat({"E": "E/E", "A": "a/a", "D": "D/nd1"})
    grey = call_horse_coat({"E": "E/e", "A": "A/a", "G": "G/n"})      # ANCHOR 2: grey epistatic (born bay)
    return (chestnut.base_color == "chestnut" and chestnut.pigment_type == "phaeomelanin"
            and bay.coat_color == "bay" and black.coat_color == "black"
            and palomino.coat_color == "palomino" and buckskin.coat_color == "buckskin"
            and cremello.coat_color.startswith("cremello")
            and grullo.coat_color.startswith("grullo")
            and grey.greying is True and grey.coat_color.startswith("grey (born bay")
            and grey.base_color == "bay")
