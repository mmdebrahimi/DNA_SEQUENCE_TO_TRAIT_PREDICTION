"""Chicken (Gallus gallus) plumage-colour decoder — the curated-catalog + EPISTASIS paradigm applied to a
BIRD (the sibling of dog/horse/cat coat; a 4th-organism cell). Its signature is Z-LINKED sex-linkage with
REVERSED hemizygosity: birds are ZW, so the FEMALE (ZW) is hemizygous — the mirror of cat's X-linked orange.

THE BIOLOGY (OMIA-curated causal variants; the standard poultry plumage model):

    E  (MC1R)     EXTENSION — the eumelanin "canvas", an allelic series by DECREASING black (Kabir 2020;
                  Ling 2003). E (extended black) > E^R (birchen) > E^Wh (dominant wheaten) > e+ (wild-type,
                  black-breasted red / partridge) > e^b (brown) > e^y (recessive wheaten). Causal E: MC1R
                  G274A (E92K) +/- M71T; buttercup H215P; rec-wheaten R213C. OMIA 000374-9031. Extended black
                  -> mostly black; wheaten -> mostly red/buff (little eumelanin for the pattern genes to act on).
    B  (CDKN2A)   SEX-LINKED BARRING (Z-linked): white bars across the eumelanin. B (barred, B1=V9D) vs b+.
                  Z-LINKED -> a MALE (ZZ) is homo/heterozygous, a FEMALE (ZW) is HEMIZYGOUS (B/W). Incomplete
                  dosage compensation. Hellstrom 2010 / Schwochow 2017. OMIA 000102-9031.
    S  (SLC45A2)  SILVER/GOLD (Z-linked): silver S restricts pheomelanin (gold -> silver/white); s+ = gold.
                  Two missense Y277C/L347M. Gunnarsson 2007. OMIA 000370-9031. Same gene as cat cream / horse cream.
    I  (PMEL17)   DOMINANT WHITE (incompletely dominant): inhibits black eumelanin (-> white; red/gold shows).
                  9-bp exon-10 insertion (Dun/Smoky are alleles). White Leghorn. Kerje 2004. OMIA 000373-9031.
    Bl (PMEL-locus) BLUE (Andalusian) — INCOMPLETELY DOMINANT eumelanin dilution: bl+/bl+ black, Bl/bl+ blue,
                  Bl/Bl splash.
    lav (MLPH)    LAVENDER / self-blue — RECESSIVE further dilution of black (lav/lav). Same gene as the mammal
                  dilute loci.
    c  (TYR)      RECESSIVE WHITE — c/c = white (tyrosinase; distinct from dominant white). Recessive.

  Resolution order: recessive-white / dominant-white (mask) -> EXTENSION canvas -> silver/gold (pheomelanin) ->
  blue/lavender (eumelanin dilution) -> barring (bars the eumelanin).

  THE EPISTASIS ANCHORS (the cases a naive has-the-allele rule mis-calls):
    (1) EXTENSION is the canvas: a wheaten (e^Wh/e^Wh) bird is mostly red/buff, so barring / blue (which act
        ON eumelanin) barely show — a naive "B -> barred" rule mis-calls a wheaten bird.
    (2) Z-LINKED barring/silver with REVERSED hemizygosity: a FEMALE (ZW) is hemizygous (one Z allele, e.g.
        B/W); a naive autosomal or mammal-style (male-hemizygous) rule mis-calls the sex-dependent zygosity.
    (3) Dominant white (I) + recessive white (c/c) INHIBIT/mask eumelanin -> white.
  `reference_integrity_ok` pins exactly these.

WHY A NON-FROZEN cell (like dog/horse/cat coat): fixed-order epistasis + Z-LINKED sex-dependent loci — not the
frozen count/OR amr_rules shape. Imports nothing from the frozen decoder surface.

HONEST SCOPE (load-bearing):
  - v0 = the loci the user named + the standard set: E/B/S/I/Bl/lav/c. Pattern loci (Columbian Co, mottling mo,
    spangling, pencilling, mille-fleur), the many E-locus sub-alleles' fine pattern, and comb/feather-structure
    genes ABSTAIN (via --present).
  - Input is PER-LOCUS allele CALLS. The Z-linked B and S loci take ONE allele for a FEMALE (ZW hemizygous, e.g.
    B=B) or TWO for a MALE (ZZ, e.g. B=B/b+); --sex may be given but is inferred from the Z-locus zygosity.
  - Calls the eumelanin CANVAS + major modifiers (barred/silver/blue/lavender/white), NOT the fine feather
    pattern, lacing, or shade. KNOWLEDGE_BASELINE: no free per-individual validation substrate.
  - Faithful-to-literature: applies the published OMIA locus/allele assignments; it is not a new model.

Pure-python, wheel-only, offline, deterministic. Benign livestock visible-trait genetics — NOT human/forensic.
"""
from __future__ import annotations

from dataclasses import dataclass, field

RULES_VERSION = "chicken-plumage-v0.1.0"


@dataclass(frozen=True)
class Locus:
    name: str
    gene: str
    alleles: tuple[str, ...]        # dominance order, most-dominant first
    source: str
    z_linked: bool = False
    note: str = ""


LOCI: dict[str, Locus] = {
    "E": Locus("E", "MC1R", ("E", "ER", "EWh", "e+", "eb", "ey"),
               source="OMIA 000374-9031 (MC1R; E extended black G274A/E92K; buttercup H215P; rec-wheaten R213C)",
               note="EXTENSION canvas by decreasing eumelanin: E extended-black > ER birchen > EWh wheaten > "
                    "e+ wild (partridge) > eb brown > ey rec-wheaten"),
    "B": Locus("B", "CDKN2A", ("B", "b+"),
               source="OMIA 000102-9031 (CDKN2A; B1 V9D barring); Hellstrom 2010 / Schwochow 2017",
               z_linked=True,
               note="Z-LINKED barring (white bars across eumelanin); FEMALE (ZW) hemizygous, MALE (ZZ) 2 alleles"),
    "S": Locus("S", "SLC45A2", ("S", "s+"),
               source="OMIA 000370-9031 (SLC45A2/MATP; Y277C/L347M silver); Gunnarsson 2007",
               z_linked=True,
               note="Z-LINKED silver: S restricts pheomelanin (gold -> silver/white); s+ = gold"),
    "I": Locus("I", "PMEL17", ("I", "i+"),
               source="OMIA 000373-9031 (PMEL17; 9-bp exon-10 insertion); Kerje 2004",
               note="DOMINANT WHITE (incompletely dominant): inhibits black eumelanin -> white (red/gold shows)"),
    "BL": Locus("BL", "PMEL-locus", ("Bl", "bl+"),
                source="Blue/Andalusian (incompletely dominant eumelanin dilution)",
                note="Bl/bl+ = blue, Bl/Bl = splash, bl+/bl+ = black (undiluted)"),
    "LAV": Locus("LAV", "MLPH", ("Lav", "lav"),
                 source="OMIA (MLPH; lavender/self-blue) — same gene as mammal dilute",
                 note="lav/lav = lavender (recessive further dilution of black)"),
    "C": Locus("C", "TYR", ("C+", "c"),
               source="OMIA (TYR; recessive white/albino series)",
               note="c/c = recessive white (distinct from dominant white)"),
}

UNMODELLED_LOCI = {
    "CO": "Columbian (restricts black to points) — pattern locus, not modelled in v0",
    "MO": "mottling (EDNRB2) — white-tipped feathers",
    "PG": "pattern/pencilling/lacing (Pg/Ml) — fine feather pattern",
    "DB": "dark brown (Db) — pheomelanin pattern modifier",
    "SP": "spangling / mille-fleur pattern",
}

_ALLELE_ALIASES = {"e+": "e+", "b+": "b+", "s+": "s+", "i+": "i+", "bl+": "bl+", "c+": "C+", "lav": "lav"}


class ChickenInputError(ValueError):
    """Unknown locus/allele or malformed genotype (never a silent wrong call)."""


@dataclass
class PlumageCall:
    plumage: str                    # composed human label
    eumelanin_canvas: str           # from E locus (extended_black / wheaten / wild_partridge / ...)
    sex_basis: str
    barred: bool
    silver: bool                    # pheomelanin restricted (gold -> silver)
    dilution: str | None            # "blue" | "splash" | "lavender" | None
    white_type: str | None          # "dominant_white" | "recessive_white" | None
    is_white_masked: bool
    confidence: str
    abstains_on: list[str]
    per_locus: dict
    rule: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "organism": "Gallus_gallus", "trait": "plumage_colour",
            "regime": "A_curated_catalog_epistatic_zlinked", "rule": self.rule,
            "plumage": self.plumage, "eumelanin_canvas": self.eumelanin_canvas, "sex_basis": self.sex_basis,
            "barred": self.barred, "silver": self.silver, "dilution": self.dilution,
            "white_type": self.white_type, "is_white_masked": self.is_white_masked,
            "confidence": self.confidence, "abstains_on": self.abstains_on, "per_locus": self.per_locus,
            "notes": self.notes,
            "scope_limit": ("v0: E/B/S/I/Bl/lav/c; calls eumelanin canvas + major modifiers (barred/silver/"
                            "blue/lavender/white), NOT fine feather pattern, lacing, pencilling, or shade"),
            "evidence_tier": "knowledge_baseline (curated OMIA catalog; no free per-individual validation substrate)",
        }


def _canon_allele(locus: str, tok: str) -> str:
    t = tok.strip()
    if t in LOCI[locus].alleles:
        return t
    if t in _ALLELE_ALIASES and _ALLELE_ALIASES[t] in LOCI[locus].alleles:
        return _ALLELE_ALIASES[t]
    raise ChickenInputError(f"unknown {locus}-locus allele {tok!r}; recognised: {list(LOCI[locus].alleles)}")


def _parse(locus: str, spec: str, z_linked: bool = False) -> tuple[str, ...]:
    parts = [p for p in spec.replace("|", "/").split("/") if p.strip()]
    if z_linked and len(parts) == 1:
        return (_canon_allele(locus, parts[0]),)
    if len(parts) != 2:
        exp = "1 (female ZW) or 2 (male ZZ) alleles" if z_linked else "a1/a2"
        raise ChickenInputError(f"{locus} genotype {spec!r} must be {exp}")
    return _canon_allele(locus, parts[0]), _canon_allele(locus, parts[1])


def _dominant(locus: str, geno: tuple[str, ...]) -> str:
    order = LOCI[locus].alleles
    return min(geno, key=lambda a: order.index(a))


_CANVAS = {
    "E": ("extended_black", "solid black"),
    "ER": ("birchen", "birchen (black with silver/gold hackles)"),
    "EWh": ("wheaten", "wheaten (mostly red/buff; little eumelanin)"),
    "e+": ("wild_partridge", "wild-type black-breasted red / partridge"),
    "eb": ("brown", "brown / partridge-brown"),
    "ey": ("recessive_wheaten", "recessive wheaten (pale buff)"),
}


def call_chicken_plumage(loci_genotypes: dict[str, str], sex: str | None = None,
                         present_loci: list[str] | None = None) -> PlumageCall:
    """Deterministic chicken plumage-colour call from per-locus allele genotypes.

    `loci_genotypes`: {locus -> 'a1/a2'} for any subset of E/B/S/I/Bl/lav/c. The Z-LINKED B and S loci take
    ONE allele for a FEMALE (ZW hemizygous, e.g. 'B') or TWO for a MALE (ZZ, e.g. 'B/b+'). `sex`
    ('male'/'female') optional — inferred from the Z-locus zygosity when omitted (REVERSED from mammals:
    1 allele -> female). `present_loci`: v0-unmodelled loci present -> the affected axis ABSTAINS.
    """
    rule = f"chicken_plumage_epistatic_v0 ({RULES_VERSION})"
    notes: list[str] = []
    per: dict = {}
    geno: dict[str, tuple[str, ...]] = {}
    for loc, spec in loci_genotypes.items():
        L = loc.strip().upper()
        if L not in LOCI:
            if L in UNMODELLED_LOCI:
                raise ChickenInputError(f"locus {loc!r} is a v0-unmodelled locus ({UNMODELLED_LOCI[L]}); pass "
                                        f"it via present_loci=[...] so the affected axis ABSTAINS")
            raise ChickenInputError(f"unknown locus {loc!r}; v0 loci: {list(LOCI)}")
        geno[L] = _parse(L, spec, z_linked=LOCI[L].z_linked)
        per[L] = "/".join(geno[L])

    abstains: list[str] = []
    for pl in (present_loci or []):
        P = pl.strip().upper()
        if P in UNMODELLED_LOCI:
            abstains.append(f"{P}: {UNMODELLED_LOCI[P]}")
        elif P not in LOCI:
            raise ChickenInputError(f"unknown present locus {pl!r}")

    # ---- sex basis (Z-linked: FEMALE hemizygous -> 1 allele; REVERSED from mammals) ----
    z_geno = geno.get("B") or geno.get("S")
    if sex:
        sex = sex.strip().lower()
        if sex not in ("male", "female"):
            raise ChickenInputError("--sex must be 'male' or 'female'")
        sex_basis = f"{sex} ({'ZW hemizygous' if sex == 'female' else 'ZZ'})"
    elif z_geno is not None:
        sex = "female" if len(z_geno) == 1 else "male"
        sex_basis = f"{sex} (inferred from Z-locus zygosity: {len(z_geno)} allele{'s' if len(z_geno) == 2 else ''}; birds are ZW)"
    else:
        sex = None
        sex_basis = "unspecified (no Z-linked B/S genotype / --sex)"

    # ---- 1. recessive white (c/c) ----
    if "C" in geno and geno["C"] == ("c", "c"):
        notes.append("c/c -> recessive white (TYR); masks colour")
        return PlumageCall("white (recessive white)", "masked", sex_basis, False, False, None,
                           "recessive_white", True, "high", abstains, per, rule, notes)

    # ---- 2. dominant white (I_) ----
    dom_white = "I" in geno and "I" in geno["I"]

    # ---- 3. Extension canvas ----
    if "E" in geno:
        edom = _dominant("E", geno["E"])
        canvas, canvas_label = _CANVAS[edom]
    else:
        canvas, canvas_label = ("wild_partridge", "wild-type (E absent -> assumed wild partridge)")
        notes.append("E (MC1R) absent -> canvas assumed wild-type partridge; pass E to resolve")

    low_eumelanin = canvas in ("wheaten", "recessive_wheaten")

    # ---- 4. silver/gold (Z-linked) ----
    silver = "S" in geno and "S" in geno["S"]

    # ---- 5. blue / lavender dilution of eumelanin ----
    dilution = None
    if "BL" in geno:
        n_bl = geno["BL"].count("Bl")
        if n_bl == 1:
            dilution = "blue"
        elif n_bl == 2:
            dilution = "splash"
    if "LAV" in geno and geno["LAV"] == ("lav", "lav"):
        dilution = "lavender" if dilution is None else f"{dilution}+lavender"

    # ---- 6. barring (Z-linked) ----
    barred = "B" in geno and "B" in geno["B"]
    if barred and low_eumelanin:
        notes.append("barring acts on EUMELANIN, but the Extension canvas is wheaten (little black) -> bars "
                     "barely show; a naive 'B -> barred' rule mis-calls this (Extension is the canvas)")

    # ---- compose ----
    if dom_white:
        eumelanin_word = "white (dominant-white-inhibited; red/gold shows through)"
    elif dilution == "splash":
        eumelanin_word = "splash (white with blue flecks)"
    elif dilution == "blue":
        eumelanin_word = "blue"
    elif dilution == "lavender":
        eumelanin_word = "lavender (self-blue)"
    elif dilution and "lavender" in dilution:
        eumelanin_word = dilution
    else:
        eumelanin_word = {"extended_black": "black", "birchen": "birchen", "wheaten": "wheaten",
                          "wild_partridge": "partridge", "brown": "brown",
                          "recessive_wheaten": "buff"}[canvas]

    pheo_word = "silver" if silver else "gold/red"
    parts_out = [eumelanin_word]
    if not dom_white and canvas not in ("extended_black",) and not low_eumelanin:
        parts_out.append(f"{pheo_word} ground")
    elif low_eumelanin:
        parts_out = [f"{pheo_word.replace('gold/red', 'red/buff')} ({canvas_label})"]
    color = ", ".join(parts_out)
    if barred and not low_eumelanin and not dom_white:
        color = f"barred {color}"
    if silver and canvas == "extended_black":
        color = f"silver {color}"

    conf = "high"
    if abstains or sex_basis.startswith("unspecified"):
        conf = "medium"
    return PlumageCall(color, canvas, sex_basis, barred, silver, dilution,
                       "dominant_white" if dom_white else None, dom_white, conf, abstains, per, rule, notes)


def reference_integrity_ok() -> bool:
    """Biology contract guard — pins known chicken genotypes -> plumage, INCLUDING the epistasis anchors a
    naive has-the-allele rule gets wrong (Extension canvas; Z-linked reversed-hemizygous barring; white masks)."""
    black = call_chicken_plumage({"E": "E/E"})                                 # extended black
    wheaten = call_chicken_plumage({"E": "EWh/EWh"})                           # mostly red/buff
    # ANCHOR 1: barring on a wheaten canvas barely shows (Extension is the canvas)
    barred_wheaten = call_chicken_plumage({"E": "EWh/EWh", "B": "B/b+"})
    # ANCHOR 2: Z-linked barring, a FEMALE is hemizygous (1 allele)
    barred_hen = call_chicken_plumage({"E": "E/E", "B": "B"})                  # 1 allele -> female
    barred_rock_male = call_chicken_plumage({"E": "E/E", "B": "B/B"})          # 2 alleles -> male
    # ANCHOR 3: white masks
    dom_white = call_chicken_plumage({"E": "E/E", "I": "I/i+"})
    rec_white = call_chicken_plumage({"E": "E/E", "C": "c/c"})
    # dilutions + silver
    blue = call_chicken_plumage({"E": "E/E", "BL": "Bl/bl+"})
    splash = call_chicken_plumage({"E": "E/E", "BL": "Bl/Bl"})
    lavender = call_chicken_plumage({"E": "E/E", "LAV": "lav/lav"})
    silver_birchen = call_chicken_plumage({"E": "ER/ER", "S": "S/S"})
    return (black.eumelanin_canvas == "extended_black" and "black" in black.plumage
            and wheaten.eumelanin_canvas == "wheaten"
            and any("naive" in n for n in barred_wheaten.notes)
            and "female" in barred_hen.sex_basis and barred_hen.barred
            and "male" in barred_rock_male.sex_basis and barred_rock_male.barred
            and dom_white.is_white_masked and rec_white.white_type == "recessive_white"
            and blue.dilution == "blue" and splash.dilution == "splash" and lavender.dilution == "lavender"
            and silver_birchen.silver)
