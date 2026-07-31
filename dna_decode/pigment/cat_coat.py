"""Cat (Felis catus) coat-colour decoder — the curated-catalog + EPISTASIS paradigm applied to a system
whose signature feature is X-LINKED mosaicism (the sibling of dog_coat / horse_coat; a 3rd-organism cell).

THE BIOLOGY (OMIA-curated causal variants; the standard feline colour model):

    W  (KIT)      DOMINANT WHITE, autosomal dominant + EPISTATIC over everything: W_ = solid white (may be
                  blue-eyed / deaf). 3-allele series W > ws > w. Causal: FERV1 endogenous-retrovirus LTR
                  (dominant white) / full FERV1 (white spotting ws). David 2014. OMIA 000209 / 001737-9685.
    O  (ARHGAP36) SEX-LINKED ORANGE — X-linked. O = orange (phaeomelanin, suppresses black-based colour,
                  EPISTATIC over B). Because X-linked, ZYGOSITY depends on SEX: a MALE (XY) is hemizygous
                  (O -> orange, o -> non-orange); a FEMALE (XX) O/O = orange, o/o = non-orange, and O/o =
                  TORTOISESHELL (a black+orange MOSAIC from random X-inactivation). Causal: 5.1-kb intron-1
                  deletion (Toh/Kaelin 2025, Current Biology). This is THE cat-specific anchor.
    A  (ASIP)     agouti: A_ = agouti (tabby banding visible), a/a = non-agouti (solid/self). Non-agouti does
                  NOT suppress the tabby pattern in ORANGE fur (red cats always show tabby). Eizirik 2003.
    B  (TYRP1)    eumelanin colour: B = black > b = chocolate > bl = cinnamon. Recessive series. Schmidt-Kuntzel 2005.
    D  (MLPH)     dilution: d/d = dilute (black->blue, chocolate->lilac, cinnamon->fawn, orange->cream). Recessive.
    C  (TYR)      albino series: C (full) > cb (Burmese sepia) = cs (Siamese POINTED, temp-sensitive) > c
                  (albino, white/blue-eyed). cs/cs = points; cb/cb = sepia; cs/cb = mink. Lyons 2005. OMIA 000202-9685.

  Resolution order (epistasis hierarchy): W dominant-white -> C albino -> Orange (sex-dependent) -> eumelanin
  base (B x D) / orange colour (D) -> agouti tabby -> colorpoint overlay -> white spotting (ws; + tortie = CALICO).

  THE EPISTASIS ANCHORS (the cases a naive has-the-allele rule mis-calls):
    (1) W/_ -> solid WHITE regardless of every colour locus (dominant-white epistasis).
    (2) a female O/o -> TORTOISESHELL (mosaic), NOT uniform orange or black (X-inactivation); a naive
        autosomal rule mis-calls. A tortie + white spotting -> CALICO.
    (3) O (orange) is EPISTATIC over B: an orange cat is red even when b/b (chocolate); the B genes show
        only in the non-orange fur.
  `reference_integrity_ok` pins exactly these.

WHY A NON-FROZEN cell (like dog/horse coat): fixed-order epistasis + an X-LINKED sex-dependent locus + a
mosaic phenotype — not the frozen count/OR amr_rules shape. Imports nothing from the frozen decoder surface.

HONEST SCOPE (load-bearing):
  - v0 = the loci the user named + brown: W/O/A/B/D/C. Tabby PATTERN sub-type (mackerel/classic/ticked/spotted;
    Taqpep/Ticked), silver/inhibitor, wideband, karpati/roan, and rare alleles ABSTAIN (via --present).
  - Input is PER-LOCUS allele CALLS (the VGL/UC-Davis report shape). The O locus takes ONE allele for a male
    (hemizygous, e.g. O=O) or TWO for a female (e.g. O=O/o); `--sex` may be given but is inferred from the O
    zygosity when omitted. A genome/VCF caller is a v0.1 follow-on.
  - Calls the COLOUR + major pattern (tortie/calico/pointed/bicolor), NOT the tabby sub-pattern, shade, or
    the exact spotting extent. KNOWLEDGE_BASELINE: no free per-individual validation substrate.
  - Faithful-to-literature: applies the published OMIA locus/allele assignments; it is not a new model.

Pure-python, wheel-only, offline, deterministic. Benign companion-animal visible-trait genetics — NOT human/forensic.
"""
from __future__ import annotations

from dataclasses import dataclass, field

RULES_VERSION = "cat-coat-colour-v0.1.0"


@dataclass(frozen=True)
class Locus:
    name: str
    gene: str
    alleles: tuple[str, ...]        # dominance order, most-dominant first
    source: str
    x_linked: bool = False
    note: str = ""


LOCI: dict[str, Locus] = {
    "W": Locus("W", "KIT", ("W", "ws", "w"),
               source="OMIA 000209-9685 (dominant white, FERV1 LTR) / 001737-9685 (white spotting, full FERV1); David 2014",
               note="W dominant-white EPISTATIC (masks all colour); ws = white spotting (bicolor/calico); w = none"),
    "O": Locus("O", "ARHGAP36", ("O", "o"),
               source="Toh 2025 / Kaelin 2025 Current Biology (X-linked; 5.1-kb ARHGAP36 intron-1 deletion)",
               x_linked=True,
               note="X-linked orange: O -> orange (EPISTATIC over B). Female O/o -> TORTOISESHELL mosaic (XCI)"),
    "A": Locus("A", "ASIP", ("A", "a"),
               source="OMIA (Agouti/ASIP; Eizirik 2003)",
               note="A_ = agouti (tabby visible), a/a = non-agouti (solid); orange fur shows tabby regardless"),
    "B": Locus("B", "TYRP1", ("B", "b", "bl"),
               source="OMIA (TYRP1; Schmidt-Kuntzel 2005) b=chocolate, bl=cinnamon (premature stop)",
               note="eumelanin colour: B black > b chocolate > bl cinnamon"),
    "D": Locus("D", "MLPH", ("D", "d"),
               source="OMIA (MLPH single-base deletion; Ishida 2006)",
               note="d/d dilute: black->blue, chocolate->lilac, cinnamon->fawn, orange->cream"),
    "C": Locus("C", "TYR", ("C", "cb", "cs", "c"),
               source="OMIA 000202-9685 (TYR; Lyons 2005) cs Siamese G302R, cb Burmese, c albino",
               note="C full > cb Burmese-sepia = cs Siamese-points > c albino; cs/cb = mink"),
}

UNMODELLED_LOCI = {
    "TA": "tabby PATTERN sub-type (mackerel/classic/spotted Taqpep; ticked Ticked) — v0 calls tabby-vs-solid, not the sub-pattern",
    "I": "silver/inhibitor (melanin inhibitor) — smoke/silver/golden series",
    "WB": "wideband — affects agouti band width (shaded/chinchilla)",
    "KA": "karpati/roan — progressive white ticking",
    "FGF5": "long hair (FGF5) — coat LENGTH, not colour",
}

_ALLELE_ALIASES = {"cs": "cs", "cb": "cb", "bl": "bl", "ws": "ws"}


class CatInputError(ValueError):
    """Unknown locus/allele or malformed genotype (never a silent wrong call)."""


@dataclass
class CatCoatCall:
    coat_color: str                 # composed human label
    base_color: str                 # underlying pigment colour (black/blue/chocolate/red/tortoiseshell/...)
    sex_basis: str                  # "male (hemizygous O)" | "female (XX O)" | "unspecified"
    is_tortoiseshell: bool
    is_epistatic_white: bool        # dominant-white or albino masks colour
    tabby: bool
    colorpoint: str | None          # "siamese_points" | "burmese_sepia" | "mink" | None
    white_pattern: str | None       # "bicolor/white-spotting" | "calico" | None
    dilute: bool
    confidence: str
    abstains_on: list[str]
    per_locus: dict
    rule: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "organism": "Felis_catus", "trait": "coat_colour",
            "regime": "A_curated_catalog_epistatic_xlinked", "rule": self.rule,
            "coat_color": self.coat_color, "base_color": self.base_color, "sex_basis": self.sex_basis,
            "is_tortoiseshell": self.is_tortoiseshell, "is_epistatic_white": self.is_epistatic_white,
            "tabby": self.tabby, "colorpoint": self.colorpoint, "white_pattern": self.white_pattern,
            "dilute": self.dilute, "confidence": self.confidence, "abstains_on": self.abstains_on,
            "per_locus": self.per_locus, "notes": self.notes,
            "scope_limit": ("v0: W/O/A/B/D/C; calls colour + major pattern (tortie/calico/pointed/bicolor), "
                            "NOT tabby sub-pattern, shade, spotting extent, or silver/wideband"),
            "evidence_tier": "knowledge_baseline (curated OMIA catalog; no free per-individual validation substrate)",
        }


def _canon_allele(locus: str, tok: str) -> str:
    t = tok.strip()
    if t in LOCI[locus].alleles:
        return t
    low = t.lower()
    if low in _ALLELE_ALIASES and _ALLELE_ALIASES[low] in LOCI[locus].alleles:
        return _ALLELE_ALIASES[low]
    raise CatInputError(f"unknown {locus}-locus allele {tok!r}; recognised: {list(LOCI[locus].alleles)}")


def _parse(locus: str, spec: str, allow_hemizygous: bool = False) -> tuple[str, ...]:
    parts = [p for p in spec.replace("|", "/").split("/") if p.strip()]
    if allow_hemizygous and len(parts) == 1:
        return (_canon_allele(locus, parts[0]),)
    if len(parts) != 2:
        exp = "1 (male) or 2 (female) alleles" if allow_hemizygous else "a1/a2"
        raise CatInputError(f"{locus} genotype {spec!r} must be {exp}")
    return _canon_allele(locus, parts[0]), _canon_allele(locus, parts[1])


def _dominant(locus: str, geno: tuple[str, ...]) -> str:
    order = LOCI[locus].alleles
    return min(geno, key=lambda a: order.index(a))


def _eumelanin(b_geno, d_geno) -> str:
    """Eumelanin colour from B (black/choc/cinnamon) x D (dilute)."""
    b = _dominant("B", b_geno) if b_geno else "B"
    dilute = bool(d_geno) and all(a == "d" for a in d_geno)
    base = {"B": "black", "b": "chocolate", "bl": "cinnamon"}[b]
    if dilute:
        return {"black": "blue", "chocolate": "lilac", "cinnamon": "fawn"}[base]
    return base


def call_cat_coat(loci_genotypes: dict[str, str], sex: str | None = None,
                  present_loci: list[str] | None = None) -> CatCoatCall:
    """Deterministic cat coat-colour call from per-locus allele genotypes.

    `loci_genotypes`: {locus -> 'a1/a2'} for any subset of W/O/A/B/D/C. The O (X-linked) locus takes ONE
    allele for a male (e.g. 'O') or TWO for a female (e.g. 'O/o'). `sex` ('male'/'female') optional — inferred
    from O zygosity when omitted. `present_loci`: v0-unmodelled loci present -> the affected axis ABSTAINS.
    """
    rule = f"cat_coat_colour_epistatic_v0 ({RULES_VERSION})"
    notes: list[str] = []
    per: dict = {}
    geno: dict[str, tuple[str, ...]] = {}
    for loc, spec in loci_genotypes.items():
        L = loc.strip().upper()
        if L not in LOCI:
            if L in UNMODELLED_LOCI:
                raise CatInputError(f"locus {loc!r} is a v0-unmodelled locus ({UNMODELLED_LOCI[L]}); pass it "
                                    f"via present_loci=[...] so the affected axis ABSTAINS")
            raise CatInputError(f"unknown locus {loc!r}; v0 loci: {list(LOCI)}")
        geno[L] = _parse(L, spec, allow_hemizygous=LOCI[L].x_linked)
        per[L] = "/".join(geno[L])

    abstains: list[str] = []
    for pl in (present_loci or []):
        P = pl.strip().upper()
        if P in UNMODELLED_LOCI:
            abstains.append(f"{P}: {UNMODELLED_LOCI[P]}")
        elif P not in LOCI:
            raise CatInputError(f"unknown present locus {pl!r}")

    # ---- sex basis (from --sex or O zygosity) ----
    o_geno = geno.get("O")
    if sex:
        sex = sex.strip().lower()
        if sex not in ("male", "female"):
            raise CatInputError("--sex must be 'male' or 'female'")
        sex_basis = f"{sex} ({'hemizygous O' if sex == 'male' else 'XX O'})"
    elif o_geno is not None:
        sex = "male" if len(o_geno) == 1 else "female"
        sex_basis = f"{sex} (inferred from O zygosity: {len(o_geno)} allele{'s' if len(o_geno) == 2 else ''})"
    else:
        sex = None
        sex_basis = "unspecified (no O genotype / --sex)"

    # ---- 1. Dominant White (W) — top epistatic ----
    if "W" in geno and "W" in geno["W"]:
        notes.append("W dominant-white MASKS all colour (epistatic); may be blue-eyed / deaf")
        return CatCoatCall("white (dominant white)", "masked", sex_basis, False, True, False, None, None, False,
                           "high", abstains, per, rule, notes)

    # ---- 2. Albino (C c/c) ----
    colorpoint = None
    if "C" in geno:
        cdom = _dominant("C", geno["C"])
        if geno["C"] == ("c", "c"):
            notes.append("C locus c/c -> albino (white, blue eyes); masks colour")
            return CatCoatCall("albino (white, blue-eyed)", "masked", sex_basis, False, True, False, None, None,
                               False, "high", abstains, per, rule, notes)
        if cdom in ("cs", "cb"):
            has_cs, has_cb = "cs" in geno["C"], "cb" in geno["C"]
            colorpoint = "mink" if has_cs and has_cb else ("siamese_points" if has_cs else "burmese_sepia")

    # ---- 3. Orange (X-linked, sex-dependent) ----
    if o_geno is None:
        orange_state = "non_orange"
        notes.append("no O genotype -> assuming non-orange (black-based); pass O (+ sex) to resolve orange/tortie")
    elif len(o_geno) == 1:                                   # male hemizygous
        orange_state = "orange" if o_geno[0] == "O" else "non_orange"
    else:                                                    # female XX
        if o_geno == ("O", "O"):
            orange_state = "orange"
        elif o_geno == ("o", "o"):
            orange_state = "non_orange"
        else:
            orange_state = "tortoiseshell"
    if orange_state == "tortoiseshell" and sex == "male":
        notes.append("O/o in a MALE is abnormal (tortie males are rare XXY / chimera / mosaic) — flagged")

    # ---- 4. base colours ----
    dilute = "D" in geno and geno["D"] == ("d", "d")
    eu = _eumelanin(geno.get("B"), geno.get("D"))            # black/blue/chocolate/lilac/cinnamon/fawn
    orange_col = "cream" if dilute else "red"
    if "B" not in geno and orange_state != "orange":
        notes.append("B (TYRP1) absent -> eumelanin colour assumes black (B/B default)")

    # ---- 5. agouti (tabby) ----
    tabby = ("A" in geno and "A" in geno["A"])
    if "A" not in geno:
        notes.append("A (agouti) absent -> tabby-vs-solid assumed SOLID (a/a default); pass A to resolve")

    # ---- 6. compose base ----
    is_tortie = orange_state == "tortoiseshell"
    if orange_state == "orange":
        base = orange_col
        tabby = True                                         # orange always shows tabby
        if "B" in geno and geno["B"] != ("B", "B"):
            notes.append("orange is EPISTATIC over B: cat is red/cream; the b/bl colour shows only in "
                         "non-orange fur (none here)")
    elif is_tortie:
        base = f"tortoiseshell ({eu} + {orange_col})"
        notes.append("X-linked O/o in a female -> TORTOISESHELL: a mosaic of black-based + orange patches "
                     "from random X-inactivation (a naive autosomal rule mis-calls this)")
    else:
        base = eu

    # ---- 7. colorpoint overlay + tabby label ----
    solid_or_tabby = "tabby" if (tabby and not is_tortie) else ("tortie-tabby (torbie)" if (tabby and is_tortie) else "solid")
    color = base
    if orange_state != "orange" and not is_tortie and tabby:
        color = f"{base} tabby"
    elif orange_state == "orange":
        color = f"{base} tabby"                              # red/cream tabby
    elif is_tortie and tabby:
        color = f"{base} (torbie/patched tabby)"

    # ---- 8. white spotting (ws) -> bicolor; tortie + ws -> calico ----
    white_pattern = None
    if "W" in geno and "ws" in geno["W"]:
        if is_tortie:
            white_pattern = "calico"
            color = f"calico ({eu} + {orange_col} + white)"
            notes.append("tortoiseshell + white spotting (ws) -> CALICO (distinct patches broken by white)")
        else:
            white_pattern = "bicolor/white-spotting"
            color = f"{color} + white (bicolor)"

    # ---- 9. colorpoint pattern note ----
    if colorpoint:
        pat = {"siamese_points": "Siamese POINTED (colour only at extremities; temp-sensitive)",
               "burmese_sepia": "Burmese sepia (reduced pigment on torso)",
               "mink": "mink (Tonkinese; cs/cb intermediate)"}[colorpoint]
        color = f"{color}, {pat}"

    conf = "high"
    if abstains or sex_basis.startswith("unspecified"):
        conf = "medium"
    return CatCoatCall(color, base, sex_basis, is_tortie, False, tabby, colorpoint, white_pattern, dilute,
                       conf, abstains, per, rule, notes)


def reference_integrity_ok() -> bool:
    """Biology contract guard — pins known cat genotypes -> colours, INCLUDING the epistasis anchors a naive
    has-the-allele rule gets wrong (dominant-white epistasis; X-linked tortoiseshell; orange epistatic over B)."""
    # ANCHOR 1: dominant white masks everything
    white = call_cat_coat({"W": "W/w", "O": "o/o", "B": "b/b"})
    # ANCHOR 2: female O/o -> tortoiseshell (mosaic)
    tortie = call_cat_coat({"O": "O/o", "B": "B/B", "D": "D/D"})
    calico = call_cat_coat({"O": "O/o", "W": "ws/w", "B": "B/B"})
    # ANCHOR 3: orange epistatic over B (a b/b orange MALE is red, not chocolate)
    red_male = call_cat_coat({"O": "O", "B": "b/b"})
    # base colours
    black = call_cat_coat({"O": "o", "A": "a/a", "B": "B/B", "D": "D/D"})   # non-orange male, solid black
    blue = call_cat_coat({"O": "o", "A": "a/a", "B": "B/B", "D": "d/d"})    # dilute -> blue
    choc = call_cat_coat({"O": "o", "A": "a/a", "B": "b/b", "D": "D/D"})    # chocolate
    siamese = call_cat_coat({"O": "o", "C": "cs/cs", "B": "B/B"})            # pointed
    return (white.is_epistatic_white and white.coat_color.startswith("white (dominant")
            and tortie.is_tortoiseshell and tortie.base_color.startswith("tortoiseshell")
            and calico.white_pattern == "calico"
            and red_male.base_color == "red" and "red" in red_male.coat_color
            and black.base_color == "black" and blue.base_color == "blue" and choc.base_color == "chocolate"
            and siamese.colorpoint == "siamese_points")
