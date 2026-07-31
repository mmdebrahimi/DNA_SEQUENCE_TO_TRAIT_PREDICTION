"""Generic mammalian coat-colour engine + per-organism OMIA catalogs — a SHARED resolver for the classic
Extension/Agouti/Brown/Dilute/Albino epistatic series, so each new organism is a CATALOG, not a bespoke module.

This is the DRY generalisation of the bespoke dog/horse/cat cells: the CORE mammalian colour epistasis is
shared (albino masks all -> extension e/e recessive-red hides agouti -> agouti pattern -> brown eumelanin
colour -> dilution -> pink-eye), and organism-specific quirks (dominant black, incomplete-dominant dilution,
dominant white / KIT, sheep's ASIP-dominant-white) are handled as locus KINDS. Ships 5 organisms:

    rabbit  A/B/C/D/E (the textbook A-E series)          cattle  E(MC1R ED>E+>e) + PMEL dilution
    mouse   A/B/C/D/E/P (foundational; +pink-eye)         pig     KIT dominant-white + E(MC1R)
    sheep   ASIP(A^Wt dominant-white/tan > a black) + E(MC1R ED dominant black)

THE EPISTASIS (encoded by locus KIND, resolved in fixed order):
    dominant_white (KIT/ASIP-white) -> masks all              albino (C c/c) -> white, masks all
    extension (E): a DOMINANT-BLACK allele (E^D/E^d) -> solid eumelanin; a homozygous RECESSIVE-RED allele
        (e/e) -> phaeomelanin (red/yellow), which HIDES agouti (the recessive-epistasis anchor)
    agouti (A): pattern (agouti/tan/self), expressed only if eumelanin-capable; sheep A^Wt -> dominant white/tan
    brown (B): eumelanin colour (black vs chocolate/brown)    dilute (D d/d): eumelanin dilution (blue/lilac)
    dilution_incomplete (PMEL): dosage dilution (cattle dun/silver)   pink_eye (P p/p): pink-eyed dilution

  ANCHORS a naive has-the-allele rule mis-calls (pinned in reference_integrity_ok per organism):
    (1) recessive-red e/e is red REGARDLESS of agouti (extension epistatic over agouti);
    (2) albino c/c (rabbit/mouse) or dominant-white (pig KIT) MASKS every colour locus;
    (3) a dominant-black E^D allele (cattle/pig/sheep) -> solid black over agouti.

OMIA-sourced causal genes (per catalog `source`); KNOWLEDGE_BASELINE (no free per-individual validation
substrate). Pure-python, wheel-only, offline, deterministic. Benign livestock/lab-animal visible-trait genetics
— NOT human/forensic. Imports nothing from the frozen decoder surface.
"""
from __future__ import annotations

from dataclasses import dataclass, field

RULES_VERSION = "mammal-colour-v0.1.0"


@dataclass(frozen=True)
class Locus:
    symbol: str
    gene: str
    kind: str                       # albino|extension|agouti|brown|dilute|dilution_incomplete|pink_eye|dominant_white
    alleles: tuple[str, ...]        # dominance order, most-dominant first
    source: str
    black_alleles: frozenset = frozenset()   # extension: alleles that give solid eumelanin (dominant black)
    red_allele: str | None = None            # extension: recessive-red allele (homozygous -> phaeomelanin)
    recessive_white_allele: str | None = None  # extension: allele whose homozygote -> WHITE (camelid e/e=white)
    albino_allele: str | None = None         # albino: allele whose homozygote -> white
    white_alleles: frozenset = frozenset()   # agouti: dominant-white/tan alleles (sheep A^Wt)
    self_allele: str | None = None           # agouti: non-agouti / self (black) allele
    recessive_allele: str | None = None      # dilute/brown/pink: the recessive (b/d/p) allele
    dilution_alleles: frozenset = frozenset()  # dilution_incomplete: PMEL dosage alleles
    dominant_allele: str | None = None       # dominant_white: the masking allele (KIT I)
    note: str = ""


@dataclass(frozen=True)
class MammalCatalog:
    organism: str                   # scientific, e.g. "Oryctolagus_cuniculus"
    common: str                     # "rabbit"
    loci: dict[str, Locus]
    anchors: tuple                  # tuple of (label, genotype-dict, predicate(call)->bool) for integrity
    tier_note: str = "KNOWLEDGE_BASELINE (curated OMIA catalog; no free per-individual validation substrate)"


class MammalInputError(ValueError):
    """Unknown locus/allele or malformed genotype (never a silent wrong call)."""


@dataclass
class MammalColorCall:
    coat_color: str
    pigment_type: str               # "eumelanin" | "phaeomelanin" | "masked"
    base_eumelanin: str | None      # black/brown/blue/lilac/... (None when phaeomelanin/masked)
    pattern: str                    # solid | agouti | tan | self | white | n/a
    dilutions: list[str]
    is_white_masked: bool
    confidence: str
    abstains_on: list[str]
    per_locus: dict
    organism: str
    rule: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "organism": self.organism, "trait": "coat_colour",
            "regime": "A_curated_catalog_epistatic_mammalian", "rule": self.rule,
            "coat_color": self.coat_color, "pigment_type": self.pigment_type,
            "base_eumelanin": self.base_eumelanin, "pattern": self.pattern, "dilutions": self.dilutions,
            "is_white_masked": self.is_white_masked, "confidence": self.confidence,
            "abstains_on": self.abstains_on, "per_locus": self.per_locus, "notes": self.notes,
            "scope_limit": "curated OMIA loci; calls colour + agouti pattern, NOT fine pattern/spotting/shade",
            "evidence_tier": "knowledge_baseline (curated OMIA catalog; no free per-individual validation substrate)",
        }


def _parse(cat: MammalCatalog, symbol: str, spec: str) -> tuple[str, str]:
    loc = cat.loci[symbol]
    parts = [p for p in spec.replace("|", "/").split("/") if p.strip()]
    if len(parts) != 2:
        raise MammalInputError(f"{symbol} genotype {spec!r} is not diploid (expected a1/a2)")
    out = []
    for p in parts:
        p = p.strip()
        if p not in loc.alleles:
            raise MammalInputError(f"unknown {symbol}-locus allele {p!r}; recognised: {list(loc.alleles)}")
        out.append(p)
    return out[0], out[1]


def _dominant(loc: Locus, geno) -> str:
    return min(geno, key=lambda a: loc.alleles.index(a))


def call_mammal_color(cat: MammalCatalog, loci_genotypes: dict[str, str]) -> MammalColorCall:
    """Deterministic mammalian coat-colour call for `cat`'s organism from per-locus allele genotypes."""
    rule = f"mammal_colour_epistatic_v0[{cat.common}] ({RULES_VERSION})"
    notes: list[str] = []
    per: dict = {}
    geno: dict[str, tuple[str, str]] = {}
    kinds: dict[str, str] = {}       # symbol -> kind
    for sym, spec in loci_genotypes.items():
        S = sym.strip()
        if S not in cat.loci:
            raise MammalInputError(f"unknown locus {sym!r} for {cat.common}; loci: {list(cat.loci)}")
        geno[S] = _parse(cat, S, spec)
        per[S] = "/".join(geno[S])
        kinds[S] = cat.loci[S].kind

    def by_kind(k):
        return [s for s in geno if kinds[s] == k]

    # ---- 1. dominant_white (KIT I) ----
    for s in by_kind("dominant_white"):
        loc = cat.loci[s]
        if loc.dominant_allele in geno[s]:
            notes.append(f"{s} dominant-white ({loc.gene}) MASKS all colour (epistatic)")
            return MammalColorCall("white (dominant white)", "masked", None, "white", [], True, "high",
                                   [], per, cat.organism, rule, notes)

    # ---- 2. albino (C c/c) ----
    for s in by_kind("albino"):
        loc = cat.loci[s]
        if loc.albino_allele and geno[s] == (loc.albino_allele, loc.albino_allele):
            notes.append(f"{s} {loc.albino_allele}/{loc.albino_allele} -> albino ({loc.gene}); masks all colour")
            return MammalColorCall("albino (white)", "masked", None, "white", [], True, "high",
                                   [], per, cat.organism, rule, notes)

    # ---- 3. extension (E) ----
    pigment = "eumelanin"
    solid_black = False
    ext = by_kind("extension")
    if ext:
        loc = cat.loci[ext[0]]
        g = geno[ext[0]]
        edom = _dominant(loc, g)
        if loc.recessive_white_allele and g == (loc.recessive_white_allele, loc.recessive_white_allele):
            notes.append(f"e/e ({loc.recessive_white_allele}) -> recessive WHITE ({loc.gene}); masks all colour "
                         "(the camelid pattern: ee is white regardless of ASIP)")
            return MammalColorCall("white (recessive white)", "masked", None, "white", [], True, "high",
                                   [], per, cat.organism, rule, notes)
        if edom in loc.black_alleles:
            solid_black = True
            notes.append(f"E dominant-black allele {edom!r} -> solid eumelanin (epistatic over agouti)")
        elif loc.red_allele and g == (loc.red_allele, loc.red_allele):
            pigment = "phaeomelanin"
            notes.append(f"e/e ({loc.red_allele}) -> phaeomelanin (red/yellow); HIDES agouti (recessive epistasis "
                         "— a naive 'agouti -> banded' rule mis-calls this)")

    # ---- 4. agouti (A) ----
    pattern = "solid"
    if pigment == "eumelanin" and not solid_black:
        for s in by_kind("agouti"):
            loc = cat.loci[s]
            adom = _dominant(loc, geno[s])
            if adom in loc.white_alleles:                # sheep A^Wt dominant white/tan
                notes.append(f"{s} dominant-white/tan allele {adom!r} ({loc.gene}) -> white/tan")
                return MammalColorCall("white/tan (dominant agouti)", "phaeomelanin", None, "white/tan", [],
                                       False, "high", [], per, cat.organism, rule, notes)
            if loc.self_allele and geno[s] == (loc.self_allele, loc.self_allele):
                pattern = "self"                          # non-agouti -> solid eumelanin
            elif adom == loc.alleles[0]:
                pattern = "agouti"
            else:
                pattern = "tan"

    # ---- 5. brown (B) ----
    eu = "black"
    for s in by_kind("brown"):
        loc = cat.loci[s]
        if loc.recessive_allele and geno[s] == (loc.recessive_allele, loc.recessive_allele):
            eu = "brown/chocolate"

    # ---- 6. dilute (D d/d) + 7. dilution_incomplete (PMEL) + 8. pink_eye ----
    dilutions: list[str] = []
    for s in by_kind("dilute"):
        loc = cat.loci[s]
        if loc.recessive_allele and geno[s] == (loc.recessive_allele, loc.recessive_allele):
            dilutions.append("dilute")
            eu = {"black": "blue", "brown/chocolate": "lilac"}.get(eu, eu)
    for s in by_kind("dilution_incomplete"):
        loc = cat.loci[s]
        n = sum(1 for a in geno[s] if a in loc.dilution_alleles)
        if n == 1:
            dilutions.append("dilution x1 (partial)")
        elif n == 2:
            dilutions.append("dilution x2 (full)")
            eu = "silver/dun" if eu == "black" else eu
    for s in by_kind("pink_eye"):
        loc = cat.loci[s]
        if loc.recessive_allele and geno[s] == (loc.recessive_allele, loc.recessive_allele):
            dilutions.append("pink-eyed dilution")

    # ---- compose ----
    if pigment == "phaeomelanin":
        core = "cream" if "dilute" in dilutions else "red/yellow"
        base_eu = None
    else:
        base_eu = eu
        if pattern == "agouti":
            core = f"agouti ({eu}-banded)"
        elif pattern == "tan":
            core = f"tan-pattern ({eu})"
        else:
            core = eu                                    # self / solid
    if dilutions and pigment == "eumelanin":
        core = f"{core} [{', '.join(dilutions)}]"
    conf = "high" if not notes or True else "medium"
    return MammalColorCall(core, pigment, base_eu, pattern, dilutions, False, "high", [],
                           per, cat.organism, rule, notes)


def reference_integrity_ok(cat: MammalCatalog) -> bool:
    """Run the catalog's pinned anchors (each a genotype -> predicate on the call)."""
    for _label, g, pred in cat.anchors:
        if not pred(call_mammal_color(cat, g)):
            return False
    return True


# ============================ per-organism OMIA catalogs ============================

def _A_series(gene, alleles, self_allele, white=frozenset(), src=""):
    return Locus("A", gene, "agouti", alleles, src, white_alleles=white, self_allele=self_allele)


RABBIT = MammalCatalog(
    "Oryctolagus_cuniculus", "rabbit",
    {
        "A": _A_series("ASIP", ("A", "at", "a"), "a", src="rabbit A locus (ASIP): A agouti > at tan > a self"),
        "B": Locus("B", "TYRP1", "brown", ("B", "b"), "rabbit B (TYRP1): B black > b chocolate", recessive_allele="b"),
        "C": Locus("C", "TYR", "albino", ("C", "cchd", "cchl", "ch", "c"),
                   "rabbit C (TYR): C full > chinchilla > Himalayan > c albino", albino_allele="c"),
        "D": Locus("D", "MLPH", "dilute", ("D", "d"), "rabbit D (MLPH): d/d dilute (blue/lilac)", recessive_allele="d"),
        "E": Locus("E", "MC1R", "extension", ("Ed", "Es", "E", "ej", "e"),
                   "rabbit E (MC1R): Ed dominant-black > Es steel > E > ej Japanese > e red",
                   black_alleles=frozenset({"Ed"}), red_allele="e"),
    },
    anchors=(
        ("albino masks", {"A": "A/a", "B": "b/b", "C": "c/c", "E": "E/e"}, lambda c: c.is_white_masked),
        ("e/e red hides agouti", {"A": "A/A", "E": "e/e"}, lambda c: c.pigment_type == "phaeomelanin"),
        ("black self", {"A": "a/a", "B": "B/B", "C": "C/C", "D": "D/D", "E": "E/E"}, lambda c: c.coat_color == "black"),
        ("blue (dilute)", {"A": "a/a", "B": "B/B", "D": "d/d", "E": "E/E"}, lambda c: c.base_eumelanin == "blue"),
        ("chocolate agouti", {"A": "A/A", "B": "b/b", "E": "E/E"}, lambda c: "brown/chocolate" in c.coat_color),
    ),
)

MOUSE = MammalCatalog(
    "Mus_musculus", "mouse",
    {
        "A": _A_series("ASIP", ("A", "at", "a"), "a", src="mouse a locus (ASIP): A agouti > at > a non-agouti"),
        "B": Locus("B", "TYRP1", "brown", ("B", "b"), "mouse b (Tyrp1): B black > b brown", recessive_allele="b"),
        "C": Locus("C", "TYR", "albino", ("C", "cch", "ch", "c"),
                   "mouse c (Tyr): C full > cch chinchilla > ch himalayan > c albino", albino_allele="c"),
        "D": Locus("D", "MYO5A", "dilute", ("D", "d"), "mouse d (Myo5a): d/d dilute", recessive_allele="d"),
        "P": Locus("P", "OCA2", "pink_eye", ("P", "p"), "mouse p (Oca2): p/p pink-eyed dilution", recessive_allele="p"),
        "E": Locus("E", "MC1R", "extension", ("E", "e"), "mouse e (Mc1r): E extension > e recessive yellow",
                   red_allele="e"),
    },
    anchors=(
        ("albino masks", {"A": "A/a", "C": "c/c", "E": "E/e"}, lambda c: c.is_white_masked),
        ("e/e yellow", {"A": "A/A", "E": "e/e"}, lambda c: c.pigment_type == "phaeomelanin"),
        ("black non-agouti", {"A": "a/a", "B": "B/B", "C": "C/C", "D": "D/D", "E": "E/E"}, lambda c: c.coat_color == "black"),
        ("brown", {"A": "a/a", "B": "b/b", "C": "C/C", "E": "E/E"}, lambda c: "brown/chocolate" in c.coat_color),
    ),
)

CATTLE = MammalCatalog(
    "Bos_taurus", "cattle",
    {
        "E": Locus("E", "MC1R", "extension", ("ED", "E+", "e"),
                   "cattle E (MC1R): ED dominant-black > E+ wild(red-brown) > e recessive-red",
                   black_alleles=frozenset({"ED"}), red_allele="e"),
        "PMEL": Locus("PMEL", "PMEL/SILV", "dilution_incomplete", ("Dc", "Dh", "n"),
                      "cattle dilution (PMEL/SILV): Dc Charolais / Dh Highland, incompletely dominant -> dun/silver",
                      dilution_alleles=frozenset({"Dc", "Dh"})),
    },
    anchors=(
        ("dominant black", {"E": "ED/e"}, lambda c: c.coat_color == "black"),
        ("recessive red", {"E": "e/e"}, lambda c: c.pigment_type == "phaeomelanin"),
        ("wild red-brown", {"E": "E+/E+"}, lambda c: c.pigment_type == "eumelanin"),
        ("dilution", {"E": "E+/E+", "PMEL": "Dc/n"}, lambda c: any("dilution" in d for d in c.dilutions)),
    ),
)

PIG = MammalCatalog(
    "Sus_scrofa", "pig",
    {
        "KIT": Locus("KIT", "KIT", "dominant_white", ("I", "i+"),
                     "pig Dominant White (KIT): I dominant-white (epistatic over E) > i+", dominant_allele="I"),
        "E": Locus("E", "MC1R", "extension", ("ED", "E+", "e"),
                   "pig E (MC1R, OMIA 001199-9823): ED dominant-black > E+ wild > e recessive-red",
                   black_alleles=frozenset({"ED"}), red_allele="e"),
    },
    anchors=(
        ("KIT dominant white masks", {"KIT": "I/i+", "E": "e/e"}, lambda c: c.is_white_masked),
        ("dominant black", {"KIT": "i+/i+", "E": "ED/e"}, lambda c: c.coat_color == "black"),
        ("recessive red", {"KIT": "i+/i+", "E": "e/e"}, lambda c: c.pigment_type == "phaeomelanin"),
    ),
)

SHEEP = MammalCatalog(
    "Ovis_aries", "sheep",
    {
        # ASIP: A^Wt dominant white/tan (190kb duplication) > a recessive black (LOF). White is DOMINANT (ASIP-up).
        "A": Locus("A", "ASIP", "agouti", ("AWt", "a"), "sheep Agouti (ASIP): A^Wt dominant white/tan > a recessive black",
                   white_alleles=frozenset({"AWt"}), self_allele="a"),
        "E": Locus("E", "MC1R", "extension", ("ED", "E+", "e"),
                   "sheep E (MC1R): ED dominant-black (M73K/D121N) > E+ > e (R67C)",
                   black_alleles=frozenset({"ED"}), red_allele="e"),
    },
    anchors=(
        ("ED dominant black overrides ASIP white", {"A": "AWt/a", "E": "ED/E+"}, lambda c: c.coat_color == "black"),
        ("dominant white/tan", {"A": "AWt/a", "E": "E+/E+"}, lambda c: c.pattern == "white/tan"),
        ("recessive black", {"A": "a/a", "E": "E+/E+"}, lambda c: c.coat_color == "black"),
    ),
)

GOAT = MammalCatalog(
    "Capra_hircus", "goat",
    {
        # ASIP is the many-pattern hub in goats (~11 alleles, CNV-driven): A^Wt white/tan dominant > a nonagouti black.
        "A": Locus("A", "ASIP", "agouti", ("AWt", "a"), "goat Agouti (ASIP, OMIA 000201-9925): A^Wt dominant "
                   "white/tan (CNV) > a recessive nonagouti black; ~11 pattern alleles collapsed to the two poles",
                   white_alleles=frozenset({"AWt"}), self_allele="a"),
        "B": Locus("B", "TYRP1", "brown", ("B", "b"), "goat B (TYRP1): b/b brown (e.g. Copperneck)", recessive_allele="b"),
    },
    anchors=(
        ("dominant white/tan", {"A": "AWt/a"}, lambda c: c.pattern == "white/tan"),
        ("recessive black", {"A": "a/a", "B": "B/B"}, lambda c: c.coat_color == "black"),
        ("brown", {"A": "a/a", "B": "b/b"}, lambda c: "brown/chocolate" in c.coat_color),
    ),
)

ALPACA = MammalCatalog(
    "Vicugna_pacos", "alpaca",
    {
        # camelid model: E/_ coloured (black if aa, fawn/agouti if A_); e/e -> WHITE regardless of ASIP.
        "E": Locus("E", "MC1R", "extension", ("E", "e"), "alpaca/llama E (MC1R): E coloured > e; e/e = recessive WHITE",
                   recessive_white_allele="e"),
        "A": Locus("A", "ASIP", "agouti", ("A", "a"), "alpaca Agouti (ASIP): A functional -> fawn/agouti (pheo) > "
                   "a loss-of-function -> black (eumelanin)", self_allele="a"),
    },
    anchors=(
        ("ee recessive white", {"E": "e/e", "A": "a/a"}, lambda c: c.is_white_masked),
        ("black (E-, aa)", {"E": "E/E", "A": "a/a"}, lambda c: c.coat_color == "black"),
        ("fawn/agouti (E-, A_)", {"E": "E/E", "A": "A/a"}, lambda c: c.pattern == "agouti"),
    ),
)

MAMMAL_CATALOGS: dict[str, MammalCatalog] = {
    "rabbit": RABBIT, "mouse": MOUSE, "cattle": CATTLE, "pig": PIG, "sheep": SHEEP,
    "goat": GOAT, "alpaca": ALPACA,
}
