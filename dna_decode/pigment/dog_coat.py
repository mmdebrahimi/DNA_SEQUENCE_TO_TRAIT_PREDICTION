"""Dog (Canis lupus familiaris) coat-colour decoder — the curated-catalog determinant->phenotype
paradigm applied to a VISIBLE PHYSICAL trait, with the twist mammalian pigmentation adds that a flat
"has-the-allele" rule does not: EPISTASIS across loci.

THE BIOLOGY (the rule this encodes — classical Little 1957 / Schmutz & Berryere, curated in OMIA):

    Coat colour is set by a small number of well-characterised loci acting in a FIXED epistatic order.
    Two pigments exist: EUMELANIN (black->brown->blue->isabella) and PHAEOMELANIN (red/yellow). The loci:

      E  (MC1R)   pigment-type SWITCH.  e/e  -> NO eumelanin in the coat -> RED/YELLOW, regardless of
                  every other colour locus. This is the top of the epistasis: an e/e dog is red even if
                  it is genetically K^B (dominant black) and b/b (brown). Causal: MC1R p.Arg306Ter (`e`).
      K  (CBD103) eumelanin DISTRIBUTION.  K^B (dominant black) -> solid eumelanin, MASKS the A locus.
                  k^y/k^y -> the A (agouti) pattern is expressed. Causal K^B: CBD103 c.67_69delGGT.
      A  (ASIP)   agouti pattern, expressed ONLY when eumelanin-capable (E-) AND k^y/k^y (no K^B):
                  A^y fawn/sable > a^w wild/agouti > a^t tan-points > a recessive-black.
      B  (TYRP1)  eumelanin COLOUR.  b/b -> brown/liver eumelanin (+ brown nose). Causal b: TYRP1 p.Gln331Ter (+ b^d/b^c).
      D  (MLPH)   eumelanin DILUTION.  d/d -> dilute (black->blue/grey, brown->isabella/lilac). Causal: MLPH c.-22G>A.

    So eumelanin colour = f(B, D):  B-/D- black | b/b brown | D-/dd blue | b/b + dd isabella/lilac.

    THE EPISTASIS ANCHOR (the case a naive "has-the-allele" rule mis-calls — the citrate/Da(1)-12 of coat
    colour): an e/e dog that is ALSO K^B and b/b is RED/YELLOW, not black and not brown. A rule that reads
    "K^B -> black" or "b/b -> brown" without the E-locus gate calls it wrong. `reference_integrity_ok`
    pins exactly this.

WHY A NON-FROZEN cell (like flowering / metabolic / TMP-SMX / TB): the rule is a fixed-order EPISTASIS
across loci with a pigment-type switch — not the frozen count/OR `amr_rules.DRUG_RULE` shape. The frozen
decoder surface is untouched (this imports nothing from it).

HONEST SCOPE (load-bearing):
  - v0 = the FIVE classic solid-colour loci (E/K/A/B/D). PATTERN loci that also change appearance —
    S/white-spotting (MITF), M/merle (PMEL/SILV), T/ticking, H/harlequin, saddle — are DELIBERATELY OUT;
    a dog carrying them ABSTAINS on the affected axis rather than a confident wrong solid-colour call.
  - Input is PER-LOCUS allele calls (the VGL/Embark report shape), e.g. `E=e/e,K=ky/ky,A=at/at,B=B/b,D=D/d`.
    A genome/VCF mode (calling these loci from the canFam4 causal-variant genotypes) is a v0.1 follow-on
    (the Darwin's Ark Dryad cohort is the validation substrate — scripts/dog_coat_darwins_ark_validate.py).
  - Calls the COLOUR (pigment type + eumelanin colour + distribution pattern) — NOT shade intensity
    (red vs cream), coat length/texture, or spotting extent. Anything not in the catalogued loci ABSTAINS.
  - Faithful-to-literature: it applies the published locus/allele assignments (OMIA); it is not a new model.

Pure-python, wheel-only, offline, deterministic. Regime-A curated catalog, NOT a learned embedding.
Scope: benign visible-trait genetics of a companion animal — NOT any human/forensic application.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# ---- curated allele catalog (each locus sourced) -------------------------------------------------------
# Per locus: the recognised alleles in DOMINANCE order (most dominant first), plus the causal gene + a
# source. Allele tokens are the standard breeder/VGL symbols. The b- and d-families collapse to a single
# recessive class for colour purposes (any b-family allele is recessive to B; any d-family recessive to D).


@dataclass(frozen=True)
class Locus:
    name: str
    gene: str
    alleles: tuple[str, ...]        # dominance order, most-dominant first
    source: str
    note: str = ""


LOCI: dict[str, Locus] = {
    "E": Locus("E", "MC1R", ("Em", "Eg", "E", "e"),
               source="OMIA 001199-9615 (MC1R p.Arg306Ter = recessive red `e`); Everts 2000; Schmutz 2003",
               note="pigment-type switch: e/e -> phaeomelanin (red/yellow), EPISTATIC over K/A/B in the coat. "
                    "Em = melanistic mask, Eg = grizzle/domino (both eumelanin-capable)"),
    "K": Locus("K", "CBD103", ("KB", "kbr", "ky"),
               source="OMIA (K locus / dominant black, CBD103 beta-defensin c.67_69delGGT); Candille 2007 Science",
               note="K^B dominant black -> solid eumelanin, MASKS A; kbr brindle; ky/ky -> A expressed"),
    "A": Locus("A", "ASIP", ("Ay", "aw", "at", "a"),
               source="OMIA (A locus ASIP); Dreger & Schmutz 2011; Bannasch 2021 (ASIP promoter)",
               note="agouti pattern, expressed only if E- and ky/ky: Ay fawn/sable > aw wild > at tan-points > a recessive-black"),
    "B": Locus("B", "TYRP1", ("B", "bs", "bd", "bc", "b"),
               source="OMIA 001249-9615 (TYRP1; bs p.Gln331Ter, bd, bc); Schmutz 2002",
               note="b-family (bs/bd/bc/b) recessive to B: b/b -> brown/liver eumelanin + brown nose"),
    "D": Locus("D", "MLPH", ("D", "d", "d2"),
               source="OMIA 000031-9615 (MLPH c.-22G>A `d`; a second variant `d2`); Drogemuller 2007; Bauer 2018",
               note="d-family (d/d2) recessive to D: d/d -> dilute (black->blue, brown->isabella/lilac)"),
}

# Pattern / spotting loci NOT modelled in v0 — presence means the appearance axis they control ABSTAINS.
UNMODELLED_LOCI = {
    "S": "white spotting (MITF) — piebald/Irish/extreme white; changes how much coat is coloured",
    "M": "merle (PMEL/SILV) — mottled dilution patches; a merle genotype is not a solid colour",
    "T": "ticking (USH2A region) — flecks of colour in white",
    "H": "harlequin (PSMB7, Great Dane) — modifies merle",
    "I": "intensity — red/yellow SHADE (deep red vs cream); v0 calls red/yellow, not the shade",
}

# tolerant allele aliases -> canonical token
_ALLELE_ALIASES = {
    "em": "Em", "eg": "Eg", "kb": "KB", "kbr": "kbr", "ay": "Ay", "aw": "aw", "at": "at",
    "bs": "bs", "bd": "bd", "bc": "bc", "d2": "d2",
}
# recessive-class membership (colour-collapsed)
_B_RECESSIVE = {"bs", "bd", "bc", "b"}
_D_RECESSIVE = {"d", "d2"}

UNSEEN_MECHANISMS = (
    "a variant that SILENTLY changes function of a locus not in the E/K/A/B/D set (this reads per-locus "
    "allele CALLS, not sequence) — e.g. a novel MC1R allele reported as `E`",
    "pattern / white-spotting extent (S/M/T/H loci) — v0 ABSTAINS on the spotting axis rather than guessing",
    "shade INTENSITY (deep red vs cream; the I locus / MFSD12 etc.) — v0 calls red/yellow, not the shade",
    "somatic / mosaic / acquired colour change (greying with age, vitiligo) — genotype calls the base coat",
)


class CoatInputError(ValueError):
    """Unknown locus/allele or malformed genotype (never a silent wrong call)."""


@dataclass
class CoatColorCall:
    coat_color: str                 # human label, e.g. "solid black" / "brown/liver" / "red/yellow" / "ABSTAIN"
    pigment_type: str               # "eumelanin" | "phaeomelanin"
    eumelanin_color: str | None     # black/brown/blue/isabella (None when phaeomelanin)
    distribution: str               # "solid" | "sable" | "agouti" | "tan_points" | "recessive_black" | "n/a_phaeomelanin"
    confidence: str                 # "high" | "medium" | "low"
    abstains_on: list[str]          # appearance axes withheld (e.g. spotting/merle present)
    per_locus: dict                 # canonical genotype used per locus, for audit
    rule: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "organism": "Canis_lupus_familiaris", "trait": "coat_colour",
            "regime": "A_curated_catalog_epistatic", "rule": self.rule,
            "coat_color": self.coat_color, "pigment_type": self.pigment_type,
            "eumelanin_color": self.eumelanin_color, "distribution": self.distribution,
            "confidence": self.confidence, "abstains_on": self.abstains_on,
            "per_locus": self.per_locus, "notes": self.notes,
            "scope_limit": ("v0: five classic solid-colour loci (E/K/A/B/D); calls colour not shade/length/"
                            "spotting extent; reads per-locus allele CALLS not sequence"),
            "undetectable_mechanisms": list(UNSEEN_MECHANISMS),
        }


def _canon_allele(locus: str, tok: str) -> str:
    # Dog allele symbols are CASE-SIGNIFICANT (E vs e, B vs b, D vs d, K^B vs k^y): case encodes dominance,
    # so matching is case-SENSITIVE. Exact match first; aliases only for tolerant multi-char tokens.
    t = tok.strip()
    if t in LOCI[locus].alleles:
        return t
    low = t.lower()
    if low in _ALLELE_ALIASES and _ALLELE_ALIASES[low] in LOCI[locus].alleles:
        return _ALLELE_ALIASES[low]
    raise CoatInputError(
        f"unknown {locus}-locus allele {tok!r}; recognised (case-sensitive): {list(LOCI[locus].alleles)}")


def parse_genotype(locus: str, spec: str) -> tuple[str, str]:
    """Parse a diploid locus genotype string ('e/e', 'KB/ky', 'B/b') into a canonical (a1, a2) pair."""
    if locus not in LOCI:
        raise CoatInputError(f"unknown locus {locus!r}; v0 loci: {list(LOCI)}")
    parts = [p for p in spec.replace("|", "/").split("/") if p.strip()]
    if len(parts) != 2:
        raise CoatInputError(f"{locus} genotype {spec!r} is not diploid (expected a1/a2, e.g. E/e)")
    return _canon_allele(locus, parts[0]), _canon_allele(locus, parts[1])


def _dominant(locus: str, geno: tuple[str, str]) -> str:
    """Most-dominant allele present at a locus (per its dominance order)."""
    order = LOCI[locus].alleles
    return min(geno, key=lambda a: order.index(a))


def _is_hom_recessive_family(geno: tuple[str, str], family: set[str]) -> bool:
    return all(a in family for a in geno)


def call_coat_color(loci_genotypes: dict[str, str], present_loci: list[str] | None = None) -> CoatColorCall:
    """Deterministic coat-colour call from per-locus allele genotypes.

    `loci_genotypes`: {locus -> 'a1/a2'} for any subset of E/K/A/B/D (E, B, D minimally drive the colour;
    K + A drive distribution). `present_loci`: optional list of ADDITIONAL loci present on the dog that v0
    does not model (e.g. ['M','S']) — the affected appearance axis is added to `abstains_on`.
    """
    rule = "dog_coat_colour_epistatic_v0"
    notes: list[str] = []
    per: dict = {}
    geno: dict[str, tuple[str, str]] = {}
    for loc, spec in loci_genotypes.items():
        L = loc.strip().upper()
        if L not in LOCI:
            if L in UNMODELLED_LOCI:
                raise CoatInputError(
                    f"locus {loc!r} is a PATTERN locus not modelled in v0 ({UNMODELLED_LOCI[L]}); pass it via "
                    f"present_loci=[...] so the affected axis ABSTAINS instead of a wrong solid-colour call")
            raise CoatInputError(f"unknown locus {loc!r}; v0 loci: {list(LOCI)}")
        geno[L] = parse_genotype(L, spec)
        per[L] = "/".join(geno[L])

    # abstention axes from unmodelled pattern loci the caller declared present
    abstains: list[str] = []
    for pl in (present_loci or []):
        P = pl.strip().upper()
        if P in UNMODELLED_LOCI:
            abstains.append(f"{P}: {UNMODELLED_LOCI[P]}")
        elif P not in LOCI:
            raise CoatInputError(f"unknown present locus {pl!r}")

    # E locus REQUIRED — it is the top of the epistasis; without it we cannot know the pigment type
    if "E" not in geno:
        raise CoatInputError("E (MC1R) genotype is required — it is the pigment-type switch at the top of "
                             "the epistasis; without it the coat colour cannot be resolved")

    # ---- 1. pigment-type switch (E locus) ----
    e_recessive_red = geno["E"] == ("e", "e")
    if e_recessive_red:
        # phaeomelanin: red/yellow REGARDLESS of K/A/B. D can dilute red -> cream (shade; v0 notes it).
        notes.append("E locus e/e -> no coat eumelanin -> RED/YELLOW (phaeomelanin); this is EPISTATIC over "
                     "K/A/B — a naive 'K^B->black' or 'b/b->brown' rule would MIS-CALL this dog")
        dil = "D" in geno and _is_hom_recessive_family(geno["D"], _D_RECESSIVE)
        if dil:
            notes.append("d/d dilutes phaeomelanin toward cream (shade axis; v0 reports red/yellow)")
        if "B" in geno and _is_hom_recessive_family(geno["B"], _B_RECESSIVE):
            notes.append("b/b present but coat is red (E-masked); b/b still lightens the NOSE to brown")
        conf = "high" if abstains == [] else "medium"
        return CoatColorCall("red/yellow", "phaeomelanin", None, "n/a_phaeomelanin", conf,
                             abstains, per, rule, notes)

    # ---- 2. eumelanin colour from B x D ----
    bb = "B" in geno and _is_hom_recessive_family(geno["B"], _B_RECESSIVE)
    dd = "D" in geno and _is_hom_recessive_family(geno["D"], _D_RECESSIVE)
    if bb and dd:
        eu = "isabella/lilac"
    elif bb:
        eu = "brown/liver"
    elif dd:
        eu = "blue/grey"
    else:
        eu = "black"
    missing_color = [x for x in ("B", "D") if x not in geno]
    if missing_color:
        notes.append(f"eumelanin colour assumes wild-type at absent locus/loci {missing_color} "
                     f"(B->black default / D->non-dilute default); confidence capped medium")

    # ---- 3. distribution from K x A ----
    dist = "solid"
    if "K" not in geno:
        notes.append("K locus absent -> distribution assumed SOLID (K^B default); pass K to resolve agouti")
        dist = "solid"
    else:
        k_dom = _dominant("K", geno["K"])
        if k_dom == "KB":
            dist = "solid"
        elif k_dom == "kbr":
            dist = "brindle"
            abstains.append("K: kbr brindle — striped eumelanin/phaeomelanin, v0 abstains on the stripe pattern")
        else:  # ky/ky -> A expressed
            if "A" not in geno:
                notes.append("ky/ky expresses the A locus but A genotype absent -> distribution UNKNOWN")
                dist = "agouti_unknown"
                abstains.append("A: agouti pattern expressed (ky/ky) but A genotype not provided")
            else:
                a_dom = _dominant("A", geno["A"])
                dist = {"Ay": "sable", "aw": "agouti", "at": "tan_points", "a": "recessive_black"}[a_dom]

    # ---- compose the human label ----
    if dist == "solid":
        coat = f"solid {eu}"
    elif dist == "recessive_black":
        coat = f"recessive-black ({eu})"
    elif dist == "sable":
        coat = f"sable/fawn ({eu}-tipped)"
    elif dist == "agouti":
        coat = f"agouti/wolf-grey ({eu} base)"
    elif dist == "tan_points":
        coat = f"tan-points ({eu} with tan)"
    elif dist == "brindle":
        coat = f"brindle ({eu} striped on red)"
    else:  # agouti_unknown
        coat = f"{eu} eumelanin, agouti pattern undetermined"

    conf = "high"
    if missing_color or dist == "agouti_unknown":
        conf = "medium"
    if abstains:
        conf = "medium" if conf == "high" else conf
    notes.append(f"eumelanin-capable (E-) -> {eu} eumelanin, {dist} distribution")
    return CoatColorCall(coat, "eumelanin", eu, dist, conf, abstains, per, rule, notes)


def reference_integrity_ok() -> bool:
    """Biology contract guard — a corrupted catalog/rule fails this. Pins known breed genotypes -> colours,
    INCLUDING the E-locus epistasis anchor a naive has-the-allele rule gets wrong."""
    # Yellow Labrador: e/e -> red/yellow (phaeomelanin), regardless of B/D
    ylab = call_coat_color({"E": "e/e", "B": "B/B", "D": "D/D"})
    # Chocolate Labrador: E- , b/b -> brown/liver
    choc = call_coat_color({"E": "E/E", "K": "KB/KB", "B": "b/b", "D": "D/D"})
    # Weimaraner: E- , b/b , d/d -> isabella/lilac
    weim = call_coat_color({"E": "E/E", "K": "KB/KB", "B": "b/b", "D": "d/d"})
    # Solid black dog: E- , K^B , B- , D-
    blk = call_coat_color({"E": "E/E", "K": "KB/KB", "B": "B/B", "D": "D/D"})
    # Blue: E-, B-, d/d
    blue = call_coat_color({"E": "E/E", "K": "KB/KB", "B": "B/B", "D": "d/d"})
    # German-Shepherd-style black-and-tan: E-, ky/ky, a^t/a^t -> tan points
    gsd = call_coat_color({"E": "E/E", "K": "ky/ky", "A": "at/at", "B": "B/B", "D": "D/D"})
    # THE EPISTASIS ANCHOR: e/e AND K^B AND b/b -> STILL red/yellow (E epistatic; naive rule mis-calls)
    anchor = call_coat_color({"E": "e/e", "K": "KB/KB", "B": "b/b", "D": "D/D"})
    return (ylab.coat_color == "red/yellow" and ylab.pigment_type == "phaeomelanin"
            and choc.eumelanin_color == "brown/liver" and choc.distribution == "solid"
            and weim.eumelanin_color == "isabella/lilac"
            and blk.coat_color == "solid black"
            and blue.eumelanin_color == "blue/grey"
            and gsd.distribution == "tan_points"
            and anchor.coat_color == "red/yellow" and anchor.pigment_type == "phaeomelanin")
