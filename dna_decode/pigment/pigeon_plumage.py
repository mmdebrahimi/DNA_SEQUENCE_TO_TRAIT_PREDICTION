"""Pigeon (Columba livia) plumage-colour decoder — a 2nd BIRD cell (sibling of chicken_plumage), and one of
the best-characterised colour systems in any organism (the Shapiro-lab rock-pigeon model, Domyan 2014 /
Vickrey 2018). Like chicken it is Z-linked (bird ZW), so the FEMALE is hemizygous.

THE BIOLOGY (molecularly-confirmed causal genes):

    B  (TYRP1)    the classical B locus, Z-LINKED, sets the BASE colour by an allelic series (Domyan 2014
                  Curr Biol, VAAST): B^A ash-red (dominant) > B+ blue/black (wild-type, ancestral — Darwin's
                  "blue rock") > b brown (recessive). Blue = mostly black eumelanin; brown = brown eumelanin;
                  ash-red = mostly pheomelanin. Z-linked -> a MALE (ZZ) is homo/heterozygous, a FEMALE (ZW)
                  is HEMIZYGOUS.
    e  (SOX10)    RECESSIVE RED — e/e makes the bird red REGARDLESS of the B locus (EPISTATIC over TYRP1).
                  Domyan 2014.
    D  (SLC45A2)  DILUTE, Z-LINKED — d/d washes the colour out (ash-red -> ash-yellow, blue -> dun, brown ->
                  khaki). Domyan 2014.
    C  (NDP)      WING PATTERN, an allelic series (Vickrey 2018 eLife): C^T T-check (T-pattern) > C checker >
                  + bar (ancestral) > c barless (a start-codon NDP mutation; homozygotes have vision defects).

  Resolution order: recessive-red (SOX10 e/e, epistatic) -> B-locus base (Z-linked) -> dilute (Z-linked) ->
  wing pattern (NDP).

  THE EPISTASIS ANCHORS (the cases a naive has-the-allele rule mis-calls):
    (1) e/e (SOX10) -> RED regardless of the B/TYRP1 base (a naive "B+ -> blue" rule mis-calls a blue-genotype
        e/e bird as blue when it is red).
    (2) Z-LINKED B/dilute with the bird ZW system: a FEMALE (ZW) is HEMIZYGOUS (one Z allele), REVERSED from
        mammals (same as chicken).
  `reference_integrity_ok` pins exactly these.

WHY A BESPOKE (non-frozen) cell: a Z-linked TYRP1 base series + SOX10-recessive-red epistasis + NDP pattern
is a distinct locus set (not the mammalian A/B/C/D/E engine, not the chicken E/B/S set). Imports nothing from
the frozen decoder surface.

HONEST SCOPE (load-bearing):
  - v0 = the 4 molecularly-confirmed loci (B/e/D/C) + optional Spread (classical). Grizzle, almond, indigo,
    recessive-opal, and the many modifier loci ABSTAIN (via --present).
  - Input is PER-LOCUS allele CALLS. The Z-linked B and D loci take ONE allele for a FEMALE (ZW, e.g. B=B+)
    or TWO for a MALE (ZZ); --sex may be given but is inferred from the Z-locus zygosity.
  - Calls the base colour + dilute + wing pattern, NOT the fine modifier phenotypes or shade. KNOWLEDGE_BASELINE.
  - Faithful-to-literature: applies the published OMIA/Shapiro-lab loci; it is not a new model.

Pure-python, wheel-only, offline, deterministic. Benign livestock/hobby visible-trait genetics — NOT human/forensic.
"""
from __future__ import annotations

from dataclasses import dataclass, field

RULES_VERSION = "pigeon-plumage-v0.1.0"


@dataclass(frozen=True)
class Locus:
    name: str
    gene: str
    alleles: tuple[str, ...]        # dominance order, most-dominant first
    source: str
    z_linked: bool = False
    note: str = ""


LOCI: dict[str, Locus] = {
    "B": Locus("B", "TYRP1", ("BA", "B+", "b"),
               source="Domyan 2014 Curr Biol (TYRP1, VAAST): B^A ash-red > B+ blue/black (wild) > b brown",
               z_linked=True,
               note="Z-LINKED base colour; FEMALE (ZW) hemizygous, MALE (ZZ) 2 alleles"),
    "E": Locus("E", "SOX10", ("E+", "e"),
               source="Domyan 2014 (SOX10): e/e recessive red, EPISTATIC over B/TYRP1",
               note="e/e -> red regardless of the B-locus base colour"),
    "D": Locus("D", "SLC45A2", ("D", "d"),
               source="Domyan 2014 (SLC45A2): d/d dilute (ash-red->ash-yellow, blue->dun, brown->khaki)",
               z_linked=True,
               note="Z-LINKED dilution"),
    "C": Locus("C", "NDP", ("CT", "C", "+", "c"),
               source="Vickrey 2018 eLife (NDP): C^T T-check > C checker > + bar (ancestral) > c barless",
               note="wing pattern series; barless c = NDP start-codon mutation (vision defects in homozygotes)"),
}

UNMODELLED_LOCI = {
    "S": "Spread — distributes pigment uniformly (solid/self); classical locus, gene not molecularly pinned in v0",
    "G": "grizzle — white-flecking modifier",
    "AL": "almond (St) — sex-linked mottling modifier",
    "IN": "indigo / recessive-opal / other dilution-family modifiers",
}

_ALLELE_ALIASES = {"ba": "BA", "b+": "B+", "e+": "E+", "ct": "CT"}


class PigeonInputError(ValueError):
    """Unknown locus/allele or malformed genotype (never a silent wrong call)."""


@dataclass
class PigeonCall:
    plumage: str                    # composed human label
    base_color: str                 # ash-red | blue/black | brown | red(recessive)
    sex_basis: str
    dilute: bool
    wing_pattern: str               # T-check | checker | bar | barless | n/a
    is_recessive_red: bool
    confidence: str
    abstains_on: list[str]
    per_locus: dict
    rule: str
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "organism": "Columba_livia", "trait": "plumage_colour",
            "regime": "A_curated_catalog_epistatic_zlinked", "rule": self.rule,
            "plumage": self.plumage, "base_color": self.base_color, "sex_basis": self.sex_basis,
            "dilute": self.dilute, "wing_pattern": self.wing_pattern, "is_recessive_red": self.is_recessive_red,
            "confidence": self.confidence, "abstains_on": self.abstains_on, "per_locus": self.per_locus,
            "notes": self.notes,
            "scope_limit": "v0: B/e/D/C molecularly-confirmed loci; calls base+dilute+wing-pattern, NOT modifiers/shade",
            "evidence_tier": "knowledge_baseline (curated OMIA/Shapiro-lab catalog; no free per-individual validation substrate)",
        }


def _canon_allele(locus: str, tok: str) -> str:
    t = tok.strip()
    if t in LOCI[locus].alleles:
        return t
    if t.lower() in _ALLELE_ALIASES and _ALLELE_ALIASES[t.lower()] in LOCI[locus].alleles:
        return _ALLELE_ALIASES[t.lower()]
    raise PigeonInputError(f"unknown {locus}-locus allele {tok!r}; recognised: {list(LOCI[locus].alleles)}")


def _parse(locus: str, spec: str, z_linked: bool = False) -> tuple[str, ...]:
    parts = [p for p in spec.replace("|", "/").split("/") if p.strip()]
    if z_linked and len(parts) == 1:
        return (_canon_allele(locus, parts[0]),)
    if len(parts) != 2:
        exp = "1 (female ZW) or 2 (male ZZ) alleles" if z_linked else "a1/a2"
        raise PigeonInputError(f"{locus} genotype {spec!r} must be {exp}")
    return _canon_allele(locus, parts[0]), _canon_allele(locus, parts[1])


def _dominant(locus: str, geno: tuple[str, ...]) -> str:
    order = LOCI[locus].alleles
    return min(geno, key=lambda a: order.index(a))


def call_pigeon_plumage(loci_genotypes: dict[str, str], sex: str | None = None,
                        present_loci: list[str] | None = None) -> PigeonCall:
    """Deterministic pigeon plumage-colour call from per-locus allele genotypes.

    `loci_genotypes`: {locus -> 'a1/a2'} for any subset of B/E/D/C. The Z-LINKED B and D loci take ONE allele
    for a FEMALE (ZW hemizygous, e.g. 'B+') or TWO for a MALE (ZZ). `sex` optional — inferred from the Z-locus
    zygosity. `present_loci`: v0-unmodelled loci present -> the affected axis ABSTAINS.
    """
    rule = f"pigeon_plumage_epistatic_v0 ({RULES_VERSION})"
    notes: list[str] = []
    per: dict = {}
    geno: dict[str, tuple[str, ...]] = {}
    for loc, spec in loci_genotypes.items():
        L = loc.strip().upper() if loc.strip().upper() in LOCI else loc.strip()
        if L not in LOCI:
            if L in UNMODELLED_LOCI:
                raise PigeonInputError(f"locus {loc!r} is a v0-unmodelled locus ({UNMODELLED_LOCI[L]}); pass it "
                                       f"via present_loci=[...] so the affected axis ABSTAINS")
            raise PigeonInputError(f"unknown locus {loc!r}; v0 loci: {list(LOCI)}")
        geno[L] = _parse(L, spec, z_linked=LOCI[L].z_linked)
        per[L] = "/".join(geno[L])

    abstains: list[str] = []
    for pl in (present_loci or []):
        P = pl.strip().upper()
        if P in UNMODELLED_LOCI:
            abstains.append(f"{P}: {UNMODELLED_LOCI[P]}")
        elif P not in LOCI:
            raise PigeonInputError(f"unknown present locus {pl!r}")

    # ---- sex basis (Z-linked: FEMALE hemizygous, REVERSED from mammals; same as chicken) ----
    z_geno = geno.get("B") or geno.get("D")
    if sex:
        sex = sex.strip().lower()
        if sex not in ("male", "female"):
            raise PigeonInputError("--sex must be 'male' or 'female'")
        sex_basis = f"{sex} ({'ZW hemizygous' if sex == 'female' else 'ZZ'})"
    elif z_geno is not None:
        sex = "female" if len(z_geno) == 1 else "male"
        sex_basis = f"{sex} (inferred from Z-locus zygosity: {len(z_geno)} allele{'s' if len(z_geno) == 2 else ''}; birds are ZW)"
    else:
        sex = None
        sex_basis = "unspecified (no Z-linked B/D genotype / --sex)"

    dilute = "D" in geno and geno["D"] == ("d",) or (("D" in geno) and geno["D"] == ("d", "d"))

    # ---- wing pattern (NDP) ----
    wing = "n/a"
    if "C" in geno:
        cdom = _dominant("C", geno["C"])
        wing = {"CT": "T-check", "C": "checker", "+": "bar", "c": "barless"}[cdom]
        if geno["C"] == ("c", "c"):
            notes.append("barless c/c (NDP start-codon mutation) — homozygotes have an increased incidence of vision defects")

    # ---- 1. recessive red (SOX10 e/e) — epistatic over B ----
    if "E" in geno and geno["E"] == ("e", "e"):
        notes.append("e/e (SOX10) -> RECESSIVE RED, EPISTATIC over the B/TYRP1 base (a naive 'B -> blue' rule "
                     "mis-calls this bird)")
        base = "red (recessive)"
        core = "recessive red" + (" (dilute -> yellow)" if dilute else "")
        return PigeonCall(f"{core}, {wing}" if wing != "n/a" else core, base, sex_basis, dilute, wing, True,
                          "high", abstains, per, rule, notes)

    # ---- 2. B-locus base (Z-linked TYRP1) ----
    if "B" in geno:
        bdom = _dominant("B", geno["B"])
        base = {"BA": "ash-red", "B+": "blue/black", "b": "brown"}[bdom]
    else:
        base = "blue/black"
        notes.append("B (TYRP1) absent -> base assumed blue/black (wild-type ancestral); pass B to resolve")

    # ---- 3. dilute (SLC45A2) ----
    core = base
    if dilute:
        core = {"ash-red": "ash-yellow", "blue/black": "dun", "brown": "khaki"}[base]

    # ---- compose ----
    plumage = f"{core}" + (f" {wing}" if wing not in ("n/a", "bar") else (" bar" if wing == "bar" else ""))
    conf = "high"
    if abstains or sex_basis.startswith("unspecified"):
        conf = "medium"
    return PigeonCall(plumage.strip(), base, sex_basis, dilute, wing, False, conf, abstains, per, rule, notes)


def reference_integrity_ok() -> bool:
    """Biology contract guard — pins known pigeon genotypes -> plumage, INCLUDING the anchors a naive rule gets
    wrong (SOX10 e/e recessive-red epistatic over TYRP1; Z-linked reversed-hemizygous B locus)."""
    # ANCHOR 1: e/e recessive red overrides a blue B-genotype
    rr = call_pigeon_plumage({"B": "B+/B+", "E": "e/e"})
    # base colours
    ash = call_pigeon_plumage({"B": "BA/B+"})               # ash-red dominant
    blue = call_pigeon_plumage({"B": "B+/b"})               # blue (brown recessive, hidden)
    brown = call_pigeon_plumage({"B": "b/b"})               # brown
    # ANCHOR 2: Z-linked B, a FEMALE is hemizygous (1 allele)
    hen = call_pigeon_plumage({"B": "BA"})                  # 1 allele -> female ash-red
    cock = call_pigeon_plumage({"B": "BA/B+"})              # 2 alleles -> male
    dun = call_pigeon_plumage({"B": "B+/B+", "D": "d/d"})   # dilute blue -> dun
    barless = call_pigeon_plumage({"B": "B+/B+", "C": "c/c"})
    return (rr.is_recessive_red and rr.base_color == "red (recessive)"
            and ash.base_color == "ash-red" and blue.base_color == "blue/black" and brown.base_color == "brown"
            and "female" in hen.sex_basis and "male" in cock.sex_basis
            and dun.dilute and "dun" in dun.plumage
            and barless.wing_pattern == "barless")
