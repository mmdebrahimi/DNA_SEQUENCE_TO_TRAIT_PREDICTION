"""Call FRI functional status from VARIANTS -- the v0.1 the flowering cell's docstring deferred.

v0 took FRI status as an INPUT (an allele name, or the source table's own call). That makes the cell a
lookup of someone else's answer. This module calls the status from the variant table itself, which is what
a decoder does -- and it makes Zhang 2020's Table S3 a VALIDATION TARGET rather than an input.

TWO RULES, and the comparison between them is the point:

  `putative`  -- the source's own heuristic: any lesion whose annotated consequence is loss-of-function
                 (frameshift / stopgain / stoploss / large indel / the Col + Ler deletions) => lof.
                 This is what Table S3's `deleterious_allele` column encodes.

  `curated`   -- putative PLUS the substitutions Zhang 2020 EXPERIMENTALLY PROVED non-functional in a
                 common null background: L276R and L294F (their Figs 5-6; transgenic + GFP localisation +
                 western blot). These are ordinary `nonsynonymous SNV` rows, so the putative rule cannot
                 see them -- and the paper's own summary column calls their carriers FUNCTIONAL, which
                 contradicts the paper's own experiments.

**We are not adding new biology.** Both rules are faithful-to-literature; `curated` is simply a more
complete reading of the SAME paper. The claim under test is only whether encoding the paper's experimental
results beats its putative-LoF shorthand on measured phenotype -- our own wrapper-vs-underlying-tool
discipline: a policy layer over a source must be validated against NAIVE use of that source.

ABSTENTION IS LOAD-BEARING. Table S1 marks a no-call as `.`. If every LoF-bearing position in an accession
is `ref` we can say functional; but if any is `.`, absence of an ALT is NOT evidence of absence of a lesion
-- that accession ABSTAINS. (Same discipline as the TB cell's uncallable-window ABSTAIN vs S-by-absence.)

Pure-python, offline, deterministic. Reads the committed CC-BY Table S1.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# --- the source's own consequence vocabulary (Table S1 `funct_consequ`) ------------------------------------
# Reconciles EXACTLY to the paper's stated 26 LoF mutations: "4 large indels, 14 frameshift indels, 7
# premature stop codons, and 1 deletion of the stop codon" == large_indel(3) + `Ler deletion`(1) = 4 large;
# frameshift deletion(9) + frameshift insertion(4) + `Col deletion`(1, a 16-bp frameshift) = 14; stopgain(7);
# stoploss(1). That reconciliation is asserted in tests -- a vocabulary drift fails loudly.
LOF_CONSEQUENCES: frozenset[str] = frozenset({
    "large_indel", "Ler deletion", "Col deletion",
    "frameshift deletion", "frameshift insertion",
    "stopgain", "stoploss",
})

# Consequences that leave the protein full-length. NOT loss-of-function *by annotation* -- but see
# VERIFIED_LOF_SUBSTITUTIONS: two of these were proved non-functional by experiment.
NON_LOF_CONSEQUENCES: frozenset[str] = frozenset({
    "nonsynonymous SNV", "synonymous SNV", "nonframeshift deletion", "nonframeshift insertion",
})

# Substitutions EXPERIMENTALLY proved to abolish FRI function (Zhang & Jimenez-Gomez 2020):
#   L276R -- prevents nuclear localisation (Fig 5: cytoplasmic-only GFP; the isolated mutation is sufficient)
#   L294F -- destabilises the protein (Fig 6: no/weak GFP + western blot; degraded)
# E302G is deliberately EXCLUDED: it co-occurs with L276R in FRI-Ba1-2, and the paper tested it in
# isolation and found it still delays flowering (i.e. FUNCTIONAL). Including it would be the naive
# "every central-domain change is bad" error the paper explicitly disproves.
VERIFIED_LOF_SUBSTITUTIONS: frozenset[str] = frozenset({"L276R", "L294F"})

RULES = ("putative", "curated")

# Table S1 genotype-matrix cell values.
_ALT, _REF, _NOCALL, _HET = "alt", "ref", ".", "het"


class FriCallerError(ValueError):
    """Raised on an unknown consequence / unknown rule -- never a silent wrong call."""


@dataclass
class FriCall:
    status: str                      # "lof" | "functional" | "unknown" (-> the cell ABSTAINs)
    rule: str
    evidence: list[str] = field(default_factory=list)   # the lesions that fired
    n_nocall_lof_positions: int = 0
    note: str = ""


@dataclass(frozen=True)
class Variant:
    """One Table S1 row, reduced to what the call needs. Frozen so it can key a genotype map."""
    cds_pos: int
    prot_change: str                 # e.g. "L294F", "M1fs", "S316X"
    consequence: str

    @property
    def is_putative_lof(self) -> bool:
        if self.consequence in LOF_CONSEQUENCES:
            return True
        if self.consequence in NON_LOF_CONSEQUENCES:
            return False
        raise FriCallerError(
            f"unknown funct_consequ {self.consequence!r} at CDS {self.cds_pos} -- refusing to guess whether "
            f"it is loss-of-function. Add it to LOF_CONSEQUENCES or NON_LOF_CONSEQUENCES.")

    def is_lof_under(self, rule: str) -> bool:
        if rule not in RULES:
            raise FriCallerError(f"unknown rule {rule!r}; expected one of {RULES}")
        if self.is_putative_lof:
            return True
        return rule == "curated" and self.prot_change in VERIFIED_LOF_SUBSTITUTIONS


def call_fri_from_variants(genotypes: dict[Variant, str], rule: str = "curated") -> FriCall:
    """Call FRI status for ONE accession. `genotypes` maps each variant -> 'alt'|'ref'|'.'|'het'.

    lof         iff the accession carries >=1 ALT at a lesion the rule counts as loss-of-function
    unknown     iff no such ALT, but >=1 counted lesion is a NO-CALL -- absence of an ALT at an uncalled
                position is not evidence of absence (the accession ABSTAINs rather than be called functional)
    functional  iff every counted lesion is confidently REF
    """
    if rule not in RULES:
        raise FriCallerError(f"unknown rule {rule!r}; expected one of {RULES}")

    hits, nocall = [], 0
    for v, gt in genotypes.items():
        if not v.is_lof_under(rule):
            continue
        if gt in (_ALT, _HET):
            hits.append(f"{v.prot_change} ({v.consequence}"
                        f"{', heterozygous' if gt == _HET else ''})")
        elif gt == _NOCALL:
            nocall += 1

    if hits:
        return FriCall("lof", rule, sorted(hits), nocall,
                       note="carries a loss-of-function lesion" +
                            (" (experimentally-verified substitution included)" if rule == "curated" else ""))
    if nocall:
        return FriCall("unknown", rule, [], nocall,
                       note=f"{nocall} loss-of-function position(s) are NO-CALL in this accession -- absence "
                            f"of an ALT there is not evidence of absence; abstaining rather than calling "
                            f"functional by default")
    return FriCall("functional", rule, [], 0, note="every loss-of-function position is confidently reference")


def load_variants(xlsx_path: Path) -> tuple[list[Variant], list[str], list[list[str]]]:
    """Parse Table S1 -> (variants, accession_ids, genotype matrix[variant][accession]).

    Needs openpyxl (an optional extra); the caller logic above is importable without it.
    """
    import openpyxl  # local import: keeps the pure caller wheel-only

    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    rows = list(wb["Table_S1"].iter_rows(values_only=True))
    header, body = rows[0], rows[1:]
    accessions = [str(c) for c in header[5:]]
    variants, matrix = [], []
    for r in body:
        variants.append(Variant(cds_pos=int(r[0]), prot_change=str(r[1]), consequence=str(r[4])))
        matrix.append([str(c) for c in r[5:]])
    return variants, accessions, matrix


def call_all_accessions(xlsx_path: Path, rule: str = "curated") -> dict[str, FriCall]:
    """Call every accession in Table S1 under `rule`."""
    variants, accessions, matrix = load_variants(xlsx_path)
    out: dict[str, FriCall] = {}
    for j, acc in enumerate(accessions):
        out[acc] = call_fri_from_variants({v: matrix[i][j] for i, v in enumerate(variants)}, rule)
    return out


def reference_integrity_ok(xlsx_path: Path | None = None) -> bool:
    """Contract guard. Pins the paper's OWN stated LoF arithmetic (4 large + 14 frameshift + 7 stopgain +
    1 stoploss = 26) and the two verified substitutions -- catches a vocabulary or catalog corruption."""
    if not VERIFIED_LOF_SUBSTITUTIONS == {"L276R", "L294F"}:
        return False
    if "E302G" in VERIFIED_LOF_SUBSTITUTIONS:      # the paper proved E302G is FUNCTIONAL in isolation
        return False
    if LOF_CONSEQUENCES & NON_LOF_CONSEQUENCES:    # a consequence cannot be both
        return False
    if xlsx_path is None or not xlsx_path.exists():
        return True
    variants, _, _ = load_variants(xlsx_path)
    return len(variants) == 171 and sum(1 for v in variants if v.is_putative_lof) == 26
