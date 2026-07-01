"""ABO blood-group O-status decoder (M4 — the serological Mendelian cell).

Extends the deterministic gene->trait decoder from pigmentation to a SEROLOGICAL trait. The O blood group is
caused by rs8176719 = a single-guanine deletion (c.261delG, Yamamoto et al. 1990 Nature 345:229) that
frameshifts the ABO glycosyltransferase -> non-functional enzyme -> O antigen. HOMOZYGOUS DELETION = blood
type O; any functional (insertion) allele = non-O (A / B / AB). This is textbook + unambiguous.

SCOPE (honest, load-bearing): this calls O vs non-O ONLY. The A-vs-B distinction among non-O is driven by
tag SNPs (rs8176746 / rs8176747) whose exact allele->A/B mapping + 23andMe strand orientation needs careful
sourcing; NOT fabricated here -> deferred (named). The O/non-O rule is the sourced, unambiguous core.

23andMe reports rs8176719 as I/D (Insertion/Deletion): D = the deletion (O allele), I = insertion
(functional). Empirically confirmed in the PGP cached files (values DD/DI/II + '--' no-calls). To avoid the
'--' no-call vs '--' homozygous-deletion collision, O is called ONLY from the unambiguous 'DD'; '--' -> INDETERMINATE.
"""
from __future__ import annotations

INDETERMINATE = "INDETERMINATE"


def call_abo_o_status(rs8176719_genotype: str) -> str:
    """O / non-O / INDETERMINATE from rs8176719 (I/D coding; also tolerates -/G ACGT-indel coding).

    'DD' (homozygous deletion) -> O. Any functional allele present ('I' or 'G') -> non-O. No-call ('--',
    single/empty) -> INDETERMINATE (never guessed)."""
    g = (rs8176719_genotype or "").strip().upper()
    if g == "DD":
        return "O"                      # homozygous c.261delG deletion
    if "I" in g or "G" in g:
        return "non-O"                  # II / DI / ID / GG / -G / G- : >=1 functional allele
    return INDETERMINATE                # '--', '', single-allele: no-call


def bin_blood_type(raw: str) -> str | None:
    """Bin a self-reported blood type ('A +', 'O -', 'AB +', "Don't know") -> 'O' / 'non-O' / None.

    None for unknown/blank (never scored)."""
    if not raw:
        return None
    s = raw.strip().lower()
    if not s or "don't know" in s or "dont know" in s or "unknown" in s or s in ("-", "n/a"):
        return None
    tok = s.replace("+", " ").replace("-", " ").split()
    t = (tok[0] if tok else s).upper()
    if t == "O":
        return "O"
    if t in ("A", "B", "AB"):
        return "non-O"
    return None
