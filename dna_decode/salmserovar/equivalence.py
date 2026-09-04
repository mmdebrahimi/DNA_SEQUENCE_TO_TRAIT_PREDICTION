"""Decide when two Salmonella serovar strings mean the same organism.

WHY THIS IS NEEDED AND WHY IT IS NOT A FUDGE. Scoring a serovar caller by exact string match punishes
NOTATION rather than biology. All of these name the same organism:

    'Typhimurium'  ==  'Typhimurium var. 5-'  ==  'I 4,[5],12:i:1,2'
    '4,[5],12:i:-' ==  'I 4,[5],12:i:-'       ==  'Typhimurium - monophasic'

while 'Johannesburg' vs 'Cubana' is a real miss. A comparison that cannot tell those two cases apart
does not measure the caller; it measures string hygiene.

THE RULE IS DELIBERATELY CONSERVATIVE. Equivalence is granted ONLY through (a) documented notation
normalisation, or (b) the committed White-Kauffmann-Le Minor formula table -- the same table the
deployed caller uses. Nothing is matched on fuzzy string similarity, because a fuzzy matcher would
quietly grant credit for near-misses like 'Newport'/'Newbrunswick'. When equivalence cannot be
established the answer is NOT_EQUIVALENT, and a miss is recorded.

SYMMETRY IS THE INTEGRITY PROPERTY. The same function scores our caller and the in-silico comparator.
Any leniency it contains is applied to both sides, so a lenient rule cannot flatter one of them.
"""
from __future__ import annotations

import re
from pathlib import Path

# 'Typhimurium var. 5-' / 'var. Copenhagen' are VARIANTS of a serovar, not different serovars: they
# denote antigen presence/absence within the same named serovar.
_VARIANT_SUFFIX = re.compile(
    r"\s*(var\.?|variant|var)\s*[a-z0-9\-\+\[\],\. ]*$", re.IGNORECASE)
# Subspecies prefixes ('I ', 'I 4,[5]...' / 'subsp. enterica') carried on a formula.
_SUBSP_PREFIX = re.compile(r"^(i{1,3}[ab]?|iv|v|vi)\s+", re.IGNORECASE)
# Monophasic Typhimurium is written a dozen ways; it IS 4,[5],12:i:-.
_MONOPHASIC = re.compile(r"typhimurium\s*[-,]?\s*(monophasic|mono)", re.IGNORECASE)


def _strip_brackets(s: str) -> str:
    """'4,[5],12:i:-' -> '4,5,12:i:-'. Brackets mark an OPTIONAL antigen; presence or absence of the
    bracket is a notation choice, not a different antigenic formula."""
    return s.replace("[", "").replace("]", "")


def canonical(s: str) -> str:
    """Notation-normalised form. Never merges two distinct serovar NAMES."""
    t = (s or "").strip().strip('"').lower()
    for cut in ("serovar ", "ser. ", "serotype "):
        if cut in t:
            t = t.split(cut, 1)[1]
    for pre in ("salmonella enterica subsp. enterica ", "salmonella enterica subspecies enterica ",
                "salmonella enterica ", "salmonella "):
        if t.startswith(pre):
            t = t[len(pre):]
    if _MONOPHASIC.search(t):
        return "4,5,12:i:-"
    t = _VARIANT_SUFFIX.sub("", t)
    t = _SUBSP_PREFIX.sub("", t)
    t = _strip_brackets(t)
    t = " ".join(t.replace("_", " ").split()).strip(" .,;")
    return t


def is_formula(s: str) -> bool:
    """An antigenic formula has the O:H1:H2 shape ('4,5,12:i:-'); a serovar name does not."""
    return s.count(":") >= 2


def load_formula_index(table_path: Path) -> dict[str, str]:
    """Map canonical antigenic formula -> canonical serovar name, from the committed W-K-L table."""
    idx: dict[str, str] = {}
    with open(table_path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:4] != ["O", "H1", "H2", "Serovar"]:
            raise ValueError(f"unexpected serovar table header: {header[:4]}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            o, h1, h2, name = f[0].strip(), f[1].strip(), f[2].strip(), f[3].strip()
            if not name:
                continue
            idx[canonical(f"{o}:{h1}:{h2}")] = canonical(name)
    return idx


def equivalent(a: str, b: str, formula_index: dict[str, str] | None = None) -> tuple[bool, str]:
    """Do these two serovar strings denote the same organism? Returns (verdict, reason).

    The reason is returned so a scored artifact can be audited row by row -- a bare boolean would hide
    WHICH rule granted the match, and a leniency nobody can see is a leniency nobody can check.
    """
    ca, cb = canonical(a), canonical(b)
    if not ca or not cb:
        return False, "empty"
    if ca == cb:
        return True, "exact-after-notation-normalisation"

    if formula_index:
        # Resolve a formula to its serovar name on either side, then re-compare.
        ra = formula_index.get(ca, ca) if is_formula(ca) else ca
        rb = formula_index.get(cb, cb) if is_formula(cb) else cb
        if ra == rb:
            return True, "formula-resolved-via-white-kauffmann-table"
    return False, "not-equivalent"
