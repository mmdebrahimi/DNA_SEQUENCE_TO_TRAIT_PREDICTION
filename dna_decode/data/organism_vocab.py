"""Canonical organism vocabulary for the closed-set identification router.

NON-FROZEN. Touches no sha256-pinned file; this is an additive translation layer and renames nothing
in `cell_registry`.

WHY THIS EXISTS
---------------
The identification router has to emit an organism token that a DOWNSTREAM consumer accepts. Two
candidate vocabularies exist in this repo and they look independent:

  * AMRFinder's `-O` allow-list   -- `Escherichia`, `Klebsiella_pneumoniae`, ... (32 values on DB
    2026-03-24.1, printed by `amrfinder -l`; that output is the AUTHORITATIVE list and is PARSED in the
    test rather than copied here, because a hand-copied allow-list drifts the moment the DB bumps).
  * the registry/rule tokens      -- `Escherichia_coli_Shigella`, `Klebsiella`, ...

and `cell_registry` is genuinely inconsistent about the second one. MEASURED over the 115 registered
cells: E. coli appears as THREE tokens (`Escherichia_coli_Shigella` on 4 AMR cells, `Escherichia_coli`
on 3 typing cells, `escherichia_coli` on 3 metabolic cells) and Klebsiella as TWO (`Klebsiella` on 6 AMR
cells, `Klebsiella_pneumoniae` on 1 typing cell).

THE PLAN'S PREMISE WAS REFUTED DURING EXECUTION, AND THIS MODULE IS THE SMALLER THING THAT SURVIVED
---------------------------------------------------------------------------------------------------
The plan asserted that "a router emitting one string satisfies neither consumer", so this module was
specified to carry an `amrfinder_organism` AND a separate `registry_organism` per organism. That is
wrong. `amr_rules.calibrated_rule_for` matches (1) case-insensitive exact, then (2) GENUS-prefix on the
first `_`-delimited token of BOTH sides -- a fallback its own docstring says exists precisely so that
"an AMRFinder `-O` value like 'Klebsiella_pneumoniae' ... resolves the genus-level registry key
('Klebsiella')".

VERIFIED EMPIRICALLY, not taken from that docstring (the meropenem rule's docstring claimed a blind spot
it does not have, measured the same day): across `Klebsiella_pneumoniae` vs `Klebsiella`, `Escherichia`
vs `Escherichia_coli_Shigella`, and `Acinetobacter_baumannii` vs `Acinetobacter`, over both cipro and
meropenem, `calibrated_rule_for` returns the IDENTICAL result in all 6 cases -- in both the
rule-resolved and the None directions.

So **ONE token suffices: the AMRFinder `-O` value.** It is what AMRFinder takes directly and what
`call_resistance(organism=...)` resolves through. `registry_labels` below is therefore descriptive
(what `cell_registry` happens to call this organism) and is NOT a second routing key. The equivalence
this relies on is PINNED BY TEST, so if anyone removes the genus fallback the dependency fails loudly
instead of silently routing to a default rule.

SCOPE
-----
Bacterial + fungal only, matching the v0 router slice. Viral (HIV / SARS-CoV-2 / HCMV) and human
clinical cells are routed by reference-CDS BLAST per cell, not by genome sketch, and are deliberately
absent -- see the plan's Open Question 3.

`routing_token` returns None for organisms AMRFinder does not cover at all (there is no
AMRFinder-for-fungi, and M. tuberculosis is absent from the `-O` list). None means "this organism is
real and supported, but it has no AMRFinder routing value" -- it is NOT "unknown", which raises.
"""
from __future__ import annotations

from dataclasses import dataclass, field


class UnknownOrganism(KeyError):
    """Raised by `canonical()` on a token this vocabulary does not know.

    Deliberately an exception and never a default. A silent default is exactly how a WRONG organism
    would reach a rule, and the organism argument selects which rule runs.
    """


@dataclass(frozen=True)
class OrganismEntry:
    """One supported organism.

    amrfinder_organism: the `-O` value, and the SINGLE routing token (see module docstring). None where
        AMRFinder has no entry for this organism at all.
    registry_labels: what `cell_registry` calls this organism. DESCRIPTIVE ONLY -- not a routing key.
        A tuple with >1 element records a real inconsistency in the registry, not an error here.
    """

    canonical: str
    amrfinder_organism: str | None
    registry_labels: tuple[str, ...] = field(default_factory=tuple)


# Every `amrfinder_organism` below is asserted a member of the live `amrfinder -l` output by
# tests/test_organism_vocab.py. Do NOT add a value here without that test passing against the image.
_ENTRIES: tuple[OrganismEntry, ...] = (
    OrganismEntry("escherichia_coli", "Escherichia",
                  ("Escherichia_coli_Shigella", "Escherichia_coli", "escherichia_coli")),
    OrganismEntry("klebsiella_pneumoniae", "Klebsiella_pneumoniae",
                  ("Klebsiella", "Klebsiella_pneumoniae")),
    OrganismEntry("campylobacter", "Campylobacter", ("Campylobacter",)),
    OrganismEntry("salmonella", "Salmonella", ("Salmonella",)),
    OrganismEntry("pseudomonas_aeruginosa", "Pseudomonas_aeruginosa", ("Pseudomonas_aeruginosa",)),
    OrganismEntry("staphylococcus_aureus", "Staphylococcus_aureus", ("Staphylococcus_aureus",)),
    OrganismEntry("acinetobacter_baumannii", "Acinetobacter_baumannii", ("Acinetobacter",)),
    OrganismEntry("streptococcus_pneumoniae", "Streptococcus_pneumoniae", ("Streptococcus_pneumoniae",)),
    OrganismEntry("neisseria_gonorrhoeae", "Neisseria_gonorrhoeae", ()),
    OrganismEntry("enterococcus_faecium", "Enterococcus_faecium", ()),
    OrganismEntry("enterobacter_cloacae", "Enterobacter_cloacae", ()),
    # AMRFinder is bacterial: there is no AMRFinder-for-fungi, so C. auris routes to the BLAST
    # ERG11/FKS1 target-site engine instead (dna_decode/data/fungal_amr.py).
    OrganismEntry("candida_auris", None, ("Candida_auris",)),
    # M. tuberculosis is absent from the `-O` allow-list; TB routes through organism_rules/tb_amr.
    OrganismEntry("mycobacterium_tuberculosis", None, ()),
)

BY_CANONICAL: dict[str, OrganismEntry] = {e.canonical: e for e in _ENTRIES}


def _alias_index() -> dict[str, str]:
    """Build the alias -> canonical index, refusing to let two organisms claim one alias.

    An alias collision would make `canonical()` answer by dict order, which is the documented
    tie-breaking trap. Raise at import rather than resolve it arbitrarily.
    """
    idx: dict[str, str] = {}
    for e in _ENTRIES:
        aliases = {e.canonical, *e.registry_labels}
        if e.amrfinder_organism:
            aliases.add(e.amrfinder_organism)
        for a in aliases:
            key = a.lower()
            prev = idx.get(key)
            if prev is not None and prev != e.canonical:
                raise ValueError(
                    f"alias {a!r} claimed by both {prev!r} and {e.canonical!r}; "
                    "resolve it explicitly rather than letting dict order decide")
            idx[key] = e.canonical
    return idx


_ALIASES: dict[str, str] = _alias_index()


def canonical(token: str) -> str:
    """Normalize any known organism spelling to its canonical key.

    Case-insensitive. Accepts a canonical key, an AMRFinder `-O` value, or any registry label.
    RAISES `UnknownOrganism` on anything else -- never defaults.
    """
    if not isinstance(token, str) or not token.strip():
        raise UnknownOrganism(f"empty organism token: {token!r}")
    got = _ALIASES.get(token.strip().lower())
    if got is None:
        raise UnknownOrganism(
            f"{token!r} is not a supported organism. Supported canonical keys: "
            f"{sorted(BY_CANONICAL)}")
    return got


def routing_token(token: str) -> str | None:
    """The SINGLE token both consumers accept: the AMRFinder `-O` value.

    None means this organism is supported but AMRFinder has no `-O` entry for it (fungal / TB), so a
    caller must route it to that organism's own engine rather than to AMRFinder.
    Raises `UnknownOrganism` for an unsupported token -- None and unknown are different answers.
    """
    return BY_CANONICAL[canonical(token)].amrfinder_organism


def supported_canonicals() -> tuple[str, ...]:
    """Every canonical key, sorted. This is the router's CLOSED SET."""
    return tuple(sorted(BY_CANONICAL))


def amrfinder_routable() -> tuple[str, ...]:
    """Canonical keys that have an AMRFinder `-O` value, sorted."""
    return tuple(sorted(c for c, e in BY_CANONICAL.items() if e.amrfinder_organism))
