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
    genus_match: genus tokens for which ANY species routes to this entry.

    GRANULARITY FOLLOWS AMRFINDER'S OWN LIST, which deliberately mixes the two levels: `Escherichia`,
    `Campylobacter` and `Salmonella` are GENUS-level `-O` values, while `Klebsiella_pneumoniae` and
    `Klebsiella_oxytoca` are SEPARATE species-level values. So a genus-level entry declares
    `genus_match` and accepts any species in it; a species-level entry leaves it empty and requires an
    exact match.

    This is not cosmetic -- it was MEASURED. Scanning 1,865 labelled local genomes, a species-strict
    reading lost 160 real in-set genomes (Campylobacter_jejuni 66, Salmonella_enterica 60,
    Campylobacter_coli 34) while correctly excluding 49 genuinely different species (Klebsiella
    aerogenes/variicola/michiganensis, Gemmata). An over-broad genus rule for Klebsiella would have
    swept K. aerogenes into K. pneumoniae's reference, which is the opposite error.
    """

    canonical: str
    amrfinder_organism: str | None
    registry_labels: tuple[str, ...] = field(default_factory=tuple)
    genus_match: tuple[str, ...] = field(default_factory=tuple)


# Every `amrfinder_organism` below is asserted a member of the live `amrfinder -l` output by
# tests/test_organism_vocab.py. Do NOT add a value here without that test passing against the image.
_ENTRIES: tuple[OrganismEntry, ...] = (
    # GENUS-level: AMRFinder's `-O Escherichia` covers E. coli AND Shigella, which is also why the
    # frozen registry token is spelled `Escherichia_coli_Shigella`.
    OrganismEntry("escherichia_coli", "Escherichia",
                  ("Escherichia_coli_Shigella", "Escherichia_coli", "escherichia_coli"),
                  genus_match=("escherichia", "shigella")),
    # GENUS-level in AMRFinder's list, and the local corpus is Campylobacter jejuni + coli.
    OrganismEntry("campylobacter", "Campylobacter", ("Campylobacter",),
                  genus_match=("campylobacter",)),
    # GENUS-level; the local corpus is Salmonella enterica.
    OrganismEntry("salmonella", "Salmonella", ("Salmonella",), genus_match=("salmonella",)),
    # SPECIES-level, and deliberately so: AMRFinder lists K. pneumoniae and K. oxytoca SEPARATELY, and
    # the local corpus also holds K. aerogenes / variicola / michiganensis which must NOT be swept in.
    OrganismEntry("klebsiella_pneumoniae", "Klebsiella_pneumoniae",
                  ("Klebsiella", "Klebsiella_pneumoniae")),
    OrganismEntry("klebsiella_oxytoca", "Klebsiella_oxytoca", ()),
    OrganismEntry("pseudomonas_aeruginosa", "Pseudomonas_aeruginosa", ("Pseudomonas_aeruginosa",)),
    OrganismEntry("staphylococcus_aureus", "Staphylococcus_aureus", ("Staphylococcus_aureus",)),
    OrganismEntry("acinetobacter_baumannii", "Acinetobacter_baumannii", ("Acinetobacter",)),
    OrganismEntry("streptococcus_pneumoniae", "Streptococcus_pneumoniae", ("Streptococcus_pneumoniae",)),
    OrganismEntry("neisseria_gonorrhoeae", "Neisseria_gonorrhoeae", ()),
    OrganismEntry("enterococcus_faecium", "Enterococcus_faecium", ()),
    OrganismEntry("enterobacter_cloacae", "Enterobacter_cloacae", ()),
    # AMRFinder is bacterial: there is no AMRFinder-for-fungi, so C. auris routes to the BLAST
    # ERG11/FKS1 target-site engine instead (dna_decode/data/fungal_amr.py).
    #
    # THREE NAMES, and the extra two are not pedantry -- they are why this organism silently fell OUT
    # of the closed set until 2026-10-04. The genus was RECLASSIFIED, so NCBI's own `ORGANISM` field on
    # the project's own C. auris assemblies reads `Candidozyma auris` (7 of 8 local genomes) and, on one
    # older record, `[Candida] auris` -- NCBI's bracket notation for a provisional genus placement.
    # The project registry still says `Candida_auris` (the fungal_amr cell), so all three must resolve.
    # SOURCED, not remembered: the strings were read off `annotations.gbk` in
    # data/raw/ar_bank_caur/genomes/ (GCA_003014415.1 -> "Candidozyma auris",
    # GCF_002759435.1 -> "Candidozyma auris B8441"). The bracketed form is listed LITERALLY rather than
    # stripped by a clever bracket-normalizer, which would over-match other provisional names.
    OrganismEntry("candida_auris", None,
                  ("Candida_auris", "Candidozyma_auris", "[Candida]_auris")),
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


def _genus_index() -> dict[str, str]:
    """genus token -> canonical, for GENUS-level entries only.

    Refuses a genus claimed by two entries: that would mean one genus routes two ways, which is a
    specification error, not something to resolve by dict order.
    """
    idx: dict[str, str] = {}
    for e in _ENTRIES:
        for g in e.genus_match:
            prev = idx.get(g.lower())
            if prev is not None and prev != e.canonical:
                raise ValueError(f"genus {g!r} claimed by both {prev!r} and {e.canonical!r}")
            idx[g.lower()] = e.canonical
    return idx


_GENERA: dict[str, str] = _genus_index()


def canonical(token: str) -> str:
    """Normalize any known organism spelling to its canonical key.

    Resolution order, and the order matters:
      1. case-insensitive EXACT alias (canonical key, AMRFinder `-O` value, or registry label)
      2. GENUS fallback, but ONLY for entries that declare `genus_match`

    Exact must win, or `Klebsiella_oxytoca` could be captured by a Klebsiella genus rule. Species in a
    genus that declares no `genus_match` stay UNKNOWN by design -- that is how K. aerogenes remains
    out-of-set rather than being silently answered as K. pneumoniae.

    RAISES `UnknownOrganism` on anything else -- never defaults.
    """
    if not isinstance(token, str) or not token.strip():
        raise UnknownOrganism(f"empty organism token: {token!r}")
    # Whitespace is normalized to "_" BEFORE the alias lookup. Without this, `canonical` answered the
    # same question two ways: "Escherichia coli" resolved (via the genus fallback, which already split
    # on whitespace) while "Candidozyma auris" raised, because its entry declares no genus_match. The
    # caller's spelling should not decide whether an alias is findable.
    t = " ".join(token.split()).strip().lower()
    got = _ALIASES.get(t) or _ALIASES.get(t.replace(" ", "_"))
    if got is not None:
        return got
    genus = t.replace(" ", "_").split("_")[0]
    got = _GENERA.get(genus)
    if got is not None:
        return got
    raise UnknownOrganism(
        f"{token!r} is not a supported organism. Supported canonical keys: "
        f"{sorted(BY_CANONICAL)}")


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
