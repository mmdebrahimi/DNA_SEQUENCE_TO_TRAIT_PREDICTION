"""Constraint registry — two orthogonal specialization axes, explicit precedence, conflict refusal.

TWO KINDS OF RULE, and the distinction is what decides transferability:

* `axiomatic` — a LAW (reading frame, mass balance, codon reachability). Universal, needs no citation,
  transfers to an organism nobody has characterised.
* `sourced` — a FACT ABOUT A CLADE (NCBI table 4 applies to Mycoplasma). Requires `source_url` +
  `verbatim_quote`; `register` REFUSES it otherwise. This mirrors `data/ecoff_catalog.py`, which raises
  rather than returning a plausible default, because a wrong constant yields a fully self-consistent
  evaluation that is entirely wrong with nothing downstream able to catch it.

TWO ORTHOGONAL AXES, not one hierarchy:

* **taxon** — molecular machinery: which genetic code, does splicing exist, is HGT available, ploidy.
* **regime** — study design: natural vs constructed population, concentrated vs distributed determinant.
  Keys are IMPORTED from `dna_decode.eval.regime`, never re-declared, so there is one source of truth.

The regime axis exists because this project MEASURED that taxonomy is not the most predictive split: a
yeast segregant cross decoded 12/12 traits while an Arabidopsis natural population failed, and the
discriminating variable was population design rather than organism complexity.

PRECEDENCE. Within the taxon axis, a deeper bucket wins (`bacteria.mycoplasma` over `bacteria`) — the
explicit-ordering convention this codebase already ratified for the genome-map feature tiers. ACROSS the
two axes there is no ordering, because none has been validated: a key claimed by both a taxon scope and a
regime scope RAISES `ConstraintConflictError` rather than silently picking one. Refusing an ambiguous
resolution is the same stance `codon_tables.table_for` takes on an unknown clade.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

REFUTED = "refuted"
SATISFIED = "satisfied"
INAPPLICABLE = "inapplicable"
VERDICTS = frozenset({REFUTED, SATISFIED, INAPPLICABLE})

AXIOMATIC = "axiomatic"
SOURCED = "sourced"

UNIVERSAL = "universal"


class UnsourcedConstraintError(ValueError):
    """Raised when a `sourced` constraint lacks source_url or verbatim_quote."""


class ConstraintConflictError(ValueError):
    """Raised when two constraints claim one key with no validated precedence between them."""


#: The taxonomic axis. Dots denote nesting; a deeper bucket inherits its ancestors' constraints.
#: NOTE: clade nuances SWAP rather than accumulate — bacteria have HGT and operons, eukaryotes have
#: splicing and chromatin. There is deliberately no "higher life form" axis, because depth does not
#: mean "more constrained".
TAXON_BUCKETS: tuple[str, ...] = (
    "virus",
    "bacteria", "bacteria.mycoplasma", "bacteria.spiroplasma",
    "archaea",
    "fungi",
    "plant", "plant.plastid",
    "eukaryote", "eukaryote.nuclear", "eukaryote.ciliate",
    "animal", "animal.mammal", "animal.bird",
    "mitochondrion", "mitochondrion.vertebrate",
)


def regime_keys() -> tuple[str, ...]:
    """The regime axis, read live from `eval.regime`. A second copy here would be the drift this repo
    has hit five times, so this is a read-through rather than a declaration."""
    from dna_decode.eval.regime import REGIMES
    return tuple(r.key for r in REGIMES)


@dataclass(frozen=True)
class ConstraintVerdict:
    """One constraint's outcome. `inapplicable` is distinct from `satisfied` on purpose — "could not
    check" and "checked and clean" must never be conflated, for the same reason the doubt layer keeps
    three states."""

    verdict: str
    detail: str = ""

    def __post_init__(self) -> None:
        if self.verdict not in VERDICTS:
            raise ValueError(f"verdict must be one of {sorted(VERDICTS)}; got {self.verdict!r}")

    @property
    def refuted(self) -> bool:
        return self.verdict == REFUTED


@dataclass(frozen=True)
class Constraint:
    key: str
    scope: str                      # UNIVERSAL | "taxon:<bucket>" | "regime:<key>"
    kind: str                       # AXIOMATIC | SOURCED
    check: Callable[[dict[str, Any]], ConstraintVerdict]
    source_url: str = ""
    verbatim_quote: str = ""
    note: str = ""
    payload: dict[str, Any] = field(default_factory=dict)

    @property
    def axis(self) -> str:
        if self.scope == UNIVERSAL:
            return UNIVERSAL
        return self.scope.split(":", 1)[0]

    @property
    def bucket(self) -> str:
        return "" if self.scope == UNIVERSAL else self.scope.split(":", 1)[1]

    @property
    def depth(self) -> int:
        """Specificity. Universal is 0; a taxon bucket is its dotted depth."""
        if self.scope == UNIVERSAL:
            return 0
        return self.bucket.count(".") + 1


_REGISTRY: dict[tuple[str, str], Constraint] = {}


def register(c: Constraint) -> Constraint:
    """Add a constraint. REFUSES an unsourced `sourced` constraint and a duplicate (key, scope)."""
    if c.kind not in (AXIOMATIC, SOURCED):
        raise ValueError(f"kind must be {AXIOMATIC!r} or {SOURCED!r}; got {c.kind!r}")
    if c.kind == SOURCED and not (c.source_url and c.verbatim_quote):
        raise UnsourcedConstraintError(
            f"constraint {c.key!r} is kind={SOURCED!r} but lacks "
            f"{'source_url' if not c.source_url else 'verbatim_quote'}. A clade fact asserted without "
            "provenance is indistinguishable from one written from memory.")
    if c.scope != UNIVERSAL:
        if c.axis == "taxon" and c.bucket not in TAXON_BUCKETS:
            raise ValueError(f"unknown taxon bucket {c.bucket!r}; declare it in TAXON_BUCKETS")
        if c.axis == "regime" and c.bucket not in regime_keys():
            raise ValueError(f"unknown regime key {c.bucket!r}; it must exist in eval.regime.REGIMES")
        if c.axis not in ("taxon", "regime"):
            raise ValueError(f"scope must be {UNIVERSAL!r}, 'taxon:<bucket>' or 'regime:<key>'")
    slot = (c.key, c.scope)
    if slot in _REGISTRY:
        raise ConstraintConflictError(f"{c.key!r} already registered for scope {c.scope!r}")
    _REGISTRY[slot] = c
    return c


def clear_registry() -> None:
    """Test-only reset."""
    _REGISTRY.clear()


def all_constraints() -> tuple[Constraint, ...]:
    return tuple(_REGISTRY.values())


def _ancestors(bucket: str) -> list[str]:
    """`a.b.c` -> ['a', 'a.b', 'a.b.c'] so a deeper bucket inherits its ancestors."""
    parts = bucket.split(".")
    return [".".join(parts[: i + 1]) for i in range(len(parts))]


def constraints_for(taxon: str | None = None, regime: str | None = None) -> dict[str, Constraint]:
    """universal ∩ taxon-specialization ∩ regime-specialization, deeper taxon winning on an equal key.

    RAISES `ConstraintConflictError` when a key is claimed by both axes for this query — there is no
    validated cross-axis precedence, so the ambiguity is surfaced rather than resolved by guesswork.
    """
    if taxon is not None and taxon not in TAXON_BUCKETS:
        raise ValueError(f"unknown taxon bucket {taxon!r}; declare it in TAXON_BUCKETS")
    if regime is not None and regime not in regime_keys():
        raise ValueError(f"unknown regime key {regime!r}")

    chain = _ancestors(taxon) if taxon else []
    chosen: dict[str, Constraint] = {}
    from_axis: dict[str, str] = {}

    for c in _REGISTRY.values():
        if c.scope == UNIVERSAL:
            applies, axis = True, UNIVERSAL
        elif c.axis == "taxon":
            applies, axis = c.bucket in chain, "taxon"
        else:
            applies, axis = (regime is not None and c.bucket == regime), "regime"
        if not applies:
            continue
        prev = chosen.get(c.key)
        if prev is None:
            chosen[c.key], from_axis[c.key] = c, axis
            continue
        prev_axis = from_axis[c.key]
        # universal is always overridable by a specialization on either axis
        if prev_axis == UNIVERSAL:
            chosen[c.key], from_axis[c.key] = c, axis
        elif axis == UNIVERSAL:
            continue
        elif prev_axis != axis:
            raise ConstraintConflictError(
                f"key {c.key!r} is claimed by both a {prev_axis} scope ({prev.scope!r}) and a {axis} "
                f"scope ({c.scope!r}). The two axes are orthogonal and no precedence between them has "
                "been validated, so this is refused rather than resolved by guesswork.")
        elif c.depth > prev.depth:
            chosen[c.key] = c
        elif c.depth == prev.depth:
            raise ConstraintConflictError(
                f"key {c.key!r} claimed twice at equal specificity: {prev.scope!r} and {c.scope!r}")
    return chosen


def detect_conflicts() -> list[str]:
    """Scan the whole registry for unresolvable claims. Returns a list of human-readable descriptions
    (empty when clean) so a caller can assert cleanliness without catching."""
    problems: list[str] = []
    by_key: dict[str, list[Constraint]] = {}
    for c in _REGISTRY.values():
        by_key.setdefault(c.key, []).append(c)
    for key, cs in by_key.items():
        axes = {c.axis for c in cs if c.scope != UNIVERSAL}
        if {"taxon", "regime"} <= axes:
            problems.append(
                f"{key!r}: claimed on BOTH axes ("
                + ", ".join(sorted(c.scope for c in cs if c.scope != UNIVERSAL))
                + ") with no validated cross-axis precedence")
        taxon_cs = [c for c in cs if c.axis == "taxon"]
        for i, a in enumerate(taxon_cs):
            for b in taxon_cs[i + 1:]:
                if a.depth == b.depth and a.bucket != b.bucket:
                    # sibling branches are only ambiguous if one bucket could be both -- it cannot,
                    # since a query resolves one ancestor chain. Recorded as informational, not a
                    # conflict, so the check does not over-fire.
                    continue
    return problems
