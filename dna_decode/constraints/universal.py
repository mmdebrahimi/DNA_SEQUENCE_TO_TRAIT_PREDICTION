"""Universal (axiomatic) constraints — Type A laws that hold for life nobody has characterised.

Each is a LAW, not a curated association, so each is `kind=axiomatic` and carries no citation: mass
balance and the reading frame are not someone's finding. They transfer to an uncharacterised organism BY
CONSTRUCTION, which is exactly what the deployed per-organism catalogs cannot do.

Each constraint can only ever REFUTE or abstain. That asymmetry is the point: a law can delete an
impossible prediction but cannot pick the right one, so the worst case for a constraint layer is that it
removes nothing. Which is why every one reports its kill count (`report.py`) — a filter that removes
nothing is not a control.

CONTEXT CONTRACT. Each `check` takes a dict and returns `inapplicable` when its preconditions are absent,
never `satisfied`. Recognised keys:

* `cds`        — nucleotide coding sequence (str)
* `clade`      — a key in `codon_tables.TABLE_IDS`; selects the genetic code. Absent -> standard.
* `codon_pos`  — 1-based codon index for a substitution
* `ref_base` / `alt_base` — single nucleotides, genome orientation already resolved
* `nt_pos`     — 1-based nucleotide index into `cds`
* `wt_aa` / `mut_aa` — single-letter amino acids
"""
from __future__ import annotations

from typing import Any

from dna_decode.constraints.codon_tables import STANDARD, table_for
from dna_decode.constraints.registry import (
    AXIOMATIC,
    INAPPLICABLE,
    REFUTED,
    SATISFIED,
    UNIVERSAL,
    Constraint,
    ConstraintConflictError,
    ConstraintVerdict,
    register,
)

_BASES = ("A", "C", "G", "T")


def _table(ctx: dict[str, Any]) -> dict[str, str]:
    """The clade's code, or the standard code when no clade is supplied."""
    clade = ctx.get("clade")
    return table_for(clade) if clade else STANDARD


def _ok(d: str = "") -> ConstraintVerdict:
    return ConstraintVerdict(SATISFIED, d)


def _no(d: str) -> ConstraintVerdict:
    return ConstraintVerdict(REFUTED, d)


def _na(d: str) -> ConstraintVerdict:
    return ConstraintVerdict(INAPPLICABLE, d)


# --------------------------------------------------------------------------- laws

def cds_length_multiple_of_three(ctx: dict[str, Any]) -> ConstraintVerdict:
    cds = ctx.get("cds")
    if not cds:
        return _na("no cds supplied")
    if len(cds) % 3 != 0:
        return _no(f"cds length {len(cds)} is not a multiple of 3")
    return _ok(f"cds length {len(cds)} is in frame")


def no_internal_stop_codon(ctx: dict[str, Any]) -> ConstraintVerdict:
    """Universal in FORM, clade-specialized in DATA — which is the whole argument for the clade tables.
    The same mid-frame TGA is a stop under table 1 and tryptophan under table 4."""
    cds = ctx.get("cds")
    if not cds:
        return _na("no cds supplied")
    if len(cds) % 3 != 0:
        return _na("cds not in frame; frame constraint must pass first")
    tab = _table(ctx)
    codons = [cds[i:i + 3].upper() for i in range(0, len(cds), 3)]
    for i, c in enumerate(codons[:-1]):          # the terminal codon is allowed to be a stop
        if tab.get(c) == "*":
            return _no(f"internal stop {c} at codon {i + 1} of {len(codons)} "
                       f"(clade={ctx.get('clade') or 'standard'})")
    return _ok(f"no internal stop across {len(codons)} codons "
               f"(clade={ctx.get('clade') or 'standard'})")


def ref_base_matches_cds(ctx: dict[str, Any]) -> ConstraintVerdict:
    """The coordinate check that catches a frame/offset error before any effect is predicted."""
    cds, pos, ref = ctx.get("cds"), ctx.get("nt_pos"), ctx.get("ref_base")
    if not cds or pos is None or not ref:
        return _na("needs cds, nt_pos and ref_base")
    if not (1 <= pos <= len(cds)):
        return _no(f"nt_pos {pos} is outside the cds (length {len(cds)})")
    actual = cds[pos - 1].upper()
    if actual != ref.upper():
        return _no(f"cds[{pos}] is {actual!r} but ref_base is {ref.upper()!r}")
    return _ok(f"cds[{pos}] == {ref.upper()!r}")


def substitution_reachable_by_single_nt(ctx: dict[str, Any]) -> ConstraintVerdict:
    """A substitution must be reachable by ONE nucleotide change.

    This is the constraint whose absence previously let a conservativeness result be over-read: DRMs
    looked chemically conservative largely because single-nucleotide neighbours are chemically similar
    by construction under the error-minimising genetic code.
    """
    tab = _table(ctx)
    wt_aa, mut_aa = ctx.get("wt_aa"), ctx.get("mut_aa")
    codon = (ctx.get("codon") or "").upper()
    if codon:
        if len(codon) != 3 or codon not in tab:
            return _na(f"codon {codon!r} is not a resolvable triplet")
        if not mut_aa:
            return _na("needs mut_aa alongside codon")
        for i in range(3):
            for b in _BASES:
                if b == codon[i]:
                    continue
                if tab.get(codon[:i] + b + codon[i + 1:]) == mut_aa:
                    return _ok(f"{codon}->{mut_aa} reachable by one nt at position {i + 1}")
        return _no(f"{codon}->{mut_aa} needs more than one nucleotide change")
    if not wt_aa or not mut_aa:
        return _na("needs (codon, mut_aa) or (wt_aa, mut_aa)")
    wt_codons = [c for c, aa in tab.items() if aa == wt_aa]
    if not wt_codons:
        return _na(f"no codon encodes wt_aa {wt_aa!r} under this code")
    for c in wt_codons:
        for i in range(3):
            for b in _BASES:
                if b == c[i]:
                    continue
                if tab.get(c[:i] + b + c[i + 1:]) == mut_aa:
                    return _ok(f"{wt_aa}->{mut_aa} reachable by one nt (via {c})")
    return _no(f"{wt_aa}->{mut_aa} needs more than one nucleotide change under any codon")


def start_codon_present(ctx: dict[str, Any]) -> ConstraintVerdict:
    """The first codon must encode methionine — but ONLY when the caller asserts a complete gene.

    REQUIRES `cds_is_complete_gene=True`. Absent or False -> `inapplicable`, never refuted.

    WHY, measured 2026-10-01: run over eight committed reference CDS this law refuted five of them, and
    all five were correct sequences — HIV RT/PR/IN/CA are in-frame extracts from the HXB2 gag/pol
    POLYPROTEIN (RT begins CCC, PR begins CCT) and SARS-CoV-2 Mpro is nsp5 cleaved from ORF1ab (begins
    AGT; Ser really is its mature N-terminus). Only the three complete C. auris ERG11 genes began ATG.
    A law cannot distinguish a legitimately truncated extract from a broken start codon by looking at
    the sequence, and a FALSE refutation is the one failure mode that makes a constraint filter
    dangerous: it would delete valid predictions. So the completeness claim belongs to the caller.

    Also scoped narrow for a second reason: NCBI's tables differ in PERMITTED ALTERNATIVE starts (that
    is the only way table 11 differs from table 1), which this layer does not model.
    """
    cds = ctx.get("cds")
    if not cds or len(cds) < 3:
        return _na("no cds supplied, or shorter than one codon")
    if not ctx.get("cds_is_complete_gene"):
        return _na("caller did not assert cds_is_complete_gene; a truncated or polyprotein-derived "
                   "extract legitimately lacks a start codon, so this abstains rather than refuting")
    tab = _table(ctx)
    first = cds[:3].upper()
    if tab.get(first) == "M":
        return _ok(f"start codon {first} encodes M")
    return _no(f"start codon {first} encodes {tab.get(first)!r}, not M")


#: The laws, in the order a caller should apply them (frame before anything frame-dependent).
UNIVERSAL_CONSTRAINTS: tuple[Constraint, ...] = (
    Constraint("cds_length_multiple_of_three", UNIVERSAL, AXIOMATIC, cds_length_multiple_of_three,
               note="reading frame; a law, not a finding"),
    Constraint("no_internal_stop_codon", UNIVERSAL, AXIOMATIC, no_internal_stop_codon,
               note="universal in form, clade-specialized in data via the genetic-code table"),
    Constraint("ref_base_matches_cds", UNIVERSAL, AXIOMATIC, ref_base_matches_cds,
               note="coordinate agreement; catches a frame/offset error before any effect prediction"),
    Constraint("substitution_reachable_by_single_nt", UNIVERSAL, AXIOMATIC,
               substitution_reachable_by_single_nt,
               note="codon accessibility under the error-minimising genetic code"),
    Constraint("start_codon_present", UNIVERSAL, AXIOMATIC, start_codon_present,
               note="canonical start only; alternative starts are not modelled here"),
)


def register_universal() -> tuple[Constraint, ...]:
    """Register the laws. Idempotent-safe for callers that may run more than once."""
    out = []
    for c in UNIVERSAL_CONSTRAINTS:
        try:
            out.append(register(c))
        except ConstraintConflictError:          # already registered -> keep the existing one
            out.append(c)
        # NOTE: deliberately NOT `except Exception`. A bare catch here also swallowed
        # UnsourcedConstraintError, so the provenance gate's own helper could not distinguish
        # "already registered" from "REFUSED for missing provenance" -- it failed OPEN.
    return tuple(out)


def evaluate_all(ctx: dict[str, Any], report=None) -> dict[str, ConstraintVerdict]:
    """Run every law over one context, recording into `report` when supplied."""
    results: dict[str, ConstraintVerdict] = {}
    for c in UNIVERSAL_CONSTRAINTS:
        v = c.check(ctx)
        results[c.key] = v
        if report is not None:
            report.record(c.key, v)
    return results
