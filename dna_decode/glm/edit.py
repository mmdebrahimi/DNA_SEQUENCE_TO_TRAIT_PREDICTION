"""`GenomeEdit` — the keystone representation every GLM edit class reduces to.

**ONE primitive.** Any substitution, deletion, insertion or block replacement is the same operation:
*replace the span `[start, start + len(ref))` on a contig with `alt`*. That is why VCF and SPDI are
shaped this way, and adopting it means the four edit classes the GLM needs are not four code paths —

    coding substitution   ref="C"      alt="T"       (len == len, == 1)
    knockout / deletion   ref=<gene>   alt=""        (len -> 0)
    heterologous insert   ref=""       alt=<cassette> (0 -> len)
    expression re-tune    ref=<prom>   alt=<new prom> (len -> len', either way)

**Two fields that must not be conflated, and the repo has paid for conflating this class of thing before:**

* `kind`   — the SHAPE. Fully DERIVABLE from (ref, alt) and therefore checkable. `derive_kind` is pure.
* `intent` — the SEMANTIC class the proposer is claiming (is this a coding substitution? an expression
  re-tune?). That is a claim about ANNOTATION, not about the sequence, so it can be WRONG and nothing in
  this module can adjudicate it. It is carried, never trusted, and `verify_intent` is where annotation
  gets to disagree. Storing one field for both would let an unverifiable claim inherit a derived field's
  credibility.

**The load-bearing safety property: `ref` must match the genome.** An edit whose `ref` does not equal
the genome's actual bases at that span was computed against a different assembly, a different strand, or
a different coordinate convention — and applying it would silently corrupt the sequence while every
downstream oracle happily returned numbers. So `apply` VERIFIES before it writes and raises
`RefMismatchError` otherwise. This is the same guard the shipped `dna-forward` CDS path already relies
on, where the REF check "is the safety net that proves the conversion landed", and the same failure mode
`constraints.universal.ref_base_matches_cds` exists to catch one level down.

**Coordinates, stated once because this repo has three recorded coordinate traps** (Ambler vs linear
numbering, minus-strand flips, a published +81 axis convention):

* `start` is **0-based, inclusive**. The span is `[start, start + len(ref))`, half-open.
* Bases are in **GENOME orientation**, never gene orientation — so an edit to a minus-strand gene carries
  the complemented bases, exactly as the shipped `dna-forward --genomic-pos` surface already requires.
  This module performs NO strand reasoning; a caller translating gene->genome coordinates must do the
  complement itself, and the `ref` check is what catches it if they forget.
* An insertion has `ref == ""` and inserts BEFORE `start`.

**Multi-edit ordering is a silent-corruption trap, so `EditSet` owns it.** Applying edits left-to-right
shifts every later coordinate by the preceding length deltas, so an edit set applied in ascending order
lands in the wrong place without raising. `EditSet.apply` therefore applies in DESCENDING start order
(later edits first, so earlier coordinates stay valid) and REFUSES overlapping edits outright rather than
silently resolving them — two edits to one span have no well-defined composition.

What this module deliberately does NOT do: no annotation lookup, no strand inference, no scoring, no
network, no GPU. It is pure so the whole keystone is testable offline.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace as _dc_replace

# ---------------------------------------------------------------------------------------------------
# Edit SHAPE — derived, checkable
# ---------------------------------------------------------------------------------------------------
SUB = "SUB"          # len(ref) == len(alt) > 0  (a 1-bp SUB is the point-mutation case)
DEL = "DEL"          # len(alt) == 0 < len(ref)
INS = "INS"          # len(ref) == 0 < len(alt)
REPLACE = "REPLACE"  # both non-empty, different lengths

EDIT_KINDS = frozenset({SUB, DEL, INS, REPLACE})

# ---------------------------------------------------------------------------------------------------
# Edit INTENT — a proposer's claim, NOT derivable from sequence
# ---------------------------------------------------------------------------------------------------
CODING_SUBSTITUTION = "coding_substitution"
KNOCKOUT = "knockout"
HETEROLOGOUS_INSERT = "heterologous_insert"
EXPRESSION_TUNE = "expression_tune"
UNKNOWN_INTENT = "unknown"

EDIT_INTENTS = frozenset({
    CODING_SUBSTITUTION, KNOCKOUT, HETEROLOGOUS_INSERT, EXPRESSION_TUNE, UNKNOWN_INTENT,
})

#: Which SHAPES are coherent with which INTENT. An intent outside its shape is a proposer bug, and
#: `verify_intent` reports it rather than this module silently relabelling the edit.
_INTENT_SHAPES: dict[str, frozenset[str]] = {
    CODING_SUBSTITUTION: frozenset({SUB}),
    KNOCKOUT: frozenset({DEL, SUB, REPLACE}),      # a stop-codon SUB is a legitimate knockout
    HETEROLOGOUS_INSERT: frozenset({INS, REPLACE}),  # REPLACE = insert-with-excision
    EXPRESSION_TUNE: frozenset({SUB, REPLACE, INS, DEL}),  # any promoter/RBS surgery
    UNKNOWN_INTENT: EDIT_KINDS,
}

_DNA = frozenset("ACGT")


class EditApplicationError(ValueError):
    """An edit could not be applied to the supplied genome."""


class RefMismatchError(EditApplicationError):
    """`ref` does not match the genome at `[start, start+len(ref))`.

    Raised rather than tolerated because the three things that cause it — wrong assembly, wrong strand,
    wrong coordinate base — all produce a plausible-looking edited sequence and a corrupted downstream
    analysis that no oracle can detect.
    """


class EditOverlapError(ValueError):
    """Two edits in one set touch the same span; their composition is undefined."""


def derive_kind(ref: str, alt: str) -> str:
    """PURE: the edit SHAPE implied by (ref, alt). Total over all string pairs except ref==alt=="".

    `ref == alt == ""` is not an edit at all and raises — a no-op edit in a proposal set is a generator
    bug worth surfacing, not a harmless identity.
    """
    if not ref and not alt:
        raise ValueError("ref and alt are both empty: that is not an edit")
    if not alt:
        return DEL
    if not ref:
        return INS
    return SUB if len(ref) == len(alt) else REPLACE


@dataclass(frozen=True)
class GenomeEdit:
    """One DNA edit: replace `[start, start+len(ref))` on `contig` with `alt`.

    Frozen on purpose — an edit is a proposal that gets passed to several oracles, and a mutable one
    could be rewritten between two scores, making a committee's disagreement unreproducible.
    """

    contig: str
    start: int
    ref: str
    alt: str
    intent: str = UNKNOWN_INTENT
    #: Free-form provenance: which proposer produced this, under what settings. Carried, never scored.
    origin: dict = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"start must be >= 0, got {self.start}")
        if not self.contig:
            raise ValueError("contig must be a non-empty identifier")
        if self.intent not in EDIT_INTENTS:
            raise ValueError(f"unknown intent {self.intent!r}; expected one of {sorted(EDIT_INTENTS)}")
        for name, seq in (("ref", self.ref), ("alt", self.alt)):
            bad = set(seq.upper()) - _DNA
            if bad:
                raise ValueError(
                    f"{name} contains non-ACGT characters {sorted(bad)}: "
                    "ambiguity codes and gaps are refused rather than coerced, because a generator "
                    "emitting N is making a claim it cannot support"
                )
        object.__setattr__(self, "ref", self.ref.upper())
        object.__setattr__(self, "alt", self.alt.upper())
        derive_kind(self.ref, self.alt)  # raises on the empty/empty no-op

    # -- derived -------------------------------------------------------------------------------------
    @property
    def kind(self) -> str:
        """The edit SHAPE. Derived, never stored, so it cannot drift from (ref, alt)."""
        return derive_kind(self.ref, self.alt)

    @property
    def end(self) -> int:
        """Exclusive end of the replaced span. `start == end` for an insertion."""
        return self.start + len(self.ref)

    @property
    def length_delta(self) -> int:
        """Net change in contig length. Negative for a net deletion."""
        return len(self.alt) - len(self.ref)

    @property
    def is_length_preserving(self) -> bool:
        """True iff the edit cannot shift any downstream coordinate — the safe-to-compose case."""
        return self.length_delta == 0

    # -- verification --------------------------------------------------------------------------------
    def verify_ref(self, genome: dict[str, str]) -> None:
        """Raise unless `ref` matches `genome[contig]` exactly at this span. The core safety check."""
        if self.contig not in genome:
            raise EditApplicationError(
                f"contig {self.contig!r} absent from genome (have: {sorted(genome)[:5]}...)"
            )
        seq = genome[self.contig]
        if self.end > len(seq):
            raise EditApplicationError(
                f"span [{self.start},{self.end}) runs past the end of {self.contig!r} "
                f"(length {len(seq)})"
            )
        observed = seq[self.start:self.end]
        if observed.upper() != self.ref:
            raise RefMismatchError(
                f"ref mismatch on {self.contig} at {self.start}: edit claims {self.ref or '(empty)'!r}, "
                f"genome has {observed or '(empty)'!r}. The usual causes are a different assembly, "
                f"gene-orientation bases where genome orientation is required, or a 1-based coordinate."
            )

    def verify_intent(self) -> tuple[bool, str]:
        """Is the claimed `intent` coherent with the derived `kind`?

        Returns `(ok, reason)`. This is the ONLY adjudication available without annotation: it catches a
        proposer claiming `coding_substitution` for a 300-bp deletion, and cannot catch a proposer
        claiming `expression_tune` for an edit that is really inside a CDS. That stronger check needs an
        annotation and belongs to the caller that has one.
        """
        allowed = _INTENT_SHAPES[self.intent]
        if self.kind in allowed:
            return True, ""
        return False, (
            f"intent {self.intent!r} is incoherent with shape {self.kind}: "
            f"that intent admits {sorted(allowed)}"
        )

    # -- application ---------------------------------------------------------------------------------
    def apply(self, genome: dict[str, str]) -> dict[str, str]:
        """Return a NEW genome dict with this edit applied. Verifies `ref` first; never mutates input."""
        self.verify_ref(genome)
        seq = genome[self.contig]
        out = dict(genome)
        out[self.contig] = seq[:self.start] + self.alt + seq[self.end:]
        return out

    def invert(self) -> GenomeEdit:
        """The edit that undoes this one — `ref` and `alt` swapped.

        Round-trip is then true by construction rather than by testing luck:
        `e.invert().apply(e.apply(g)) == g`. The inverse of an insertion is a deletion at the same start,
        which is why the single-primitive representation is what makes reversibility free.
        """
        return _dc_replace(self, ref=self.alt, alt=self.ref)

    # -- interop -------------------------------------------------------------------------------------
    def constraint_ctx(self, **extra) -> dict:
        """Build the `ctx` dict `constraints.universal.evaluate_all` consumes.

        Only the fields this type can know are filled; anything requiring annotation (`cds`, `protein`,
        `codon_table`) must be supplied by the caller via `**extra`. Deliberately NOT guessed — a wrong
        codon table silently turns 71.6% of real *M. genitalium* CDS into false nonsense, which is the
        measured reason `constraints/` exists.
        """
        ctx = {
            "contig": self.contig,
            "pos": self.start,
            "ref": self.ref,
            "alt": self.alt,
            "edit_kind": self.kind,
            "edit_intent": self.intent,
        }
        ctx.update(extra)
        return ctx

    def as_dict(self) -> dict:
        ok, reason = self.verify_intent()
        return {
            "record": "glm-genome-edit-v1",
            "contig": self.contig,
            "start": self.start,
            "end": self.end,
            "ref": self.ref,
            "alt": self.alt,
            "kind": self.kind,
            "intent": self.intent,
            "intent_coherent": ok,
            "intent_note": reason,
            "length_delta": self.length_delta,
            "origin": dict(self.origin),
        }

    def hgvs_like(self) -> str:
        """A compact human-readable form. NOT standards-compliant HGVS — named `_like` so nobody
        pipes it into a tool expecting the real grammar."""
        if self.kind == INS:
            return f"{self.contig}:{self.start}ins{self.alt}"
        if self.kind == DEL:
            return f"{self.contig}:{self.start}_{self.end}del"
        if self.kind == SUB and len(self.ref) == 1:
            return f"{self.contig}:{self.start}{self.ref}>{self.alt}"
        return f"{self.contig}:{self.start}_{self.end}delins{self.alt}"


@dataclass(frozen=True)
class EditSet:
    """An ordered set of edits applied as one atomic change, with the ordering trap handled.

    Why this exists rather than `for e in edits: g = e.apply(g)`: each applied edit shifts every
    subsequent coordinate by its `length_delta`, so the naive loop silently misplaces every edit after
    the first length-changing one. `apply` solves that by assembling the result in ONE pass instead of
    splicing repeatedly (see its docstring — descending-order splicing was tried first and is NOT
    sufficient). Overlapping edits are refused outright, because two edits to one span have no defined
    composition.
    """

    edits: tuple[GenomeEdit, ...]
    label: str = ""

    def __post_init__(self) -> None:
        if not self.edits:
            raise ValueError("an EditSet needs at least one edit")
        self._assert_disjoint()

    def _assert_disjoint(self) -> None:
        by_contig: dict[str, list[GenomeEdit]] = {}
        for e in self.edits:
            by_contig.setdefault(e.contig, []).append(e)
        for contig, group in by_contig.items():
            ordered = sorted(group, key=lambda e: (e.start, e.end))
            for a, b in zip(ordered, ordered[1:]):
                # Half-open spans touch iff the next start is strictly inside the previous span. Two
                # insertions at the SAME point also conflict: their relative order would decide the
                # result, so an arbitrary choice would be silent.
                if b.start < a.end or (a.start == b.start and a.kind == INS and b.kind == INS):
                    raise EditOverlapError(
                        f"overlapping edits on {contig}: {a.hgvs_like()} and {b.hgvs_like()}. "
                        "Composition of two edits to one span is undefined; merge them into a single "
                        "REPLACE if that is what was meant."
                    )

    @property
    def contigs(self) -> tuple[str, ...]:
        return tuple(sorted({e.contig for e in self.edits}))

    @property
    def n_edits(self) -> int:
        return len(self.edits)

    def apply(self, genome: dict[str, str]) -> dict[str, str]:
        """Apply every edit atomically. ALL refs are verified BEFORE any is written.

        Verifying up-front matters: a half-applied set would leave a genome that is neither the original
        nor the intended design, and the `ref` checks of later edits are only meaningful against the
        unedited sequence.

        **SINGLE PASS, no sequential splicing** — and the reason is a bug this very test suite caught.
        The obvious fix for coordinate drift is to splice in DESCENDING start order so a later edit cannot
        move an earlier one's coordinates. That is not sufficient: inverting a set can legitimately put an
        INSERTION and a DELETION at the SAME start (undoing a deletion inserts at p, while undoing an
        insertion deletes at p, and they coincide precisely because the deletion removed what separated
        them). Two splices at one coordinate have no stable order, and the second then lands on bases the
        first just wrote — measured: it deleted `CCCC` instead of `ACGT`, silently.

        Walking the contig once and emitting segments removes the ambiguity entirely rather than ordering
        around it: each edit contributes `seq[pos:start] + alt` and advances `pos` to `end`, so no edit
        ever reads a region another edit has written. A zero-width insertion sorts before a deletion at
        the same coordinate via `(start, end)`, which is exactly "insert before" semantics.
        """
        for e in self.edits:
            e.verify_ref(genome)
        out = dict(genome)
        by_contig: dict[str, list[GenomeEdit]] = {}
        for e in self.edits:
            by_contig.setdefault(e.contig, []).append(e)
        for contig, group in by_contig.items():
            seq = out[contig]
            pieces: list[str] = []
            pos = 0
            for e in sorted(group, key=lambda e: (e.start, e.end)):
                pieces.append(seq[pos:e.start])
                pieces.append(e.alt)
                pos = e.end
            pieces.append(seq[pos:])
            out[contig] = "".join(pieces)
        return out

    def invert(self) -> EditSet:
        """The set that undoes this one.

        The inverted edits' coordinates are those of the EDITED genome, so this is only a valid round-trip
        against `self.apply(g)` — which is exactly the use. Recomputing shifted coordinates is done here
        rather than left to the caller, because getting it wrong is the silent trap this class exists for.
        """
        shifted: list[GenomeEdit] = []
        for e in self.edits:
            prior = sum(
                o.length_delta for o in self.edits
                if o.contig == e.contig and o.start < e.start
            )
            shifted.append(_dc_replace(e, start=e.start + prior, ref=e.alt, alt=e.ref))
        return EditSet(edits=tuple(shifted), label=f"invert({self.label})" if self.label else "")

    @property
    def total_length_delta(self) -> int:
        return sum(e.length_delta for e in self.edits)

    def as_dict(self) -> dict:
        return {
            "record": "glm-edit-set-v1",
            "label": self.label,
            "n_edits": self.n_edits,
            "contigs": list(self.contigs),
            "total_length_delta": self.total_length_delta,
            "edits": [e.as_dict() for e in self.edits],
        }
