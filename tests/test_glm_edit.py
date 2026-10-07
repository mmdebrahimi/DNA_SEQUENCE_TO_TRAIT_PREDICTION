"""Guards for the GLM keystone edit representation.

The expensive failures this type exists to prevent are SILENT: a mis-coordinated edit, a wrong-strand
edit and a mis-ordered multi-edit all produce a plausible sequence that every downstream oracle happily
scores. So these tests concentrate on (a) the ref check actually refusing, (b) the multi-edit ordering
being demonstrably load-bearing, and (c) round-trip holding for all four edit classes.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.glm.edit import (
    DEL,
    INS,
    REPLACE,
    SUB,
    EditApplicationError,
    EditOverlapError,
    EditSet,
    GenomeEdit,
    RefMismatchError,
    derive_kind,
)

ROOT = Path(__file__).resolve().parents[1]

# A toy contig with distinct 10-base blocks so a misplaced edit is visible by eye.
G = {"c1": "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"}


# ---------------------------------------------------------------------------------------------------
# shape derivation — the DERIVABLE half
# ---------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(("ref", "alt", "kind"), [
    ("C", "T", SUB),
    ("CCC", "TTT", SUB),
    ("CCC", "", DEL),
    ("", "GGG", INS),
    ("C", "TTTT", REPLACE),
    ("CCCC", "T", REPLACE),
])
def test_derive_kind_covers_all_four_shapes(ref, alt, kind):
    assert derive_kind(ref, alt) == kind
    assert GenomeEdit("c1", 0, ref, alt).kind == kind


def test_an_empty_empty_edit_is_refused_not_treated_as_a_harmless_noop():
    """A no-op in a proposal set is a GENERATOR bug; silently accepting it would hide it."""
    with pytest.raises(ValueError, match="not an edit"):
        derive_kind("", "")
    with pytest.raises(ValueError, match="not an edit"):
        GenomeEdit("c1", 0, "", "")


def test_kind_is_derived_so_it_cannot_drift_from_ref_alt():
    """`kind` is a property, never a stored field — a stored shape could contradict the sequences."""
    e = GenomeEdit("c1", 10, "C", "T")
    assert e.kind == SUB
    assert "kind" not in e.__dataclass_fields__, "kind must stay derived, not stored"


# ---------------------------------------------------------------------------------------------------
# the core safety property: ref must match
# ---------------------------------------------------------------------------------------------------
def test_ref_mismatch_raises_rather_than_corrupting_the_sequence():
    bad = GenomeEdit("c1", 10, "A", "T")      # position 10 is 'C', not 'A'
    with pytest.raises(RefMismatchError, match="ref mismatch"):
        bad.apply(G)


def test_the_ref_check_is_NON_VACUOUS_the_correct_edit_at_the_same_site_passes():
    """Pins that the guard discriminates: same contig, same position, right base -> applies."""
    good = GenomeEdit("c1", 10, "C", "T")
    out = good.apply(G)
    assert out["c1"][10] == "T"
    assert len(out["c1"]) == len(G["c1"])


def test_a_one_based_coordinate_is_caught_by_the_ref_check():
    """The documented coordinate trap: a caller passing 1-based coords lands one base early.

    At 0-based 20 the base is 'G'; a 1-based caller meaning 'the 20th base' passes 19, which is 'C'.
    The ref check is the only thing standing between that and a silently wrong analysis.
    """
    one_based = GenomeEdit("c1", 19, "G", "A")   # meant 0-based 20
    with pytest.raises(RefMismatchError):
        one_based.apply(G)
    assert GenomeEdit("c1", 20, "G", "A").apply(G)["c1"][20] == "A"


def test_gene_orientation_bases_on_a_minus_strand_edit_are_caught():
    """Bases are GENOME orientation. A caller supplying the complement gets refused, not accepted.

    This is the shipped `dna-forward --genomic-pos` convention; without the check a minus-strand edit
    would apply the wrong base and still return a confident effect prediction.
    """
    complemented = GenomeEdit("c1", 10, "G", "A")   # genome has C; G is its complement
    with pytest.raises(RefMismatchError):
        complemented.apply(G)


def test_a_span_past_the_end_of_the_contig_is_refused():
    with pytest.raises(EditApplicationError, match="past the end"):
        GenomeEdit("c1", 38, "TTTT", "A").apply(G)


def test_an_absent_contig_is_refused_and_names_what_is_present():
    with pytest.raises(EditApplicationError, match="absent from genome"):
        GenomeEdit("nope", 0, "A", "T").apply(G)


def test_apply_never_mutates_the_input_genome():
    before = dict(G)
    GenomeEdit("c1", 10, "C", "T").apply(G)
    assert G == before


# ---------------------------------------------------------------------------------------------------
# round-trip, all four classes
# ---------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(("start", "ref", "alt"), [
    (10, "C", "T"),                 # SUB  — coding substitution
    (10, "CCCCCCCCCC", ""),         # DEL  — knockout
    (20, "", "ACGTACGT"),           # INS  — heterologous insert
    (10, "CCCCCCCCCC", "ACGT"),     # REPLACE — expression re-tune
])
def test_round_trip_holds_for_every_edit_class(start, ref, alt):
    e = GenomeEdit("c1", start, ref, alt)
    edited = e.apply(G)
    assert edited != G, "the edit must actually change something"
    back = e.invert().apply(edited)
    assert back == G, f"{e.hgvs_like()} did not round-trip"


def test_invert_of_an_insertion_is_a_deletion_at_the_same_start():
    ins = GenomeEdit("c1", 20, "", "ACGT")
    inv = ins.invert()
    assert inv.kind == DEL and inv.start == 20 and inv.ref == "ACGT"


# ---------------------------------------------------------------------------------------------------
# multi-edit ordering — the silent-corruption trap, proven load-bearing
# ---------------------------------------------------------------------------------------------------
def test_edit_set_ordering_is_NON_VACUOUS_naive_ascending_order_gives_a_DIFFERENT_wrong_answer():
    """THE test for this class. A length-changing edit before a later one shifts its coordinates.

    Two edits: delete the 10-base C block at 10, and substitute the G at 20. Applied DESCENDING (correct)
    the G at 20 becomes A and then the C block vanishes. Applied ASCENDING (naive) the deletion happens
    first, so index 20 now holds a T — and the second edit's ref check would fail, or without a ref check
    would silently edit the wrong base. Asserting the two orders DISAGREE is what proves EditSet is doing
    real work rather than decoration.
    """
    d = GenomeEdit("c1", 10, "CCCCCCCCCC", "")
    s = GenomeEdit("c1", 20, "G", "A")
    correct = EditSet((d, s)).apply(G)

    # the naive ascending loop, written out to show exactly what EditSet avoids
    naive = dict(G)
    for e in sorted((d, s), key=lambda e: e.start):
        seq = naive[e.contig]
        naive[e.contig] = seq[:e.start] + e.alt + seq[e.end:]

    assert correct != naive, "if these agreed, the descending-order fix would be untested"
    assert correct["c1"] == "AAAAAAAAAA" + "A" + "GGGGGGGGG" + "TTTTTTTTTT"
    # and the correct result is exactly what applying them one at a time in the safe order gives
    assert correct == d.apply(s.apply(G))


def test_all_refs_are_verified_before_anything_is_written():
    """A half-applied set is neither the original nor the design. One bad ref must abort the whole set."""
    good = GenomeEdit("c1", 10, "C", "T")
    bad = GenomeEdit("c1", 30, "A", "C")       # position 30 is 'T'
    before = dict(G)
    with pytest.raises(RefMismatchError):
        EditSet((good, bad)).apply(G)
    assert G == before


def test_overlapping_edits_are_refused():
    a = GenomeEdit("c1", 10, "CCCCC", "A")
    b = GenomeEdit("c1", 12, "C", "T")          # inside a's span
    with pytest.raises(EditOverlapError, match="overlapping edits"):
        EditSet((a, b))


def test_two_insertions_at_the_same_point_are_refused_because_their_order_would_decide_the_result():
    a = GenomeEdit("c1", 20, "", "AAA")
    b = GenomeEdit("c1", 20, "", "TTT")
    with pytest.raises(EditOverlapError):
        EditSet((a, b))


def test_abutting_but_non_overlapping_edits_are_allowed():
    """Half-open spans: [10,15) and [15,20) touch without overlapping and compose fine."""
    a = GenomeEdit("c1", 10, "CCCCC", "A")
    b = GenomeEdit("c1", 15, "CCCCC", "T")
    out = EditSet((a, b)).apply(G)
    assert out["c1"] == "AAAAAAAAAA" + "AT" + "GGGGGGGGGG" + "TTTTTTTTTT"


def test_an_insertion_and_a_deletion_at_the_SAME_start_compose_correctly():
    """Regression: the bug that killed descending-order splicing, found by the round-trip test.

    Undoing a deletion inserts at p; undoing an insertion deletes at p. After inversion those coincide,
    because the deletion had removed what separated them. Two SPLICES at one coordinate have no stable
    order, and the second lands on bases the first just wrote — measured, it removed `CCCC` instead of
    `ACGT`. The single-pass assembly in `EditSet.apply` has no such ordering, and a zero-width insertion
    sorts before a same-start deletion, which is "insert before" semantics.
    """
    g = {"c1": "AAAAAAAAAA" + "ACGT" + "GGGG"}
    ins = GenomeEdit("c1", 10, "", "CCCC")        # insert BEFORE the ACGT
    dele = GenomeEdit("c1", 10, "ACGT", "")       # and remove the ACGT
    out = EditSet((ins, dele)).apply(g)
    assert out["c1"] == "AAAAAAAAAA" + "CCCC" + "GGGG", (
        "the insertion must land before the deleted span, and the deletion must remove ACGT — "
        "not the bases the insertion just wrote"
    )
    # Non-vacuity, stated precisely. Sequential splicing is ORDER-DEPENDENT here, and that is the defect:
    # (delete, insert) happens to be right, (insert, delete) is wrong, and NO sort key on `start` can tell
    # them apart because both start at 10 — a stable sort yields input order, i.e. the wrong one. So the
    # claim is not "every order fails", it is "the answer depends on an order nothing can determine".
    def _splice(order):
        naive = dict(g)
        for e in order:
            s = naive["c1"]
            naive["c1"] = s[:e.start] + e.alt + s[e.end:]
        return naive["c1"]

    assert _splice((ins, dele)) != _splice((dele, ins)), (
        "if splice order did not matter here, this regression test would be vacuous"
    )
    assert _splice((ins, dele)) != out["c1"], "the input-order splice is the wrong one"
    assert _splice((dele, ins)) == out["c1"], (
        "and the other order is right — which is the point: correctness hinged on an order that "
        "sorting by start cannot recover"
    )


def test_edit_set_round_trips_through_its_own_inverse():
    es = EditSet((
        GenomeEdit("c1", 10, "CCCCCCCCCC", ""),
        GenomeEdit("c1", 20, "", "ACGT"),
        GenomeEdit("c1", 30, "T", "A"),
    ), label="mixed")
    edited = es.apply(G)
    assert es.invert().apply(edited) == G


def test_an_empty_edit_set_is_refused():
    with pytest.raises(ValueError, match="at least one edit"):
        EditSet(())


def test_edits_on_different_contigs_never_interact():
    g = {"c1": "AAAACCCC", "c2": "GGGGTTTT"}
    es = EditSet((GenomeEdit("c1", 0, "AAAA", ""), GenomeEdit("c2", 4, "TTTT", "A")))
    out = es.apply(g)
    assert out == {"c1": "CCCC", "c2": "GGGGA"}


# ---------------------------------------------------------------------------------------------------
# intent — the CLAIM half, carried but never trusted
# ---------------------------------------------------------------------------------------------------
def test_intent_incoherent_with_shape_is_reported_not_silently_relabelled():
    e = GenomeEdit("c1", 10, "CCCCCCCCCC", "", intent="coding_substitution")
    ok, reason = e.verify_intent()
    assert not ok and "incoherent" in reason
    assert e.intent == "coding_substitution", "the claim is carried verbatim, not corrected"
    assert e.as_dict()["intent_coherent"] is False


def test_a_coherent_intent_verifies():
    assert GenomeEdit("c1", 10, "C", "T", intent="coding_substitution").verify_intent()[0]
    assert GenomeEdit("c1", 20, "", "ACGT", intent="heterologous_insert").verify_intent()[0]
    assert GenomeEdit("c1", 10, "CCCCCCCCCC", "", intent="knockout").verify_intent()[0]


def test_a_stop_codon_substitution_counts_as_a_knockout_intent():
    """Shape-level coherence must not force a knockout to be a deletion — a nonsense SUB is a knockout."""
    assert GenomeEdit("c1", 10, "C", "T", intent="knockout").verify_intent()[0]


def test_an_unknown_intent_string_is_refused_at_construction():
    with pytest.raises(ValueError, match="unknown intent"):
        GenomeEdit("c1", 10, "C", "T", intent="make_it_better")


# ---------------------------------------------------------------------------------------------------
# input hygiene
# ---------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("seq", ["N", "ACGTN", "ACGT-", "acgtx"])
def test_non_acgt_is_refused_rather_than_coerced(seq):
    with pytest.raises(ValueError, match="non-ACGT"):
        GenomeEdit("c1", 0, "A", seq)


def test_lowercase_input_is_normalised_not_rejected():
    e = GenomeEdit("c1", 10, "c", "t")
    assert (e.ref, e.alt) == ("C", "T")
    assert e.apply(G)["c1"][10] == "T"


def test_a_negative_start_is_refused():
    with pytest.raises(ValueError, match="start must be"):
        GenomeEdit("c1", -1, "A", "T")


# ---------------------------------------------------------------------------------------------------
# interop
# ---------------------------------------------------------------------------------------------------
def test_constraint_ctx_fills_only_what_this_type_can_know():
    """A guessed codon table silently turns 71.6% of real M. genitalium CDS into false nonsense, so
    annotation-dependent fields must be absent unless the caller supplies them."""
    ctx = GenomeEdit("c1", 10, "C", "T").constraint_ctx()
    assert ctx["ref"] == "C" and ctx["pos"] == 10 and ctx["edit_kind"] == SUB
    for guessed in ("cds", "protein", "codon_table"):
        assert guessed not in ctx, f"{guessed} must not be invented"
    assert GenomeEdit("c1", 10, "C", "T").constraint_ctx(cds="ATG")["cds"] == "ATG"


def test_length_delta_and_is_length_preserving():
    assert GenomeEdit("c1", 10, "C", "T").is_length_preserving
    assert GenomeEdit("c1", 10, "CCC", "").length_delta == -3
    assert GenomeEdit("c1", 20, "", "ACGT").length_delta == 4
    assert EditSet((
        GenomeEdit("c1", 10, "CCC", ""), GenomeEdit("c1", 20, "", "ACGTAC"),
    )).total_length_delta == 3


def test_hgvs_like_is_named_like_because_it_is_not_real_hgvs():
    assert GenomeEdit("c1", 10, "C", "T").hgvs_like() == "c1:10C>T"
    assert GenomeEdit("c1", 20, "", "ACGT").hgvs_like() == "c1:20insACGT"
    assert GenomeEdit("c1", 10, "CCC", "").hgvs_like() == "c1:10_13del"
    assert GenomeEdit("c1", 10, "CCC", "A").hgvs_like() == "c1:10_13delinsA"


def test_the_edit_is_frozen_so_a_committee_cannot_score_a_moving_target():
    e = GenomeEdit("c1", 10, "C", "T")
    with pytest.raises(Exception):
        e.start = 11  # type: ignore[misc]


def test_as_dict_is_a_stable_versioned_record():
    d = GenomeEdit("c1", 10, "C", "T", intent="coding_substitution", origin={"by": "test"}).as_dict()
    assert d["record"] == "glm-genome-edit-v1"
    assert d["start"] == 10 and d["end"] == 11 and d["kind"] == SUB
    assert d["origin"] == {"by": "test"}


def test_origin_does_not_affect_equality_so_the_same_edit_from_two_proposers_is_one_edit():
    a = GenomeEdit("c1", 10, "C", "T", origin={"by": "a"})
    b = GenomeEdit("c1", 10, "C", "T", origin={"by": "b"})
    assert a == b and len({a, b}) == 1


# ---------------------------------------------------------------------------------------------------
# the real reference — a genuine minus-strand coding edit
# ---------------------------------------------------------------------------------------------------
def _mg1655():
    """The committed MG1655 reference, if present on this host (it lives under the D: junction)."""
    import glob
    for pat in ("data/genomes/*MG1655*.fna", "D:/dna_decode_cache/refseq/GCF_000005845.2/genome.fna",
                "data/**/GCF_000005845*/genome.fna"):
        for hit in glob.glob(str(ROOT / pat) if not pat.startswith("D:") else pat, recursive=True):
            seqs, name, buf = {}, None, []
            for line in Path(hit).read_text().splitlines():
                if line.startswith(">"):
                    if name:
                        seqs[name] = "".join(buf)
                    name, buf = line[1:].split()[0], []
                else:
                    buf.append(line.strip())
            if name:
                seqs[name] = "".join(buf)
            return seqs
    return None


def test_a_real_minus_strand_coding_edit_on_the_real_reference():
    """gyrA S83L on MG1655 is genomic NC_000913.3:2339173 G>A (0-based 2339172).

    gyrA is MINUS strand, so the genome-orientation base is G where the CDS reads C — exactly the case
    where a caller who forgot to complement gets caught. Skips when the reference is not on this host.
    """
    g = _mg1655()
    if not g:
        pytest.skip("MG1655 reference not available on this host")
    contig = next(c for c in g if c.startswith("NC_000913"))
    e = GenomeEdit(contig, 2339172, "G", "A", intent="coding_substitution")
    e.verify_ref(g)                                  # raises if the coordinate convention is wrong
    out = e.apply(g)
    assert out[contig][2339172] == "A"
    assert len(out[contig]) == len(g[contig])
    assert e.invert().apply(out) == g
    # the complemented base must be refused at the same site
    with pytest.raises(RefMismatchError):
        GenomeEdit(contig, 2339172, "C", "T").verify_ref(g)
