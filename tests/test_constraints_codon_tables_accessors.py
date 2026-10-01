"""Gap coverage for `constraints/codon_tables.py` accessors and traps (test-epilogue, 2026-10-01).

The plan's `test_constraints_codon_tables.py` pins the DATA (derivation from the verbatim NCBI strings,
the stated differences, the refusal on an unknown clade). What it does not reach:

* the accessors' input handling -- `None`, whitespace padding, and `table_meta`'s own refusal path
  (only `table_for`'s is covered);
* two of `_decode`'s three guards: the empty-`AAs` arm of `verified`, and the duplicate-codon check;
* `TABLE_IDS` holds six keys that are NOT registry `TAXON_BUCKETS`, so a clade can be translatable but
  not constraint-resolvable. That asymmetry is real and currently undocumented anywhere;
* **`STOPS` is table-1 only.** Its name carries no clade, it is exported from the package `__init__`,
  and `TGA` is in it while table 4 reads `TGA` as tryptophan. A clade-aware caller reaching for `STOPS`
  gets the standard answer silently -- exactly the latent defect the module was created to remove,
  surviving in the one derived constant;
* `STANDARD` is a plain mutable dict aliased by eleven consumers, whereas `table_for` returns a FRESH
  dict per call. The two have opposite mutation safety and nothing said so.
"""

from __future__ import annotations

import dataclasses

import pytest

from dna_decode.constraints import codon_tables as T
from dna_decode.constraints import registry as R

# --------------------------------------------------------------------- accessor input handling

@pytest.mark.parametrize("clade", [None, "", "   ", "nonsense_clade", "Bacteria", "bacteria."])
def test_table_for_refuses_every_unresolvable_clade_including_None(clade):
    """`(clade or "").strip()` means `None` and whitespace both arrive as `""`. Matching is EXACT and
    case-sensitive, so `"Bacteria"` must refuse rather than resolve by a lenient match."""
    with pytest.raises(T.UnknownCladeError):
        T.table_for(clade)


@pytest.mark.parametrize("clade", [None, "", "   ", "nonsense_clade"])
def test_table_meta_refuses_the_same_inputs_as_table_for(clade):
    """Only `table_for`'s refusal was covered. `table_meta` is a separate public entry point and a
    caller reading provenance must hit the same wall, not a `KeyError` from a raw dict lookup."""
    with pytest.raises(T.UnknownCladeError):
        T.table_meta(clade)


@pytest.mark.parametrize("padded", ["  bacteria  ", "\tbacteria.mycoplasma\n", " eukaryote.ciliate"])
def test_surrounding_whitespace_is_stripped_rather_than_refused(padded):
    """A clade read from a TSV cell or a CLI argument commonly arrives padded; the strip is deliberate,
    so pin it in both accessors rather than leaving it as incidental behaviour."""
    bare = padded.strip()
    assert T.table_for(padded) == T.table_for(bare)
    assert T.table_meta(padded) is T.table_meta(bare)


def test_every_declared_clade_decodes_to_a_complete_64_codon_map():
    """Includes the six clades with no registry bucket (mitochondrion.mold/protozoan/coelenterate,
    eukaryote.dasycladacean/hexamita, standard), which no existing test decodes."""
    for clade in T.TABLE_IDS:
        table = T.table_for(clade)
        assert len(table) == 64, clade
        assert set(table) == set(T.STANDARD), f"{clade} decodes a different codon set"
        assert all(len(aa) == 1 for aa in table.values()), clade


def test_every_TABLE_ID_points_at_a_declared_TABLE():
    """A clade mapped to an id with no `TranslationTable` would raise `KeyError`, not the module's own
    `UnknownCladeError`, so the refusal contract would break for a declared clade."""
    missing = sorted({i for i in T.TABLE_IDS.values()} - set(T.TABLES))
    assert not missing, f"TABLE_IDS references undeclared table ids {missing}"
    for table_id, table in T.TABLES.items():
        assert table.table_id == table_id, f"TABLES[{table_id}] self-reports {table.table_id}"


# --------------------------------------------------------------------- translatable != resolvable

def test_a_clade_can_be_TRANSLATABLE_without_being_CONSTRAINT_RESOLVABLE():
    """Two vocabularies, deliberately different sizes, and the gap is undocumented until now:
    `TABLE_IDS` answers "which genetic code", `TAXON_BUCKETS` answers "which constraints apply".
    `table_for('mitochondrion.mold')` works while `constraints_for(taxon='mitochondrion.mold')`
    REFUSES. That is coherent (a code is known, a constraint scope is not declared) but a caller
    passing one vocabulary's key to the other needs to know which direction fails."""
    translatable_only = sorted(set(T.TABLE_IDS) - set(R.TAXON_BUCKETS))
    assert translatable_only == [
        "eukaryote.dasycladacean", "eukaryote.hexamita", "mitochondrion.coelenterate",
        "mitochondrion.mold", "mitochondrion.protozoan", "standard",
    ], translatable_only

    assert T.table_for("mitochondrion.mold")["TGA"] == "W"
    with pytest.raises(ValueError, match="unknown taxon bucket"):
        R.constraints_for(taxon="mitochondrion.mold")


def test_every_registry_bucket_with_a_shipped_code_specialization_is_translatable():
    """The reverse direction must hold for anything the specializations actually ship, or building the
    specialization would raise at import time."""
    from dna_decode.constraints.specializations import _CODE_CLADES
    for clade in _CODE_CLADES:
        assert clade in T.TABLE_IDS, clade
        assert clade in R.TAXON_BUCKETS, clade


# --------------------------------------------------------------------- _decode's remaining guards

def test_an_EMPTY_AAs_line_is_refused_by_the_verified_gate():
    """`verified` requires BOTH `source_url` and `aas`. Only the source_url arm was covered; an empty
    AAs line with a real URL would otherwise decode 64 codons to empty strings."""
    good = T.TABLES[1]
    blank = dataclasses.replace(good, table_id=97, aas="")
    assert blank.verified is False
    with pytest.raises(T.UnverifiedTableError):
        T._decode(blank)


def test_duplicate_base_strings_are_refused_rather_than_silently_collapsing():
    """The three Base lines are what make the 64 codons distinct. A degenerate Base3 collapses them to
    16 keys, each silently taking the LAST amino acid written -- a fully self-consistent wrong table.
    Non-vacuous: the real table decodes 64."""
    good = T.TABLES[1]
    assert len(T._decode(good)) == 64
    degenerate = dataclasses.replace(good, table_id=96, base3="T" * 64)
    with pytest.raises(ValueError, match="decoded 16 distinct codons"):
        T._decode(degenerate)


def test_the_verbatim_quote_reproduces_all_four_ncbi_lines():
    """It is the provenance a `sourced` constraint carries, so it must contain the actual strings and
    not a summary of them."""
    for table in T.TABLES.values():
        q = table.verbatim_quote
        for label, line in (("AAs", table.aas), ("Base1", table.base1),
                            ("Base2", table.base2), ("Base3", table.base3)):
            assert f"{label}  = {line}" in q, (table.table_id, label)


# --------------------------------------------------------------------- STOPS is not clade-aware

def test_STANDARD_STOPS_names_its_own_scope_and_the_clade_blind_name_is_GONE():
    """FIXED 2026-10-01. The constant was exported as a clade-free `STOPS` while holding the STANDARD
    code's stops -- the trap this module exists to remove, surviving in its one derived constant. It is
    now `STANDARD_STOPS`, whose name carries the scope, and `stops_for(clade)` is the clade-aware way
    to ask. The per-clade facts that made the old name dangerous still hold and are pinned here:

      * `TGA` is a standard stop but is TRYPTOPHAN under table 4 (Mycoplasma/Spiroplasma);
      * `TAA`/`TAG` are standard stops but are GLUTAMINE under table 6 (ciliates);
      * `AGA`/`AGG` are NOT standard stops but ARE stops under table 2 (vertebrate mitochondria).
    """
    assert T.STANDARD_STOPS == frozenset({"TAA", "TAG", "TGA"})
    assert not hasattr(T, "STOPS"), "the clade-blind name is back"
    import dna_decode.constraints as pkg
    assert "STOPS" not in pkg.__all__ and "STANDARD_STOPS" in pkg.__all__
    assert "stops_for" in pkg.__all__

    assert "TGA" in T.STANDARD_STOPS and T.table_for("bacteria.mycoplasma")["TGA"] == "W"
    ciliate = T.table_for("eukaryote.ciliate")
    assert ciliate["TAA"] == "Q" and ciliate["TAG"] == "Q"
    mito = T.table_for("mitochondrion.vertebrate")
    assert {"AGA", "AGG"}.isdisjoint(T.STANDARD_STOPS) and mito["AGA"] == mito["AGG"] == "*"


def test_stops_for_gives_a_DIFFERENT_answer_per_clade():
    """The replacement idiom is now SHIPPED as `stops_for`, not a local helper in a test."""
    assert T.stops_for("standard") == T.STANDARD_STOPS
    assert T.stops_for("bacteria") == T.STANDARD_STOPS                 # table 11 matches internally
    assert T.stops_for("bacteria.mycoplasma") == frozenset({"TAA", "TAG"})
    assert T.stops_for("eukaryote.ciliate") == frozenset({"TGA"})
    assert T.stops_for("mitochondrion.vertebrate") == frozenset({"TAA", "TAG", "AGA", "AGG"})


def test_stops_for_refuses_an_unknown_clade_like_table_for():
    import pytest as _pytest
    with _pytest.raises(T.UnknownCladeError):
        T.stops_for("not_a_clade")


# --------------------------------------------------------------------- mutation safety

def test_table_for_returns_a_FRESH_dict_so_a_caller_cannot_corrupt_the_source():
    """Each call re-derives from the verbatim strings, so mutating a returned table is local."""
    a, b = T.table_for("bacteria"), T.table_for("bacteria")
    assert a == b and a is not b
    a["TGA"] = "CORRUPTED"
    assert T.table_for("bacteria")["TGA"] == "*"
    assert T.STANDARD["TGA"] == "*"


def test_STANDARD_is_a_SHARED_MUTABLE_dict_and_that_is_the_refactors_one_new_risk():
    """Opposite safety from `table_for`. Before consolidation each consumer owned a private copy, so a
    stray write was local; now eleven aliases point at ONE dict and a write would reach all of them AND
    diverge from `table_for('standard')`, which re-derives. Documented here; the static guard that no
    consumer actually writes to its alias lives in `test_codon_consumers_behavioural.py`."""
    import dna_decode.forward.genome_edit as GE
    import dna_decode.forward.inverse as INV
    import dna_decode.typing.codon_map as CM

    assert isinstance(T.STANDARD, dict), "not a MappingProxy, so writes are not blocked"
    assert GE._CODON is INV.CODON_TABLE is CM.CODON is T.STANDARD

    snapshot = dict(T.STANDARD)
    T.STANDARD["TGA"] = "W"                       # simulate one consumer's stray write
    try:
        assert GE._CODON["TGA"] == "W", "the write did not propagate; shared-alias premise is wrong"
        assert T.table_for("standard")["TGA"] == "*", "table_for must stay derived, not aliased"
    finally:
        T.STANDARD.clear()
        T.STANDARD.update(snapshot)
    assert T.STANDARD == snapshot and T.STANDARD["TGA"] == "*"


def test_the_module_level_clade_aliases_match_their_clade_lookup():
    """`MYCOPLASMA_SPIROPLASMA` and friends are separate `_decode` calls; if one were wired to the wrong
    table id the probe's clade demo would still look coherent."""
    for alias, clade in ((T.VERTEBRATE_MITOCHONDRIAL, "mitochondrion.vertebrate"),
                         (T.MYCOPLASMA_SPIROPLASMA, "bacteria.mycoplasma"),
                         (T.CILIATE_NUCLEAR, "eukaryote.ciliate"),
                         (T.BACTERIAL_PLASTID, "bacteria"),
                         (T.STANDARD, "standard")):
        assert alias == T.table_for(clade), clade
    for alias in (T.VERTEBRATE_MITOCHONDRIAL, T.MYCOPLASMA_SPIROPLASMA, T.CILIATE_NUCLEAR):
        assert alias is not T.STANDARD and alias != T.STANDARD
