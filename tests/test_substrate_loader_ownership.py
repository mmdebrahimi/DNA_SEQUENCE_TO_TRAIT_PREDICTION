"""ONE loader per substrate — the guard that catches a second loader being written.

WHY THIS EXISTS, measured 2026-10-13. `dna_decode/glm/expression.py` had, for weeks:

  * the hyphen-vs-underscore join fix (`_norm_key`, with the mutation arrow protected),
  * a `MIN_JOIN_RATE = 0.95` floor that RAISES on a degraded join,
  * the trap written out in its first 25 lines,
  * and `tests/test_glm_expression.py` pinning `rep.naive_join == 3821` as the counterfactual.

A second loader (`dna_decode/glm/constructed.py`) was then written beside it and performed exactly the
naive join that file exists to prevent, publishing 3,821 as the substrate size for a day.

**THE PROTECTION WAS INSIDE THE FUNCTION THAT GOT BYPASSED.** `MIN_JOIN_RATE` guards callers of
`load_pairs`. Writing a second reader of the same files is the one action that defeats every guard
attached to the first one — and it is also the action that feels like making progress.

So this guard is attached to the DATA instead: a substrate's accession filenames may be named in
exactly ONE module. A second reader has to name the files to open them, so it trips this test instead
of silently re-deriving a solved trap.

**Measured context for why "I couldn't find it" is not the explanation:** of ten plausible search
terms, NINE find `expression.py` (`sigma70`, `GSE108535`, `hyphen`, `underscore`, `norm_key`,
`barcode_mapping`, `join_rate`, `naive_join`, `loader`). The only miss was `constructed` — the word
chosen for the duplicate module. The search was conducted in the vocabulary of the CONCEPT being
built, not of the DATA being loaded.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parent.parent / "dna_decode"

#: A GEO series accession, bare. **Deliberately NOT `GSE\d+_<filename>`.** The first draft of this
#: guard required the underscore-plus-filename form, and measuring it showed `acc = "GSE108535"`
#: followed by an f-string evaded detection entirely — a module can name the accession once and build
#: every filename from it. Matching the bare accession closes that; the cost is that a PROSE mention
#: in a third module also trips, which `ALLOWED_MENTIONS` handles as a visible, deliberate exemption
#: rather than a silent hole.
#: NOTE the missing trailing `\b` — it is wrong here and that was measured, not reasoned. `\bGSE\d+\b`
#: fails to match `GSE108535_sigma70_variant_data.txt.gz` because `_` is a WORD character, so there is
#: no boundary after the digits: the anchored form matched the bare accession and missed every actual
#: filename, i.e. exactly backwards from what is needed.
_ACCESSION = re.compile(r"\bGSE\d+")

#: Modules permitted to NAME an accession without owning it — prose, cross-references, or this test's
#: own registry. Adding a row here is the intended friction: it is a one-line admission that a module
#: talks about a substrate it does not read, and it is reviewable in a diff.
#: Each row is a module that TALKS ABOUT a substrate it does not read. Verified by reading the hits,
#: not by assuming: both are module-level prose cross-references, and the prose is load-bearing (the
#: genomewide docstring is where the 0.7925 sequence-model ceiling and "GSE108535 could not give us
#: one" are recorded). A guard that fired on these would be the documented trap where four guards
#: fired on the comment explaining what they banned.
ALLOWED_MENTIONS: dict[str, frozenset[str]] = {
    # 6 docstring mentions contrasting the two substrates' ceilings and traps; reads no GSE108535 file
    "GSE108535": frozenset({"glm/genomewide.py"}),
    # 1 comment recording the measured <=150 bp substrate bound; reads no GSE144621 file
    # glm/lentimpra.py: 1 comment (line 47) saying its human cache is namespace-separate "so a human
    # table can never be mistaken for a GSE144621 one". VERIFIED BY READING THE HITS, per this file's
    # own standard: exactly one `GSE144621` occurrence, inside a `#:` comment; zero `GSE\d+_<file>`
    # patterns; no import of genomewide; no `.txt.gz`; it opens only
    # `Table S3/S6 - ....xlsx` from `D:/dna_decode_cache/mpra_human/`.
    # Worth keeping in the record: the comment that tripped this detector is the comment explaining
    # why the two substrates CANNOT be confused -- the documented shape where a guard fires on the
    # prose that prevents the thing it guards against.
    # RESIDUAL, named rather than hidden: this exemption is per-MODULE, so if GSE144621-reading code
    # is ever added to lentimpra.py the detector will no longer catch it. Same property as the
    # encoder.py row above. The exemption buys legibility, not immunity.
    "GSE144621": frozenset({"glm/encoder.py", "glm/lentimpra.py"}),
}

#: Accession -> (owning module relative to dna_decode/, the trap that owner's loader encodes).
#: The REASON is load-bearing: when this test fails, the message must tell the next session what the
#: owning loader already solves, or they will solve it again -- which is the whole failure mode.
OWNERS: dict[str, tuple[str, str]] = {
    "GSE108535": (
        "glm/expression.py",
        "the UP-element name is hyphenated in the expression table and underscored in the barcode "
        "map, so a naive join silently returns 3,821 of 10,898 rows and looks healthy. "
        "`_norm_key` fixes it (protecting the `->` mutation arrow) and `load_pairs` enforces "
        "MIN_JOIN_RATE = 0.95. Call `expression.load_pairs`; do not re-read these files.",
    ),
    "GSE144621": (
        "glm/genomewide.py",
        "the fragment files are SPACE-delimited and the peak-tile file is TAB-delimited, so a global "
        "split parses one of them to zero rows -- which reads as 'no signal', not 'no data'. "
        "`sniff_delimiter` decides per file. Call `genomewide.load_fragments` / `load_peak_tiles`.",
    ),
}


def _accession_references() -> dict[str, set[str]]:
    """accession -> set of module paths (relative to dna_decode/) that NAME one of its files."""
    found: dict[str, set[str]] = {}
    for py in sorted(PKG.rglob("*.py")):
        text = py.read_text(encoding="utf-8", errors="replace")
        for hit in _ACCESSION_REFS(text):
            found.setdefault(hit, set()).add(py.relative_to(PKG).as_posix())
    return found


def _ACCESSION_REFS(text: str) -> set[str]:
    return set(_ACCESSION.findall(text))


# --------------------------------------------------------------------------------------------------
# the guard
# --------------------------------------------------------------------------------------------------
def test_each_substrate_accession_is_read_by_exactly_ONE_module():
    """THE load-bearing test. A second module naming the same accession's files is a second loader,
    and a second loader does not inherit the first one's guards."""
    refs = _accession_references()
    for accession, modules in sorted(refs.items()):
        if accession not in OWNERS:
            continue
        owner, trap = OWNERS[accession]
        extra = sorted(modules - {owner} - ALLOWED_MENTIONS.get(accession, frozenset()))
        assert not extra, (
            f"{accession} is named in {extra} as well as its owner {owner}.\n"
            f"A second reader of this substrate does NOT inherit the owner's guards, and the one "
            f"it will re-derive is:\n  {trap}\n"
            f"If a new consumer needs this data, CALL the owner rather than opening the files.\n"
            f"If one of these only MENTIONS the accession in prose and reads nothing, add it to "
            f"ALLOWED_MENTIONS['{accession}'] in this file — deliberately a visible one-line "
            f"exemption, not a silent pass."
        )


def test_every_declared_owner_actually_names_its_accession():
    """Guards the guard. If an owner is renamed or the constants move, the entry above silently stops
    covering anything -- a registry pointing at nothing passes the test above vacuously."""
    refs = _accession_references()
    for accession, (owner, _) in sorted(OWNERS.items()):
        assert accession in refs, f"{accession} is no longer named anywhere in the package"
        assert owner in refs[accession], (
            f"declared owner {owner} does not name {accession}; the registry is stale and this "
            f"guard is covering nothing for that substrate"
        )


def test_the_detector_ACTUALLY_DETECTS_a_second_reader():
    """Non-vacuity, run against synthetic text rather than by planting a real file.

    A guard that cannot see the thing it bans is worse than no guard: it reads as a clean bill."""
    assert _ACCESSION_REFS('VARIANT_FILE = "GSE108535_sigma70_variant_data.txt.gz"') == {"GSE108535"}
    assert _ACCESSION_REFS('f = "GSE144621_peak_tile_expression_formatted_std.txt.gz"') == {"GSE144621"}
    # the real duplicate that shipped named BOTH files
    dup = 'VARIANT_FILE = "GSE108535_sigma70_variant_data.txt.gz"\nMAPPING_FILE = "GSE108535_barcode_mapping.txt.gz"'
    assert _ACCESSION_REFS(dup) == {"GSE108535"}
    # substrate prose carrying no accession must NOT trip it
    assert _ACCESSION_REFS("the sigma70 designed library is a near-complete factorial") == set()


def test_the_failure_message_names_the_owner_AND_the_trap():
    """When this fires, the next session must get the ANSWER, not just 'duplicate found'. The
    duplicate was written because the trap's solution was unknown to the writer, so a message that
    only reports a collision leaves them to re-derive it."""
    for accession, (owner, trap) in OWNERS.items():
        assert owner.endswith(".py") and "/" in owner, accession
        assert len(trap) > 80, f"{accession}: the trap description must be actionable, not a label"
        # it must name a callable entry point, or 'call the owner' is not actionable
        assert "load_" in trap, f"{accession}: name the function a new consumer should call"


# --------------------------------------------------------------------------------------------------
# honest limits, pinned so they are not forgotten
# --------------------------------------------------------------------------------------------------
def test_what_this_guard_catches_is_MEASURED_not_asserted():
    """The first draft of this file asserted that a bare `"GSE108535"` string was caught. It was NOT
    -- the regex then required `GSE\\d+_<filename>`, so a module naming the accession once and
    building filenames from it evaded the guard completely. That claim was written without running
    it, which is the same error this whole guard exists to prevent, recurring inside the fix.

    The regex is now the bare accession, and these cases are each measured."""
    # a full literal filename -- the measured real-world duplicate
    assert _ACCESSION_REFS('"GSE108535_sigma70_variant_data.txt.gz"') == {"GSE108535"}
    # a bare accession plus a dynamic filename -- the hole the first draft missed, now CLOSED
    assert _ACCESSION_REFS('acc = "GSE108535"; p = f"{acc}_sigma70_variant_data.txt.gz"') == {
        "GSE108535"
    }
    # a prose mention is ALSO caught -- deliberate; ALLOWED_MENTIONS is the visible exemption
    assert _ACCESSION_REFS("# substrate is GSE144621, see the census") == {"GSE144621"}


def test_the_REMAINING_hole_is_recorded_not_claimed_closed():
    """Still evades: a module that never writes the accession digits at all -- importing the owner's
    filename constant, or composing `f"GSE{series}_..."` from a number held elsewhere.

    Recorded rather than papered over. Closing it needs import-graph or AST analysis of what gets
    handed to `gzip.open`, which is a different and much more brittle instrument; the measured case
    is a literal copied between modules, and that is now caught twice over. Claiming more would be
    the over-claim this repo keeps finding in its own `--help` strings."""
    assert _ACCESSION_REFS('p = f"GSE{series}_{body}.txt.gz"') == set()
    assert _ACCESSION_REFS("from dna_decode.glm.expression import EXPR_FILE") == set()


@pytest.mark.parametrize("accession", sorted(OWNERS))
def test_owner_module_exists_on_disk(accession):
    owner, _ = OWNERS[accession]
    assert (PKG / owner).is_file(), f"{owner} (owner of {accession}) is missing"
