"""Guards for the derived data inventory.

The tool's whole value is answering "has anything already used this dataset?" correctly. The tests that
matter are therefore the ones that pin the ways that answer can go quietly WRONG -- above all the
self-contamination invariant, which broke FOUR times during the build and never once announced itself:
every run printed a clean success while the orphan signal drained away.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Imported as a package member, not via spec_from_file_location: a dynamically loaded module is not in
# sys.modules, and @dataclass resolves annotations through sys.modules[cls.__module__] -- so the
# file-loader form dies at import with a bare AttributeError on NoneType.
from scripts import data_inventory as dinv  # noqa: E402

ARTIFACT = ROOT / "wiki" / "data_inventory.json"


# --- the self-contamination invariant --------------------------------------------------------

def test_describing_a_dataset_must_not_make_it_look_consumed(tmp_path, monkeypatch):
    """THE invariant. It broke four times, each worse than the last:

      1. the emitted artifacts land in an artifact root and name every dataset -> run #2 showed
         everything as referenced (orphans 59 -> 8);
      2. the curated notes module lives in a code root, so writing "this is unused" marked it used;
      3. the scanner itself, because a comment explaining (2) NAMED an orphan while doing so;
      4. this test file, because the guard against naming datasets had to name two of them.

    All four share one shape: prose ABOUT a dataset is lexically identical to a reference TO it.
    """
    ds_root = tmp_path / "data" / "raw"
    (ds_root / "forgotten_cohort").mkdir(parents=True)
    (ds_root / "forgotten_cohort" / "x.tsv").write_text("a\n", encoding="utf-8")

    # A file that only DESCRIBES the dataset, living in a scanned root.
    describer = tmp_path / "wiki" / "data_inventory.md"
    describer.parent.mkdir(parents=True, exist_ok=True)
    describer.write_text("| `forgotten_cohort` | nothing uses this |\n", encoding="utf-8")

    monkeypatch.setattr(dinv, "ROOT", tmp_path)
    monkeypatch.setattr(dinv, "SCAN_ROOTS", [("raw", ds_root, "cohort")])
    monkeypatch.setattr(dinv, "CODE_ROOTS", [])
    monkeypatch.setattr(dinv, "ARTIFACT_ROOTS", [tmp_path / "wiki"])

    guarded = dinv.attach_uses(dinv.scan_datasets())
    assert guarded[0].orphan is True, "a dataset mentioned ONLY by the inventory's own output must stay orphan"

    # Non-vacuity: with the guard removed the same input produces the WRONG answer.
    monkeypatch.setattr(dinv, "SELF_REFERENCE_PATTERNS", ())
    unguarded = dinv.attach_uses(dinv.scan_datasets())
    assert unguarded[0].orphan is False, (
        "if this passes, the fixture never exercised the contamination path and the test above proves "
        "nothing")


def test_every_inventory_owned_file_is_excluded_including_ones_not_yet_written():
    """The exclusion is DERIVED, not enumerated. Four separate regressions came from a hand-listed set
    that a newly-added inventory file was missing from, so the test that matters is that a file which
    does not exist yet would already be covered."""
    for existing in ("wiki/data_inventory.md", "wiki/data_inventory.json",
                     "dna_decode/data/inventory_notes.py", "scripts/data_inventory.py",
                     "tests/test_data_inventory.py"):
        assert dinv.is_self_reference(existing), existing
    # The whole point: a future file nobody remembered to register.
    assert dinv.is_self_reference("scripts/data_inventory_v2.py")
    assert dinv.is_self_reference("wiki/data_inventory_2027-01-01.json")
    # ...without swallowing unrelated files, which would hide real consumers.
    assert not dinv.is_self_reference("scripts/salmserovar_validate.py")
    assert not dinv.is_self_reference("wiki/decoder_validation_report_card.json")


def test_the_scanner_module_never_names_a_real_dataset_key():
    """Belt-and-braces for defect (3). The scanner is excluded by pattern now, so this cannot corrupt
    the index today -- it pins the HABIT, so that if the exclusion is ever narrowed the blast radius is
    already zero. Keys are read from the curated notes rather than written literally here, because
    hardcoding them would reintroduce defect (4) in the very test that guards against it."""
    from dna_decode.data.inventory_notes import NOTES
    src = (ROOT / "scripts" / "data_inventory.py").read_text(encoding="utf-8")
    # Tokenised with the index's OWN regex, not a substring scan. A guard stricter than the mechanism
    # it protects raises false alarms and gets switched off: `embeddings` inside `embeddings_cache` is
    # not a reference, and the index already knows that.
    toks = set(dinv._TOKEN.findall(src))
    named = sorted(k for k in NOTES if k in toks)
    assert not named, f"the scanner module names real dataset keys as standalone tokens: {named}"


def test_no_INDEX_SCANNED_module_names_a_real_dataset_key_IN_PROSE():
    """The sibling of the test above, for the rest of the scanned tree -- PROSE only.

    `build_reverse_index` scans every CODE_ROOT, so any module whose COMMENTS or
    DOCSTRINGS use a dataset key as a bare token manufactures a false consumer for
    that dataset. Found live 2026-10-19: a new `dna_decode/data/external_root.py`
    used, in a comment, the bare word that is also a dataset key.

    EXECUTABLE literals are deliberately EXEMPT, and the distinction is the whole
    point. The first version of this guard also flagged `refseq` in that same module
    -- but there it is `CACHE_REQUIRED_SUBDIRS = ("refseq",)`, a value that MUST equal
    the real directory name and which makes the module a genuine (if shallow) reader
    of that dataset. Recording that is TRUE. Flagging it would have pushed me to
    rename a load-bearing string to satisfy a test, which is the failure mode of a
    guard stricter than the thing it protects. Same reasoning the test above gives for
    `embeddings` inside `embeddings_cache`.

    Scoped to the modules it actually reads, so it states a fact it checks rather than
    implying tree-wide coverage.
    """
    import ast
    import io
    import tokenize

    from dna_decode.data.inventory_notes import NOTES

    def prose_of(src: str) -> str:
        """Comments + docstrings only -- never an executable string literal."""
        chunks: list[str] = []
        try:
            for tok in tokenize.generate_tokens(io.StringIO(src).readline):
                if tok.type == tokenize.COMMENT:
                    chunks.append(tok.string)
        except (tokenize.TokenError, IndentationError):  # pragma: no cover
            pass
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                doc = ast.get_docstring(node)
                if doc:
                    chunks.append(doc)
        return chr(10).join(chunks)

    offenders = {}
    for rel in ("dna_decode/data/external_root.py", "scripts/data_inventory.py"):
        path = ROOT / rel
        if not path.is_file():
            continue
        toks = set(dinv._TOKEN.findall(prose_of(path.read_text(encoding="utf-8"))))
        named = sorted(k for k in NOTES if k in toks)
        if named:
            offenders[rel] = named
    assert not offenders, (
        f"modules naming real dataset keys as bare tokens IN PROSE: {offenders}. "
        "Rephrase the prose -- do not rename an executable literal to satisfy this."
    )


def test_the_prose_guard_is_not_vacuous():
    """It must still catch the real 2026-10-19 defect: a dataset key bare in a comment."""
    import ast
    import io
    import tokenize

    from dna_decode.data.inventory_notes import NOTES

    key = next(iter(sorted(NOTES)))
    planted = f"# a comment mentioning {key} as a bare word" + chr(10) + "X = 1" + chr(10)
    chunks = [t.string for t in tokenize.generate_tokens(io.StringIO(planted).readline)
              if t.type == tokenize.COMMENT]
    ast.parse(planted)
    toks = set(dinv._TOKEN.findall(chr(10).join(chunks)))
    assert key in toks, "the prose extractor stopped seeing comments; the guard above is now vacuous"
