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


# --- the staleness guard ----------------------------------------------------------------------

def test_a_note_for_a_dataset_that_does_not_exist_fails_self_check():
    """The curated layer is the only hand-written part, so it is the only part that can rot. A note
    naming a path that is gone is the signal, and it must be a GATE, not a warning."""
    rep = {"stale_notes": ["a_dataset_that_was_deleted"], "counts": {"datasets_present": 5}}
    ok, problems = dinv.self_check(rep)
    assert ok is False
    assert any("do not exist" in p for p in problems)


def test_an_empty_scan_fails_self_check():
    """A scan that found nothing would render an inventory claiming the project has no data -- the
    most confidently wrong output this tool could produce."""
    ok, problems = dinv.self_check({"stale_notes": [], "counts": {"datasets_present": 0}})
    assert ok is False
    assert any("vacuous" in p for p in problems)


def test_a_clean_report_passes_self_check():
    ok, problems = dinv.self_check({"stale_notes": [], "counts": {"datasets_present": 5}})
    assert ok is True and problems == []


# --- index correctness ------------------------------------------------------------------------

def test_the_index_tokenises_and_does_not_substring_match(tmp_path):
    """Substring matching would make `cache` match `embeddings_cache`, `refs`, `torch_cache`... and
    every dataset would look universally used, which is worse than no signal at all."""
    src = tmp_path / "s.py"
    src.write_text('P = "data/processed/embeddings_cache/x.h5"\n', encoding="utf-8")
    idx = dinv.build_reverse_index({"cache", "embeddings_cache"}, [tmp_path], {".py"})
    assert idx["embeddings_cache"], "the exact directory name must match"
    assert not idx["cache"], "a substring of another dataset's name must NOT count as a reference"


def test_a_hardcoded_posix_path_resolves_to_its_directory_key(tmp_path):
    """The dominant real form: paths appear as string literals with forward slashes."""
    src = tmp_path / "s.py"
    src.write_text('ap.add_argument("--x", default="data/raw/ar_bank_caur_extval_micafungin/sel.tsv")\n',
                   encoding="utf-8")
    idx = dinv.build_reverse_index({"ar_bank_caur_extval_micafungin"}, [tmp_path], {".py"})
    assert idx["ar_bank_caur_extval_micafungin"]


# --- missing roots are information, not failure -----------------------------------------------

def test_an_unreachable_root_is_recorded_not_silently_dropped(tmp_path, monkeypatch):
    """D: is an external drive. If it is unplugged, the honest report is 'this root was unreachable',
    NOT an inventory that quietly claims those datasets do not exist."""
    monkeypatch.setattr(dinv, "SCAN_ROOTS", [("cache", tmp_path / "nope", "external_cache")])
    ds = dinv.scan_datasets()
    assert len(ds) == 1 and ds[0].exists is False
    rep = dinv.build_report(ds, [])
    assert rep["counts"]["roots_missing"] == 1
    assert rep["roots_missing"][0]["root"] == "cache"


# --- the committed artifact -------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not generated")
def test_the_committed_artifact_reports_its_own_coverage():
    """An inventory that hid its gaps would be the failure it exists to prevent."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    c = d["counts"]
    assert c["documented"] + c["undocumented"] == c["datasets_present"]
    assert c["undocumented"] > 0, "if everything is documented, verify that rather than assuming it"
    assert d["derived_fields"] and d["curated_fields"], "the split must be explicit to a reader"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not generated")
def test_the_orphan_signal_is_alive_in_the_committed_artifact():
    """Guards the regression directly: if a future change re-contaminates the index, orphans collapse
    toward zero and this fails. 59 of 248 was the measured state at first generation."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    n_orphan = d["counts"]["orphans_no_reference_anywhere"]
    assert n_orphan > 20, (
        f"only {n_orphan} orphans; the contamination bug collapsed this from 59 to 8 while still "
        f"printing a successful run -- check SELF_REFERENCES before accepting this number")
