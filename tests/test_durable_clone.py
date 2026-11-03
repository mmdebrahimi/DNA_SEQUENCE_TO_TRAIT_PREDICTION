"""Tests for the durable-clone mirror.

The point of this script is not copying -- `shutil` does that -- it is being able to SHOW the copy
matches. So the load-bearing tests are the ones that prove a bad mirror is detected:

* a corrupted destination byte must be caught, even when the file SIZE still matches, because the
  cheap pre-filter is size-based and hashing is what closes it;
* `verification_summary` must require EVERY file, since a percentage would let a partial mirror read
  as a qualified success and for a checkpoint that is just a failure;
* and the mirror must be USABLE, not merely byte-equal -- `load_atlas` is pointed at it and the
  parameter count re-checked, which is the only assertion that proves the thing the arms need.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from dna_decode.glm.atlas_embed import DEFAULT_PA_REPO, EXPECTED_PARAMS
from scripts import durable_clone as DC

REPO = pathlib.Path(__file__).resolve().parents[1]
ARTIFACT = REPO / "wiki/durable_clone_2026-11-03.json"
SRC = pathlib.Path(DEFAULT_PA_REPO)
_SRC_PRESENT = (SRC / DC.BASE_CHECKPOINT).is_file() and (SRC / DC.REQUIRED_SOURCE).is_file()
needs_src = pytest.mark.skipif(not _SRC_PRESENT,
                               reason="source clone absent (volatile temp dir has gone)")


# --------------------------------------------------------------------------------------------------
# synthetic tree
# --------------------------------------------------------------------------------------------------
def _tree(root: pathlib.Path, *, bad_pt: bool = False) -> pathlib.Path:
    """A miniature clone with the same SHAPE as the real one.

    The `.pt` files are REAL torch archives, not placeholder bytes, so the static pickle audit
    genuinely passes on the happy path. That matters: with placeholder bytes every audit failed, so
    a test asserting "exit 3 on corruption" would have returned 3 whatever it was given.

    `bad_pt=True` adds one deliberately unreadable checkpoint for the audit-failure guard.
    """
    torch = pytest.importorskip("torch", reason="the synthetic tree writes real torch archives")
    (root / "trained_weights/base_model").mkdir(parents=True)
    (root / "src/promoter_atlas/models").mkdir(parents=True)
    torch.save({"model_state_dict": {"w": torch.zeros(3)}},
               root / "trained_weights/base_model/promoteratlas-base.pt")
    torch.save({"model_state_dict": {"w": torch.ones(2)}}, root / "trained_weights/head.pt")
    if bad_pt:
        (root / "trained_weights/truncated.pt").write_bytes(b"NOT-AN-ARCHIVE")
    (root / DC.REQUIRED_SOURCE).write_text("import torch\n", encoding="utf-8")
    # noise that must NOT be mirrored
    (root / "src/promoter_atlas/models/__pycache__").mkdir()
    (root / "src/promoter_atlas/models/__pycache__/x.cpython-311.pyc").write_bytes(b"\x00")
    (root / "src/promoter_atlas/stale.pyc").write_bytes(b"\x00")
    return root


# --------------------------------------------------------------------------------------------------
# copy_manifest -- pure
# --------------------------------------------------------------------------------------------------
def test_the_manifest_is_DETERMINISTIC_and_puts_trained_weights_first(tmp_path):
    """Ordering is per-subtree, NOT globally lexicographic, and that is deliberate.

    `trained_weights` comes first because it is the part whose integrity the arms' numbers depend
    on. I first asserted global sortedness here and it failed — the looser word "sorted" in the
    docstring is what misled me, and that wording is now fixed too.
    """
    root = _tree(tmp_path / "src")
    a = DC.copy_manifest(root)
    b = DC.copy_manifest(root)
    assert a == b, "the manifest must be deterministic or the mirror is not reproducible"

    rels = [r["relpath"] for r in a]
    assert "trained_weights/base_model/promoteratlas-base.pt" in rels
    assert DC.REQUIRED_SOURCE in rels

    # sorted WITHIN each subtree
    for sub in DC.SUBTREES:
        inside = [r for r in rels if r.startswith(sub + "/")]
        assert inside == sorted(inside), "%s is not internally sorted" % sub
    # and the weights subtree precedes the source subtree
    assert rels.index("trained_weights/base_model/promoteratlas-base.pt") < rels.index(DC.REQUIRED_SOURCE)


def test_the_manifest_EXCLUDES_bytecode_and_pycache(tmp_path):
    """Mirroring compiled bytecode from another interpreter is noise at best."""
    rels = [r["relpath"] for r in DC.copy_manifest(_tree(tmp_path / "src"))]
    assert not any(r.endswith(".pyc") for r in rels)
    assert not any("__pycache__" in r for r in rels)


def test_the_manifest_carries_real_sizes(tmp_path):
    root = _tree(tmp_path / "src")
    rows = {r["relpath"]: r["bytes"] for r in DC.copy_manifest(root)}
    for rel, size in rows.items():
        assert size == (root / rel).stat().st_size, rel
    assert rows["trained_weights/base_model/promoteratlas-base.pt"] > 0


def test_an_absent_subtree_is_REFUSED_rather_than_silently_shortening_the_manifest(tmp_path):
    """A short manifest mirrors a partial clone and reports success."""
    root = tmp_path / "src"
    (root / "trained_weights").mkdir(parents=True)
    (root / "trained_weights/x.pt").write_bytes(b"x")
    with pytest.raises(DC.DurableCloneError, match="source subtree absent"):
        DC.copy_manifest(root)


def test_an_empty_manifest_is_refused(tmp_path):
    root = tmp_path / "src"
    for sub in DC.SUBTREES:
        (root / sub).mkdir(parents=True)
    with pytest.raises(DC.DurableCloneError, match="manifest is empty"):
        DC.copy_manifest(root)


# --------------------------------------------------------------------------------------------------
# verification_summary -- pure, and it must demand EVERY file
# --------------------------------------------------------------------------------------------------
def test_the_summary_is_ok_only_when_every_file_matched():
    rows = [{"relpath": "a", "bytes": 1, "digest_matches": True},
            {"relpath": "b", "bytes": 2, "digest_matches": True}]
    s = DC.verification_summary(rows)
    assert s["ok"] is True and s["n_matched"] == 2 and s["total_bytes"] == 3


def test_ONE_mismatch_makes_the_whole_summary_NOT_ok():
    rows = [{"relpath": "a", "digest_matches": True},
            {"relpath": "b", "digest_matches": False}]
    s = DC.verification_summary(rows)
    assert s["ok"] is False
    assert s["mismatched"] == ["b"]
    assert "percentage" in s["rule"], "the reason a percentage is refused must travel with the verdict"


def test_an_UNCHECKED_file_also_blocks_ok_rather_than_being_assumed_fine():
    rows = [{"relpath": "a", "digest_matches": True}, {"relpath": "b", "digest_matches": None}]
    s = DC.verification_summary(rows)
    assert s["ok"] is False and s["unchecked"] == ["b"]


def test_an_empty_summary_is_not_ok():
    assert DC.verification_summary([])["ok"] is False


# --------------------------------------------------------------------------------------------------
# needs_copy -- the cheap pre-filter, and the gap hashing closes
# --------------------------------------------------------------------------------------------------
def test_needs_copy_is_true_for_an_absent_or_differently_sized_destination(tmp_path):
    s, d = tmp_path / "s", tmp_path / "d"
    s.write_bytes(b"hello")
    assert DC.needs_copy(s, d) is True
    d.write_bytes(b"hi")
    assert DC.needs_copy(s, d) is True
    d.write_bytes(b"hello")
    assert DC.needs_copy(s, d) is False


def test_needs_copy_CANNOT_see_a_same_size_corruption_which_is_why_hashing_follows(tmp_path):
    """Pins the pre-filter's known blind spot, so nobody later 'optimises away' the hash pass."""
    s, d = tmp_path / "s", tmp_path / "d"
    s.write_bytes(b"AAAA")
    d.write_bytes(b"BBBB")
    assert DC.needs_copy(s, d) is False, "same size -- the cheap filter is blind here by design"


# --------------------------------------------------------------------------------------------------
# mirror -- end to end on the synthetic tree
# --------------------------------------------------------------------------------------------------
def test_mirror_copies_then_verifies_every_file(tmp_path):
    src = _tree(tmp_path / "src")
    rec = DC.mirror(src_root=src, dest_root=tmp_path / "dest", verbose=False)
    assert rec["verification"]["ok"] is True
    assert rec["n_copied"] == rec["verification"]["n_files"] > 0
    assert rec["required_source_present"] is True
    assert rec["PA_REPO"] == str(tmp_path / "dest")


def test_a_SECOND_mirror_copies_nothing_and_still_verifies(tmp_path):
    src = _tree(tmp_path / "src")
    DC.mirror(src_root=src, dest_root=tmp_path / "dest", verbose=False)
    again = DC.mirror(src_root=src, dest_root=tmp_path / "dest", verbose=False)
    assert again["n_copied"] == 0
    assert again["idempotent"] is True
    assert again["verification"]["ok"] is True, "idempotent must not mean unverified"


def test_a_CORRUPTED_destination_byte_is_DETECTED_even_at_the_same_size(tmp_path):
    """The whole reason the mirror hashes both sides instead of trusting sizes."""
    src = _tree(tmp_path / "src")
    dest = tmp_path / "dest"
    DC.mirror(src_root=src, dest_root=dest, verbose=False)

    victim = dest / "trained_weights/base_model/promoteratlas-base.pt"
    original = victim.read_bytes()
    # flip ONE byte in place, preserving length, so the size-based pre-filter stays blind and only
    # the hash pass can catch it
    buf = bytearray(original)
    mid = len(buf) // 2
    buf[mid] ^= 0xFF
    victim.write_bytes(bytes(buf))
    assert len(victim.read_bytes()) == len(original)
    assert victim.read_bytes() != original
    assert DC.needs_copy(src / "trained_weights/base_model/promoteratlas-base.pt", victim) is False, (
        "the cheap pre-filter must be blind here, or this test is not exercising the hash pass"
    )

    rec = DC.mirror(src_root=src, dest_root=dest, verbose=False)
    assert rec["verification"]["ok"] is False
    assert "trained_weights/base_model/promoteratlas-base.pt" in rec["verification"]["mismatched"]


def test_the_corruption_check_is_non_vacuous_the_same_mirror_passes_when_restored(tmp_path):
    src = _tree(tmp_path / "src")
    dest = tmp_path / "dest"
    DC.mirror(src_root=src, dest_root=dest, verbose=False)
    assert DC.mirror(src_root=src, dest_root=dest, verbose=False)["verification"]["ok"] is True


def test_an_UNREADABLE_checkpoint_is_a_FAILED_AUDIT_not_an_escaping_exception(tmp_path):
    """Regression guard for a defect found while writing these tests.

    `audit_checkpoint` opens a `.pt` as a zip archive, so a truncated or non-archive file raises
    `zipfile.BadZipFile` — NOT `AtlasEmbedError`. The first version caught only the latter, so a
    corrupt checkpoint crashed the mirror with a zipfile traceback instead of recording that the
    audit failed. This project has already been bitten by files of the right NAME and the wrong
    CONTENTS, so an unreadable checkpoint must be a verdict, never a crash.
    """
    src = _tree(tmp_path / "src", bad_pt=True)
    rec = DC.mirror(src_root=src, dest_root=tmp_path / "dest", verbose=False)

    audits = rec["checkpoint_audits"]
    assert len(audits) == 3, "two real archives plus the deliberately unreadable one"
    bad = {k: v for k, v in audits.items() if v["passed"] is False}
    assert list(bad) == ["trained_weights/truncated.pt"], (
        "only the unreadable file may fail; if the real archives also failed, this test would pass "
        "for the wrong reason"
    )
    assert bad["trained_weights/truncated.pt"]["kind"] == "unreadable_or_not_an_archive"
    assert "BadZipFile" in bad["trained_weights/truncated.pt"]["error"]

    # the mirror's overall verdict reflects it rather than reporting success
    assert rec["all_checkpoints_audited"] is False

    # ...while the BYTE-level verification is INDEPENDENT and still passes: the unreadable file was
    # copied faithfully. A failed audit and a failed copy are different claims.
    assert rec["verification"]["ok"] is True


def test_mirror_refuses_an_absent_source_clone(tmp_path):
    with pytest.raises(DC.DurableCloneError, match="source clone absent"):
        DC.mirror(src_root=tmp_path / "nothing", dest_root=tmp_path / "dest")


def test_mirror_refuses_a_source_missing_the_one_file_load_atlas_needs(tmp_path):
    src = _tree(tmp_path / "src")
    (src / DC.REQUIRED_SOURCE).unlink()
    with pytest.raises(DC.DurableCloneError, match="missing"):
        DC.mirror(src_root=src, dest_root=tmp_path / "dest")


def test_main_returns_nonzero_and_writes_nothing_on_an_absent_source(tmp_path):
    out = tmp_path / "x.json"
    assert DC.main(["--src", str(tmp_path / "nope"), "--out", str(out)]) == 2
    assert not out.exists()


def test_main_returns_3_on_a_verification_failure_and_writes_nothing(tmp_path):
    """A failed mirror must not leave an artifact that reads like a success.

    **Non-vacuity note.** A synthetic tree can never reach exit 0, because the base checkpoint
    cannot match the REAL pinned sha256 — so "exit 3" alone would be uninformative here. What this
    isolates is the corruption's own effect: byte-verification passes before and fails after, with
    the mismatch naming the file. Exit 0 is covered against the real mirror, which does match the
    pin, by the committed-artifact test below.
    """
    src = _tree(tmp_path / "src")
    dest = tmp_path / "dest"
    before = DC.mirror(src_root=src, dest_root=dest, verbose=False)
    assert before["verification"]["ok"] is True

    victim = dest / "trained_weights/head.pt"
    victim.write_bytes(b"X" * len(victim.read_bytes()))         # same length, corrupted content

    after = DC.mirror(src_root=src, dest_root=dest, verbose=False)
    assert after["verification"]["ok"] is False
    assert after["verification"]["mismatched"] == ["trained_weights/head.pt"]

    out = tmp_path / "x.json"
    assert DC.main(["--src", str(src), "--dest", str(dest), "--out", str(out)]) == 3
    assert not out.exists(), "a refused mirror must leave no artifact behind"


# --------------------------------------------------------------------------------------------------
# the real mirror -- the assertion that proves it is USABLE, not merely byte-equal
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not ARTIFACT.is_file(), reason="mirror has not been made")
def test_the_committed_artifact_records_a_fully_verified_mirror():
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert rec["ok"] is True
    assert rec["verification"]["ok"] is True
    assert rec["verification"]["mismatched"] == [] and rec["verification"]["unchecked"] == []
    assert rec["pinned_base_checkpoint"]["matches_pin"] is True
    assert rec["pinned_base_checkpoint"]["sha256"].startswith("bb413bc8602f05b6")
    assert rec["pinned_base_checkpoint"]["bytes"] == 22_313_814
    assert rec["all_checkpoints_audited"] is True
    assert rec["required_source_present"] is True


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="mirror has not been made")
def test_the_artifact_states_the_C_drive_inversion_with_its_measured_evidence():
    inv = json.loads(ARTIFACT.read_text(encoding="utf-8"))["target_is_C_drive_deliberately"]
    assert "D:" in inv["inverts"]
    assert "0xc0000006" in inv["why"], "the death signature must be a figure, not an adjective"
    assert "rc=0" in inv["why"], "the one survivor is the positive half of the evidence"
    why_ok = inv["why_the_rule_is_not_violated"]
    assert "HEAVY artifacts" in why_ok, "the rule's PURPOSE is what the exception turns on"
    assert "bulk caches stay on D:" in why_ok, "the rule must be shown still in force where it aims"
    assert "NOT-proven" in inv["confidence"], "this must not read as a guarantee"


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="mirror has not been made")
def test_the_MIRROR_IS_USABLE_load_atlas_loads_from_it_and_the_param_count_holds():
    """Byte-equality is not the property the arms need; loadability is.

    `load_atlas` refuses on a parameter-count mismatch, so pointing it at the mirror and getting a
    model back is a stronger statement than any digest comparison.
    """
    pytest.importorskip("torch", reason="loading the encoder needs torch")
    from dna_decode.glm.atlas_embed import load_atlas

    dest = pathlib.Path(json.loads(ARTIFACT.read_text(encoding="utf-8"))["PA_REPO"])
    if not (dest / DC.BASE_CHECKPOINT).is_file():
        pytest.skip("mirror directory has gone")
    model = load_atlas(dest)
    assert sum(p.numel() for p in model.parameters()) == EXPECTED_PARAMS


@needs_src
def test_the_real_source_clone_contains_what_the_manifest_requires():
    rels = {r["relpath"] for r in DC.copy_manifest(SRC)}
    assert DC.BASE_CHECKPOINT in rels
    assert DC.REQUIRED_SOURCE in rels
    assert sum(1 for r in rels if r.endswith(".pt")) >= 6, (
        "base + five expression heads + segmentation were measured present"
    )
