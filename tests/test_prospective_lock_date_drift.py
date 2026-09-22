"""The prospective sweep's eligibility cutoff must DERIVE from the active lock, never restate it.

FOUND 2026-09-22 while pre-flighting the accrual sweep: `scripts/fetch_prospective_cohort.py` carried
`LOCK_DATE = "2026-06-13"` as a hardcoded module constant and used it as the `--lock-date` default. The
v2 gentamicin lock had moved the cutoff to 2026-08-31 three weeks earlier and the constant did not
follow, so a sweep would have collected isolates dated after the RETIRED v1 cutoff and labelled them
prospective. For the deployed rule those isolates are PRE-lock -- a false prospective claim, which is
precisely the leakage the lock exists to prevent.

This is the same class as the salmserovar CLI that ran a coverage threshold replaced five days earlier
("a shipped default is a promise"), and it survived for the same reason: nothing compared the CLI's
default against the constant it mirrors.

The resolution rule is DEFINITIONAL, not conventional: a manifest pins sha256 hashes of the frozen
surface, so the manifest that still VERIFIES describes the live decoder and one that no longer verifies
describes a retired one. Retiring a lock IS the act of making its manifest stop verifying, so nothing
needs hand-syncing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.eval.prospective_lock import (  # noqa: E402
    NoActiveLock, active_lock_date, resolve_active_lock, verify_lock,
)


def test_exactly_one_committed_manifest_describes_the_live_surface():
    """The v1 manifest must NOT verify -- the gentamicin rescue changed amr_rules.py, which is what
    retiring a lock means. If both verified, the cutoff would be ambiguous."""
    manifests = sorted((ROOT / "wiki").glob("prospective_lock_manifest_*.json"))
    assert len(manifests) >= 2, "need the retired v1 and the live v2 to make this test meaningful"
    ok = [p.name for p in manifests
          if verify_lock(json.loads(p.read_text(encoding="utf-8"))).ok]
    assert ok == ["prospective_lock_manifest_2026-08-31.json"], ok


def test_the_active_lock_date_is_the_v2_cutoff_not_the_retired_one():
    assert active_lock_date() == "2026-08-31"
    assert active_lock_date() != "2026-06-13", "this is the retired v1 cutoff"


def test_the_sweep_default_is_derived_not_hardcoded():
    """THE guard. The source must not carry a literal cutoff to drift -- equality today drifts tomorrow,
    which is exactly what happened."""
    src = (ROOT / "scripts" / "fetch_prospective_cohort.py").read_text(encoding="utf-8")
    assert 'LOCK_DATE = "' not in src, "the cutoff is hardcoded again"
    assert "active_lock_date" in src, "the sweep no longer derives its cutoff"
    # The DOCSTRING is the --help text, so a stale date asserted there is the same defect one layer up
    # (it told a reader the cutoff was 2026-06-13 for three weeks after it moved). Scoped to the
    # docstring rather than the whole file on purpose: banning the string everywhere would also flag
    # the comment EXPLAINING the history and an unrelated lexicographic-comparison example, and a rule
    # that fires on correct code gets suppressed.
    doc = src.split('"""')[1]
    assert "2026-06-13" not in doc, "the help text asserts the retired v1 cutoff"


def test_the_sweep_arg_default_is_none_so_derivation_can_happen():
    """A non-None argparse default would shadow the derivation silently -- the salmserovar failure mode
    exactly, where main() passed the CLI default over the validated constant on every call."""
    import argparse
    from unittest.mock import patch
    import scripts.fetch_prospective_cohort as m

    captured = {}
    real = argparse.ArgumentParser.add_argument

    def spy(self, *a, **kw):
        if a and a[0] == "--lock-date":
            captured["default"] = kw.get("default", "MISSING")
        return real(self, *a, **kw)

    with patch.object(argparse.ArgumentParser, "add_argument", spy):
        with patch.object(argparse.ArgumentParser, "parse_args", side_effect=SystemExit):
            with pytest.raises(SystemExit):
                m.main([])
    assert captured.get("default") is None, captured


def test_resolution_fails_closed_when_nothing_verifies(tmp_path):
    """Falling back to the newest file, or to any default, is how a stale cutoff silently becomes a
    prospective claim. With no manifest at all the resolver must REFUSE."""
    (tmp_path / "wiki").mkdir()
    with pytest.raises(NoActiveLock):
        resolve_active_lock(repo=tmp_path)


def test_the_resolver_returns_the_manifest_it_used():
    """The caller must be able to SAY which lock it swept under; a bare date is unauditable."""
    path, manifest = resolve_active_lock()
    assert path.name == "prospective_lock_manifest_2026-08-31.json"
    assert manifest["lock_date"] == "2026-08-31"


# --- a manifest that pins NOTHING must not verify (found in review, 2026-09-22) -------------------

def test_a_manifest_pinning_nothing_does_not_verify():
    """REGRESSION. `verify_lock` iterated `surface_sha256.items()`; an EMPTY dict finds no drift and no
    missing file, so it returned ok=True -- a lock that checks nothing, reported as a lock that holds.

    Harmless while every caller named ONE trusted manifest. NOT harmless once `resolve_active_lock`
    globs every committed manifest and selects a verifying one: a vacuous manifest would be picked
    EXACTLY when the real lock stops verifying, i.e. it would defeat fail-closed at the only moment
    fail-closed matters. Same shape as the vacuous filters this project has caught three times."""
    from dna_decode.eval.prospective_lock import verify_lock
    assert verify_lock({"lock_date": "2099-01-01", "surface_sha256": {}}).ok is False
    assert verify_lock({"lock_date": "2099-01-01"}).ok is False


def test_a_manifest_pinning_only_a_subset_does_not_verify():
    """A partial pin is the dangerous middle case: every hash it DOES carry is correct, so the old
    loop found no drift. The decoder could have changed in any unpinned file."""
    import hashlib
    from dna_decode.eval.prospective_lock import FROZEN_SURFACE_FILES, verify_lock
    one = FROZEN_SURFACE_FILES[2]
    digest = hashlib.sha256((ROOT / one).read_bytes()).hexdigest()
    res = verify_lock({"surface_sha256": {one: digest}})
    assert res.ok is False
    assert res.incomplete_pin, "the reason must be reported, not just the refusal"


def test_incomplete_pin_is_reported_separately_from_missing():
    """`missing` means 'pinned but absent from disk'; `incomplete_pin` means the commitment was never
    made. Collapsing them would report a manifest that pins nothing as a manifest whose files vanished."""
    from dna_decode.eval.prospective_lock import verify_lock
    res = verify_lock({"surface_sha256": {}})
    assert res.incomplete_pin and not res.drifted


def test_the_real_active_manifest_still_verifies_after_the_tightening():
    """NON-VACUITY for the tightening itself: it must reject the vacuous forms WITHOUT rejecting the
    genuine lock. A guard that refuses everything is not fail-closed, it is broken."""
    assert active_lock_date() == "2026-08-31"


# --- ambiguity REFUSES rather than silently picking a cutoff (corrected 2026-09-22) ---------------

def _write_manifest(wiki, name, lock_date, hashes):
    (wiki / name).write_text(json.dumps({"schema": "prospective-lock-manifest-v1",
                                         "lock_date": lock_date, "surface_sha256": hashes}),
                             encoding="utf-8")


def _stage(tmp_path):
    """A tmp repo carrying real copies of the frozen surface, so manifests can genuinely verify."""
    from dna_decode.eval.prospective_lock import FROZEN_SURFACE_FILES, surface_hashes
    (tmp_path / "wiki").mkdir()
    for rel in FROZEN_SURFACE_FILES:
        dst = tmp_path / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes((ROOT / rel).read_bytes())
    return surface_hashes(tmp_path)


def test_two_verifying_manifests_with_DIFFERENT_cutoffs_refuse(tmp_path):
    """THE corrected rule. Taking max(lock_date) looked conservative -- a later cutoff admits fewer
    isolates -- but two manifests can only BOTH verify if they pin identical hashes, i.e. they describe
    the SAME decoder while claiming DIFFERENT cutoffs. That is a protocol contradiction, and silently
    taking the later one would let an arbitrary cutoff DISCARD already-accrued unfavourable evidence
    while looking cautious."""
    h = _stage(tmp_path)
    _write_manifest(tmp_path / "wiki", "prospective_lock_manifest_2026-01-01.json", "2026-01-01", h)
    _write_manifest(tmp_path / "wiki", "prospective_lock_manifest_2026-08-31.json", "2026-08-31", h)
    with pytest.raises(NoActiveLock) as e:
        resolve_active_lock(repo=tmp_path)
    assert "AMBIGUOUS" in str(e.value)


def test_two_verifying_manifests_with_the_SAME_cutoff_are_not_ambiguous(tmp_path):
    """Identical cutoffs make the same commitment, so either answers the question. Refusing here would
    be a false alarm -- and a guard that cries wolf gets suppressed."""
    h = _stage(tmp_path)
    _write_manifest(tmp_path / "wiki", "prospective_lock_manifest_a.json", "2026-08-31", h)
    _write_manifest(tmp_path / "wiki", "prospective_lock_manifest_b.json", "2026-08-31", h)
    path, m = resolve_active_lock(repo=tmp_path)
    assert m["lock_date"] == "2026-08-31" and path.name == "prospective_lock_manifest_a.json"


def test_a_non_verifying_manifest_never_creates_ambiguity(tmp_path):
    """A retired lock coexists with the live one by design -- that is what retiring MEANS. It must not
    trip the ambiguity refusal, or every future lock would break the resolver."""
    h = _stage(tmp_path)
    _write_manifest(tmp_path / "wiki", "prospective_lock_manifest_live.json", "2026-08-31", h)
    _write_manifest(tmp_path / "wiki", "prospective_lock_manifest_retired.json", "2026-01-01",
                    {**h, sorted(h)[0]: "0" * 64})
    assert resolve_active_lock(repo=tmp_path)[1]["lock_date"] == "2026-08-31"


# --- both CLIs derive their manifest, and an override is marked as not-a-live-claim ---------------

def test_the_validator_default_is_derived_not_hardcoded():
    src = (ROOT / "scripts" / "prospective_lock_validate.py").read_text(encoding="utf-8")
    assert "prospective_lock_manifest_2026-08-31.json" not in src, "manifest filename hardcoded again"
    assert "resolve_active_lock" in src


def test_the_fetcher_stamps_manifest_IDENTITY_not_just_the_date():
    """A cohort recording a cutoff but not WHICH lock produced it can be re-interpreted later under a
    different manifest claiming the same date -- the same class as the stale-constant bug."""
    src = (ROOT / "scripts" / "fetch_prospective_cohort.py").read_text(encoding="utf-8")
    assert "lock_provenance" in src and "resolve_active_lock" in src
    for field in ("manifest_path", "frozen_commit", "surface_sha256"):
        assert field in src, field


def test_an_overridden_cutoff_is_marked_as_not_a_live_prospective_claim():
    """An override is not backed by a verified manifest. The console says so, but the console scrolls
    away; the durable artifact has to carry it."""
    src = (ROOT / "scripts" / "fetch_prospective_cohort.py").read_text(encoding="utf-8")
    assert "historical_reproduction" in src
    assert "overridden_on_command_line" in src


# --- hash agreement alone does not make a file a LOCK (found in review, 2026-09-22) ---------------

def test_a_manifest_with_correct_hashes_but_no_schema_is_not_selectable(tmp_path):
    """VERIFIED before fixing: a manifest carrying the right hashes, no schema and
    `lock_date: "not-a-date"` WAS selected as the active lock -- and that bogus cutoff would have gone
    straight into the eligibility filter. Verifying asks "does this describe the live decoder";
    well-formedness asks "is this a commitment at all"."""
    h = _stage(tmp_path)
    (tmp_path / "wiki" / "prospective_lock_manifest_bogus.json").write_text(
        json.dumps({"lock_date": "not-a-date", "surface_sha256": h}), encoding="utf-8")
    with pytest.raises(NoActiveLock) as e:
        resolve_active_lock(repo=tmp_path)
    assert "Malformed" in str(e.value), "the reason must be reported, not silently skipped"


def test_malformed_manifests_are_named_not_silently_skipped(tmp_path):
    """One sitting unnoticed in wiki/ is exactly how a bad cutoff enters later."""
    from dna_decode.eval.prospective_lock import is_well_formed_manifest
    problems = is_well_formed_manifest({"lock_date": "not-a-date", "surface_sha256": {}})
    assert any("schema" in p for p in problems)
    assert any("lock_date" in p for p in problems)


def test_a_non_hex_digest_is_rejected(tmp_path):
    from dna_decode.eval.prospective_lock import FROZEN_SURFACE_FILES, is_well_formed_manifest
    bad = {rel: "zz" * 32 for rel in FROZEN_SURFACE_FILES}
    assert any("64-hex" in p for p in is_well_formed_manifest(
        {"schema": "prospective-lock-manifest-v1", "lock_date": "2026-08-31", "surface_sha256": bad}))


def test_the_real_committed_manifests_are_BOTH_well_formed():
    """NON-VACUITY, and it protects the retired one too: the v1 manifest must fail on HASH DRIFT (it is
    retired), never on being malformed -- otherwise the refusal reason would be wrong and a future
    reader would think the v1 lock was junk rather than superseded."""
    from dna_decode.eval.prospective_lock import is_well_formed_manifest
    for p in sorted((ROOT / "wiki").glob("prospective_lock_manifest_*.json")):
        assert is_well_formed_manifest(json.loads(p.read_text(encoding="utf-8"))) == [], p.name


def test_filename_date_need_NOT_match_lock_date():
    """Pinned because it is a tempting invariant that is FALSE here: the v1 manifest is named
    ..._2026-06-22.json and legitimately carries lock_date 2026-06-13, because a manifest RECORDS a
    freeze that already happened rather than creating one. Enforcing agreement would reject the repo's
    own real manifest."""
    m = json.loads((ROOT / "wiki" / "prospective_lock_manifest_2026-06-22.json").read_text(
        encoding="utf-8"))
    assert m["lock_date"] == "2026-06-13"
