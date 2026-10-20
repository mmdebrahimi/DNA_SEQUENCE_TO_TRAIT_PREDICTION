"""Tests for the external-root probe.

All offline: a "drive" is just a tmp_path directory, so the CRLF-free, letter-free
cases are testable on a host that has only one real disk.

THE LOAD-BEARING TEST is `test_an_unfit_candidate_is_never_substituted`. The whole
point of this module is that a directory of the right NAME is not evidence of the
right CONTENTS -- on the host that motivated it, `E:` carried a `dna_decode_cache`
holding only `amr_portal` and `mpra`, with no `refseq`. Substituting it would have
turned a loud missing-ROOT into a quiet missing-DATASET across 389 call sites.
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.data.external_root import (
    CACHE_REQUIRED_SUBDIRS,
    SIGNATURE_DIRS,
    Candidate,
    candidate_drive_roots,
    describe,
    probe,
    resolve_external_root,
)


def _make(root: Path, rel: str, *entries: str) -> Path:
    d = root / rel
    d.mkdir(parents=True, exist_ok=True)
    for e in entries:
        (d / e).mkdir(exist_ok=True)
    return d


def test_probe_finds_a_signature_dir_at_a_drive_root(tmp_path):
    _make(tmp_path, "dna_decode_cache", "refseq")
    found = probe(drives=[tmp_path])
    assert [c.signature for c in found] == ["dna_decode_cache"]
    assert found[0].entries == ("refseq",)


def test_probe_finds_a_nested_layout(tmp_path):
    """The MEASURED layout -- a root-only scan found nothing on the real host."""
    _make(tmp_path, "PythonProjects/dna_decode_cache", "amr_portal", "mpra")
    found = probe(drives=[tmp_path])
    assert len(found) == 1
    assert found[0].path.name == "dna_decode_cache"
    assert found[0].entries == ("amr_portal", "mpra")


def test_probe_reports_an_empty_candidate_rather_than_skipping_it(tmp_path):
    _make(tmp_path, "dna_decode_cache")
    found = probe(drives=[tmp_path])
    assert len(found) == 1 and found[0].entries == ()


def test_probe_does_not_double_report_the_same_directory(tmp_path):
    _make(tmp_path, "dna_decode_cache", "refseq")
    found = probe(drives=[tmp_path, tmp_path])
    assert len(found) == 1


def test_probe_finds_nothing_when_there_is_nothing(tmp_path):
    assert probe(drives=[tmp_path]) == []


def test_both_signature_names_are_recognised(tmp_path):
    for sig in SIGNATURE_DIRS:
        _make(tmp_path, sig)
    found = {c.signature for c in probe(drives=[tmp_path])}
    assert found == set(SIGNATURE_DIRS)


def test_candidate_has_and_missing_are_case_insensitive():
    c = Candidate(path=Path("x"), signature="dna_decode_cache", entries=("RefSeq", "mpra"))
    assert c.has(("refseq",))
    assert c.missing(("refseq",)) == ()
    assert c.missing(("refseq", "embeddings")) == ("embeddings",)


def test_configured_root_is_always_preferred_when_it_exists(tmp_path):
    conf = _make(tmp_path, "configured_cache", "refseq")
    _make(tmp_path, "dna_decode_cache", "refseq")
    path, status = resolve_external_root(conf, allow_discovered=True,
                                         required=CACHE_REQUIRED_SUBDIRS)
    assert status == "configured" and path == conf


def test_absent_configured_with_nothing_found_returns_the_configured_path(tmp_path):
    path, status = resolve_external_root(tmp_path / "nope", candidates=[])
    assert status == "absent" and path == tmp_path / "nope"


def test_default_does_not_substitute_even_a_FIT_candidate(tmp_path):
    """Substitution is opt-in, so wiring this in changes no existing behaviour."""
    fit = Candidate(path=tmp_path / "found", signature="dna_decode_cache", entries=("refseq",))
    path, status = resolve_external_root(tmp_path / "nope", candidates=[fit],
                                         required=CACHE_REQUIRED_SUBDIRS)
    assert path == tmp_path / "nope", "the default must never silently repoint the project"
    assert status != "discovered"


def test_opt_in_substitutes_a_fit_candidate(tmp_path):
    fit = Candidate(path=tmp_path / "found", signature="dna_decode_cache", entries=("refseq",))
    path, status = resolve_external_root(tmp_path / "nope", candidates=[fit],
                                         required=CACHE_REQUIRED_SUBDIRS,
                                         allow_discovered=True)
    assert status == "discovered" and path == tmp_path / "found"


def test_an_unfit_candidate_is_never_substituted(tmp_path):
    """THE load-bearing case: the real E: drive, name-matching but refseq-less."""
    unfit = Candidate(path=tmp_path / "E_PythonProjects_cache", signature="dna_decode_cache",
                      entries=("amr_portal", "mpra"))
    path, status = resolve_external_root(tmp_path / "nope", candidates=[unfit],
                                         required=CACHE_REQUIRED_SUBDIRS,
                                         allow_discovered=True)
    assert status == "absent-candidates-unfit"
    assert path == tmp_path / "nope", "an unfit root must not be adopted even on opt-in"


def test_describe_names_the_missing_requirement(tmp_path):
    """HERMETIC as of 2026-10-20. It was not, and that is why it broke.

    This used to call `describe` with no candidates, so it scanned the REAL host and its
    UNFIT/refseq assertion only held while some attached volume carried a refseq-less
    `dna_decode_cache`. It passed for as long as `E:` did, then failed the hour that volume
    detached -- a unit test whose verdict tracked the operator's hardware. The candidate is
    now injected, exactly as the sibling `resolve_external_root` tests already do.
    """
    unfit = Candidate(path=tmp_path / "E_PythonProjects_cache", signature="dna_decode_cache",
                      entries=("amr_portal", "mpra"))
    text = describe(tmp_path / "nope", required=CACHE_REQUIRED_SUBDIRS, candidates=[unfit])
    assert "ABSENT" in text
    assert "UNFIT" in text and "refseq" in text
    assert "do not point the project here" in text


def test_describe_with_injected_candidates_does_NOT_touch_the_real_host(monkeypatch):
    """Non-vacuity for the seam: if it still called probe(), this would raise.

    Without this, the fix above could pass while `describe` quietly ignored the argument.
    """
    import dna_decode.data.external_root as mod

    def _boom(*_a, **_k):
        raise AssertionError("describe() probed the real host despite injected candidates")

    monkeypatch.setattr(mod, "probe", _boom)
    fit = Candidate(path=Path("X:/cache"), signature="dna_decode_cache", entries=("refseq",))
    text = mod.describe("Q:/nope", required=CACHE_REQUIRED_SUBDIRS, candidates=[fit])
    assert "fit" in text


def test_describe_with_NO_candidates_still_reports_the_live_host(monkeypatch):
    """The live path must be unchanged -- the seam is additive, not a behaviour switch."""
    import dna_decode.data.external_root as mod

    monkeypatch.setattr(mod, "probe", lambda *_a, **_k: [])
    assert "none on any resolving drive" in mod.describe("Q:/nope",
                                                         required=CACHE_REQUIRED_SUBDIRS)


def test_describe_reports_a_clean_configured_root(tmp_path):
    conf = _make(tmp_path, "dna_decode_cache", "refseq")
    assert "EXISTS" in describe(conf, required=CACHE_REQUIRED_SUBDIRS)


def test_candidate_drive_roots_sees_at_least_one_real_drive():
    """Real-surface smoke: a drive with no media must not appear, but C: must."""
    roots = {str(p).upper()[0] for p in candidate_drive_roots()}
    assert roots, "no drive resolved at all, which cannot be true of a running host"
