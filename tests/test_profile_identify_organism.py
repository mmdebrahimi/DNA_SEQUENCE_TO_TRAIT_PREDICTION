"""`profile --identify-organism`: replace the E. coli ASSUMPTION with a MEASUREMENT.

Offline and deterministic -- the identify query is stubbed. What is pinned here is the PROVENANCE
contract, because the whole point is that a reader can tell these four states apart:

    explicit          the user said so (and that always wins over inference)
    identified        dna-identify named a supported organism -> organism_assumed becomes FALSE
    assumed_identify_abstained    a correct biological answer: not a supported organism
    assumed_identify_unavailable  an infrastructure FAULT

Collapsing the last two is the failure this guards: a wedged Docker mount must never read as
"this genome is not a supported organism".

WHAT THIS IS NOT: a defect fix. `profile` already DISCLOSED the assumption (`organism_assumed=true`
plus a printed note). This upgrades a disclosed assumption to a measurement where one is available.
"""
from __future__ import annotations

import pytest

import dna_decode.profile.cli as pcli
from dna_decode.identify.core import AbstainReason, Hit, IdentifyCall
from dna_decode.identify.runner import ReferenceUnavailable


def _call(organism="klebsiella_pneumoniae", token="Klebsiella_pneumoniae", d=0.001):
    return IdentifyCall(organism, token, False, nearest_distance=d)


def _abstain(reason=AbstainReason.ABOVE_MAX_DISTANCE, d=0.2):
    return IdentifyCall(None, None, True, reason, nearest_distance=d)


@pytest.fixture
def stub(monkeypatch, tmp_path):
    """Stub the identify seam at the two points the resolver imports it from."""
    man_dir = tmp_path / "wiki"
    man_dir.mkdir()
    (man_dir / "identify_reference_manifest_2026-01-01.json").write_text(
        '{"accessions_by_organism": {"klebsiella_pneumoniae": ["GCA_1.1"]}}', encoding="utf-8")

    def _install(call=None, exc=None):
        import dna_decode.identify.core as core
        import dna_decode.identify.runner as runner
        monkeypatch.setattr(pcli, "__file__", str(tmp_path / "dna_decode" / "profile" / "cli.py"))
        (tmp_path / "dna_decode" / "profile").mkdir(parents=True, exist_ok=True)
        if exc is not None:
            monkeypatch.setattr(runner, "query", lambda *a, **k: (_ for _ in ()).throw(exc))
        else:
            monkeypatch.setattr(runner, "query",
                                lambda *a, **k: type("R", (), {"hits": (Hit("/x/GCA_1.1/g.fna", "q",
                                                                           0.001, 0.0, "1/1"),),
                                                               "stderr": ""})())
            monkeypatch.setattr(core, "decide", lambda *a, **k: call)
    return _install


# --------------------------------------------------------------------------- user authority

def test_an_explicit_organism_ALWAYS_WINS_over_identification(stub):
    """User authority outranks inference, and an explicit value is also the escape hatch when the
    router abstains. Identification must not even be attempted."""
    import dna_decode.identify.runner as runner
    def boom(*a, **k):
        raise AssertionError("identification ran despite an explicit --amr-organism")
    runner.query = boom
    r = pcli._resolve_amr_organism("g.fna", "Salmonella", True)
    assert r["organism_used"] == "Salmonella"
    assert r["organism_assumed"] is False
    assert r["organism_source"] == "explicit"


def test_identification_is_OPT_IN(stub):
    """`profile`'s auto-run path (used by `dna-decode decode`) is offline-safe. Identification needs
    Docker, so making it automatic would break that guarantee."""
    import dna_decode.identify.runner as runner
    def boom(*a, **k):
        raise AssertionError("identification ran without --identify-organism")
    runner.query = boom
    r = pcli._resolve_amr_organism("g.fna", None, False)
    assert r["organism_used"] == "Escherichia"
    assert r["organism_assumed"] is True
    assert r["organism_source"] == "assumed_default"


# --------------------------------------------------------------------------- the measurement

def test_a_successful_identification_makes_organism_assumed_FALSE(stub):
    """The field name means what it says: once measured it is no longer an assumption."""
    stub(call=_call())
    r = pcli._resolve_amr_organism("g.fna", None, True)
    assert r["organism_used"] == "Klebsiella_pneumoniae"
    assert r["organism_assumed"] is False
    assert r["organism_source"] == "identified"
    assert r["organism_identify"]["organism"] == "klebsiella_pneumoniae"


def test_a_supported_organism_with_NO_amrfinder_value_does_not_fake_one(stub):
    """C. auris / M. tuberculosis are supported but have no AMRFinder `-O` value -- they route to
    their own engines. The AMR section must not pretend otherwise, so the assumption STANDS and says
    why."""
    stub(call=IdentifyCall("candida_auris", None, False, nearest_distance=0.001))
    r = pcli._resolve_amr_organism("g.fna", None, True)
    assert r["organism_used"] == "Escherichia"
    assert r["organism_assumed"] is True
    assert r["organism_source"] == "assumed_identified_but_not_amrfinder_routable"
    assert r["organism_identify"]["routing_token"] is None


# --------------------------------------------------------------------------- the two failure classes

def test_an_ABSTENTION_keeps_the_disclosed_assumption_and_records_why(stub):
    """An abstention is a CORRECT answer ('not a supported organism'), so the documented E. coli
    assumption stands -- but it must not look identical to never having tried."""
    stub(call=_abstain())
    r = pcli._resolve_amr_organism("g.fna", None, True)
    assert r["organism_used"] == "Escherichia"
    assert r["organism_assumed"] is True
    assert r["organism_source"] == "assumed_identify_abstained"
    assert r["organism_identify"]["reason"] == "above_max_distance"


def test_a_FAULT_is_recorded_SEPARATELY_from_an_abstention(stub):
    """THE DISTINCTION THIS FILE EXISTS FOR. A wedged Docker D: mount is a KNOWN failure on this host;
    if it collapsed into the abstention bucket it would read as a biological finding."""
    stub(exc=ReferenceUnavailable("sketch missing"))
    r = pcli._resolve_amr_organism("g.fna", None, True)
    assert r["organism_source"] == "assumed_identify_unavailable"
    assert r["organism_assumed"] is True
    assert r["organism_identify"]["status"] == "reference_unavailable"
    assert r["organism_source"] != "assumed_identify_abstained"


def test_an_unexpected_exception_is_also_a_FAULT_not_an_abstention(stub):
    """Fail-safe: any other error still reports as a fault, never as 'organism unsupported', and never
    propagates out to abort a whole profile run."""
    stub(exc=RuntimeError("docker daemon gone"))
    r = pcli._resolve_amr_organism("g.fna", None, True)
    assert r["organism_source"] == "assumed_identify_failed"
    assert r["organism_assumed"] is True
    assert "RuntimeError" in r["organism_identify"]["detail"]


# --------------------------------------------------------------------------- surface

def test_the_flag_is_ACCEPTED_by_the_real_parser_and_a_bogus_one_is_not(tmp_path, capsys):
    """Advertised flags are promises, so resolve this one through the REAL parser.

    Discriminating, not just "it ran": `--identify-organism` must parse and reach the
    file-not-found path (rc 2), while an unknown flag must be REJECTED by argparse (SystemExit).
    Without the second half the test would pass for a flag argparse silently ignored.

    (`main(["--help"])` returns None here rather than raising, so help-text scraping is not a valid
    probe for this CLI -- measured, not assumed.)
    """
    missing = tmp_path / "nope.fna"
    assert pcli.main([str(missing), "--identify-organism"]) == 2
    assert "not found" in capsys.readouterr().err.lower()

    with pytest.raises(SystemExit):
        pcli.main([str(missing), "--definitely-not-a-flag"])


def test_the_flag_is_documented_in_its_own_help_string():
    """The help text is the user-facing promise; it must name the tool and the fallback behaviour."""
    import inspect
    src = inspect.getsource(pcli.main)
    assert "--identify-organism" in src
    assert "dna-identify" in src
    assert "wins" in src, "the help must say --amr-organism still wins over identification"


def test_every_source_value_is_distinguishable():
    """Non-vacuity: the four states must be four VALUES, not one flag plus prose."""
    seen = set()
    import dna_decode.identify.runner as runner
    runner.query = lambda *a, **k: (_ for _ in ()).throw(ReferenceUnavailable("x"))
    seen.add(pcli._resolve_amr_organism("g.fna", "Salmonella", True)["organism_source"])
    seen.add(pcli._resolve_amr_organism("g.fna", None, False)["organism_source"])
    seen.add(pcli._resolve_amr_organism("g.fna", None, True)["organism_source"])
    assert len(seen) == 3, seen
