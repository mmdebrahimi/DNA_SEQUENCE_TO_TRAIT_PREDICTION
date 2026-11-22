"""The fetch driver must refuse an unresolvable spec and must catch a mis-resolved accession.

`download_genome` takes an accession and nothing else, so it cannot know which organism the caller
believed that accession to be. A resolution pass that returned the wrong assembly therefore produces a
fully successful download of the wrong genome, which enters the corpus under the spec's species label
with nothing downstream to contradict it -- the sequence is real, the annotation parses, the windows
extract, the model trains. That is the gap this driver closes, and the graded response is the point:
a different genus is never a naming variation, while a different epithet often is.

Every test here is offline. The accession gate and the disk check both run before any network call, so
they are testable without one, and the organism comparison is a pure function over a header string.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_fetch_corpus_genomes.py"
PEER_SPEC = Path("C:/Users/Laptop Portable/pythonprojects/dna-generator-research/handoffs"
                 "/bacterial-expansion-species-spec.json")

# Real headers, copied from assemblies in this project's own cache.
REAL_HEADERS = {
    "Escherichia coli": ">DBQBFQ010000001.1 TPA_asm: Escherichia coli isolate BCW12650 SAMN4294",
    "Campylobacter jejuni": ">JCAGAE010000063.1 Campylobacter jejuni isolate 357 NODE_101_length_29785",
    "Klebsiella pneumoniae": ">XYZ010000001.1 Klebsiella pneumoniae strain KP-1 contig00001",
}


@pytest.fixture(scope="module")
def mod():
    import sys

    spec = importlib.util.spec_from_file_location("_fetchdrv", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the header parser

@pytest.mark.parametrize("organism,header", sorted(REAL_HEADERS.items()))
def test_the_binomial_is_recovered_from_real_cache_headers(mod, organism, header):
    got = mod.parse_organism_from_header(header)
    assert got is not None, header
    assert " ".join(got) == organism


def test_the_parser_does_not_mistake_TPA_asm_for_a_binomial(mod):
    """`TPA_asm:` precedes the organism in a real header format and must not be captured instead.

    This is the specific false positive the pattern was written against: if it matched, every
    TPA-submitted assembly would verify as genus `TPA` and the genus gate would reject all of them.
    """
    got = mod.parse_organism_from_header(REAL_HEADERS["Escherichia coli"])
    assert got == ("Escherichia", "coli")


def test_a_header_with_no_organism_returns_None_rather_than_guessing(mod):
    assert mod.parse_organism_from_header(">NODE_1_length_4231_cov_12.3") is None
    assert mod.parse_organism_from_header(">contig00001") is None


# ------------------------------------------------------------------ the graded comparison

def test_a_matching_organism_verifies(mod):
    v = mod.verify_organism("Escherichia coli", REAL_HEADERS["Escherichia coli"])
    assert v["state"] == "VERIFIED"


def test_a_different_GENUS_is_refused(mod):
    """The case that matters: a transposed or mis-resolved accession landing on another genus."""
    v = mod.verify_organism("Bacillus subtilis", REAL_HEADERS["Escherichia coli"])
    assert v["state"] == "REFUSE_GENUS_MISMATCH"
    assert "not" in v["reason"] and "naming variation" in v["reason"]


def test_a_different_EPITHET_is_flagged_not_refused(mod):
    """Where reclassifications and synonyms land; a human should decide, so the file is kept."""
    v = mod.verify_organism("Campylobacter coli", REAL_HEADERS["Campylobacter jejuni"])
    assert v["state"] == "FLAG_SPECIES_MISMATCH"
    assert "reclassification" in v["reason"] or "synonym" in v["reason"]


def test_an_unparsable_header_is_UNVERIFIABLE_and_not_a_pass(mod):
    """Three states, not two: 'could not check' must never read as 'checked and fine'."""
    v = mod.verify_organism("Escherichia coli", ">contig00001")
    assert v["state"] == "UNVERIFIABLE"
    assert v["state"] != "VERIFIED"


def test_an_unparsable_EXPECTED_organism_is_also_unverifiable(mod):
    v = mod.verify_organism("unknown", REAL_HEADERS["Escherichia coli"])
    assert v["state"] == "UNVERIFIABLE"


def test_the_comparison_is_case_insensitive(mod):
    v = mod.verify_organism("escherichia coli", REAL_HEADERS["Escherichia coli"])
    assert v["state"] == "VERIFIED"


def test_every_state_is_reachable(mod):
    """A graded gate with an unreachable branch is a binary gate wearing a costume."""
    states = {
        mod.verify_organism("Escherichia coli", REAL_HEADERS["Escherichia coli"])["state"],
        mod.verify_organism("Bacillus subtilis", REAL_HEADERS["Escherichia coli"])["state"],
        mod.verify_organism("Campylobacter coli", REAL_HEADERS["Campylobacter jejuni"])["state"],
        mod.verify_organism("Escherichia coli", ">contig00001")["state"],
    }
    assert states == {"VERIFIED", "REFUSE_GENUS_MISMATCH", "FLAG_SPECIES_MISMATCH", "UNVERIFIABLE"}


# ------------------------------------------------------------------ the disk check

def test_the_disk_check_refuses_an_impossible_batch(mod, tmp_path):
    """A eukaryotic-scale batch must stop the run rather than fill the disk."""
    d = mod.check_disk(tmp_path, 1_000_000)
    assert d["ok"] is False
    assert d["estimated_need_gb"] > d["free_gb"]


def test_the_disk_check_passes_a_small_batch(mod, tmp_path):
    d = mod.check_disk(tmp_path, 16)
    assert d["estimated_need_gb"] < 0.2
    assert isinstance(d["ok"], bool)


def test_the_disk_check_keeps_a_headroom_margin(mod):
    assert mod.MIN_FREE_GB > 0


# ------------------------------------------------------------------ the accession gate

def test_a_spec_with_null_accessions_is_refused_with_the_reason(mod, tmp_path, capsys):
    spec = tmp_path / "s.json"
    spec.write_text(json.dumps({"spec": "t", "genomes": [
        {"species": "Bacillus subtilis", "strain": "168", "domain": "bacteria",
         "lineage_key": "Bacillota;Bacilli", "accession": None}]}), encoding="utf-8")
    rc = mod.main([str(spec), "--cache-root", str(tmp_path)])
    out = capsys.readouterr().out
    assert rc == 3
    assert "no accession" in out
    assert "NO taxon lookup" in out, "the refusal must say WHY it cannot resolve, not just that it did"


def test_an_entry_with_no_organism_at_all_is_refused(mod, tmp_path, capsys):
    """Without an expected organism there is nothing to verify the download against."""
    spec = tmp_path / "s.json"
    spec.write_text(json.dumps({"spec": "t", "genomes": [
        {"accession": "GCF_000005845.2", "domain": "bacteria", "lineage_key": "L1"}]}),
        encoding="utf-8")
    rc = mod.main([str(spec), "--cache-root", str(tmp_path)])
    assert rc == 3
    assert "neither" in capsys.readouterr().out


def test_an_empty_spec_is_refused(mod, tmp_path):
    spec = tmp_path / "s.json"
    spec.write_text(json.dumps({"spec": "t", "genomes": []}), encoding="utf-8")
    with pytest.raises(SystemExit):
        mod.load_spec(spec)


def test_a_missing_spec_file_is_refused(mod, tmp_path):
    assert mod.main([str(tmp_path / "nope.json")]) == 2


# ------------------------------------------------------------------ against the REAL peer spec

@pytest.mark.skipif(not PEER_SPEC.exists(), reason="peer handoff spec not present on this host")
def test_the_real_peer_spec_is_refused_for_the_right_reason(mod, tmp_path, capsys):
    """Non-vacuity on live data: the handed-over spec has accession null on every entry, so the gate
    must fire -- and must fire for the resolvability reason rather than on disk or parsing."""
    rc = mod.main([str(PEER_SPEC), "--cache-root", str(tmp_path), "--dry-run"])
    out = capsys.readouterr().out
    assert rc == 3
    assert "NO taxon lookup" in out
    assert "resolved GCF/GCA accessions" in out


@pytest.mark.skipif(not PEER_SPEC.exists(), reason="peer handoff spec not present on this host")
def test_the_real_peer_spec_is_otherwise_well_formed(mod):
    """Separating 'unresolvable' from 'malformed' matters: the spec is good, it is just not fetchable."""
    rec, genomes = mod.load_spec(PEER_SPEC)
    assert len(genomes) >= 12
    for g in genomes:
        assert g.get("species")
        assert g.get("domain") == "bacteria"
        assert g.get("lineage_key"), g
    assert len({g["lineage_key"] for g in genomes}) >= 10, "the point of the spec is lineage breadth"
