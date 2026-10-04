"""Tests for the identification reference-sketch builder (Step 2).

Offline and deterministic: every test builds a synthetic genome tree under tmp_path. Nothing here
needs Docker or the D: drive. The one real-surface fact the builder depends on -- that every genome
resolves under D: so a SINGLE mount makes one `mash sketch` call possible -- is enforced in code by
`to_container_path` and pinned by `test_a_path_outside_the_mounted_drive_RAISES`.
"""
from __future__ import annotations

import json

import pytest

from scripts.build_identify_reference import (
    MIN_ORGANISMS,
    ReferenceBuildError,
    SKETCH_PATH,
    build_manifest,
    enumerate_labelled_genomes,
    genus_species,
    main,
    organism_from_genbank,
    to_container_path,
)

GBK = """LOCUS       SYNTH               1000 bp    DNA     circular CON
DEFINITION  {name}, complete genome.
ACCESSION   SYNTH
VERSION     SYNTH.1
KEYWORDS    .
SOURCE      {name}
  ORGANISM  {name}
            Bacteria; Proteobacteria.
"""


def _genome(root, accession: str, organism: str | None, *, fasta: bool = True, gbk: bool = True):
    d = root / accession
    d.mkdir(parents=True, exist_ok=True)
    if fasta:
        (d / "genome.fna").write_text(">c1\nACGT\n", encoding="utf-8")
    if gbk and organism is not None:
        (d / "annotations.gbk").write_text(GBK.format(name=organism), encoding="utf-8")
    return d


# --------------------------------------------------------------------------- label parsing

def test_organism_is_read_from_the_genbank_header(tmp_path):
    d = _genome(tmp_path, "GCA_1.1", "Escherichia coli 'BL21-Gold(DE3)pLysS AG'")
    assert organism_from_genbank(d / "annotations.gbk") == "Escherichia coli 'BL21-Gold(DE3)pLysS AG'"


def test_a_missing_or_unreadable_gbk_returns_None_rather_than_raising(tmp_path):
    """One bad genome must not abort a 2,000-genome enumeration -- but the caller counts the Nones,
    because a silently-shrinking corpus is this script's failure mode."""
    assert organism_from_genbank(tmp_path / "nope.gbk") is None
    (tmp_path / "empty.gbk").write_text("", encoding="utf-8")
    assert organism_from_genbank(tmp_path / "empty.gbk") is None


def test_an_ORGANISM_line_past_the_header_window_is_not_read(tmp_path):
    """Bounded read: the field belongs in the header, and scanning a whole 40 MB GenBank file per
    genome would dominate the enumeration."""
    f = tmp_path / "deep.gbk"
    f.write_text("x\n" * 200 + "  ORGANISM  Escherichia coli\n", encoding="utf-8")
    assert organism_from_genbank(f, max_header_lines=60) is None


@pytest.mark.parametrize("name,expect", [
    ("Escherichia coli 'BL21-Gold(DE3)pLysS AG'", "Escherichia_coli"),
    ("Klebsiella pneumoniae subsp. pneumoniae HS11286", "Klebsiella_pneumoniae"),
    ("Campylobacter jejuni", "Campylobacter_jejuni"),
    ("Gemmata", "Gemmata"),
    ("   ", ""),
])
def test_genus_species_drops_the_strain_suffix(name, expect):
    assert genus_species(name) == expect


# --------------------------------------------------------------------------- the mount guard

def test_a_path_outside_the_mounted_drive_RAISES(tmp_path):
    """FAIL-CLOSED, and this is the test that matters for sketch integrity.

    A genome outside the mounted drive simply does not exist inside the container, so `mash` skips it
    and the sketch silently covers fewer genomes than the manifest claims. Raising before any
    container starts is the only way that stays visible.
    """
    with pytest.raises(ReferenceBuildError):
        to_container_path(tmp_path / "genome.fna", drive=type(tmp_path)("Z:/nonexistent-root"))


def test_a_path_under_the_mounted_drive_is_rewritten(tmp_path):
    got = to_container_path(tmp_path / "sub" / "genome.fna", drive=tmp_path, mount="/d")
    assert got == "/d/sub/genome.fna", got


def test_the_input_list_is_written_with_LF_endings(tmp_path, monkeypatch):
    """A REAL failure, not a hypothetical. On Windows `write_text` turns "\\n" into "\\r\\n", so every
    line of the mash file-of-filenames ended with a stray CR; mash then tried to open
    `genome.fna\\r` and reported `could not open ... for reading` for a file that was demonstrably
    present and readable inside the container. The CR additionally returned the terminal cursor and
    mangled mash's own error text, which sent the diagnosis toward the mount.

    Asserts on the BYTES, because the bug is invisible in any text-mode read.
    """
    import scripts.build_identify_reference as mod

    sketch_dir = tmp_path / "identify"
    monkeypatch.setattr(mod, "SKETCH_DIR", sketch_dir)
    monkeypatch.setattr(mod, "HOST_DRIVE", tmp_path)

    g1 = _genome(tmp_path / "refseq", "GCA_1.1", "Escherichia coli") / "genome.fna"
    g2 = _genome(tmp_path / "refseq", "GCA_2.1", "Campylobacter jejuni") / "genome.fna"

    captured = {}

    def fake_run(image, args, **kw):
        captured["args"] = args
        raw = (sketch_dir / "reference_inputs.txt").read_bytes()
        assert b"\r" not in raw, f"CRLF leaked into the mash input list: {raw[:80]!r}"
        (sketch_dir / "reference.msh").write_bytes(b"stub")
        class R:  # minimal CompletedProcess stand-in
            stdout = stderr = ""
        return R()

    import tools.docker_runner as dr
    monkeypatch.setattr(dr, "run", fake_run)
    mod._mash_sketch([g1, g2], sketch_dir / "reference.msh")

    raw = (sketch_dir / "reference_inputs.txt").read_bytes()
    assert b"\r" not in raw
    assert raw.decode().splitlines() == ["/d/refseq/GCA_1.1/genome.fna",
                                         "/d/refseq/GCA_2.1/genome.fna"]
    assert "-l" in captured["args"] and "-k" in captured["args"]


def test_mash_reporting_success_without_producing_a_sketch_RAISES(tmp_path, monkeypatch):
    """exit 0 is not evidence. This run's own harness reported the sketch build as "exit code 0"
    while the script had exited 1, so the file-existence check is the thing that actually decides."""
    import scripts.build_identify_reference as mod
    import tools.docker_runner as dr

    sketch_dir = tmp_path / "identify"
    monkeypatch.setattr(mod, "SKETCH_DIR", sketch_dir)
    monkeypatch.setattr(mod, "HOST_DRIVE", tmp_path)
    g = _genome(tmp_path / "refseq", "GCA_1.1", "Escherichia coli") / "genome.fna"

    class R:
        stdout = stderr = "looked fine"
    monkeypatch.setattr(dr, "run", lambda *a, **k: R())   # produces NO .msh

    with pytest.raises(ReferenceBuildError, match="does not exist"):
        mod._mash_sketch([g], sketch_dir / "reference.msh")


def test_the_sketch_lives_on_D_and_never_inside_the_repo():
    """Standing host rule: heavy artifacts go on D:, never C:. The sketch is regenerable and is
    deliberately not committed -- only its manifest is."""
    s = SKETCH_PATH.as_posix().lower()
    assert s.startswith("d:/"), SKETCH_PATH
    assert "pythonprojects" not in s, SKETCH_PATH


# --------------------------------------------------------------------------- enumeration

def test_enumeration_labels_dedups_and_counts_out_of_set(tmp_path):
    r1, r2 = tmp_path / "r1", tmp_path / "r2"
    _genome(r1, "GCA_1.1", "Escherichia coli K-12")
    _genome(r1, "GCA_2.1", "Campylobacter jejuni")          # genus-level -> in-set
    _genome(r1, "GCA_3.1", "Klebsiella aerogenes")          # congener -> OUT of set
    _genome(r1, "GCA_4.1", "Escherichia coli", fasta=False)  # no fasta -> skipped
    _genome(r1, "GCA_5.1", None, gbk=False)                 # no organism -> counted
    _genome(r2, "GCA_1.1", "Escherichia coli K-12")          # SAME accession in another root

    per, stats = enumerate_labelled_genomes([r1, r2])

    assert sorted(per) == ["campylobacter", "escherichia_coli"]
    assert len(per["escherichia_coli"]) == 1, "the same accession must not enter the reference twice"
    assert stats["n_duplicate_accession"] == 1
    assert stats["n_no_fasta"] == 1
    assert stats["n_no_organism"] == 1
    assert stats["out_of_set_counts"] == {"Klebsiella_aerogenes": 1}
    assert stats["n_out_of_set"] == 1


# --------------------------------------------------------------------------- manifest

def test_manifest_reports_imbalance_and_unscorable_organisms(tmp_path):
    per = {
        "escherichia_coli": [(f"E{i}", tmp_path / f"E{i}") for i in range(7)],
        "campylobacter": [("C1", tmp_path / "C1"), ("C2", tmp_path / "C2")],
        "klebsiella_oxytoca": [("K1", tmp_path / "K1")],      # n=1
    }
    m = build_manifest(per, {"n_out_of_set": 3, "out_of_set_counts": {"Klebsiella_aerogenes": 3}})

    assert m["n_genomes"] == 10 and m["n_organisms"] == 3
    assert m["largest_organism_share"] == pytest.approx(0.7)
    assert m["unscorable_under_leave_one_out"] == ["klebsiella_oxytoca"], (
        "an organism with one genome cannot be leave-one-out scored -- removing it leaves no "
        "reference for that organism, so it must be named, not silently averaged in")
    assert m["label_is_wet_lab"] is False
    assert m["sketch_committed"] is False
    assert "ORGANISM" in m["label_source"]
    assert m["out_of_set"]["n_genomes"] == 3
    assert json.loads(json.dumps(m)) == m, "manifest must be JSON-round-trippable"


def test_per_organism_is_ordered_largest_first(tmp_path):
    per = {"a": [("1", tmp_path)], "b": [("2", tmp_path), ("3", tmp_path)]}
    assert list(build_manifest(per, {})["per_organism"]) == ["b", "a"]


# --------------------------------------------------------------------------- the refusal

def test_main_REFUSES_a_reference_covering_fewer_than_two_organisms(tmp_path, monkeypatch):
    """A one-organism reference cannot discriminate between organisms, which is the only thing it
    exists to do. Refuse rather than emit something that answers every query with its single label."""
    import scripts.build_identify_reference as mod

    root = tmp_path / "refseq"
    _genome(root, "GCA_1.1", "Escherichia coli")
    _genome(root, "GCA_2.1", "Escherichia coli")
    monkeypatch.setattr(mod, "refseq_roots", lambda *a, **k: [root])

    rc = main(["--dry-run", "--manifest-out", str(tmp_path / "m.json")])
    assert rc == 3, f"expected the MIN_ORGANISMS refusal (exit 3), got {rc}"
    assert not (tmp_path / "m.json").exists(), "a refused build must not leave a manifest behind"
    assert MIN_ORGANISMS == 2


def test_main_writes_a_manifest_and_skips_mash_on_dry_run(tmp_path, monkeypatch):
    import scripts.build_identify_reference as mod

    root = tmp_path / "refseq"
    _genome(root, "GCA_1.1", "Escherichia coli")
    _genome(root, "GCA_2.1", "Campylobacter jejuni")
    monkeypatch.setattr(mod, "refseq_roots", lambda *a, **k: [root])

    def boom(*a, **k):  # pragma: no cover - must not run
        raise AssertionError("--dry-run must not invoke mash")
    monkeypatch.setattr(mod, "_mash_sketch", boom)

    out = tmp_path / "m.json"
    assert main(["--dry-run", "--manifest-out", str(out)]) == 0
    m = json.loads(out.read_text(encoding="utf-8"))
    assert m["n_organisms"] == 2
    assert m["schema"] == "identify-reference-manifest-v1"
