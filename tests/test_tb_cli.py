"""Tests for the `dna-tb` CLI surface + its registration (2026-10-04).

Offline. The WHO catalogue master is gitignored, so every test that needs a real determinant match
SKIPS when it is absent rather than asserting on a fixture that pretends to be it -- a catalogue-shaped
fake would make the refusal path untestable, and the refusal path is the point of this cell.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_decode.organism_rules import tb_cli

REPO = Path(__file__).resolve().parent.parent
CAT = REPO / "data" / "raw" / "who_tb_catalogue"
VCF_DIR = REPO / "data" / "raw" / "tb_indep" / "vcf"


def _catalogue_ready() -> bool:
    try:
        from dna_decode.data import tb_who_catalogue as cat
        return cat.catalogue_available(CAT) and all(cat.verify_pins(CAT).values())
    except Exception:
        return False


def _a_real_vcf() -> Path | None:
    got = sorted(VCF_DIR.glob("*.vcf")) if VCF_DIR.is_dir() else []
    return got[0] if got else None


# --------------------------------------------------------------- argparse surface

def test_the_advertised_examples_all_parse():
    """Every command the module docstring tells a user to run must resolve through the REAL parser.
    `tests/test_advertised_commands.py` enforces this repo-wide; this is the local pin."""
    p = tb_cli.build_parser()
    p.parse_args(["--vcf", "x.vcf", "--drug", "rifampicin"])
    p.parse_args(["--vcf", "x.vcf", "--drug", "isoniazid", "--json"])
    p.parse_args(["--list-drugs"])


def test_an_unknown_drug_is_refused_not_guessed(tmp_path, capsys):
    v = tmp_path / "g.vcf"; v.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    assert tb_cli.main(["--vcf", str(v), "--drug", "not-a-drug"]) == 2
    assert "not in the WHO catalogue" in capsys.readouterr().err


def test_a_missing_vcf_is_refused(tmp_path, capsys):
    assert tb_cli.main(["--vcf", str(tmp_path / "nope.vcf"), "--drug", "rifampicin"]) == 2
    assert "VCF not found" in capsys.readouterr().err


# --------------------------------------------------------------- the FILTER-vocabulary defect

def test_filter_dot_is_normalized_to_PASS():
    """THE DEFECT THIS CELL SHIPPED WITH FOR ONE RUN, found on real data.

    `parse_masked_calls` requires FILTER==PASS (right for CRyPTIC minos VCFs). The minimap2/paftools.js
    VCFs the AMR-Portal arm actually produced write FILTER `.`, so without this the caller parsed ZERO
    calls from a VCF holding 1,382 genuine GT=1/1 rows and returned a confident S for every drug.
    """
    raw = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
           "NC_000962.3\t761155\t.\tC\tT\t60\t.\tQNAME=x\tGT\t1/1\n")
    out, rewrote = tb_cli.normalize_filter_column(raw)
    assert rewrote is True
    assert "\tPASS\t" in out
    from dna_decode.organism_rules import tb_vcf
    assert tb_vcf.parse_masked_calls(raw) == {}, "precondition: the raw form parses to nothing"
    assert 761155 in tb_vcf.parse_masked_calls(out)


def test_a_real_PASS_vocabulary_is_left_alone():
    """A caller that DOES filter marks failures with a NAME, never `.`, so only `.` may be promoted --
    otherwise this would smuggle a genuinely-failing record through as a call."""
    raw = ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
           "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1/1\n"
           "NC_000962.3\t761156\t.\tA\tG\t60\tMIN_DP\t.\tGT\t1/1\n")
    out, rewrote = tb_cli.normalize_filter_column(raw)
    assert rewrote is False
    assert out.splitlines() == raw.splitlines()
    assert "MIN_DP" in out, "a FAILING record must stay failed"


def test_a_vcf_with_rows_that_parses_to_nothing_REFUSES(tmp_path, capsys):
    """An unparseable VCF and a clean genome would both print `PREDICTION: S`. That is the
    indistinguishable-failure shape, and it is exactly how the FILTER bug presented."""
    if not _catalogue_ready():
        pytest.skip("gitignored WHO catalogue absent")
    v = tmp_path / "weird.vcf"
    # data rows present, but GT is reference -- nothing is a call, so this is a PARSE/format mismatch
    v.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                 "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t0/0\n", encoding="utf-8")
    rc = tb_cli.main(["--vcf", str(v), "--drug", "rifampicin", "--catalogue-dir", str(CAT)])
    assert rc == 4
    err = capsys.readouterr().err
    assert "PARSE FAILURE" in err and "not a susceptible genome" in err


def test_an_empty_vcf_does_NOT_trip_the_parse_guard(tmp_path):
    """Zero data rows is a legitimately empty VCF, not a parse failure -- the guard keys on
    rows-present-but-nothing-parsed, so it must not fire here."""
    if not _catalogue_ready():
        pytest.skip("gitignored WHO catalogue absent")
    v = tmp_path / "empty.vcf"
    v.write_text("##fileformat=VCFv4.2\n"
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n", encoding="utf-8")
    assert tb_cli.count_data_rows(v.read_text(encoding="utf-8")) == 0
    assert tb_cli.main(["--vcf", str(v), "--drug", "rifampicin",
                        "--catalogue-dir", str(CAT)]) == 0


# --------------------------------------------------------------- catalogue refusal

def test_an_absent_catalogue_REFUSES_and_never_returns_a_call(tmp_path, capsys):
    """'no determinant found' and 'no catalogue loaded' must not both print S."""
    v = tmp_path / "g.vcf"
    v.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                 "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1/1\n", encoding="utf-8")
    rc = tb_cli.main(["--vcf", str(v), "--drug", "rifampicin",
                      "--catalogue-dir", str(tmp_path / "no-catalogue")])
    cap = capsys.readouterr()
    assert rc == 3
    assert "CATALOGUE UNUSABLE" in cap.err
    assert "PREDICTION" not in cap.out, "a refusal must not emit a call"
    assert "S" not in cap.out


def test_a_TAMPERED_catalogue_refuses_too_not_just_an_absent_one(tmp_path, capsys):
    """An absent catalogue and a DRIFTED one are different failures and both must refuse.

    Verified against a real tamper: copying the catalogue and appending ONE byte to the pinned
    coordinates file is caught, named by filename with both hashes, exit 3. A pin that is only checked
    for existence would pass a catalogue whose contents no longer match what the published numbers were
    measured against.
    """
    if not _catalogue_ready():
        pytest.skip("gitignored WHO catalogue absent")
    import shutil
    dst = tmp_path / "cat"
    shutil.copytree(CAT, dst)
    victim = next(p for p in sorted(dst.iterdir())
                  if p.is_file() and "checksum" not in p.name.lower())
    with victim.open("ab") as fh:
        fh.write(b"x")
    v = tmp_path / "g.vcf"
    v.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                 "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1/1\n", encoding="utf-8")
    rc = tb_cli.main(["--vcf", str(v), "--drug", "rifampicin", "--catalogue-dir", str(dst)])
    cap = capsys.readouterr()
    assert rc == 3
    assert "CATALOGUE UNUSABLE" in cap.err
    assert victim.name in cap.err, "the refusal must name WHICH file drifted"
    assert "PREDICTION" not in cap.out


def test_the_refusal_fires_BEFORE_any_scoring(tmp_path, monkeypatch):
    """Order matters: scoring first and refusing after would still have loaded a determinant set."""
    import dna_decode.organism_rules.tb_amr as tb_amr

    def boom(*a, **k):
        raise AssertionError("score_drug must not be reached when the catalogue is unusable")
    monkeypatch.setattr(tb_amr, "score_drug", boom)
    v = tmp_path / "g.vcf"; v.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    assert tb_cli.main(["--vcf", str(v), "--drug", "rifampicin",
                        "--catalogue-dir", str(tmp_path / "gone")]) == 3


# --------------------------------------------------------------- real data

def test_a_real_resistant_isolate_calls_R_with_the_canonical_determinant(capsys):
    """END-TO-END on a real VCF. Pinned on rpoB S450L specifically -- the single most common
    rifampicin determinant worldwide and the one the coordinate-alignment probe verified at
    761155 C>T -- so this fails loudly on a frame/coordinate error rather than drifting."""
    if not _catalogue_ready():
        pytest.skip("gitignored WHO catalogue absent")
    vcf = VCF_DIR / "SAMEA1015921.vcf"
    if not vcf.exists():
        pytest.skip("real TB VCF not on this host")
    rc = tb_cli.main(["--vcf", str(vcf), "--drug", "rifampicin",
                      "--catalogue-dir", str(CAT), "--json"])
    assert rc == 0
    rec = json.loads(capsys.readouterr().out)
    assert rec["prediction"] == "R"
    assert any("450" in d for d in rec["matched_determinants"]), rec["matched_determinants"]
    assert rec["n_nonref_calls"] > 0
    assert rec["filter_column_normalized"] is True
    assert "no_caller_quality_floor" in rec, "the missing quality floor must be DISCLOSED"


def test_the_determinant_name_is_not_doubled(capsys):
    """`variant` already carries the gene prefix, so a naive f'{gene}_{variant}' gives
    `rpoB_rpoB_p.Ser450Leu`."""
    if not _catalogue_ready():
        pytest.skip("gitignored WHO catalogue absent")
    vcf = _a_real_vcf()
    if vcf is None:
        pytest.skip("no real TB VCF on this host")
    tb_cli.main(["--vcf", str(vcf), "--drug", "isoniazid", "--catalogue-dir", str(CAT), "--json"])
    for d in json.loads(capsys.readouterr().out)["matched_determinants"]:
        gene = d.split("_")[0]
        assert not d.startswith(f"{gene}_{gene}_"), f"doubled gene name: {d}"


# --------------------------------------------------------------- the two-tier honesty rail

def test_only_rif_and_inh_claim_an_independent_number():
    assert set(tb_cli.INDEPENDENT_DRUGS) == {"rifampicin", "isoniazid"}


def test_the_independent_numbers_match_the_artifact():
    """Written FROM the artifact, never from memory. If the artifact is re-scored this must be
    re-read, not adjusted to agree."""
    art = REPO / "wiki" / "tb_independent_amr_portal_lineage_collapsed.json"
    if not art.exists():
        pytest.skip("artifact absent")
    d = json.loads(art.read_text(encoding="utf-8"))["drugs"]
    for drug, ours in tb_cli.INDEPENDENT_DRUGS.items():
        assert ours["lineage_sens"] == d[drug]["lineage_collapsed"]["sens"]
        assert ours["lineage_spec"] == d[drug]["lineage_collapsed"]["spec"]
        assert ours["raw_sens"] == d[drug]["raw"]["sens"]
        assert ours["raw_spec"] == d[drug]["raw"]["spec"]


def test_the_lineage_number_is_the_one_presented_as_the_headline(capsys):
    """TB is monomorphic and heavily clonal: raw counts one vote per isolate over ~2,845 isolates that
    collapse to ~67 lineages, so quoting raw would inflate sensitivity roughly 2x."""
    if not _catalogue_ready():
        pytest.skip("gitignored WHO catalogue absent")
    vcf = _a_real_vcf()
    if vcf is None:
        pytest.skip("no real TB VCF on this host")
    tb_cli.main(["--vcf", str(vcf), "--drug", "rifampicin", "--catalogue-dir", str(CAT)])
    out = capsys.readouterr().out
    assert "LINEAGE-COLLAPSED" in out and "0.444" in out
    assert "clonality-INFLATED" in out, "the raw figure must be labelled where it appears"


def test_an_in_distribution_drug_does_NOT_borrow_the_independent_numbers(capsys):
    if not _catalogue_ready():
        pytest.skip("gitignored WHO catalogue absent")
    vcf = _a_real_vcf()
    if vcf is None:
        pytest.skip("no real TB VCF on this host")
    tb_cli.main(["--vcf", str(vcf), "--drug", "ethambutol", "--catalogue-dir", str(CAT)])
    out = capsys.readouterr().out
    assert "IN-DISTRIBUTION" in out and "NOT independent validation" in out
    assert "0.444" not in out and "0.321" not in out, "another drug's number must not leak in"


def test_list_drugs_tags_EVERY_drug_with_its_OWN_tier(capsys):
    """Asserted PER LINE, not by document order: the listing is alphabetical, so an in-distribution
    drug (amikacin) prints first and a split-on-first-marker test would pass for the wrong reason."""
    assert tb_cli.main(["--list-drugs"]) == 0
    lines = {ln.split()[0]: ln for ln in capsys.readouterr().out.splitlines()
             if ln.startswith("   ")}
    assert len(lines) == 12, lines
    for drug, ln in lines.items():
        if drug in tb_cli.INDEPENDENT_DRUGS:
            assert "INDEPENDENT" in ln and "IN-DISTRIBUTION" not in ln, ln
        else:
            assert "IN-DISTRIBUTION" in ln and "NOT validation" in ln, ln
            assert "0.444" not in ln and "0.321" not in ln, f"borrowed number: {ln}"


# --------------------------------------------------------------- registration

def test_the_route_declares_a_regime():
    from dna_decode.data.cell_regime import regime_for_route
    assert regime_for_route("dna-tb") == "curated_catalog_exists"


def test_every_routable_tb_drug_has_exactly_one_contract():
    """COVERAGE, and NON-VACUOUS BY CONSTRUCTION: the routable set is read from the CATALOGUE (what
    `--drug` validates against) while the contracts are written by hand, so the two sides are
    independent. Deriving the routable set from the contracts would compare a set to itself and pass
    however many contracts were missing -- which is precisely how HCMV shipped five decoders with none."""
    from dna_decode.data.cell_registry import cells
    from dna_decode.data.routable_drugs import all_routable_tb_drugs
    tb = [c for c in cells() if c.route == "dna-tb"]
    covered = {c.target for c in tb}
    routable = all_routable_tb_drugs()
    assert covered == routable, f"extra={covered - routable} missing={routable - covered}"
    assert len(tb) == len(covered), "duplicate tb cell"


def test_the_tb_track_is_separate_from_the_FROZEN_amr_projection():
    """The `amr` track is projected verbatim from the frozen `shipped_decoder_surface` and a guard
    asserts set equality, so a TB cell on that track would break it. TB is NOT in the frozen surface by
    design -- adding it would edit a sha256-pinned file and retire the active prospective lock."""
    from dna_decode.data import shipped_decoder_surface as sds
    from dna_decode.data.cell_registry import cells
    tb = [c for c in cells() if c.route == "dna-tb"]
    assert tb and all(c.track == "tb" for c in tb)
    surface_organisms = {k[0] for k in sds.surface_index()}
    assert not any("tubercul" in o.lower() for o in surface_organisms), \
        "TB must stay out of the frozen surface"


def test_the_two_tiers_are_not_collapsed():
    """Ten drugs have only an in-distribution number; two have a provenance-disjoint one. One tier for
    all twelve would over-claim ten or under-claim two, and under-claiming is as much a trust-surface
    falsehood as over-claiming."""
    from dna_decode.data.cell_registry import EvidenceTier, cells
    tb = {c.target: c for c in cells() if c.route == "dna-tb"}
    assert tb["rifampicin"].evidence_tier is EvidenceTier.NEAR_INDEPENDENT
    assert tb["isoniazid"].evidence_tier is EvidenceTier.NEAR_INDEPENDENT
    for d in ("ethambutol", "bedaquiline", "linezolid", "amikacin"):
        assert tb[d].evidence_tier is EvidenceTier.KNOWLEDGE_BASELINE, d
        assert "0.444" not in tb[d].validation_slice, f"{d} must not borrow rifampicin's number"


def test_the_console_entry_exists_in_pyproject():
    src = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    assert 'dna-tb = "dna_decode.organism_rules.tb_cli:main"' in src
