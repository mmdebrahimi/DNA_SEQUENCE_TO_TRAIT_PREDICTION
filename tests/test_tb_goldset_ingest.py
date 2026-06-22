"""Offline tests for the TB gold-set independence check + manifest ingester (deliverable-b helper)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.organism_rules import tb_goldset  # noqa: E402
from scripts.build_tb_goldset_manifest import manifest_rows, read_candidates  # noqa: E402


def test_cryptic_accessions_unions_columns(tmp_path):
    p = tmp_path / "reuse.csv"
    p.write_text("ENA_RUN,UNIQUEID,ENA_SAMPLE\nERR1,site.01,ERS9\nERR2,site.02,\n", encoding="utf-8")
    acc = tb_goldset.cryptic_accessions(p)
    assert acc == {"ERR1", "ERR2", "SITE.01", "SITE.02", "ERS9"}   # union, upper-cased; blanks skipped


def test_assert_independent_partitions_and_is_case_insensitive():
    cryptic = {"ERR1", "ERS9"}
    rep = tb_goldset.assert_independent(["err1", "ERRX", "ERS9", "SAMEA_NEW"], cryptic)
    assert set(rep.leaked) == {"err1", "ERS9"}        # err1 matches ERR1 case-insensitively
    assert set(rep.clean) == {"ERRX", "SAMEA_NEW"}
    assert rep.n_checked == 4


def test_manifest_rows_filters_to_clean_RS_per_drug():
    cand = [
        {"strain_id": "g1", "masked_vcf": "g1.vcf", "regeno_vcf": "g1.rg.vcf", "rif_label": "R", "inh_label": "S"},
        {"strain_id": "g2", "masked_vcf": "g2.vcf", "regeno_vcf": "", "rif_label": "", "inh_label": "R"},
        {"strain_id": "g3", "masked_vcf": "g3.vcf", "regeno_vcf": "", "rif_label": "U", "inh_label": "I"},  # both ambiguous
    ]
    rif = manifest_rows(cand, "RIF")
    assert [r["strain_id"] for r in rif] == ["g1"]                 # g2 blank, g3 'U' excluded
    assert rif[0] == {"strain_id": "g1", "masked_vcf": "g1.vcf", "regeno_vcf": "g1.rg.vcf", "label": "R"}
    inh = manifest_rows(cand, "INH")
    assert [r["strain_id"] for r in inh] == ["g1", "g2"]           # g3 'I' excluded
    assert inh[1]["regeno_vcf"] is None                            # blank regeno -> None


def test_read_candidates_tsv(tmp_path):
    p = tmp_path / "cand.tsv"
    p.write_text("strain_id\tena_accession\tmasked_vcf\tregeno_vcf\trif_label\tinh_label\n"
                 "s1\tERRX\ts1.vcf\t\tR\tS\n", encoding="utf-8")
    rows = read_candidates(p)
    assert rows == [{"strain_id": "s1", "ena_accession": "ERRX", "masked_vcf": "s1.vcf",
                     "regeno_vcf": "", "rif_label": "R", "inh_label": "S"}]


def test_assert_independent_aliased_catches_any_alias():
    # CRyPTIC carries ENA-style; a candidate is SRA-style but shares its biosample alias -> leaked.
    cryptic = {"ERR1", "SAMEA999"}
    aliases = {
        "indep1": ["SRR_NEW", "SRS_NEW", "SAMN_NEW"],      # no overlap -> clean
        "leak_by_biosample": ["SRR_X", "SAMEA999"],         # biosample alias in CRyPTIC -> leaked
        "leak_by_run": ["err1", "SAMN_Y"],                  # run alias matches ERR1 (case-insensitive) -> leaked
    }
    rep = tb_goldset.assert_independent_aliased(aliases, cryptic)
    assert set(rep.clean) == {"indep1"}
    assert set(rep.leaked) == {"leak_by_biosample", "leak_by_run"}
    assert rep.n_checked == 3


def test_validate_candidates_ok_and_balance():
    from scripts.validate_tb_goldset_candidates import validate
    header = ["strain_id", "run_accession", "masked_vcf", "regeno_vcf", "rif_label", "inh_label"]
    rows = [
        {"strain_id": "a", "run_accession": "SRR1", "masked_vcf": "", "regeno_vcf": "", "rif_label": "S", "inh_label": "S"},
        {"strain_id": "b", "run_accession": "SRR2", "masked_vcf": "", "regeno_vcf": "", "rif_label": "R", "inh_label": "R"},
    ]
    ok, errors, summary = validate(rows, header, require_vcf=False)   # no --require-vcf -> blank VCF ok
    assert ok and not errors
    assert summary["n_usable"] == 2
    assert summary["rif_inh_buckets"]["RIF=S/INH=S"] == 1 and summary["rif_inh_buckets"]["RIF=R/INH=R"] == 1


def test_validate_candidates_failures():
    from scripts.validate_tb_goldset_candidates import validate
    base_header = ["strain_id", "run_accession", "masked_vcf", "regeno_vcf", "rif_label", "inh_label"]
    # no accession-alias column at all -> structural fail
    ok, errs, _ = validate([], ["strain_id", "masked_vcf", "regeno_vcf", "rif_label", "inh_label"])
    assert not ok and any("accession-alias" in e for e in errs)
    # bad label value
    ok, errs, _ = validate([{"strain_id": "a", "run_accession": "SRR1", "masked_vcf": "", "regeno_vcf": "",
                             "rif_label": "RESISTANT", "inh_label": "S"}], base_header)
    assert not ok and any("bad labels in rif_label" in e for e in errs)
    # duplicate strain_id
    dup = [{"strain_id": "x", "run_accession": "SRR1", "masked_vcf": "", "regeno_vcf": "", "rif_label": "R", "inh_label": "S"},
           {"strain_id": "x", "run_accession": "SRR2", "masked_vcf": "", "regeno_vcf": "", "rif_label": "S", "inh_label": "S"}]
    ok, errs, _ = validate(dup, base_header)
    assert not ok and any("duplicate strain_id" in e for e in errs)


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
