"""Pins the GeT-RM concordance harness (scripts/pgx_getrm_concordance.py).

Pure normalization/equivalence logic always runs. The full concordance number (72/72 core) re-runs only
when the gitignored 1000G VCF is present (Docker-fetched) -> it pins the published report-card number +
the *38==*1 phenotype-equivalence + the genuine-residual bucketing.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from pgx_getrm_concordance import REF_EQUIV, _norm  # noqa: E402

_VCF = REPO / "data" / "pgx_1000g" / "cyp2c19_1000g.vcf.gz"
_TRUTH = REPO / "tests" / "data" / "pgx_getrm" / "star-allele-comparison_common.tsv"


def test_norm_handles_slash_and_pipe_and_sorts():
    assert _norm("*17|*1") == ("*1", "*17")
    assert _norm("*1/*1") == ("*1", "*1")
    assert _norm("*2/*35") == ("*2", "*35")
    assert _norm("") == ()


def test_star38_is_reference_equivalent():
    # *38 is the true variant-free reference -> phenotype-identical to *1
    assert REF_EQUIV.get("*38") == "*1"


def test_truth_set_committed_and_has_cyp2c19_column():
    import csv
    with open(_TRUTH, encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        assert "CYP2C19_getrm_ngs" in rdr.fieldnames
        rows = list(rdr)
    assert len(rows) >= 80   # ~87 GeT-RM x 1000G overlap samples
    # the documented NA19122 *2/*35 case is present
    assert any(r["Coriell"] == "NA19122" and "*35" in r["CYP2C19_getrm_ngs"] for r in rows)


@pytest.mark.skipif(not _VCF.exists(), reason="1000G VCF not present (Docker-fetched, gitignored)")
def test_getrm_core_concordance_is_72_of_72():
    from pgx_getrm_concordance import main
    rc = main(["--gene", "cyp2c19"])
    assert rc == 0
    import json
    rep = json.loads((REPO / "wiki" / "pgx_getrm_concordance_2026-06-25.json").read_text(encoding="utf-8"))
    assert rep["core_diplotype_hits"] == "72/72"          # caller perfect on the alleles it covers
    assert rep["core_diplotype_concordance"] == 1.0
    assert rep["noncore_correctly_withheld"] == 2          # *4 + *35 sentinels fired (incl. NA19122)
    assert rep["genuine_silent_miscall"] == 6              # *8/*13/*15/*39 residual (honest blind spot)
    assert rep["caller_is_independent_of_consensus_tools"] is True


_C9_VCF = REPO / "data" / "pgx_1000g" / "cyp2c9_1000g.vcf.gz"


@pytest.mark.skipif(not _C9_VCF.exists(), reason="CYP2C9 1000G region VCF not present (Docker-fetched)")
def test_cyp2c9_getrm_core_concordance_is_73_of_73():
    from pgx_getrm_concordance import main
    rc = main(["--gene", "cyp2c9"])
    assert rc == 0
    import json
    rep = json.loads((REPO / "wiki" / "pgx_getrm_concordance_cyp2c9_2026-06-25.json").read_text(encoding="utf-8"))
    assert rep["gene"] == "CYP2C9"
    assert rep["core_diplotype_hits"] == "73/73"          # caller perfect on *1/*2/*3
    # v0.1 sentinels (*5/*8/*9/*11) now WITHHOLD non-core SNP alleles -> residual is the indel/uncatalogued tail
    assert rep["noncore_correctly_withheld"] == 10
    assert rep["genuine_silent_miscall"] == 4             # *6 (indel) + *61 + undefined -- documented residual


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
