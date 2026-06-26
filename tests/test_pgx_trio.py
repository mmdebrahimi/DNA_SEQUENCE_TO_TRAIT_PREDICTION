"""Pins the PGx trio Mendelian co-segregation QC (Unit B). Pure logic always runs; the real-data run
asserts the published consistency number when the 1000G region VCFs are present (Docker-fetched)."""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from pgx_trio_concordance import _mendelian_ok  # noqa: E402

_C19_VCF = REPO / "data" / "pgx_1000g" / "cyp2c19_1000g.vcf.gz"


@pytest.mark.parametrize("child,father,mother,ok", [
    (["*1", "*2"], ["*1", "*1"], ["*1", "*2"], True),    # *2 from mother, *1 from father
    (["*1", "*17"], ["*1", "*1"], ["*17", "*17"], True),
    (["*2", "*2"], ["*1", "*1"], ["*1", "*1"], False),    # neither parent carries *2 -> violation
    (["*2", "*3"], ["*1", "*2"], ["*1", "*3"], True),     # *2 from f, *3 from m
    (["*2", "*3"], ["*1", "*1"], ["*1", "*3"], False),    # *2 from neither -> violation
    (["*1", "*1"], ["*1", "*2"], ["*1", "*17"], True),
])
def test_mendelian_ok(child, father, mother, ok):
    assert _mendelian_ok(child, father, mother) is ok


@pytest.mark.skipif(not _C19_VCF.exists(), reason="1000G region VCF not present (Docker-fetched)")
def test_trio_run_produces_consistency_numbers():
    from pgx_trio_concordance import main
    rc = main()
    assert rc == 0
    import json
    rep = json.loads((REPO / "wiki" / "pgx_trio_mendelian_2026-06-25.json").read_text(encoding="utf-8"))
    c19 = rep["genes"]["CYP2C19"]
    assert c19["status"] == "ok"
    assert c19["n_callable"] > 50                      # ~600 trios in the panel
    assert c19["mendelian_consistency"] >= 0.95        # a correct caller should be highly consistent


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
