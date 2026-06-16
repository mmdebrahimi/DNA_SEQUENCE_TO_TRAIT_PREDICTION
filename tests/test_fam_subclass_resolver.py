"""Pin the deployed fam.tsv per-gene Subclass resolution that the gent+cef cohort scoring depends on.

The resolver (scripts/fam_subclass_resolver.py) reproduces AMRFinder v4.2.7 FAMILY-level Class/Subclass
from the deployed DB so a gene-presence table (e.g. Sci234) can be scored by the frozen decoder rule.
These tests pin BOTH the faithful resolutions (the E. coli drivers) AND the documented scope-limits
(rare ESBL/carbapenemase variants on a generic BETA-LACTAM family node) so a future fam.tsv bump that
silently shifts curation is caught.
"""
import pytest

from scripts.fam_subclass_resolver import FamResolver, FAM_TSV_DEPLOYED

pytestmark = pytest.mark.skipif(not FAM_TSV_DEPLOYED.exists(),
                                reason="deployed fam.tsv not present on this host")


@pytest.fixture(scope="module")
def R():
    return FamResolver()


@pytest.mark.parametrize("gene, sub_must_contain", [
    # cef drivers in E. coli -> EXTENDED-spectrum subclass (counted by the ceftriaxone rule)
    ("blaCTX-M-15", "CEPHALOSPORIN"),
    ("blaCMY-2", "CEPHALOSPORIN"),
    ("blaDHA-1", "CEPHALOSPORIN"),
    ("blaACT-1", "CEPHALOSPORIN"),
    ("blaNDM-5", "CARBAPENEM"),
    ("blaKPC-2", "CARBAPENEM"),
    # gent determinants -> GENTAMICIN subclass (counted by the gentamicin rule)
    ("aac(3)-IIa", "GENTAMICIN"),
    ("aac(3)-IId", "GENTAMICIN"),
    ("aac(3)-IV", "GENTAMICIN"),
    ("armA", "GENTAMICIN"),
    ("ant(2'')-Ia", "GENTAMICIN"),
])
def test_driver_resolutions(R, gene, sub_must_contain):
    cls, sub, kind = R.resolve(gene)
    assert sub is not None and sub_must_contain in sub, (gene, cls, sub, kind)


@pytest.mark.parametrize("gene, sub_must_not_contain", [
    # narrow/intrinsic beta-lactamases: ampicillin-R but NOT ceftriaxone-R -> NOT extended-spectrum
    ("blaTEM-1B", "CEPHALOSPORIN"),
    ("blaSHV-1", "CEPHALOSPORIN"),
    ("blaOXA-1", "CEPHALOSPORIN"),
    # plain aac(6')-Ib confers amikacin/tobramycin, NOT gentamicin (family node = generic AMINOGLYCOSIDE)
    ("aac(6')-Ib", "GENTAMICIN"),
])
def test_non_driver_resolutions(R, gene, sub_must_not_contain):
    cls, sub, kind = R.resolve(gene)
    assert sub is not None and sub_must_not_contain not in sub, (gene, cls, sub, kind)


@pytest.mark.parametrize("gene", ["blaTEM-52C", "blaSHV-12", "blaOXA-48"])
def test_documented_scope_limit_esbl_carbapenemase_on_generic_family(R, gene):
    """ESBL/carbapenemase variants whose fam.tsv FAMILY node is generic BETA-LACTAM resolve to
    BETA-LACTAM (a conservative under-call for cef). This is the SURFACED scope-limit, not a bug —
    pinned so a future curation change that fixes it is noticed."""
    cls, sub, kind = R.resolve(gene)
    assert sub == "BETA-LACTAM", (gene, cls, sub, kind)


def test_unresolved_returns_none(R):
    cls, sub, kind = R.resolve("not_a_real_gene_xyz")
    assert (cls, sub, kind) == (None, None, "unresolved")


def test_mdfA_intrinsic_efflux_unresolved(R):
    """mdf(A) is the intrinsic E. coli efflux pump (present in ~all E. coli); it has no curated
    gent/cef Subclass and is correctly dropped (unresolved) rather than miscounted."""
    cls, sub, kind = R.resolve("mdf(A)")
    assert cls is None
