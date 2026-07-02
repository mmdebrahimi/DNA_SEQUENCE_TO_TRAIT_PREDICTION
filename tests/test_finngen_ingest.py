"""Pin the FinnGen ingest clump-to-lead-per-gene + GW-sig filter on a synthetic gz (no network/big-file)."""
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.finngen_ingest import ingest  # noqa: E402

_HEADER = "#chrom\tpos\tref\talt\trsids\tnearest_genes\tpval\tmlogp\tbeta\tsebeta\taf_alt\taf_alt_cases\taf_alt_controls"
_ROWS = [
    # gene HLA: two GW-sig variants -> keep the stronger (mlogp 8.5); a risk beta
    "6\t100\tA\tG\trs1\tHLA\t3e-9\t8.5\t0.4\t0.06\t0.1\t0.12\t0.1",
    "6\t200\tC\tT\trs2\tHLA\t1e-8\t8.0\t0.3\t0.06\t0.1\t0.11\t0.1",
    # gene FOO: one GW-sig protective variant
    "1\t50\tG\tA\trs3\tFOO\t2e-8\t7.7\t-0.5\t0.08\t0.05\t0.04\t0.05",
    # gene BAR: NOT significant (mlogp 4) -> excluded
    "2\t99\tT\tC\trs4\tBAR\t1e-4\t4.0\t0.2\t0.05\t0.2\t0.21\t0.2",
]


def _write_gz(tmp: Path) -> Path:
    p = tmp / "syn.gz"
    with gzip.open(p, "wt") as f:
        f.write(_HEADER + "\n" + "\n".join(_ROWS) + "\n")
    return p


def test_gw_sig_filter_and_clump(tmp_path):
    res = ingest(_write_gz(tmp_path), "SYN", n_cases=100, n_controls=900)
    assert res["n_variants"] == 4
    assert res["n_genome_wide_sig"] == 3          # HLAx2 + FOO; BAR excluded
    assert res["n_gene_loci"] == 2                # clumped to HLA + FOO
    by_gene = {r["gene"]: r for r in res["lead_loci_by_gene"]}
    assert by_gene["HLA"]["mlogp"] == 8.5          # stronger of the two HLA variants
    assert by_gene["HLA"]["direction"] == "risk"
    assert by_gene["FOO"]["direction"] == "protective"


def test_ordered_by_mlogp(tmp_path):
    res = ingest(_write_gz(tmp_path), "SYN")
    mlogps = [r["mlogp"] for r in res["lead_loci_by_gene"]]
    assert mlogps == sorted(mlogps, reverse=True)


def test_scope_note_present(tmp_path):
    # honesty rail: the artifact must carry the host-genetics / not-pathogen-decoder scope note
    res = ingest(_write_gz(tmp_path), "SYN")
    assert "HOST" in res["scope_note"] and "does NOT feed the" in res["scope_note"]
