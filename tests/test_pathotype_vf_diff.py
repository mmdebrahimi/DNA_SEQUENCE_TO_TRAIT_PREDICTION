"""Contract tests for the v0.1 canonical-VirulenceFinder side-by-side diff (vf_diff).

Pins the ledger-promised "side-by-side comparison against CGE VirulenceFinder gene-call
output" that completes the v0 spec. The diff compares the fast k-mer-seed resolver
(detect.py) against a canonical blastn caller over the SAME VF DB — an AUDIT of the fast
caller, NOT an independent baseline (interrogation Q2, 2026-06-04).

Offline-safe: structural + honesty + degradation assertions run WITHOUT blastn (so the
test_pathotype_vf_diff endpoint is green in CI on a host with no binary). The real
agreement assertion runs only when blastn is installed. Runnable via pytest OR standalone
(`python tests/test_pathotype_vf_diff.py`) so it satisfies a test-exit-0 gate without pytest.
"""
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.pathotype.cli import DEFAULT_DB, analyze
from dna_decode.pathotype.vf_runner import (
    NON_INDEPENDENCE_CAVEAT, build_vf_diff, find_blastn, run_canonical_vf,
)

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / DEFAULT_DB

# CI / fresh checkout: the VirulenceFinder DB lives under gitignored data/. Skip the module
# when absent (repo convention for data-gated tests) instead of hard-failing.
pytestmark = pytest.mark.skipif(
    not DB.exists(),
    reason="gitignored VirulenceFinder DB absent (data/virulencefinder_db/virulence_ecoli.fsa)",
)


def _first_allele(db_path: Path, prefix: str) -> tuple[str, str]:
    """Return (header_id, sequence) of the first DB allele whose id starts with `prefix`."""
    name, buf, want = None, [], None
    for line in db_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if want is not None:
                return want, "".join(buf)
            hid = line[1:].split(":")[0]
            if hid.lower().startswith(prefix.lower()):
                want, name, buf = line[1:].split()[0], hid, []
            else:
                name, buf = None, []
        elif name is not None or want is not None:
            buf.append(line.strip())
    if want is not None:
        return want, "".join(buf)
    raise AssertionError(f"no allele with prefix {prefix} in {db_path}")


def _synthetic_stx2_assembly(tmp: Path) -> Path:
    """Build a deterministic assembly FASTA embedding a real stx2 allele (so BOTH callers
    detect STX2). Flanked by fixed filler so it reads as a contig, not a bare gene."""
    _, seq = _first_allele(DB, "stx2")
    flank = ("ACGTACGTAC" * 60)
    contig = flank + seq + flank
    fa = tmp / "synthetic_stx2.fna"
    fa.write_text(f">synthetic_stx2_contig\n{contig}\n", encoding="utf-8")
    return fa


def test_vf_diff_section_shape():
    """vf_diff exists with per-gene + per-cluster + honesty-flag + concordance shape."""
    with tempfile.TemporaryDirectory() as td:
        fa = _synthetic_stx2_assembly(Path(td))
        rec = analyze(str(fa), str(DB), "shape-test", vf_diff=True)
    diff = rec["vf_diff"]
    assert "vf_diff" in rec
    assert diff["status"] in ("ok", "unavailable")
    assert "per_gene" in diff and "per_cluster" in diff
    assert "caveat" in diff and diff["caveat"] == NON_INDEPENDENCE_CAVEAT
    if diff["status"] == "ok":
        assert "concordance" in diff and diff["concordance"]["n_clusters"] >= 1
        # per-cluster rows carry both calls + an agreement flag
        row = diff["per_cluster"][0]
        for k in ("cluster", "resolver_called", "canonical_called", "agreement"):
            assert k in row


def test_caller_is_not_independent_baseline():
    """The non-independence honesty flag is carried into the diff section (circularity)."""
    with tempfile.TemporaryDirectory() as td:
        fa = _synthetic_stx2_assembly(Path(td))
        rec = analyze(str(fa), str(DB), "honesty-test", vf_diff=True)
    assert rec["vf_diff"]["caller_is_independent_baseline"] is False
    # the top-level caller block already records this; the diff must agree.
    assert rec["caller"]["caller_is_independent_baseline"] is False


def test_run_canonical_vf_additive_contract_keys():
    """v2 adds per_hit + db_sha to run_canonical_vf WITHOUT disturbing the per_gene/
    per_cluster best-hit contract build_vf_diff reads (regression pin)."""
    res = run_canonical_vf(str(DB), str(DB), blastn_bin="/nonexistent/blastn")
    for k in ("status", "per_gene", "per_cluster", "per_hit", "db_sha"):
        assert k in res
    assert isinstance(res["per_hit"], list)
    # build_vf_diff still consumes the unavailable result unchanged
    diff = build_vf_diff({"STX2": True}, res, resolver_marker_hits=[])
    assert diff["status"] == "unavailable"


def test_offline_safe_degrades_when_blastn_absent():
    """When canonical VF is unavailable the diff section is RETAINED as status=unavailable
    with a reason — never silently dropped (CI-on-no-binary contract)."""
    unavailable = {"status": "unavailable", "reason": "blastn not found", "per_gene": {},
                   "per_cluster": {}}
    diff = build_vf_diff({"STX2": True}, unavailable, resolver_marker_hits=[])
    assert diff["status"] == "unavailable"
    assert diff["reason"]
    assert diff["caller_is_independent_baseline"] is False
    assert diff["caveat"] == NON_INDEPENDENCE_CAVEAT
    # also via a bogus blastn path through the real runner
    res = run_canonical_vf(str(DB), str(DB), blastn_bin="/nonexistent/blastn")
    assert res["status"] == "unavailable" and res["reason"]


def test_canonical_and_resolver_agree_on_stx2():
    """On a known stx2-positive assembly, canonical VF and the resolver AGREE that STX2 is
    present. Resolver assertion runs always; canonical assertion runs only with blastn."""
    with tempfile.TemporaryDirectory() as td:
        fa = _synthetic_stx2_assembly(Path(td))
        rec = analyze(str(fa), str(DB), "stx2-agree", vf_diff=True)
    # resolver (k-mer) must see STX2 — its full allele is embedded in the contig.
    assert rec["cluster_profile"].get("STX2") is True
    diff = rec["vf_diff"]
    if find_blastn():
        assert diff["status"] == "ok"
        stx2 = next(r for r in diff["per_cluster"] if r["cluster"] == "STX2")
        assert stx2["canonical_called"] is True
        assert stx2["resolver_called"] is True
        assert stx2["agreement"] == "both_present"
    else:
        # no binary on this host → degraded but contract intact
        assert diff["status"] == "unavailable"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
