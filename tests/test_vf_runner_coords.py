"""Step 1 tests — coord-retaining VF caller (per_hit + the no-parse_seqids pin).

Two layers:
  - PURE parse tests over synthetic blastn `-outfmt 6 qseqid sseqid sstart send pident
    length qlen` stdout (run WITHOUT blastn): minus-strand normalization, tandem-copy
    retention, overlap dedup, unclustered alleles kept with cluster=None, sub-threshold
    HSPs excluded from per_hit, and per_gene/per_cluster byte-identical to the best-hit
    contract `build_vf_diff` depends on.
  - a LIVE integration test (skipped when blastn is absent) that pins, on real blastn,
    that `sseqid` equals the exact FASTA first-token and that `-parse_seqids` is never
    passed (the empirically-verified contig-name-map contract).
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import pytest

from dna_decode.pathotype.cli import DEFAULT_DB
from dna_decode.pathotype.vf_runner import (
    VF_COVERAGE_THRESHOLD,
    VF_IDENTITY_THRESHOLD,
    _interval_dedup,
    find_blastn,
    parse_blastn_outfmt6,
    run_canonical_vf,
)

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / DEFAULT_DB

ID = VF_IDENTITY_THRESHOLD
COV = VF_COVERAGE_THRESHOLD


def _row(qseqid, sseqid, sstart, send, pident, length, qlen) -> str:
    return f"{qseqid}\t{sseqid}\t{sstart}\t{send}\t{pident}\t{length}\t{qlen}"


def _synthetic_stdout() -> str:
    return "\n".join([
        # clustered, plus strand, full coverage -> called
        _row("stx2A:acc1", "contigX", 100, 1341, 99.5, 1242, 1242),
        # clustered, MINUS strand (sstart > send) -> normalized + strand "-"
        _row("stx1A:acc2", "contigX", 2000, 1001, 98.0, 1000, 1000),
        # UNCLUSTERED allele (no CLUSTER_MARKERS prefix), two TANDEM copies at distinct coords
        _row("blaTEM-1:acc3", "contigY", 5000, 5860, 100.0, 861, 861),
        _row("blaTEM-1:acc3", "contigY", 7000, 7860, 100.0, 861, 861),
        # an OVERLAPPING HSP of the first copy -> must collapse into one (dedup)
        _row("blaTEM-1:acc3", "contigY", 5010, 5870, 99.0, 860, 861),
        # clustered but SUB-THRESHOLD coverage (300/900 = 33%) -> excluded from per_hit
        _row("eae:acc4", "contigZ", 10, 500, 99.0, 300, 900),
    ])


def test_per_hit_normalizes_strand_and_keeps_tandem_and_unclustered():
    calls = parse_blastn_outfmt6(_synthetic_stdout(), ID, COV)
    per_hit = calls["per_hit"]
    by_gene = {}
    for h in per_hit:
        by_gene.setdefault(h["vf_gene"], []).append(h)

    # stx2 — plus strand, normalized coords
    stx2 = by_gene["stx2A"][0]
    assert stx2["cluster"] == "STX2"
    assert (stx2["start"], stx2["stop"], stx2["strand"]) == (100, 1341, "+")

    # stx1 — minus strand normalized to start<=stop with strand "-"
    stx1 = by_gene["stx1A"][0]
    assert stx1["cluster"] == "STX1"
    assert (stx1["start"], stx1["stop"], stx1["strand"]) == (1001, 2000, "-")

    # blaTEM — UNCLUSTERED (cluster=None) but RETAINED (C2); two distinct copies after dedup
    tem = sorted(by_gene["blaTEM-1"], key=lambda h: h["start"])
    assert all(h["cluster"] is None for h in tem)
    assert len(tem) == 2, "tandem copies retained; the overlapping HSP collapsed into one"
    assert tem[0]["start"] == 5000 and tem[0]["stop"] == 5870  # merged 5000-5860 + 5010-5870
    assert tem[1]["start"] == 7000 and tem[1]["stop"] == 7860

    # eae — sub-threshold coverage -> NOT a called per_hit
    assert "eae" not in by_gene
    assert all(h["called"] is True for h in per_hit)


def test_per_gene_per_cluster_byte_identical_best_hit_contract():
    """per_gene/per_cluster stay best-hit + cluster-scoped (unclustered alleles excluded)
    so build_vf_diff is unchanged. Pins the exact pre-change shape."""
    calls = parse_blastn_outfmt6(_synthetic_stdout(), ID, COV)
    per_gene, per_cluster = calls["per_gene"], calls["per_cluster"]

    # only clustered alleles appear in per_gene; blaTEM (unclustered) excluded
    assert set(per_gene) == {"stx2A:acc1", "stx1A:acc2", "eae:acc4"}
    assert per_gene["stx2A:acc1"] == {
        "cluster": "STX2", "percent_identity": 99.5, "percent_coverage": 100.0, "called": True}
    assert per_gene["eae:acc4"]["called"] is False  # sub-threshold

    # per_cluster carries every catalog cluster; STX2 called, LEE present-but-not-called
    assert per_cluster["STX2"]["called"] is True
    assert per_cluster["STX2"]["best_gene"] == "stx2A:acc1"
    assert per_cluster["LEE"]["called"] is False
    # a cluster with zero allele hits is still present as a not-called row
    assert per_cluster["ST"]["called"] is False and per_cluster["ST"]["best_gene"] is None


@pytest.mark.skipif(not DB.exists(), reason="gitignored VirulenceFinder DB absent (db_sha needs a readable DB)")
def test_run_canonical_vf_unavailable_carries_per_hit_and_db_sha():
    """The additive contract holds even on the offline-degrade path (no blastn)."""
    res = run_canonical_vf(str(DB), str(DB), blastn_bin="/nonexistent/blastn")
    assert res["status"] == "unavailable"
    assert res["per_hit"] == []
    assert res["db_sha"] and len(res["db_sha"]) == 16  # DB is readable -> sha present


def test_parse_blastn_outfmt6_skips_malformed_and_zero_qlen():
    """Robust to a non-tab-7-column line, a non-numeric field, and qlen=0 (no ZeroDivision).
    None of the junk rows produce a call; the one good row still parses."""
    stdout = "\n".join([
        "too\tfew\tcols",                                   # < 7 fields -> skipped
        _row("stx2A:bad", "c", "notanint", 900, 99.0, 100, 100),  # ValueError -> skipped
        _row("stx2A:zero", "c", 1, 900, 99.0, 100, 0),      # qlen=0 -> cov 0 -> not called
        _row("stx2A:good", "contigX", 100, 1341, 99.5, 1242, 1242),  # the only real call
    ])
    calls = parse_blastn_outfmt6(stdout, ID, COV)
    assert [h["allele_id"] for h in calls["per_hit"]] == ["stx2A:good"]
    # the bad/zero rows never created a per_gene entry
    assert "stx2A:bad" not in calls["per_gene"]
    # zero-qlen row recorded a best-hit but at 0% coverage -> not called (no crash)
    assert calls["per_gene"].get("stx2A:zero", {}).get("called") in (None, False)


def test_interval_dedup_merge_keeps_higher_coverage_hsp_attrs():
    """When two overlapping HSPs of one copy merge, the merged interval inherits the
    HIGHER-coverage HSP's identity/coverage/strand and the widest span; a disjoint HSP
    stays a separate copy (the tandem-retention contract)."""
    hsps = [
        {"start": 100, "stop": 200, "strand": "+", "percent_identity": 90.0, "percent_coverage": 70.0},
        {"start": 150, "stop": 260, "strand": "-", "percent_identity": 99.0, "percent_coverage": 95.0},
        {"start": 5000, "stop": 5100, "strand": "+", "percent_identity": 98.0, "percent_coverage": 88.0},
    ]
    out = _interval_dedup(hsps)
    assert len(out) == 2  # the first two overlap -> one copy; the third is disjoint
    merged = out[0]
    assert (merged["start"], merged["stop"]) == (100, 260)            # widest span
    assert merged["percent_identity"] == 99.0 and merged["percent_coverage"] == 95.0  # higher-cov wins
    assert merged["strand"] == "-"
    assert (out[1]["start"], out[1]["stop"]) == (5000, 5100)


# ---- LIVE integration (real blastn) ----

def _minus_strand_stx2_assembly(tmp: Path) -> Path:
    """Embed the reverse complement of a real stx2 allele so blastn reports it on the
    MINUS strand (sstart>send) — exercises coord normalization end-to-end."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    name, buf, want = None, [], None
    for line in DB.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if want is not None:
                break
            hid = line[1:].split(":")[0]
            if hid.lower().startswith("stx2"):
                want, buf = line[1:].split()[0], []
            else:
                want = None
        elif want is not None:
            buf.append(line.strip())
    seq = "".join(buf)
    revcomp = seq.translate(comp)[::-1]
    flank = "ACGTACGTAC" * 60
    fa = tmp / "synthetic_stx2_minus.fna"
    fa.write_text(f">synthetic_stx2_contig\n{flank}{revcomp}{flank}\n", encoding="utf-8")
    return fa


def test_live_sseqid_is_fasta_first_token_and_no_parse_seqids(monkeypatch):
    if not find_blastn():
        pytest.skip("blastn not installed on this host")
    with tempfile.TemporaryDirectory() as td:
        fa = _minus_strand_stx2_assembly(Path(td))

        captured: list[list] = []
        real_run = subprocess.run

        def spy(cmd, *a, **k):
            captured.append([str(x) for x in cmd])
            return real_run(cmd, *a, **k)

        monkeypatch.setattr(subprocess, "run", spy)
        res = run_canonical_vf(str(fa), str(DB), all_hits=True)

    assert res["status"] == "ok"
    # the no-parse_seqids pin: never passed to makeblastdb OR blastn
    assert all("-parse_seqids" not in cmd for cmd in captured)
    # the coord outfmt is used
    assert any("6 qseqid sseqid sstart send pident length qlen" in cmd for cmd in captured)
    # sseqid is the EXACT FASTA first-token (the shared contig-name-map contract)
    stx2 = [h for h in res["per_hit"] if h["vf_gene"].lower().startswith("stx2")]
    assert stx2, "stx2 should be called on the embedded allele"
    h = stx2[0]
    assert h["sseqid"] == "synthetic_stx2_contig"
    assert h["start"] < h["stop"]
    assert h["strand"] == "-"  # embedded as reverse complement
    assert res["db_sha"] and len(res["db_sha"]) == 16


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q"]))
