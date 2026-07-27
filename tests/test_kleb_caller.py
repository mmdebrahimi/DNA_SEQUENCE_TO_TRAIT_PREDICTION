"""Offline tests for the Klebsiella depolymerase->KL-type caller (fetch-only cell).

Pure caller logic (top-K ranking, abstain), the reference loader, LOO accounting, the resolve_reference
fetch-only degrade, and the CLI's no-reference INDETERMINATE path. No network, no bundled data.
"""
from __future__ import annotations

import json

import dna_decode.kleb.depolymerase_caller as dc
import dna_decode.kleb.cli as kcli
from dna_decode.cli import main as unified_main


def _ref(tmp_path, entries):
    p = tmp_path / "ref.faa"
    with open(p, "w", encoding="utf-8") as fh:
        fh.write("; test ref\n")
        for i, (kl, seq) in enumerate(entries):
            fh.write(f">{kl}|{i}\n{seq}\n")
    return p


def test_load_reference_parses_kltype(tmp_path):
    p = _ref(tmp_path, [("KL1", "MKVLAAAAWWWW"), ("KL2", "QQQQYYYYZZZZ")])
    km, kl = dc.load_reference(p, k=4)
    assert len(km) == 2 and set(kl.values()) == {"KL1", "KL2"}


def test_call_ranks_topk_and_abstains(tmp_path):
    km, kl = dc.load_reference(_ref(tmp_path, [
        ("KL1", "MKVLAAAAWWWWTTTT"), ("KL1", "MKVLAAAAWWWWTTTS"),
        ("KL2", "QQQQYYYYZZZZGGGG")]), k=4)
    # query identical to a KL1 domain -> KL1 ranked first
    call = dc.call_kltype("MKVLAAAAWWWWTTTT", km, kl, top_k=3)
    assert call.status == "CALLED" and call.ranked_kltypes[0] == "KL1"
    # a query with no shared k-mers -> INDETERMINATE (abstain, not a fabricated call)
    ab = dc.call_kltype("PPPPCCCCDDDDEEEE", km, kl, min_similarity=0.05)
    assert ab.status == "INDETERMINATE" and ab.ranked_kltypes == ()


def test_leave_one_out_topk(tmp_path):
    km, kl = dc.load_reference(_ref(tmp_path, [
        ("KL1", "MKVLAAAAWWWWTTTT"), ("KL1", "MKVLAAAAWWWWTTTS"),
        ("KL2", "QQQQYYYYZZZZGGGG"), ("KL2", "QQQQYYYYZZZZGGGH")]), k=4)
    r = dc.leave_one_out(km, kl, top_k=5)
    assert r.n == 4 and r.called == 4 and r.top1 == 4          # each finds its same-KL twin
    assert r.per_kl["KL1"] == [2, 2] and r.per_kl["KL2"] == [2, 2]


def test_resolve_reference_fetch_only_returns_none(tmp_path, monkeypatch):
    # no env, no known local file -> None (fetch-only; nothing bundled)
    monkeypatch.delenv("DPO_KLEB_REFERENCE", raising=False)
    monkeypatch.setattr(dc, "_DEFAULT_REF_CANDIDATES", (str(tmp_path / "nope.faa"),))
    assert dc.resolve_reference(None) is None
    # explicit path that exists -> returned
    p = _ref(tmp_path, [("KL1", "MKVL")])
    assert dc.resolve_reference(str(p)) == p


def test_cli_degrades_without_reference(tmp_path, monkeypatch, capsys):
    monkeypatch.delenv("DPO_KLEB_REFERENCE", raising=False)
    monkeypatch.setattr(dc, "_DEFAULT_REF_CANDIDATES", (str(tmp_path / "nope.faa"),))
    q = tmp_path / "q.faa"; q.write_text(">q\nMKVLAAAA\n", encoding="utf-8")
    rc = kcli.main(["--depolymerase-fasta", str(q), "--json"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 0 and out["status"] == "INDETERMINATE"
    assert "fetch_dpotropisearch" in out["reason"]          # actionable, license-noticed


def test_cli_calls_with_reference_and_routes_through_unified(tmp_path, monkeypatch, capsys):
    p = _ref(tmp_path, [("KL5", "MKVLAAAAWWWWTTTT"), ("KL5", "MKVLAAAAWWWWTTTS"), ("KL9", "QQQQYYYYZZZZ")])
    q = tmp_path / "q.faa"; q.write_text(">q\nMKVLAAAAWWWWTTTT\n", encoding="utf-8")
    rc = unified_main(["kleb", "--depolymerase-fasta", str(q), "--reference", str(p), "--json"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 0 and out["status"] == "CALLED" and out["ranked_kltypes"][0] == "KL5"
