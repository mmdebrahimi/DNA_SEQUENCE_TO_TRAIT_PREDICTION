"""Unified genome profile — offline-safe composition (each section degrades independently)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.profile.cli import main  # noqa: E402


def test_cli_missing_fasta_exit2():
    assert main(["/nonexistent.fna"]) == 2


def test_profile_all_unavailable_still_succeeds(tmp_path):
    # point every DB at a nonexistent path -> every section 'unavailable', profile still returns 0
    g = tmp_path / "g.fna"; g.write_text(">c\nACGTACGTACGT\n", encoding="utf-8")
    out = tmp_path / "p.json"
    rc = main([str(g), "--pathotype-db", str(tmp_path / "no.fsa"),
               "--plasmid-db", str(tmp_path / "no.fsa"), "--serotype-db", str(tmp_path / "no.fsa"),
               "--resfinder-db-dir", str(tmp_path / "nodir"),
               "--pointfinder-db-dir", str(tmp_path / "nodir"), "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["schema"] == "genome-profile-v0"
    assert set(rec["decoders"]) == {"pathotype", "serotype", "plasmid", "resfinder", "pointfinder", "amr"}
    assert all(d["status"] != "ok" for d in rec["decoders"].values())   # all unavailable, none crashed
    assert rec["decoders_ok"] == 0 and rec["decoders_total"] == 6


def test_pathotype_section_self_guards_on_bad_db(tmp_path):
    # _pathotype must never raise — a present-but-garbage DB returns an error/unavailable status, not a crash
    from dna_decode.profile.cli import _pathotype
    bad = tmp_path / "bad.fsa"; bad.write_text("not a fasta\n", encoding="utf-8")
    g = tmp_path / "g.fna"; g.write_text(">c\nACGTACGT\n", encoding="utf-8")
    res = _pathotype(str(g), str(bad))
    assert res["status"] in ("error", "ok", "unavailable")   # never raised


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
