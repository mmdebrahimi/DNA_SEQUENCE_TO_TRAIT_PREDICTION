"""Offline tests for the novel-protein MSA fetcher (dna_decode/forward/msa_fetch.py).

No network: the HTTP calls are injected. Pins the cache-first contract (etiquette), the submit->poll->
download->extract protocol, the a3m extraction from a ColabFold-shaped tarball, and the input guard.
"""
import io
import sys
import tarfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward.msa_fetch import (  # noqa: E402
    MsaFetchError, _extract_uniref_a3m, _seq_key, evolution_table_for_sequence, fetch_msa,
)

UBIQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def _fake_tarball(a3m_text: str) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        data = a3m_text.encode()
        info = tarfile.TarInfo(name="uniref.a3m")
        info.size = len(data)
        tf.addfile(info, io.BytesIO(data))
    return buf.getvalue()


A3M = ">query\nMKAY\n>h1\nMKAY\n>h2\nLKAF\n>h3\nMRGY\n"


def test_seq_key_is_stable_and_mode_sensitive():
    assert _seq_key("MK", "env") == _seq_key("mk ", "env")     # normalized
    assert _seq_key("MK", "env") != _seq_key("MK", "")          # mode participates


def test_extract_uniref_a3m_from_tarball():
    assert _extract_uniref_a3m(_fake_tarball(A3M)) == A3M


def test_bad_sequence_is_refused(tmp_path):
    with pytest.raises(MsaFetchError, match="amino-acid"):
        fetch_msa(">header not stripped", tmp_path)
    with pytest.raises(MsaFetchError, match="amino-acid"):
        fetch_msa("MK-AY", tmp_path)                            # gaps not allowed in a query


def test_offline_ok_raises_on_cache_miss(tmp_path):
    with pytest.raises(MsaFetchError, match="not cached"):
        fetch_msa("MKAY", tmp_path, offline_ok=True)


def test_fetch_is_cache_first_no_network_on_hit(tmp_path):
    # pre-seed the cache; the injected HTTP MUST NOT be called
    key = _seq_key("MKAY", "env")
    (tmp_path / f"{key}.a3m").write_text(A3M, encoding="utf-8")

    def boom(*a, **k):
        raise AssertionError("network called on a cache hit")

    p = fetch_msa("MKAY", tmp_path, _post=boom, _get_text=boom, _get_bytes=boom)
    assert p.read_text(encoding="utf-8") == A3M


def test_fetch_submit_poll_download_protocol(tmp_path):
    calls = {"post": 0, "poll": 0, "dl": 0}

    def post(url, data, timeout=60):
        calls["post"] += 1
        assert data["q"].startswith(">query\n") and data["mode"] == "env"
        return '{"id": "TID", "status": "PENDING"}'

    def get_text(url, timeout=60):
        calls["poll"] += 1
        assert url.endswith("/ticket/TID")
        return '{"status": "COMPLETE"}'          # completes on first poll

    def get_bytes(url, timeout=120):
        calls["dl"] += 1
        assert url.endswith("/result/download/TID")
        return _fake_tarball(A3M)

    p = fetch_msa(UBIQ, tmp_path, poll_interval=0, _post=post, _get_text=get_text, _get_bytes=get_bytes)
    assert p.exists() and p.read_text(encoding="utf-8") == A3M
    assert calls == {"post": 1, "poll": 1, "dl": 1}
    # a re-fetch is now a cache hit -> no further network
    fetch_msa(UBIQ, tmp_path, _post=post, _get_text=get_text, _get_bytes=get_bytes)
    assert calls["post"] == 1


def test_job_error_raises(tmp_path):
    def post(url, data, timeout=60):
        return '{"id": "T", "status": "ERROR"}'
    with pytest.raises(MsaFetchError, match="did not complete"):
        fetch_msa(UBIQ, tmp_path, poll_interval=0, _post=post,
                  _get_text=lambda *a, **k: '{"status":"ERROR"}', _get_bytes=lambda *a, **k: b"")


def test_evolution_table_for_sequence_end_to_end_with_stubbed_fetch(tmp_path):
    # seed cache so the fetch is offline; the table then flows through site_independent_table
    key = _seq_key("MKAY", "env")
    (tmp_path / f"{key}.a3m").write_text(A3M, encoding="utf-8")
    tab = evolution_table_for_sequence("MKAY", tmp_path, offline_ok=True)
    assert tab and all(len(m) >= 3 for m in tab)
    assert any(m.startswith("M1") for m in tab)      # position 1 scored


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
