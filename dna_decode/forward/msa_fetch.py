"""Novel-protein MSA search: sequence -> a homolog MSA, the SEARCH half of the run-time evolution pipeline.

`msa_evolution.py` turns an MSA into an evolution-score table; this module GETS the MSA for an arbitrary
protein sequence (one that has no precomputed MSA on disk). On a disk-tight host a local UniRef/BFD database
(tens of GB) is out, so the deployable route is a free network search: the **ColabFold MMseqs2 API**
(`api.colabfold.com`) returns a ready a3m alignment. `evolution_table_for_sequence` chains
fetch -> `site_independent_table`, so the whole modality-hybrid evolution component runs on any sequence.

**Etiquette (load-bearing):** the ColabFold API is a FREE SHARED resource for the community. This module is
**cache-first** (one search per unique sequence, cached to disk forever) and single-query -- NEVER batch or
loop it against the API. If you need many MSAs, run MMseqs2 locally. `fetch_msa(..., offline_ok=True)` on a
cache miss raises rather than hitting the network, for CI / airgapped use.

a3m note: the returned alignment's query row is all-uppercase (defines the match columns); homolog inserts
are lowercase. `msa_evolution.parse_a2m` extracts match columns by keeping uppercase + '-', so a3m and a2m
both parse correctly.
"""
from __future__ import annotations

import hashlib
import io
import json
import tarfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

COLABFOLD_API = "https://api.colabfold.com"
DEFAULT_CACHE = Path("D:/dna_decode_cache/msa")


class MsaFetchError(RuntimeError):
    pass


def _http_post(url: str, data: dict, timeout: int = 60) -> str:
    body = urllib.parse.urlencode(data).encode()
    req = urllib.request.Request(url, data=body, method="POST",
                                 headers={"User-Agent": "dna_decode/forward msa_fetch (single-query, cached)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode()


def _http_get_text(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "dna_decode/forward msa_fetch"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode()


def _http_get_bytes(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "dna_decode/forward msa_fetch"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _seq_key(sequence: str, mode: str) -> str:
    return hashlib.sha256(f"{mode}:{sequence.strip().upper()}".encode()).hexdigest()[:16]


def _extract_uniref_a3m(tgz_bytes: bytes) -> str:
    """Pull the uniref.a3m member out of the ColabFold result tarball (prefer uniref; fall back to any a3m)."""
    with tarfile.open(fileobj=io.BytesIO(tgz_bytes), mode="r:gz") as tf:
        names = tf.getnames()
        pick = next((n for n in names if n.endswith("uniref.a3m")), None) or \
            next((n for n in names if n.endswith(".a3m")), None)
        if pick is None:
            raise MsaFetchError(f"no .a3m in ColabFold result (members: {names})")
        f = tf.extractfile(pick)
        if f is None:
            raise MsaFetchError(f"could not read {pick}")
        return f.read().decode()


def fetch_msa(sequence: str, cache_dir: str | Path = DEFAULT_CACHE, *, mode: str = "env",
              api: str = COLABFOLD_API, poll_interval: float = 5.0, max_wait: float = 900.0,
              offline_ok: bool = False,
              _post=None, _get_text=None, _get_bytes=None) -> Path:
    """Return a path to an a3m MSA for `sequence`. CACHE-FIRST (one network search per unique sequence).

    On a cache miss: submit to the ColabFold MMseqs2 API, poll the ticket to COMPLETE, download + extract the
    uniref a3m, cache it. `offline_ok=True` raises on a cache miss instead of touching the network. The
    `_post`/`_get_*` hooks are injection points for offline tests (never used in production)."""
    seq = "".join(sequence.split()).upper()
    if not seq or any(c not in "ACDEFGHIKLMNPQRSTVWYXBZUO" for c in seq):
        raise MsaFetchError("sequence must be plain amino-acid letters (no header, no gaps)")
    cache_dir = Path(cache_dir)
    cache = cache_dir / f"{_seq_key(seq, mode)}.a3m"
    if cache.exists():
        return cache
    if offline_ok:
        raise MsaFetchError(f"MSA not cached for this sequence and offline_ok=True ({cache})")

    post = _post or _http_post
    get_text = _get_text or _http_get_text
    get_bytes = _get_bytes or _http_get_bytes

    resp = json.loads(post(f"{api}/ticket/msa", {"q": f">query\n{seq}", "mode": mode}))
    tid, status = resp.get("id"), resp.get("status")
    if not tid:
        raise MsaFetchError(f"no ticket id in submit response: {resp}")
    waited = 0.0
    while status not in ("COMPLETE", "ERROR") and waited < max_wait:
        time.sleep(poll_interval)
        waited += poll_interval
        status = json.loads(get_text(f"{api}/ticket/{tid}")).get("status")
    if status != "COMPLETE":
        raise MsaFetchError(f"MSA job did not complete (status={status}, waited={waited}s)")

    a3m = _extract_uniref_a3m(get_bytes(f"{api}/result/download/{tid}"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache.write_text(a3m, encoding="utf-8")
    return cache


def evolution_table_for_sequence(sequence: str, cache_dir: str | Path = DEFAULT_CACHE, *,
                                 mode: str = "env", offline_ok: bool = False, **fetch_kw) -> dict[str, float]:
    """The deployable novel-protein evolution component: sequence -> MSA (fetched/cached) -> site-independent
    score table ({mutation: score}, higher=preserved), ready for `variant_effect.rank_average_hybrid`.

    NOTE the measured tier (`wiki/forward_modality_hybrid_2026-07-17.md`): site-independent is the FLOOR and
    does not lift ESM2 on its own -- for the real hybrid lift, pass the fetched MSA to a coevolution-grade
    model (GEMME / MSA-Transformer) and feed its table through `msa_evolution.evolution_table_from_scores`."""
    from .msa_evolution import site_independent_table
    msa = fetch_msa(sequence, cache_dir, mode=mode, offline_ok=offline_ok, **fetch_kw)
    _depth = sum(1 for line in Path(msa).read_text(encoding="utf-8").splitlines() if line.startswith(">"))
    return site_independent_table(msa, weights=[1.0] * _depth)
