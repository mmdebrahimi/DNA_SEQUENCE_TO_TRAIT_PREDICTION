"""Klebsiella phage depolymerase -> capsule (KL-type) caller — the cross-organism cell, FETCH-ONLY.

Predicts which Klebsiella capsule type(s) a phage depolymerase (its enzymatic tail-spike domain) targets, by
protein k-mer nearest-neighbour transfer against a REFERENCE the user builds locally from the DpoTropiSearch
dataset (Concha-Eloko et al., Nat Commun 2025; Zenodo 10.5281/zenodo.14065540). Depolymerases are promiscuous,
so the caller returns a RANKED top-K KL-type shortlist (the phage-therapy-useful form).

LICENSE / FETCH-ONLY (load-bearing): this package ships ONLY this caller code (mine) + a fetch helper — it
BUNDLES NONE of the DpoTropiSearch data. That data carries a Zenodo CC-BY-4.0 record AND a repo "Decapsulate
Non-Commercial License v1.1"; the conflict is the user's to resolve for THEIR use. By never redistributing it
(the user fetches + accepts the license via `scripts/fetch_dpotropisearch.py`), the shipped package is
unaffected either way. Same pattern as the AMR cells (AMRFinder DB) + phage genome-mode (genome fetch).

Validation (in-distribution, prophage-LCA labels — NOT independent wet-lab): clonality-corrected LOO top-1
~0.45 / top-5 ~0.60 (+0.49 over null), the paradigm GENERALIZES cross-organism (wiki/klebsiella_crossorganism*).

Pure functions here are offline-testable; the reference load + real transfer need the fetched reference.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

from dna_decode.phage.rbp_caller import kmer_similarity, protein_kmers

# Default local reference locations searched (in order); NONE ship with the package (fetch-only).
# The fetch helper writes here; $DPO_KLEB_REFERENCE overrides.
_DEFAULT_REF_CANDIDATES = (
    "data/kleb_ref/depolymerase_kltype_reference.faa",           # user-built, gitignored
    "D:/dna_decode_cache/kleb/depolymerase_kltype_reference.faa",  # this host's build location
)


def resolve_reference(path: str | None = None) -> Path | None:
    """$DPO_KLEB_REFERENCE -> explicit arg -> known local build locations. None if unbuilt (fetch-only)."""
    env = os.environ.get("DPO_KLEB_REFERENCE")
    for cand in ([path, env] if (path or env) else []) + list(_DEFAULT_REF_CANDIDATES):
        if cand and Path(cand).exists():
            return Path(cand)
    return None


def load_reference(path: str | Path, k: int = 4) -> tuple[dict[str, frozenset[str]], dict[str, str]]:
    """Load a `>{KL_type}|{id}` depolymerase-domain reference. Returns (kmers{id->set}, kltype{id->KL})."""
    kmers: dict[str, frozenset[str]] = {}
    kltype: dict[str, str] = {}
    label = kl = None
    parts: list[str] = []

    def _flush():
        if label and parts:
            kmers[label] = protein_kmers("".join(parts), k=k)
            kltype[label] = kl
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        idx = 0
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                _flush()
                kl = line[1:].partition("|")[0].strip()
                label = f"{kl}|{idx}"; idx += 1
                parts = []
            else:
                parts.append(line.strip())
    _flush()
    return kmers, kltype


@dataclass(frozen=True)
class KLCall:
    status: str                       # "CALLED" | "INDETERMINATE"
    ranked_kltypes: tuple[str, ...]   # top-K predicted capsule types, most-similar first
    top_similarity: float | None
    reason: str = ""
    method: str = "depolymerase_domain_kmer_transfer_v0"


def call_kltype(depolymerase_protein: str, ref_kmers: dict[str, frozenset[str]],
                ref_kltype: dict[str, str], *, k: int = 4, top_k: int = 5,
                min_similarity: float = 0.05) -> KLCall:
    """Rank capsule KL-types by nearest depolymerase-domain k-mer neighbours. INDETERMINATE (abstain, never a
    fabricated call) when no reference domain clears `min_similarity`."""
    q = protein_kmers(depolymerase_protein, k=k)
    sims = []
    for label, ks in ref_kmers.items():
        s = kmer_similarity(q, ks)
        if s >= min_similarity:
            sims.append((s, ref_kltype[label]))
    if not sims:
        return KLCall("INDETERMINATE", (), None,
                      reason=f"no reference depolymerase cleared min_similarity={min_similarity}")
    sims.sort(key=lambda x: x[0], reverse=True)
    ranked: list[str] = []
    for _, kl in sims:
        if kl not in ranked:
            ranked.append(kl)
        if len(ranked) >= top_k:
            break
    return KLCall("CALLED", tuple(ranked), sims[0][0])


@dataclass
class KLLoo:
    n: int = 0
    called: int = 0
    top1: int = 0
    topk: int = 0
    per_kl: dict[str, list[int]] = field(default_factory=dict)


def leave_one_out(ref_kmers, ref_kltype, *, top_k: int = 5, min_similarity: float = 0.05) -> KLLoo:
    """LOO over the reference (each domain predicted from the others). Reused by the validation harness."""
    ids = list(ref_kmers)
    res = KLLoo()
    for i in ids:
        res.n += 1
        sims = []
        for j in ids:
            if j == i:
                continue
            s = kmer_similarity(ref_kmers[i], ref_kmers[j])
            if s >= min_similarity:
                sims.append((s, ref_kltype[j]))
        if not sims:
            continue
        res.called += 1
        sims.sort(key=lambda x: x[0], reverse=True)
        ranked = []
        for _, kl in sims:
            if kl not in ranked:
                ranked.append(kl)
        true = ref_kltype[i]
        b = res.per_kl.setdefault(true, [0, 0]); b[1] += 1
        if ranked and ranked[0] == true:
            res.top1 += 1; b[0] += 1
        if true in ranked[:top_k]:
            res.topk += 1
    return res
