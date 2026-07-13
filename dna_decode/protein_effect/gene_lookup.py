"""Gene -> protein-sequence lookup for the rung-2 mutation-effect predictor (2026-07-13).

Convenience layer over the predictor: resolve an E. coli gene NAME to its canonical protein sequence so a
caller can say `--gene gyrA` instead of pasting 875 residues. Fetches from UniProt (free, no API key) and
caches locally.

The design review's deferral reason was STRAIN / NUMBERING ambiguity — handled head-on:
  * PINNED reference: E. coli K-12 (UniProt organism_id 83333). One strain => deterministic 1-based residue
    numbering over the canonical sequence. (Verified: gyrA -> P0AES4, residue 83 = S, matching the literature
    QRDR gyrA S83.)
  * Prefer REVIEWED (Swiss-Prot) canonical; fall back to unreviewed with a flag.
  * FAIL-CLOSED on ambiguity: >1 match => GeneLookupAmbiguous (caller must pass an accession), never a guess.
  * Offline-safe: a cache hit needs no network; a network failure with no cache raises GeneLookupUnavailable
    (never fabricates a sequence).
  * Every result stamps provenance: the UniProt accession + "canonical numbering may differ from PDB /
    mature-protein conventions" so a mutation position is interpreted against the right frame.
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

ECOLI_K12 = {"organism_id": "83333", "label": "Escherichia coli K-12"}
_BASE = "https://rest.uniprot.org/uniprotkb/search"
_FIELDS = "accession,gene_names,organism_name,protein_name,sequence,reviewed"
_DEFAULT_CACHE = Path(__file__).resolve().parent.parent.parent / "data" / "processed" / "gene_lookup_cache"


class GeneLookupError(RuntimeError):
    """Base for gene-lookup failures."""


class GeneLookupAmbiguous(GeneLookupError):
    """>1 UniProt match for the gene — the caller must disambiguate by accession."""


class GeneLookupNotFound(GeneLookupError):
    """No UniProt match for the gene in the pinned organism."""


class GeneLookupUnavailable(GeneLookupError):
    """Network unreachable and no cached result — refuse to fabricate a sequence."""


def _cache_path(gene: str, organism_id: str, cache_dir: Path) -> Path:
    return Path(cache_dir) / f"{organism_id}_{gene.lower()}.json"


def _query(gene: str, organism_id: str, reviewed: bool) -> list[dict]:
    q = f"(gene_exact:{gene}) AND (organism_id:{organism_id})" + (" AND (reviewed:true)" if reviewed else "")
    url = _BASE + "?" + urllib.parse.urlencode(
        {"query": q, "fields": _FIELDS, "format": "json", "size": "10"})
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.load(r).get("results", [])


def _record(gene: str, organism_id: str, hits: list[dict], reviewed: bool) -> dict:
    if not hits:
        raise GeneLookupNotFound(f"no UniProt match for gene {gene!r} in organism {organism_id}")
    if len(hits) > 1:
        accs = [h["primaryAccession"] for h in hits]
        raise GeneLookupAmbiguous(f"{len(hits)} matches for {gene!r} ({accs}); pass an accession to disambiguate")
    h = hits[0]
    seq = h["sequence"]["value"]
    return {
        "gene": gene, "accession": h["primaryAccession"], "sequence": seq, "length": len(seq),
        "organism": h["organism"]["scientificName"],
        "reviewed": h.get("entryType", "").startswith("UniProtKB reviewed") or reviewed,
        "protein_name": (h.get("proteinDescription", {}).get("recommendedName", {})
                         .get("fullName", {}).get("value")),
        "provenance": {
            "source": "UniProt", "organism_id": organism_id,
            "numbering": ("1-based over the UniProt CANONICAL sequence for the pinned strain; may differ from "
                          "PDB / mature-protein (signal-peptide-cleaved) numbering — interpret the mutation "
                          "position against this frame"),
            "reviewed_preferred": reviewed,
        },
    }


def fetch_protein_sequence(gene: str, organism_id: str = ECOLI_K12["organism_id"],
                           cache_dir: str | Path = _DEFAULT_CACHE, force: bool = False) -> dict:
    """Resolve `gene` -> canonical protein record for the pinned organism. Cached; offline-safe; fail-closed.

    Raises GeneLookupAmbiguous / GeneLookupNotFound / GeneLookupUnavailable (never fabricates)."""
    gene = gene.strip()
    cp = _cache_path(gene, organism_id, cache_dir)
    if cp.exists() and not force:
        return json.loads(cp.read_text(encoding="utf-8"))
    try:
        hits = _query(gene, organism_id, reviewed=True)
        used_reviewed = True
        if not hits:                                   # fall back to unreviewed (TrEMBL) with a flag
            hits = _query(gene, organism_id, reviewed=False)
            used_reviewed = False
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise GeneLookupUnavailable(f"UniProt unreachable + no cache for {gene!r}: {type(e).__name__}") from e
    rec = _record(gene, organism_id, hits, reviewed=used_reviewed)
    cp.parent.mkdir(parents=True, exist_ok=True)
    cp.write_text(json.dumps(rec, indent=1), encoding="utf-8")
    return rec
