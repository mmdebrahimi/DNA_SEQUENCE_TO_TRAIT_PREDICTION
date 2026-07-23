"""Build the held-out ESM2+ProSST hybrid manifest — local structure quantization (the Kaggle-prep step).

The full leakage-free hybrid at scale needs, per held-out MaveDB assay (gene NOT in ProteinGym): the assay
sequence + DMS (fetched on Kaggle from MaveDB) AND the ProSST structure tokens (which need torch_geometric +
the ProSST repo — the risky dep). This keeps the PROVEN local Windows quantizer here and ships ONLY the
pre-quantized structure tokens to Kaggle, so the Kaggle notebook runs ESM2 + ProSST(transformer-only) with NO
torch_geometric.

Per gene: fetch AlphaFold PDB (v6) for the UniProt, quantize_structure -> full-protein tokens, then SLICE to
the assay region `tokens[offset:offset+L]` (MaveDB assays are often domain/construct-numbered; the offset is
shipped in the score-set metadata). The DMS mutation keys stay ASSAY-LOCAL (ESM2/ProSST score the assay seq
directly); the offset is used ONLY to align the full-protein structure tokens to the assay region.

Writes wiki/holdout_hybrid_manifest.json = {urn: {uniprot, offset, seq_len, n_tokens, structure_tokens}} —
the compact payload the Kaggle notebook reads (uploaded as a dataset).

  HF_HOME=D:/hf_cache PROSST_REPO=D:/prosst_repo uv run python scripts/build_holdout_hybrid_manifest.py --limit 45

Frozen AMR surface byte-unchanged (READ-only).
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.clinical_gene_landscape_census import enumerate_human_mavedb  # noqa: E402
from scripts.clinical_am_hybrid_auroc import af_pdb, struct_tokens  # noqa: E402
from scripts.mavedb_prospective_holdout import proteingym_gene_symbols  # noqa: E402

API = "https://api.mavedb.org/api/v1"
MANIFEST = Path("wiki/holdout_hybrid_manifest.json")
MAXLEN = 1022       # ESM2 context
MAX_STRUCT = 1400   # cap the FULL-protein AlphaFold structure size -- a domain assay (offset>0) of a GIANT
                    # protein (OBSCN ~7900aa, HECTD1 ~2600aa) makes ProSST quantize wedge; skip those.


MIN_STRUCT_IDENTITY = 0.95   # assay seq must MATCH the AlphaFold structure at the offset, not merely fit by length

_A3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def _pdb_ca_sequence(pdb_path: Path) -> str:
    """Per-residue 1-letter sequence from an AlphaFold PDB's CA atoms (order = residue order)."""
    out = []
    with open(pdb_path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
                out.append(_A3.get(ln[17:20].strip(), "X"))
    return "".join(out)


def _pdb_residue_count(pdb_path: Path) -> int:
    """Count CA atoms (= residues) in an AlphaFold PDB -- a cheap size guard before the slow quantize."""
    return len(_pdb_ca_sequence(pdb_path))


def struct_identity(assay_seq_str: str, pdb_path: Path, offset: int) -> float:
    """Fraction of positions where the assay sequence matches the structure's residues at [offset:offset+L].
    A LENGTH match is not an alignment -- run-3's manifest carried 7 assays that fit by length but encoded a
    different region (PSD95 0.06, OSTF1 0.05, AB42 0.12), which feeds ProSST a mismatched structure. PURE."""
    ps = _pdb_ca_sequence(pdb_path)
    region = ps[offset:offset + len(assay_seq_str)]
    if len(region) != len(assay_seq_str) or not region:
        return 0.0
    return sum(1 for a, b in zip(assay_seq_str, region) if a == b) / len(region)

_AA3TO1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
           "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
           "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
_BASES = "TCAG"
_CODON = dict(zip([a + b + c for a in _BASES for b in _BASES for c in _BASES],
                  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def _translate(dna: str) -> str:
    dna = dna.upper().replace("U", "T")
    out = []
    for i in range(0, len(dna) - 2, 3):
        aa = _CODON.get(dna[i:i + 3], "X")
        if aa == "*":
            break
        out.append(aa)
    return "".join(out)


def _get(url: str) -> str:
    with urllib.request.urlopen(url, timeout=90) as r:
        return r.read().decode("utf-8", "replace")


def assay_seq(urn: str) -> str | None:
    """The assay's protein target sequence (translate DNA), assay-local numbering."""
    rec = json.loads(_get(f"{API}/score-sets/{urn}"))
    tgs = rec.get("targetGenes") or []
    if not tgs:
        return None
    ts = tgs[0].get("targetSequence") or {}
    seq = (ts.get("sequence") or "").strip().upper()
    if not seq:
        return None
    if (ts.get("sequenceType") or "").lower() == "dna":
        seq = _translate(seq)
    return seq or None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=45, help="max held-out genes to quantize")
    a = ap.parse_args()

    pg = {s.upper() for s in proteingym_gene_symbols()}
    landscape = enumerate_human_mavedb(use_cache=True)
    held = [(g, m) for g, m in landscape.items()
            if m.get("uniprot") and g.upper() not in pg and (m.get("n_variants") or 0) >= 500]
    held.sort(key=lambda gm: -(gm[1].get("n_variants") or 0))

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8")) if MANIFEST.exists() else {}
    done_urns = set(manifest)
    built = len([u for u in manifest])
    print(f"held-out human candidates: {len(held)}; already in manifest: {built}", flush=True)

    for g, m in held:
        if len([u for u in manifest]) >= a.limit:
            break
        urn, up, off = m["urn"], m["uniprot"], int(m.get("offset", 0) or 0)
        if urn in done_urns:
            continue
        try:
            seq = assay_seq(urn)
            if not seq or len(seq) > MAXLEN:
                print(f"  {g:12s} SKIP seq (len={len(seq) if seq else 0})", flush=True)
                continue
            pdb = af_pdb(up)
            if pdb is None:
                print(f"  {g:12s} SKIP no AlphaFold ({up})", flush=True)
                continue
            nres = _pdb_residue_count(pdb)
            if nres > MAX_STRUCT:
                print(f"  {g:12s} SKIP giant structure ({up}, {nres} residues > {MAX_STRUCT})", flush=True)
                continue
            # SEQUENCE-IDENTITY GATE (before the slow quantize): a length fit is NOT an alignment.
            ident = struct_identity(seq, pdb, off)
            if ident < MIN_STRUCT_IDENTITY:
                print(f"  {g:12s} SKIP misaligned (identity {ident:.2f} < {MIN_STRUCT_IDENTITY})", flush=True)
                continue
            toks_full = struct_tokens(up)   # cached per uniprot on D:
            if not toks_full:
                print(f"  {g:12s} SKIP quantize failed ({up})", flush=True)
                continue
            sliced = toks_full[off:off + len(seq)]
            if len(sliced) != len(seq):
                print(f"  {g:12s} SKIP align: struct[{off}:{off+len(seq)}]={len(sliced)} != seq {len(seq)} "
                      f"(full={len(toks_full)})", flush=True)
                continue
            manifest[urn] = {"gene": g, "uniprot": up, "offset": off, "seq_len": len(seq),
                             "n_tokens": len(sliced), "structure_tokens": sliced,
                             "struct_identity": round(ident, 3)}
            MANIFEST.write_text(json.dumps(manifest), encoding="utf-8")   # checkpoint every gene (restartable)
            print(f"  {g:12s} OK  len={len(seq)} off={off} tokens={len(sliced)}  [{len(manifest)} total]", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"  {g:12s} ERR {type(e).__name__}: {str(e)[:70]}", flush=True)

    print(f"\nmanifest: {MANIFEST} ({len(manifest)} genes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
