"""Determinant overlay — DeterminantHit + the hard join-quality gate (Step 4, catch C3).

The determinant->feature join is the INTEGRITY CRUX (the documented gene-symbol
0%-overlap trap). `amr_rules.cipro_determinants_from_main` / `call_resistance`
return a determinant's SYMBOL/CLASS, not its coordinates — joining a feature by
`gene_symbol` is the trap (locus tags never equal gene symbols; the gene_symbol
column is populated for only ~11% of CDS). So this module parses the RAW
AMRFinder `main.tsv` into a `DeterminantHit` that RETAINS the protein-id / contig
/ coordinates, and joins by an explicit hierarchy:

    protein-id (exact)  -> join_confidence = "protein_id"   (HIGH)
    coordinate overlap  -> join_confidence = "coord"        (HIGH)
    gene-symbol match   -> join_confidence = "symbol_fallback" (LOW — the trap)
    (none)              -> unjoined

**Symbol-fallback joins are VISIBLE but do NOT earn the `determinant-phenotype`
primary tier and do NOT count for G1.** The spike NO-GOs if ALL of a genome's
determinant joins are symbol-fallback.

Contig reconciliation: AMRFinder runs `-n` on the ORIGINAL FASTA (NCBI contig
headers) while Bakta RENAMES contigs (`contig_1`…). A coordinate join therefore
needs a `contig_name_map` (AMRFinder contig -> Bakta seqid). `build_contig_name_map`
reconciles by matching unique contig LENGTHS (Bakta's `##sequence-region`
lengths vs the FASTA contig lengths). When names already match, no map is needed.

This module is READ-ONLY w.r.t. the frozen AMR surface.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import csv

# Column-name candidates (AMRFinder versions differ). Resolved per-file from the
# actual header; mirrors scripts/genome_map_tool_surface column detection.
_PROTEIN_ID_COLS = ("Protein identifier", "Protein id")
_CONTIG_COLS = ("Contig id", "Contig")
_START_COLS = ("Start",)
_STOP_COLS = ("Stop", "End")
_SYMBOL_COLS = ("Element symbol", "Gene symbol")
_NAME_COLS = ("Element name", "Sequence name")
_CLASS_COLS = ("Class",)
_SUBCLASS_COLS = ("Subclass",)
_METHOD_COLS = ("Method",)

# join_confidence values that earn the determinant-phenotype tier + G1 credit.
HIGH_CONFIDENCE_JOINS = frozenset({"protein_id", "coord"})
SYMBOL_FALLBACK = "symbol_fallback"

# AMRFinder's sentinel for "no value" in a coordinate/protein column.
_NA_VALUES = {"", "NA", "na", "N/A"}

# Bakta feature types that are NOT genes and must be EXCLUDED from the coordinate
# join — chiefly `region`, the whole-contig (1..length) metadata feature that
# would otherwise win every coordinate overlap and swallow all determinants.
_NON_GENE_COORD_TYPES = frozenset({"region", "oriC", "oriT", "regulatory_region"})


@dataclass
class DeterminantHit:
    """One AMRFinder main.tsv determinant row, retaining join coordinates.

    `protein_id` / `contig` / `start` / `stop` are None when the column is
    absent or holds AMRFinder's NA sentinel (the signal that drives a fallback).
    """
    symbol: str
    name: str
    cls: str
    subclass: str
    method: str
    protein_id: str | None
    contig: str | None
    start: int | None
    stop: int | None
    raw: dict = field(default_factory=dict)


@dataclass
class JoinedHit:
    """A DeterminantHit joined (or not) to a feature row index."""
    hit: DeterminantHit
    feature_index: int | None
    join_confidence: str | None  # "protein_id" | "coord" | "symbol_fallback" | None
    join_key: str = ""

    @property
    def is_high_confidence(self) -> bool:
        return self.join_confidence in HIGH_CONFIDENCE_JOINS


def _first_present(headers: list[str], candidates: tuple[str, ...]) -> str | None:
    hset = set(headers)
    for c in candidates:
        if c in hset:
            return c
    return None


def _to_int(val: str | None) -> int | None:
    if val is None or str(val).strip() in _NA_VALUES:
        return None
    try:
        return int(str(val).strip())
    except ValueError:
        return None


def _clean(val: str | None) -> str | None:
    if val is None:
        return None
    s = str(val).strip()
    return None if s in _NA_VALUES else s


def parse_determinant_hits(main_tsv: Path | str) -> list[DeterminantHit]:
    """Parse a raw AMRFinder main.tsv into DeterminantHits (all rows, coords retained).

    Resolves the symbol/coordinate column names from the actual header (version
    tolerant). A missing/empty file returns []. Coordinate/protein fields that
    are absent or NA become None (the fallback signal).
    """
    p = Path(main_tsv)
    if not p.exists() or p.stat().st_size == 0:
        return []
    with p.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames or []
        protein_c = _first_present(headers, _PROTEIN_ID_COLS)
        contig_c = _first_present(headers, _CONTIG_COLS)
        start_c = _first_present(headers, _START_COLS)
        stop_c = _first_present(headers, _STOP_COLS)
        symbol_c = _first_present(headers, _SYMBOL_COLS)
        name_c = _first_present(headers, _NAME_COLS)
        class_c = _first_present(headers, _CLASS_COLS)
        subclass_c = _first_present(headers, _SUBCLASS_COLS)
        method_c = _first_present(headers, _METHOD_COLS)
        out: list[DeterminantHit] = []
        for row in reader:
            out.append(
                DeterminantHit(
                    symbol=(row.get(symbol_c, "") if symbol_c else "").strip(),
                    name=(row.get(name_c, "") if name_c else "").strip(),
                    cls=(row.get(class_c, "") if class_c else "").strip(),
                    subclass=(row.get(subclass_c, "") if subclass_c else "").strip(),
                    method=(row.get(method_c, "") if method_c else "").strip(),
                    protein_id=_clean(row.get(protein_c) if protein_c else None),
                    contig=_clean(row.get(contig_c) if contig_c else None),
                    start=_to_int(row.get(start_c) if start_c else None),
                    stop=_to_int(row.get(stop_c) if stop_c else None),
                    raw=dict(row),
                )
            )
    return out


def build_contig_name_map(
    fasta_contig_lengths: dict[str, int],
    bakta_contig_lengths: dict[str, int],
) -> dict[str, str]:
    """Reconcile AMRFinder contig names (FASTA headers) -> Bakta seqids by UNIQUE length.

    Bakta renames + may reorder contigs, but each contig's LENGTH is preserved.
    Matches only lengths that are unique on BOTH sides (ambiguous lengths are left
    unmapped -> those hits fall to symbol-fallback, surfaced honestly).
    """
    # lengths unique on each side
    def _unique_by_len(d: dict[str, int]) -> dict[int, str]:
        by_len: dict[int, list[str]] = {}
        for name, ln in d.items():
            by_len.setdefault(ln, []).append(name)
        return {ln: names[0] for ln, names in by_len.items() if len(names) == 1}

    fasta_u = _unique_by_len(fasta_contig_lengths)
    bakta_u = _unique_by_len(bakta_contig_lengths)
    out: dict[str, str] = {}
    for ln, fasta_name in fasta_u.items():
        if ln in bakta_u:
            out[fasta_name] = bakta_u[ln]
    return out


def _coord_overlap(a_start: int, a_stop: int, b_start: int, b_stop: int) -> int:
    """Length of the overlap between two inclusive ranges (0 if disjoint)."""
    lo = max(min(a_start, a_stop), min(b_start, b_stop))
    hi = min(max(a_start, a_stop), max(b_start, b_stop))
    return max(0, hi - lo + 1)


def join_hits(
    features,
    hits: list[DeterminantHit],
    *,
    contig_name_map: dict[str, str] | None = None,
) -> tuple[list[JoinedHit], dict]:
    """Join each DeterminantHit to a feature row by protein-id > coord > symbol-fallback.

    `features` is the AnnotationTable (a DataFrame with seqid/start/end/gene_id/
    locus_tag/gene_symbol). `contig_name_map` maps an AMRFinder contig name to the
    Bakta seqid (from build_contig_name_map); when None, contig names are compared
    directly.

    Returns (joined_hits, counts) where counts =
    {n_main_rows, n_high_confidence_join, n_symbol_fallback, n_unjoined}.
    """
    # Pre-index features for the three join modes.
    locus_to_idx: dict[str, int] = {}
    symbol_to_idx: dict[str, int] = {}
    feats: list[dict] = []
    for i, (_, r) in enumerate(features.iterrows()):
        gid = str(r.get("gene_id") or "")
        ltag = str(r.get("locus_tag") or "")
        sym = str(r.get("gene_symbol") or "")
        seqid = str(r.get("seqid") or "")
        ftype = str(r.get("type") or "")
        start = int(r.get("start") or 0)
        end = int(r.get("end") or 0)
        feats.append({"seqid": seqid, "start": start, "end": end, "type": ftype,
                      "gene_id": gid, "locus_tag": ltag, "gene_symbol": sym})
        for key in (gid, ltag):
            if key:
                locus_to_idx.setdefault(key, i)
        if sym and sym not in symbol_to_idx:
            symbol_to_idx[sym] = i

    cmap = contig_name_map or {}
    joined: list[JoinedHit] = []
    n_high = n_symbol = n_unjoined = 0

    for h in hits:
        idx: int | None = None
        conf: str | None = None
        key = ""

        # 1. protein-id exact join (against gene_id / locus_tag).
        if h.protein_id and h.protein_id in locus_to_idx:
            idx = locus_to_idx[h.protein_id]
            conf = "protein_id"
            key = h.protein_id

        # 2. coordinate-overlap join (contig reconciled + ranges overlap).
        #    Skips whole-contig metadata features (Bakta emits a `region` row per
        #    contig spanning 1..len that would otherwise win every overlap); on a
        #    tie prefers the SMALLEST (most specific) gene, so a determinant maps
        #    to its CDS, not an enclosing feature.
        if idx is None and h.contig and h.start is not None and h.stop is not None:
            bakta_seqid = cmap.get(h.contig, h.contig)
            best_i, best_ov, best_len = None, 0, None
            for i, fe in enumerate(feats):
                if fe["seqid"] != bakta_seqid or fe["type"] in _NON_GENE_COORD_TYPES:
                    continue
                ov = _coord_overlap(h.start, h.stop, fe["start"], fe["end"])
                if ov <= 0:
                    continue
                flen = abs(fe["end"] - fe["start"]) + 1
                if ov > best_ov or (ov == best_ov and (best_len is None or flen < best_len)):
                    best_ov, best_i, best_len = ov, i, flen
            if best_i is not None:
                idx = best_i
                conf = "coord"
                key = f"{bakta_seqid}:{h.start}-{h.stop}"

        # 3. gene-symbol fallback (the trap — visible, low-confidence).
        if idx is None and h.symbol:
            # AMRFinder symbols can carry a point-mutation suffix (gyrA_S83L);
            # match the bare gene token too.
            bare = h.symbol.split("_", 1)[0]
            if h.symbol in symbol_to_idx:
                idx, conf, key = symbol_to_idx[h.symbol], SYMBOL_FALLBACK, h.symbol
            elif bare in symbol_to_idx:
                idx, conf, key = symbol_to_idx[bare], SYMBOL_FALLBACK, bare

        jh = JoinedHit(hit=h, feature_index=idx, join_confidence=conf, join_key=key)
        joined.append(jh)
        if jh.is_high_confidence:
            n_high += 1
        elif conf == SYMBOL_FALLBACK:
            n_symbol += 1
        else:
            n_unjoined += 1

    counts = {
        "n_main_rows": len(hits),
        "n_high_confidence_join": n_high,
        "n_symbol_fallback": n_symbol,
        "n_unjoined": n_unjoined,
    }
    return joined, counts


def all_joins_symbol_fallback(counts: dict) -> bool:
    """True iff there are determinant rows AND none cleared a high-confidence join.

    The gene-symbol-trap guard: a genome whose determinant joins are ALL
    symbol-fallback (or unjoined) cannot honestly surface determinant-phenotype
    features -> the spike NO-GOs on it.
    """
    return counts.get("n_main_rows", 0) > 0 and counts.get("n_high_confidence_join", 0) == 0


def determinant_phenotype_field(joined_hit: JoinedHit, drug: str, verdict: dict | None) -> dict | None:
    """Build the phenotype annotation for a HIGH-confidence determinant feature.

    Returns None for a non-high-confidence join (the phenotype wall — symbol-
    fallback never carries a phenotype). When the genome-level determinant cell
    ABSTAINs/SUSPENDs (`verdict["prediction"] in {ABSTAIN, INDETERMINATE}`), the
    feature shows ABSTAIN, never a forced call (AC8).
    """
    if not joined_hit.is_high_confidence:
        return None
    h = joined_hit.hit
    pred = (verdict or {}).get("prediction")
    abstain = pred in {"ABSTAIN", "INDETERMINATE", "SUSPEND"}
    return {
        "drug": drug,
        "determinant_symbol": h.symbol,
        "amrfinder_class": h.cls,
        "amrfinder_subclass": h.subclass,
        "method": h.method,
        "join_confidence": joined_hit.join_confidence,
        "phenotype": "ABSTAIN" if abstain else (verdict or {}).get("prediction"),
        "provenance": (verdict or {}).get("rule", "amrfinder_curated_determinant"),
        "abstain": abstain,
    }
