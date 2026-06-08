"""Pure concordance core — gene-family normalization + set comparison. No I/O, fully unit-testable.

AMRFinder and ResFinder name acquired genes differently at the allele level (blaNDM-19 vs blaNDM-1), so a
fair cross-tool comparison normalizes to the gene FAMILY: lowercase + strip a trailing '-<digits>'
allele-variant suffix. Conservatively KEEPS trailing bare digits that are part of the gene name
(sul1 != sul2, qnrS1) and non-numeric variant tags (aac(6')-Ib-cr) — only the hyphen-number variant is
dropped, which is what distinguishes blaNDM-1/blaNDM-19 (same family) from sul1/sul2 (different genes).
"""
from __future__ import annotations

import re
from collections.abc import Iterable

_VARIANT_SUFFIX = re.compile(r"-\d+$")


def family_normalize(gene: str) -> str:
    """'blaNDM-19' -> 'blandm'; 'blaCTX-M-15' -> 'blactx-m'; 'sul1' -> 'sul1'; 'tet(A)' -> 'tet(a)';
    \"aac(6')-Ib-cr\" -> \"aac(6')-ib-cr\". Lowercase + strip ONE trailing -<digits> variant suffix."""
    return _VARIANT_SUFFIX.sub("", gene.strip().lower())


def compare(amr_genes: Iterable[str], resfinder_genes: Iterable[str]) -> dict:
    """Family-level concordance between two acquired-gene call sets.

    Returns both / amr_only / resfinder_only (sorted family names) + counts + Jaccard agreement.
    `representatives` keeps one original gene name per family per side for human-readable output.
    """
    amr_list = [g for g in amr_genes if g and g.strip()]
    res_list = [g for g in resfinder_genes if g and g.strip()]
    amr_fam = {family_normalize(g) for g in amr_list}
    res_fam = {family_normalize(g) for g in res_list}
    both = amr_fam & res_fam
    amr_only = amr_fam - res_fam
    res_only = res_fam - amr_fam
    union = amr_fam | res_fam
    reps_amr = _reps(amr_list)
    reps_res = _reps(res_list)
    return {
        "both": sorted(both),
        "amr_only": sorted(amr_only),
        "resfinder_only": sorted(res_only),
        "n_both": len(both), "n_amr_only": len(amr_only), "n_resfinder_only": len(res_only),
        "n_amr_total": len(amr_fam), "n_resfinder_total": len(res_fam),
        "agreement": round(len(both) / len(union), 3) if union else None,   # Jaccard over families
        "representatives": {"amr": reps_amr, "resfinder": reps_res},
    }


def _reps(genes: list[str]) -> dict[str, str]:
    """family -> first original gene name seen (for readable reporting)."""
    out: dict[str, str] = {}
    for g in genes:
        out.setdefault(family_normalize(g), g)
    return out


def amr_acquired_genes_from_main(main_tsv: str) -> set[str]:
    """Parse an AMRFinder main.tsv -> set of ACQUIRED AMR gene symbols (excludes POINT mutations).

    AMRFinder main.tsv columns include 'Gene symbol', 'Element type' (AMR/...), 'Method' (EXACTX/BLASTX/
    POINTX/...). Acquired genes = Element type AMR AND Method NOT starting 'POINT'. Tolerant of column order.
    """
    from pathlib import Path
    lines = Path(main_tsv).read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        return set()
    header = lines[0].split("\t")
    def idx(name):
        return header.index(name) if name in header else None
    gi, ei, mi = idx("Gene symbol"), idx("Element type"), idx("Method")
    if gi is None:
        return set()
    out = set()
    for ln in lines[1:]:
        f = ln.split("\t")
        if gi >= len(f):
            continue
        etype = f[ei] if (ei is not None and ei < len(f)) else "AMR"
        method = f[mi] if (mi is not None and mi < len(f)) else ""
        if etype.upper() == "AMR" and not method.upper().startswith("POINT"):
            sym = f[gi].strip()
            if sym:
                out.add(sym)
    return out
