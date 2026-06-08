"""Pure MLST core — profiles-table parse + profile->ST lookup + allele-id parsing. No I/O.

PubMLST profiles TSV: header `ST\t<locus1>\t...\t<locusN>[\tclonal_complex]`, one ST per row. The profile
key is the ordered tuple of per-locus allele numbers. allele FASTA headers are `<locus>_<number>` (e.g.
'adk_4'). Fully unit-testable without network or BLAST.
"""
from __future__ import annotations


def parse_profiles(tsv_text: str) -> tuple[list[str], dict[tuple[int, ...], str], dict[tuple[int, ...], str]]:
    """Parse a PubMLST profiles TSV. Returns (loci_order, profile->ST, profile->clonal_complex).

    loci_order = the locus columns between 'ST' and any 'clonal_complex'/'CC' trailing column, in file order.
    """
    lines = [ln for ln in tsv_text.splitlines() if ln.strip()]
    if not lines:
        return [], {}, {}
    header = lines[0].split("\t")
    # locus columns = everything after ST up to a trailing clonal_complex/CC column
    cc_idx = next((i for i, h in enumerate(header) if h.strip().lower() in ("clonal_complex", "cc")), None)
    end = cc_idx if cc_idx is not None else len(header)
    loci = [h.strip() for h in header[1:end]]
    st_of: dict[tuple[int, ...], str] = {}
    cc_of: dict[tuple[int, ...], str] = {}
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) <= len(loci):
            continue
        try:
            prof = tuple(int(f[1 + i]) for i in range(len(loci)))
        except ValueError:
            continue
        st = f[0].strip()
        st_of[prof] = st
        if cc_idx is not None and cc_idx < len(f):
            cc_of[prof] = f[cc_idx].strip()
    return loci, st_of, cc_of


def allele_number(allele_id: str) -> tuple[str, int] | None:
    """'adk_4' -> ('adk', 4). None if not parseable."""
    if "_" not in allele_id:
        return None
    locus, _, num = allele_id.rpartition("_")
    try:
        return locus, int(num)
    except ValueError:
        return None


def lookup_st(profile: dict[str, int], loci_order: list[str],
              st_of: dict[tuple[int, ...], str]) -> str | None:
    """profile {locus: allele_number} + loci order -> ST, or None (novel/incomplete profile)."""
    if any(profile.get(loc) is None for loc in loci_order):
        return None
    key = tuple(profile[loc] for loc in loci_order)
    return st_of.get(key)
