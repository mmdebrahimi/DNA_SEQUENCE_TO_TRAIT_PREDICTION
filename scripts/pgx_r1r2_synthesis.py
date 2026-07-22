"""R1 x R2 pharmacogene synthesis — does the free molecular DMS (R2) corroborate the star-allele clinical
function (R1) for the SAME gene?

The project has BOTH regimes for CYP2C9: an R1 deterministic catalog cell (star-allele -> CPIC metabolizer
phenotype, GeT-RM 1.0 concordance, `dna_decode/pgx/cyp2c9_catalog.py`) AND today's R2 molecular substrate
(MaveDB CYP2C9 deep-mutational-scanning). This joins them: for each star-allele-defining MISSENSE variant with
a known CPIC function, look up its MEASURED DMS score and ask whether the CPIC-reduced-function alleles land in
the damaging tail.

FINDING (2026-07-22): they do NOT corroborate for the clinically-important alleles — because ALL free CYP2C9
MaveDB DMS is an ABUNDANCE assay (VAMP-seq, per the score-set titles), while CYP2C9 *2 (R144C) and *3 (I359L)
are CATALYTIC-ACTIVITY defects — stable proteins that are simply less active. An abundance assay cannot see an
activity defect. Only *11 (R335W, a destabilizing variant) lands in the abundance-damaging tail (1/5).

LESSON: R1 x R2 corroboration is PROPERTY-SPECIFIC. The R2 molecular assay validates the R1 clinical call ONLY
when it measures the SAME molecular property the phenotype depends on (activity vs abundance vs binding). The
free CYP2C9 DMS measures the wrong property for the metabolizer phenotype. This is the g2p regime boundary
(fitness-alignment) at finer grain, and it echoes the HIV-DRM finding (conservative substitutions at conserved
sites evade abundance/likelihood scorers). Honest negative — kept, not tuned.

CYP2C19 is OUT of this synthesis: its core alleles are non-missense (*2 splice, *3 stop, *17 promoter), so
neither DMS nor ESM sees them.

  uv run python scripts/pgx_r1r2_synthesis.py     # real MaveDB fetch; writes wiki/pgx_r1r2_synthesis_<date>.json
"""
from __future__ import annotations

import bisect
import csv
import io
import json
import statistics as st
import sys
import urllib.request
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

API = "https://api.mavedb.org/api/v1"
# CYP2C9 abundance DMS assays (VAMP-seq); the metabolizer phenotype depends on ACTIVITY, not abundance.
CYP2C9_ABUNDANCE_URNS = ["urn:mavedb:00000095-b-1", "urn:mavedb:00000095-a-1"]
# star-allele-defining MISSENSE variants + CPIC function (from dna_decode/pgx/cyp2c9_catalog.py + CPIC).
CYP2C9_STAR_MISSENSE = {
    "R144C": ("*2", "decreased (activity 0.5)"),
    "I359L": ("*3", "no function (activity 0.0)"),
    "D360E": ("*5", "no function"),
    "R150H": ("*8", "decreased"),
    "H251R": ("*9", "decreased"),
    "R335W": ("*11", "decreased/no function"),
}
_A3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}


def hgvs_pro_short(h: str) -> str | None:
    """p.Arg144Cys -> 'R144C'; p.Trp212Ter -> 'W212*'; None for non-single-substitution. PURE."""
    if not h or not h.startswith("p."):
        return None
    b = h[2:].strip()
    if b in ("=", "?") or "[" in b or "fs" in b or "del" in b or "ins" in b:
        return None
    if b[:3] not in _A3:
        return None
    i = 3
    while i < len(b) and b[i].isdigit():
        i += 1
    if i == 3 or i >= len(b):
        return None
    alt = b[i:]
    if alt == "Ter":
        return _A3[b[:3]] + b[3:i] + "*"
    return _A3[b[:3]] + b[3:i] + _A3[alt] if alt in _A3 else None


def _get(url: str, timeout: int = 90) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def score_synthesis(urn: str) -> dict:
    """Fetch one CYP2C9 abundance DMS; join the star-allele missense variants; return the concordance record."""
    rec = json.loads(_get(f"{API}/score-sets/{urn}"))
    title = rec.get("title", "")
    rows: dict[str, float] = {}
    for r in csv.DictReader(io.StringIO(_get(f"{API}/score-sets/{urn}/scores"))):
        sh = hgvs_pro_short(r.get("hgvs_pro", ""))
        raw = r.get("score", "")
        if sh and raw not in ("", "NA", None):
            try:
                rows[sh] = float(raw)
            except ValueError:
                continue
    vals = sorted(rows.values())
    nonsense = [v for k, v in rows.items() if k.endswith("*")]
    lower_is_damaging = (st.median(nonsense) < st.median(list(rows.values()))) if nonsense else True
    per_allele, hits, scored = [], 0, 0
    for v, (star, fn) in CYP2C9_STAR_MISSENSE.items():
        if v not in rows:
            per_allele.append({"star": star, "variant": v, "cpic_function": fn, "in_assay": False})
            continue
        pctile = round(bisect.bisect_left(vals, rows[v]) / len(vals), 3)
        damaging = (pctile <= 0.25) if lower_is_damaging else (pctile >= 0.75)
        hits += damaging
        scored += 1
        per_allele.append({"star": star, "variant": v, "cpic_function": fn, "in_assay": True,
                           "dms_score": round(rows[v], 3), "percentile": pctile, "in_damaging_tail": damaging})
    return {"urn": urn, "title": title, "readout": "abundance (VAMP-seq)", "n_missense": len(rows),
            "lower_is_damaging": lower_is_damaging, "n_star_scored": scored, "n_star_damaging": hits,
            "per_allele": per_allele}


def main() -> int:
    results = [score_synthesis(u) for u in CYP2C9_ABUNDANCE_URNS]
    primary = results[0]
    print(f"CYP2C9 R1xR2 synthesis (primary {primary['urn']}, {primary['readout']}):")
    for a in primary["per_allele"]:
        if a["in_assay"]:
            print(f"  {a['star']:4s} {a['variant']:6s} CPIC={a['cpic_function']:26s} "
                  f"DMS pctile={a['percentile']:.3f} damaging={'YES' if a['in_damaging_tail'] else 'no'}")
        else:
            print(f"  {a['star']:4s} {a['variant']:6s} CPIC={a['cpic_function']:26s} (not in assay)")
    print(f"  -> {primary['n_star_damaging']}/{primary['n_star_scored']} CPIC-reduced alleles in the "
          f"abundance-damaging tail")
    art = {"_schema": "pgx-r1r2-synthesis-v1", "date": _date.today().isoformat(), "gene": "CYP2C9",
           "finding": "PROPERTY-MISMATCH: free CYP2C9 DMS is ABUNDANCE; *2/*3 are activity defects (stable) -> "
                      "abundance R2 does NOT corroborate the activity-based R1 clinical function. Only *11 "
                      "(destabilizing) concords. R1xR2 corroboration is property-specific.",
           "cyp2c19_excluded": "core alleles non-missense (*2 splice / *3 stop / *17 promoter)",
           "assays": results}
    out = Path(f"wiki/pgx_r1r2_synthesis_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"artifact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
