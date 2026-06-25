"""Salmonella enterica serovar caller — antigen allele DB via the shared blastn engine.

The Kauffmann-White-Le Minor scheme defines a serovar by an ANTIGENIC FORMULA: O-antigen group : H1
(phase-1 flagellin = fliC) : H2 (phase-2 flagellin = fljB). E.g. Typhimurium = `1,4,[5],12 : i : 1,2`.
SeqSero2 determines this from WGS by detecting the O-group + fliC + fljB allele sequences, then looks the
formula up in the White-Kauffmann table. This is the SAME shape as the E. coli serotype caller (best O
antigen + best H antigen) with a serovar-name lookup on top (the way ktype adds a wzi->KL map).

DB layout (a directory):
  * `salmonella_antigens.fasta` -- allele sequences, headers `<axis>__<antigen>__<id>` where
    axis in {O, H1, H2}, e.g. `O__9__01`, `H1__d__01`, `H2__1,2__01`. (The shared engine keys on the
    full header; this runner parses axis + antigen from it.)
  * `serovar_table.tsv` -- `O<tab>H1<tab>H2<tab>Serovar`, the White-Kauffmann-Le Minor formula table.

HONESTY (load-bearing): faithful to the SeqSero2 / Kauffmann-White method (blastn over the antigen DB +
formula lookup); NOT an independent baseline. Many formulas are shared by >1 serovar / are phase-incomplete
-> the call is the FORMULA + a serovar IFF the formula resolves uniquely (else serovar=None, like O?/H?).
Salmonella enterica only. Offline-safe via the shared engine (missing blastn/DB -> status 'unavailable').
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.typing.blast_caller import call_alleles

# Salmonella antigen alleles are near-exact; SeqSero2 uses high-identity matching.
SEROVAR_IDENTITY_THRESHOLD = 90.0
SEROVAR_COVERAGE_THRESHOLD = 80.0


def parse_axis_antigen(allele_id: str) -> tuple[str, str] | None:
    """`O__9__01` -> ('O', '9'); `H1__d__01` -> ('H1', 'd'); `H2__1,2__01` -> ('H2', '1,2'). None if malformed."""
    parts = allele_id.split("__")
    if len(parts) < 2 or parts[0] not in ("O", "H1", "H2"):
        return None
    return parts[0], parts[1]


def load_serovar_table(tsv: str | Path) -> dict[tuple[str, str, str], str]:
    """Parse `O<tab>H1<tab>H2<tab>Serovar` -> {(O, H1, H2): serovar}. Header row (`O`/`Serovar`) skipped."""
    table: dict[tuple[str, str, str], str] = {}
    for line in Path(tsv).read_text(encoding="utf-8").splitlines():
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 4 or parts[0].lower() == "o" or parts[3].lower() == "serovar":
            continue
        o, h1, h2, serovar = parts[0], parts[1], parts[2], parts[3]
        if o and serovar:
            table[(o, h1, h2)] = serovar
    return table


def _best_per_axis(per_allele: dict) -> dict[str, dict]:
    """For each axis (O/H1/H2), the best-coverage CALLED antigen."""
    axis_best: dict[str, dict] = {}
    for allele_id, hit in per_allele.items():
        if not hit["called"]:
            continue
        pa = parse_axis_antigen(allele_id)
        if pa is None:
            continue
        axis, antigen = pa
        cov = hit["percent_coverage"]
        cur = axis_best.get(axis)
        if cur is None or cov > cur["percent_coverage"]:
            axis_best[axis] = {"axis": axis, "antigen": antigen, "best_allele": allele_id,
                               "percent_identity": hit["percent_identity"], "percent_coverage": cov}
    return axis_best


def call_serovar(fasta: str | Path, db_dir: str | Path, *,
                 identity_threshold: float = SEROVAR_IDENTITY_THRESHOLD,
                 coverage_threshold: float = SEROVAR_COVERAGE_THRESHOLD,
                 blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn the Salmonella antigen DB vs `fasta`; assemble the O:H1:H2 formula + look up the serovar."""
    db_dir = Path(db_dir)
    antigens_fasta, serovar_tsv = db_dir / "salmonella_antigens.fasta", db_dir / "serovar_table.tsv"
    if not antigens_fasta.exists() or not serovar_tsv.exists():
        return {"status": "unavailable", "serovar": None, "antigenic_formula": None,
                "reason": f"Salmonella antigen DB not found in {db_dir} "
                          "(need salmonella_antigens.fasta + serovar_table.tsv)"}
    res = call_alleles(fasta, antigens_fasta, identity_threshold=identity_threshold,
                       coverage_threshold=coverage_threshold, blastn_bin=blastn_bin, timeout=timeout)
    if res["status"] != "ok":
        return {"status": "unavailable", "tool": res.get("tool"), "serovar": None,
                "antigenic_formula": None, "reason": res.get("reason")}

    axis_best = _best_per_axis(res["per_allele"])
    o = axis_best.get("O", {}).get("antigen")
    h1 = axis_best.get("H1", {}).get("antigen")
    h2 = axis_best.get("H2", {}).get("antigen")
    formula = f"{o or 'O?'}:{h1 or 'H?'}:{h2 or '-'}" if (o or h1) else None

    serovar = None
    if o and h1:
        table = load_serovar_table(serovar_tsv)
        serovar = table.get((o, h1, h2 or "-")) or table.get((o, h1, h2 or ""))
        # phase-incomplete fallback: match O+H1 ignoring H2 IFF that resolves uniquely
        if serovar is None:
            cands = {sv for (to, th1, _), sv in table.items() if to == o and th1 == h1}
            serovar = next(iter(cands)) if len(cands) == 1 else None

    base = {"status": "ok", "tool": "blastn", "method": "seqsero2_style_antigen_blastn_v0",
            "parameters": {"identity_threshold": identity_threshold, "coverage_threshold": coverage_threshold},
            "o_antigen": o, "h1_antigen": h1, "h2_antigen": h2,
            "antigenic_formula": formula, "serovar": serovar,
            "antigens": sorted(axis_best.values(), key=lambda v: v["axis"])}
    return base
