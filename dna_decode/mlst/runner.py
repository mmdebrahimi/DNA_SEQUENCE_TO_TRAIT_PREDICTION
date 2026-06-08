"""MLST caller — exact per-locus allele via the shared engine + profile->ST lookup.

For each locus FASTA: blastn its alleles vs the genome at 100% identity / 100% coverage (an MLST allele
call is an EXACT match); the called allele's number is the locus call. Assemble the profile -> look up the
ST in the PubMLST profiles table. Novel/incomplete profile -> ST None (reported honestly, never guessed).
Offline-safe via the engine.
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.mlst.core import allele_number, lookup_st, parse_profiles
from dna_decode.typing.blast_caller import call_alleles

MLST_IDENTITY = 100.0   # exact allele match
MLST_COVERAGE = 100.0


def call_mlst(genome_fasta: str | Path, locus_fastas: dict[str, str | Path], profiles_tsv: str | Path, *,
              blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """Type a genome: exact allele per locus -> profile -> ST. `locus_fastas` = {locus: allele .fasta}."""
    if not Path(profiles_tsv).exists():
        return {"status": "unavailable", "reason": f"profiles table not found at {profiles_tsv}",
                "st": None, "profile": {}}
    loci_order, st_of, cc_of = parse_profiles(Path(profiles_tsv).read_text(encoding="utf-8", errors="replace"))

    profile: dict[str, int | None] = {}
    for locus, fa in locus_fastas.items():
        if not Path(fa).exists():
            profile[locus] = None
            continue
        res = call_alleles(genome_fasta, fa, identity_threshold=MLST_IDENTITY,
                           coverage_threshold=MLST_COVERAGE, blastn_bin=blastn_bin, timeout=timeout)
        if res["status"] != "ok":
            return {"status": "unavailable", "tool": res.get("tool"), "reason": res.get("reason"),
                    "st": None, "profile": {}}
        # exact (100/100) called alleles for this locus; pick the one (highest coverage, then lowest number)
        exact = []
        for aid, hit in res["per_allele"].items():
            if hit["called"]:
                parsed = allele_number(aid)
                if parsed and parsed[0] == locus:
                    exact.append(parsed[1])
        profile[locus] = min(exact) if exact else None   # min: lowest allele number on a (rare) tie

    order = loci_order or sorted(locus_fastas)
    st = lookup_st({k: v for k, v in profile.items()}, order, st_of)
    complete = all(profile.get(loc) is not None for loc in order)
    key = tuple(profile.get(loc) for loc in order) if complete else None
    return {
        "status": "ok", "tool": "blastn", "method": "pubmlst_blastn_exact_v0",
        "scheme_loci": order,
        "profile": {loc: profile.get(loc) for loc in order},
        "st": st,
        "clonal_complex": cc_of.get(key) if (key and st) else None,
        "complete": complete,
        "novel": complete and st is None,    # full profile but not in the table = novel ST
    }
