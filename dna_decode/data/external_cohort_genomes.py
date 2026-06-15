"""External-cohort genome resolution — BioSample -> downloadable GCA accession.

The cohort isolate identity is the BioSample (every ENA run + every NCBI assembly
carries one; the MIC table is keyed by, or resolvable to, BioSample). This module
maps each labeled BioSample to a downloadable assembly accession:

  - FREE              : BioSamples with >=1 linked GCA/GCF -> {biosample: gca}
  - ASSEMBLY_REQUIRED : BioSamples with NO linked assembly (reads-only) -> excluded
                        from the free pilot, COUNTED + REPORTED (never silently dropped)

IMPORTANT (brainstorm C2): this returns ACCESSIONS, not FASTA paths. The downstream
scorer feeds each GCA accession to `organism_drug_validate.ensure_run`, which performs
the download itself (it does NOT accept a FASTA path). No download happens here.
"""
from __future__ import annotations

from typing import Iterable, Protocol


class _Resolver(Protocol):
    def biosample_to_assemblies(self, biosample: str) -> list[str]: ...


def pick_assembly(gcas: Iterable[str]) -> str | None:
    """Deterministically pick ONE assembly accession for a BioSample.

    Preference: a RefSeq `GCF_` accession (annotated), else a GenBank `GCA_`, else
    the lexicographically smallest. None if the list is empty (reads-only BioSample).
    Deterministic so re-runs select the same genome.
    """
    items = sorted({g for g in gcas if g})
    if not items:
        return None
    gcf = [g for g in items if g.startswith("GCF_")]
    if gcf:
        return gcf[0]
    gca = [g for g in items if g.startswith("GCA_")]
    if gca:
        return gca[0]
    return items[0]


def resolve_cohort_genomes(biosamples: Iterable[str], resolver: _Resolver) -> dict:
    """Resolve labeled BioSamples to GCA accessions.

    Returns {"free": {biosample: gca}, "assembly_required": [biosample,...],
    "n_free": int, "n_assembly_required": int}. `assembly_required` is the
    reads-only subset — reported, not dropped.
    """
    free: dict[str, str] = {}
    assembly_required: list[str] = []
    for bs in sorted(set(biosamples)):
        gca = pick_assembly(resolver.biosample_to_assemblies(bs))
        if gca:
            free[bs] = gca
        else:
            assembly_required.append(bs)
    return {
        "free": free,
        "assembly_required": assembly_required,
        "n_free": len(free),
        "n_assembly_required": len(assembly_required),
    }
