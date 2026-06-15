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

# BioSample accession prefixes: SAMN (NCBI), SAMEA (EBI/ENA), SAMD (DDBJ).
BIOSAMPLE_PREFIXES = ("SAMN", "SAMEA", "SAMD")


class ExternalKeyError(ValueError):
    """Cohort isolate keys are not BioSample accessions (the join-chain contract).

    The cohort identity MUST be a BioSample so BioSample->assembly resolution works.
    A MIC table keyed by sample_alias / secondary_sample_accession / run_accession /
    isolate name silently collapses to ASSEMBLY_REQUIRED — caught here instead. Supply
    a `crosswalk` (alias -> BioSample) from the MIC-ingester to resolve such keys.
    """


class _Resolver(Protocol):
    def biosample_to_assemblies(self, biosample: str) -> list[str]: ...


def is_biosample_key(key: str) -> bool:
    return isinstance(key, str) and key.strip().upper().startswith(BIOSAMPLE_PREFIXES)


def validate_biosample_keys(keys: Iterable[str]) -> dict:
    """Partition keys into BioSample-shaped vs not. Advisory shape-check (not resolution)."""
    bad = sorted({k for k in keys if not is_biosample_key(k)})
    return {"ok": not bad, "bad_keys": bad, "n_bad": len(bad)}


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


def resolve_cohort_genomes(biosamples: Iterable[str], resolver: _Resolver, *,
                           crosswalk: dict[str, str] | None = None,
                           require_biosample_keys: bool = True) -> dict:
    """Resolve labeled BioSamples to GCA accessions.

    Cohort keys MUST be BioSample accessions (the join-chain contract). A `crosswalk`
    (alias -> BioSample) from the MIC-ingester is applied first; with
    `require_biosample_keys` (default True), any key that is still not BioSample-shaped
    raises `ExternalKeyError` rather than silently collapsing to ASSEMBLY_REQUIRED.

    Returns {"free": {biosample: gca}, "assembly_required": [biosample,...],
    "n_free": int, "n_assembly_required": int}. `assembly_required` is the reads-only
    subset — reported, not dropped.
    """
    crosswalk = crosswalk or {}
    resolved_keys = sorted({crosswalk.get(bs, bs) for bs in biosamples})
    if require_biosample_keys:
        check = validate_biosample_keys(resolved_keys)
        if not check["ok"]:
            raise ExternalKeyError(
                f"{check['n_bad']} cohort key(s) are not BioSample accessions: {check['bad_keys'][:5]}"
                f"{'...' if check['n_bad'] > 5 else ''}. Supply a crosswalk (alias->BioSample) from "
                f"the MIC-ingester, or pass require_biosample_keys=False to bypass (NOT recommended)."
            )
    free: dict[str, str] = {}
    assembly_required: list[str] = []
    for bs in resolved_keys:
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
