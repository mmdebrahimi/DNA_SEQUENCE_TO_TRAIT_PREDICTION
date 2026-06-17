"""TB lineage assignment from the pinned barcode -> the `clusters` map for clonality collapse.

Ratified F (brainstorm C2b): replaces de-novo SNP-distance clustering. Each isolate is assigned its
DEEPEST supported sublineage by matching the barcode's derived alleles against its masked-VCF calls
(Step 1). The resulting `strain_id -> cluster_id` dict is exactly the shape
`clonality.cluster_weighted_confusion` consumes.

UNASSIGNED isolates (no barcode allele carried — e.g. a near-reference lineage-4 genome) each get their
OWN singleton cluster: they are reported, NEVER silently merged into one giant bucket (which would
fabricate a huge fake clone and crush effective_lineage_n).
"""
from __future__ import annotations

from dna_decode.data.tb_lineage_barcode import BarcodeSNP
from dna_decode.organism_rules.tb_vcf import VariantCall

UNASSIGNED = "UNASSIGNED"


def assign_lineage(masked_calls: dict[int, VariantCall], barcode: list[BarcodeSNP]) -> str:
    """Deepest supported sublineage, or UNASSIGNED. Carried = isolate's alt at a barcode pos == its allele."""
    support: dict[str, int] = {}
    for snp in barcode:
        c = masked_calls.get(snp.pos)
        if c is not None and c.alt == snp.allele:
            support[snp.lineage] = support.get(snp.lineage, 0) + 1
    if not support:
        return UNASSIGNED
    # deepest = most dotted; tie-break by support count then name (deterministic)
    return max(support, key=lambda L: (L.count("."), support[L], L))


def lineage_clusters(calls_by_strain: dict[str, dict[int, VariantCall]],
                     barcode: list[BarcodeSNP]) -> dict[str, int]:
    """{strain_id: cluster_id}. Same-sublineage strains share an id; each UNASSIGNED is a singleton."""
    assigned = {sid: assign_lineage(calls, barcode) for sid, calls in calls_by_strain.items()}
    out: dict[str, int] = {}
    lin_to_id: dict[str, int] = {}
    nxt = 0
    for sid in sorted(assigned):
        lin = assigned[sid]
        if lin == UNASSIGNED:
            out[sid] = nxt
            nxt += 1
        else:
            if lin not in lin_to_id:
                lin_to_id[lin] = nxt
                nxt += 1
            out[sid] = lin_to_id[lin]
    return out


def lineage_assignments(calls_by_strain: dict[str, dict[int, VariantCall]],
                        barcode: list[BarcodeSNP]) -> dict[str, str]:
    """{strain_id: lineage_string} for reporting (incl. UNASSIGNED count)."""
    return {sid: assign_lineage(calls, barcode) for sid, calls in calls_by_strain.items()}
