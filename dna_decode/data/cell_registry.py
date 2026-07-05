"""Evidence-Contract Registry (v0.1) — one checked-in, test-enforced contract per shipped cell.

A `CellContract` declares, for each deployed decoder cell, WHAT it claims, at what HONEST evidence tier,
on what validation SLICE, with what label provenance, what abstention vocabulary it speaks, and (declared,
not executed) its falsifier / incoming-data gate / demotion rule. The validation report card reads its AMR
grid from here so a shipped decoder cannot ship invisibly and abstention has ONE vocabulary.

v0.1 SCOPE (full CLI-routable surface — the brainstorm C1 per-route manifest):
  - `amr`   (route dna-amr): the frozen `shipped_decoder_surface` (bacterial + fungal + antimalarial +
            influenza-antiviral) projected verbatim.
  - `viral` (route dna-amr): HIV-1 + SARS-CoV-2 drugs — CLI-routable via `dna-amr --drug` but NOT in the
            AMR surface (their own report cards). track=viral, route=amr (brainstorm C2's track/route split).
  - `pgx`   (route dna-pgx): CYP2C19 / CYP2C9 / VKORC1.
  - `typing`/`finder` (route dna-<trait>): the 10 whole-tool typing + determinant-finder decoders.

INTEGRITY RAILS (load-bearing):
- `cell_id` is a DISPLAY string ONLY (brainstorm C2). The AMR join key is `cell_key.canonical_cell_key`
  (organism, drug); `amr_projection_keys()` returns exactly that set and the consistency test asserts it
  EQUALS the frozen surface's keys. Never join AMR cells by raw `cell_id` string.
- AMR contracts are a PROJECTION of `shipped_decoder_surface.shipped_decoder_rows()` built programmatically,
  so the projection == surface BY CONSTRUCTION. `surface_index()` re-exports the surface-shaped dict FROM
  the registry so the report card reads its grid from here (equal by construction → 0 behavior change).
- `cli_routable_manifest()` is derived from the LIVE CLI catalogs (the `dna-amr --drug` union, `dna-pgx
  --gene`, `dna_decode.cli.TRAITS`) so the coverage test cannot silently drift as the CLI grows.
- NO numeric confidence field exists anywhere (anti-"trust-layer-theater" guardrail). `evidence_tier` is a
  categorical honesty label; `claim_status` carries the structural status separately (brainstorm M1).
- Imports `shipped_decoder_surface` READ-ONLY and touches NO frozen file -> `test_tb_leak_guard.py` green.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from dna_decode.data.cell_key import canonical_cell_key
from dna_decode.data.cell_registry_vocab import AbstentionVocab
from dna_decode.data.shipped_decoder_surface import shipped_decoder_rows


class EvidenceTier(str, Enum):
    """Honest evidence tier per cell. Categorical, NOT numeric (guardrail)."""

    INDEPENDENT_MEASURED = "independent_measured"   # free INDEPENDENT isolate-level wet-lab label (HIV PhenoSense)
    NEAR_INDEPENDENT = "near_independent"           # provenance-disjoint stress test / consensus panel (NCBI-PD, GeT-RM)
    FAITHFUL_TO_TOOL = "faithful_to_tool"           # faithful to a reference tool/DB/guideline, not an independent label
    KNOWLEDGE_BASELINE = "knowledge_baseline"       # literature/catalogue assignment, in-distribution
    NO_FREE_SOURCE = "no_free_source"               # no free isolate-level phenotype source exists
    NOT_CENSUSED = "not_censused"                    # CLI-routable but never validated (no report-card row)


HIV_UNDERPOWERED_N = 50  # report-card n below this -> measured-but-UNDERPOWERED rather than full SCORED


# Tracks (v0.1). amr/viral both route through `dna-amr`; route ≠ track (brainstorm C2 split).
TRACKS = ("amr", "viral", "pgx", "typing", "finder")


@dataclass(frozen=True)
class CellContract:
    """One shipped decoder cell's evidence contract. Frozen; NO numeric confidence field by design."""

    cell_id: str               # DISPLAY string only ("track:organism:target"); NOT the join key
    track: str                 # one of TRACKS
    route: str                 # the CLI entrypoint: "dna-amr" | "dna-pgx" | "dna-<trait>"
    organism: str
    target: str                # drug (amr/viral) | gene (pgx) | scheme/tool name (typing/finder)
    claim: str                 # one-line plain claim the cell makes
    evidence_tier: EvidenceTier
    claim_status: str          # structural status (phenotype_source_status for amr; calling-status etc) — M1 split
    validation_slice: str      # the slice the tier was earned on
    label_provenance: str      # where the labels came from
    abstention_vocab: AbstentionVocab  # this cell's abstention KIND, collapsed to the controlled vocab
    native_abstention: str     # the cell's own raw in-tree abstention term
    falsifier_ref: str         # path to a falsifier script, or "none" (DECLARED, not executed)
    incoming_data_gate: str    # which of the 8 rejection gates apply, or "n/a" (DECLARED)
    demotion_rule: str         # free-text v0: the trigger that would demote this cell's tier
    # AMR-only surface fields (carried so surface_index() re-exports the surface-shaped dict from the registry):
    engine: str | None = None
    organism_scope: str | None = None
    census_group: str | None = None


# --- AMR phenotype_source_status -> (evidence_tier, abstention_vocab, native) ---
_AMR_STATUS_MAP: dict[str, tuple[EvidenceTier, AbstentionVocab, str]] = {
    "ncbi_pd":          (EvidenceTier.NEAR_INDEPENDENT, AbstentionVocab.SCORED, "SCORED"),
    "label_confounded": (EvidenceTier.FAITHFUL_TO_TOOL, AbstentionVocab.LABEL_CONFOUNDED, "LABEL_CONFOUNDED"),
    "no_free_source":   (EvidenceTier.NO_FREE_SOURCE,   AbstentionVocab.NO_FREE_SOURCE, "NO_FREE_PHENOTYPE"),
}


def _amr_contracts() -> list[CellContract]:
    """Project every frozen `shipped_decoder_surface` row to an AMR CellContract (== surface by construction)."""
    out: list[CellContract] = []
    for r in shipped_decoder_rows():
        org, drug, status = r["organism"], r["drug"], r["phenotype_source_status"]
        tier, vocab, native = _AMR_STATUS_MAP[status]
        scoreable = status == "ncbi_pd"
        out.append(CellContract(
            cell_id=f"amr:{org}:{drug}", track="amr", route="dna-amr", organism=org, target=drug,
            claim=f"{r['engine']} R/S call for {org} x {drug}",
            evidence_tier=tier, claim_status=status,
            validation_slice=("NCBI-PD provenance-disjoint stress test (lineage-disclosed)" if scoreable
                              else "label-confounded surrogate (cefoxitin is the CLSI surrogate)"
                              if status == "label_confounded" else "no free isolate-level phenotype source"),
            label_provenance=("NCBI Pathogen Detection AST_phenotypes" if scoreable else "none (structural non-cell)"),
            abstention_vocab=vocab, native_abstention=native,
            falsifier_ref="scripts/provenance_disjoint_validate.py" if scoreable else "none",
            incoming_data_gate="G1,G7,G8" if scoreable else "n/a",
            demotion_rule=("SCORED -> UNDERPOWERED below the powering floor; lineage-collapse can demote the "
                           "disclosed metric" if scoreable else "n/a (no free label to demote against)"),
            engine=r["engine"], organism_scope=r["organism_scope"], census_group=r["census_group"],
        ))
    return out


# --- Viral cells: HIV-1 + SARS-CoV-2 drugs route via `dna-amr --drug` but are NOT in the AMR surface. ---
def _hiv_card_drugs() -> dict[str, dict]:
    """{drug: card_row} from the PACKAGED HIV report card (trust_surface loader; wheel-safe). {} if absent."""
    from dna_decode.data import trust_surface
    card = trust_surface._load("hiv_decoder_report_card.json")
    if not card:
        return {}
    return {c["drug"]: c for c in card.get("cells", []) if c.get("drug")}


def _viral_contracts() -> list[CellContract]:
    from dna_decode.data.hiv_amr import all_supported_hiv_drugs
    from dna_decode.data.sarscov2_amr import all_supported_sarscov2_drugs
    out: list[CellContract] = []
    card = _hiv_card_drugs()  # data-drive HIV tiers from the report card (brainstorm C1: no overclaim)
    for d in sorted(all_supported_hiv_drugs()):
        row = card.get(d)
        if row is None:
            # CLI-routable but NOT in the validation report card (e.g. delavirdine) -> NOT_CENSUSED, never measured
            tier, status, vocab, native = (EvidenceTier.NOT_CENSUSED, "cli_routable_not_validated",
                                           AbstentionVocab.NOT_CENSUSED, "NOT_CENSUSED")
            vslice = "CLI-routable; NOT in the HIV report card (uncensused)"
            falsifier = "none"
        elif (row.get("n") or 0) >= HIV_UNDERPOWERED_N:
            tier, status, vocab, native = (EvidenceTier.INDEPENDENT_MEASURED, "independent_wetlab_validated",
                                           AbstentionVocab.SCORED, "SCORED")
            vslice = f"Stanford HIVDB PhenoSense fold-change (n={row.get('n')}, class={row.get('drug_class')})"
            falsifier = "scripts/hiv_targetsite_validate.py"
        else:
            tier, status, vocab, native = (EvidenceTier.INDEPENDENT_MEASURED, "independent_wetlab_underpowered",
                                           AbstentionVocab.UNDERPOWERED, "UNDERPOWERED")
            vslice = f"Stanford HIVDB PhenoSense fold-change, UNDERPOWERED (n={row.get('n')})"
            falsifier = "scripts/hiv_targetsite_validate.py"
        out.append(CellContract(
            cell_id=f"viral:HIV-1:{d}", track="viral", route="dna-amr", organism="HIV-1", target=d,
            claim=f"HIV-1 RT/PR/IN/CA target-site resistance call for {d}",
            evidence_tier=tier, claim_status=status, validation_slice=vslice,
            label_provenance=("Stanford HIVDB PhenoSense (Rhee 2003); catalog from the HIVDB dataset page"
                              if row is not None else "none (CLI-routable, not yet validated)"),
            abstention_vocab=vocab, native_abstention=native,
            falsifier_ref=falsifier, incoming_data_gate="n/a",
            demotion_rule="re-tiers from the HIV report card on revalidation",
        ))
    for d in sorted(all_supported_sarscov2_drugs()):
        out.append(CellContract(
            cell_id=f"viral:SARS-CoV-2:{d}", track="viral", route="dna-amr", organism="SARS-CoV-2", target=d,
            claim=f"SARS-CoV-2 Mpro inhibitor resistance call for {d}",
            evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="cov_rdb_in_distribution",
            validation_slice="CoV-RDB invitro_selection fold-change (in-distribution knowledge baseline)",
            label_provenance="Stanford CoV-RDB (covid-drdb-payload); UNDERPOWERED / TN-starved",
            abstention_vocab=AbstentionVocab.UNDERPOWERED, native_abstention="UNDERPOWERED",
            falsifier_ref="scripts/sarscov2_mpro_validate.py", incoming_data_gate="n/a",
            demotion_rule="held-out CoV-RDB or clinical fold would re-tier toward independent",
        ))
    return out


# --- PGx cells (dna-pgx --gene). ---
_PGX_CONTRACTS: list[CellContract] = [
    CellContract(
        cell_id="pgx:human:cyp2c19", track="pgx", route="dna-pgx", organism="human", target="cyp2c19",
        claim="CYP2C19 star-allele diplotype + CPIC metabolizer phenotype from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="calling_validated_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM consensus core-diplotype concordance on real 1000G + trio Mendelian QC",
        label_provenance="GeT-RM consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) join 1000G",
        abstention_vocab=AbstentionVocab.WITHHELD_NONCORE, native_abstention="phenotype_withheld",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="a non-core *4/*35 sentinel hit -> phenotype withheld rather than mis-called"),
    CellContract(
        cell_id="pgx:human:cyp2c9", track="pgx", route="dna-pgx", organism="human", target="cyp2c9",
        claim="CYP2C9 star-allele diplotype + CPIC activity-score phenotype from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="calling_validated_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM consensus core-diplotype concordance 73/73 on real 1000G + trio Mendelian QC",
        label_provenance="GeT-RM consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) join 1000G",
        abstention_vocab=AbstentionVocab.WITHHELD_NONCORE, native_abstention="phenotype_withheld",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="a non-core *5/*8/*9/*11 sentinel hit -> phenotype withheld"),
    CellContract(
        cell_id="pgx:human:cyp2c8", track="pgx", route="dna-pgx", organism="human", target="cyp2c8",
        claim="CYP2C8 star-allele diplotype (*2/*3/*4) from a phased VCF — CALLING only, NO CPIC phenotype",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="calling_validated_no_cpic_phenotype_substrate_dependent",
        validation_slice="GeT-RM CYP2C8_getrm_ngs core-diplotype concordance 82/82 on real 1000G (Docker-free tabix-HTTP region fetch)",
        label_provenance="GeT-RM consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) CYP2C8_getrm_ngs join 1000G",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="rare non-core CYP2C8 allele mis-called *1 (no sentinel layer in v0); function is substrate-dependent so NO PM/IM/NM is ever emitted"),
    CellContract(
        cell_id="pgx:human:vkorc1", track="pgx", route="dna-pgx", organism="human", target="vkorc1",
        claim="VKORC1 -1639G>A (rs9923231) warfarin-sensitivity genotype from a phased VCF",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="single_snp_genotype_to_sensitivity",
        validation_slice="direct genotype readout (minus-strand encoded); not a star/diplotype system",
        label_provenance="literature sensitivity assignment (no independent panel validation in-repo)",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="n/a", demotion_rule="n/a (deterministic single-SNP readout)"),
]

# --- Typing + determinant-finder whole-tool cells (route dna-<trait>). Faithful-to-tool curated-DB callers. ---
# (track, trait, organism-scope, one-line claim, native-abstention)
_TYPING_FINDER: list[tuple[str, str, str, str, str]] = [
    ("typing", "pathotype", "Escherichia_coli", "E. coli pathotype compatibility call + abstention (VirulenceFinder resolver)", "ABSTAIN"),
    ("typing", "serotype", "Escherichia_coli", "E. coli O:H serotype (SerotypeFinder allele DB)", "O?"),
    ("typing", "mlst", "bacteria", "multi-locus sequence type (PubMLST allele->profile->ST)", "ABSTAIN"),
    ("typing", "ktype", "Klebsiella", "Klebsiella K/O capsule type (Kaptive)", "ABSTAIN"),
    ("typing", "salmserovar", "Salmonella", "Salmonella serovar (antigenic-formula resolver)", "ABSTAIN"),
    ("typing", "pneumoserotype", "Streptococcus_pneumoniae", "pneumococcal capsular serotype", "ABSTAIN"),
    ("finder", "plasmid", "bacteria", "plasmid Inc-replicon typing (PlasmidFinder allele DB)", "ABSTAIN"),
    ("finder", "resfinder", "bacteria", "acquired AMR genes (ResFinder allele DB) — independent cross-tool check vs amr", "ABSTAIN"),
    ("finder", "pointfinder", "Escherichia_coli", "chromosomal AMR point mutations (PointFinder) — independent vs amr POINT", "ABSTAIN"),
    ("finder", "disinfinder", "bacteria", "biocide/disinfectant resistance genes (DisinFinder)", "ABSTAIN"),
]


def _typing_finder_contracts() -> list[CellContract]:
    from dna_decode.data.cell_registry_vocab import to_vocab
    out: list[CellContract] = []
    for track, trait, scope, claim, native in _TYPING_FINDER:
        out.append(CellContract(
            cell_id=f"{track}:{scope}:{trait}", track=track, route=f"dna-{trait}", organism=scope, target=trait,
            claim=claim, evidence_tier=EvidenceTier.FAITHFUL_TO_TOOL,
            claim_status="curated_db_caller_faithful_to_tool",
            validation_slice="deterministic curated-allele-DB caller; faithful-to-tool, not an independent baseline",
            label_provenance="the tool's own reference allele DB",
            abstention_vocab=to_vocab(native), native_abstention=native,
            falsifier_ref="none", incoming_data_gate="n/a",
            demotion_rule="an independent-baseline comparison on disjoint data would re-tier"))
    return out


def cells() -> list[CellContract]:
    """Every v0.1 cell contract (AMR projection + viral + PGx + typing/finder)."""
    return _amr_contracts() + _viral_contracts() + list(_PGX_CONTRACTS) + _typing_finder_contracts()


def by_cell_id() -> dict[str, CellContract]:
    return {c.cell_id: c for c in cells()}


def amr_cells() -> list[CellContract]:
    return [c for c in cells() if c.track == "amr"]


def pgx_cells() -> list[CellContract]:
    return [c for c in cells() if c.track == "pgx"]


def amr_projection_keys() -> set[tuple[str, str]]:
    """AMR cells' canonical (organism, drug) join keys — for the surface-consistency test (NOT cell_id)."""
    return {canonical_cell_key(c.organism, c.target) for c in amr_cells()}


def surface_index() -> dict[tuple[str, str], dict]:
    """(organism.lower, drug.lower) -> surface-shaped row dict, re-exported FROM the registry's AMR cells.

    The validation report card reads its grid from here (== the frozen surface_index by construction).
    """
    out: dict[tuple[str, str], dict] = {}
    for c in _amr_contracts():  # project directly from AMR cells (narrow import surface; M2)
        out[canonical_cell_key(c.organism, c.target)] = {
            "organism": c.organism, "drug": c.target, "engine": c.engine,
            "organism_scope": c.organism_scope, "phenotype_source_status": c.claim_status,
            "census_group": c.census_group,
        }
    return out


def cli_routable_manifest() -> dict[str, set[str]]:
    """The authoritative v0.1 CLI-routable set, derived LIVE from the CLI catalogs (drift-proof)."""
    from dna_decode.cli import TRAITS
    from dna_decode.pgx import PGX_GENES
    from dna_decode.data.antimalarial_amr import supported_antimalarial_drugs
    from dna_decode.data.antiviral_amr import supported_antiviral_drugs
    from dna_decode.data.fungal_amr import supported_fungal_drugs
    from dna_decode.data.hiv_amr import all_supported_hiv_drugs
    from dna_decode.data.mic_tiers import supported_drugs
    from dna_decode.data.sarscov2_amr import all_supported_sarscov2_drugs
    amr_drugs = (set(supported_drugs()) | set(supported_fungal_drugs())
                 | set(supported_antimalarial_drugs()) | set(supported_antiviral_drugs())
                 | set(all_supported_hiv_drugs()) | set(all_supported_sarscov2_drugs()))
    return {
        "dna-amr": {d.lower() for d in amr_drugs},
        "dna-pgx": set(PGX_GENES),
        "traits": set(TRAITS) - {"amr", "pgx"},  # the typing/finder whole-tool traits
    }


def cli_routable_cell_ids() -> set[str]:
    return set(by_cell_id())
