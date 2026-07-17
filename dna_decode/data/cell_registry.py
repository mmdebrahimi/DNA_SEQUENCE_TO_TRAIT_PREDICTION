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
TRACKS = ("amr", "viral", "pgx", "hla", "mendelian", "typing", "finder")


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
        cell_id="pgx:human:cyp3a5", track="pgx", route="dna-pgx", organism="human", target="cyp3a5",
        claim="CYP3A5 star-allele diplotype (*3/*6/*7) + CPIC expressor/non-expressor phenotype (tacrolimus) from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="calling_validated_underpowered_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM CDC multi-lab consensus (CYP3A5_getrm_cons) 8/8 core-diplotype on 1000G-overlap; UNDERPOWERED n=8; covers *1/*3/*6/*7 incl. *7 insertion",
        label_provenance="GeT-RM CDC reference-material multi-lab consensus (CYP3A4/CYP3A5 J Mol Diagn 2023 table) join 1000G",
        abstention_vocab=AbstentionVocab.UNDERPOWERED, native_abstention="UNDERPOWERED",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="only ~8 GeT-RM CYP3A5 samples overlap 1000G (UNDERPOWERED); rare non-core alleles mis-called *1 (no sentinel layer v0)"),
    CellContract(
        cell_id="pgx:human:tpmt", track="pgx", route="dna-pgx", organism="human", target="tpmt",
        claim="TPMT COMPOUND star-allele diplotype (*3A=*3B+*3C) + CPIC thiopurine metabolizer phenotype from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="compound_calling_validated_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM CDC consolidated consensus 85/85 core-comparable on 1000G-overlap (truth *1/*3A/*3B/*3C); compound *3A path exercised (6 *3A + 8 *3C samples)",
        label_provenance="GeT-RM CDC consolidated 363-sample PGx consensus (TPMT/NUDT15 J Mol Diagn 2022) join 1000G",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="first true compound-allele cell (>=2 SNPs in cis -> *3A); rare non-core alleles (*2/*8/*16) mis-called *1 (no sentinel layer v0)"),
    CellContract(
        cell_id="pgx:human:cyp2b6", track="pgx", route="dna-pgx", organism="human", target="cyp2b6",
        claim="CYP2B6 *6-proxy (516G>T signal) + CPIC efavirenz metabolizer phenotype from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="single_snp_proxy_calling_validated_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM CDC consolidated consensus 62/62 on clean *1/*6 truth on 1000G-overlap; SINGLE-SNP *6-proxy (516G>T) — cannot split *6/*9 (rs2279343/785A>G absent from 1000G 30x panel)",
        label_provenance="GeT-RM CDC consolidated 363-sample PGx consensus (CYP2B6) join 1000G",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="single-SNP *6-proxy: rs2279343 (785A>G) absent from the 1000G 30x panel so *6 can't be split from *9; a callset with 785 upgrades to the compound resolver (v0.1)"),
    CellContract(
        cell_id="pgx:human:cyp2d6", track="pgx", route="dna-pgx", organism="human", target="cyp2d6",
        claim="CYP2D6 SNP-surface star-allele diplotype (core {*2,*3,*4,*6,*9,*10,*17,*29,*35,*41}) + CPIC activity-score phenotype from a phased VCF — structural alleles UNASSESSED",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="snp_surface_calling_validated_structural_unassessed_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM CDC/ursaPGx consensus (CYP2D6_getrm_cons) core-comparable SNP-diplotype concordance on the SNP-decodable 1000G-overlap subset; structural alleles (*5/*13/*36/*68/*xN; ~28/87) BAM-required and EXCLUDED (cnv_hybrid_unassessed)",
        label_provenance="GeT-RM consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) CYP2D6_getrm_cons join 1000G",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="SNP surface: structural alleles NOT withheld (may be SILENTLY mis-called). Structural surface off a BAM/CRAM resolves COPY NUMBER (*5/*xN, 26/26), HYBRID PRESENCE (CYP2D7 depth, sens 0.62/spec 1.0), and HYBRID IDENTITY via read-level PSV D6-fraction (cyp2d6_hybrid_identity; Cyrius 117-PSV method; full-N GO, spec 1.0, *68 4/4 / *36 6/8); subtle *36 conversions + *13 (n=1 unpowered) abstain; non-core SNP alleles (*14/*15/*21/*40/*46) mis-called (no sentinel v0)"),
    CellContract(
        cell_id="pgx:human:dpyd", track="pgx", route="dna-pgx", organism="human", target="dpyd",
        claim="DPYD fluoropyrimidine-toxicity phenotype: CPIC activity-score over the 4 actionable DPD-deficiency haplotypes (*2A/*13 no-function, c.2846A>T/HapB3 decreased) from a phased VCF",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="cpic_activity_score_deployment_validated_no_getrm_concordance_yet",
        validation_slice="v0 DEPLOYMENT tier: decoded end-to-end on 5 real PGP-UK humans (all *1/*1 NM, no false-positive deficiency call); the 4 haplotype coords are Ensembl-GRCh38-verified. GeT-RM DPYD concordance = v0.1 (CDC characterized DPYD in the 2016/2019 rounds -> fetch+join)",
        label_provenance="CPIC DPYD guideline (Amstutz 2018) allele-functionality + PharmVar DPYD; deployment on PGP-UK PRJEB17529",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_decode_pgp_uk.py", incoming_data_gate="n/a",
        demotion_rule="all-SNP, NO structural blind spot (unlike CYP2D6); NO sentinel layer -> rarer uncertain-function DPYD alleles called *1 (CPIC's own non-actionable posture — only the 4 actionable haplotypes change fluoropyrimidine dosing)"),
    CellContract(
        cell_id="pgx:human:nudt15", track="pgx", route="dna-pgx", organism="human", target="nudt15",
        claim="NUDT15 thiopurine-toxicity phenotype: CPIC activity-score over the dominant no-function *3 (rs116855232) from a phased VCF",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="cpic_activity_score_deployment_validated_no_getrm_concordance_yet",
        validation_slice="v0 DEPLOYMENT tier: caller runs end-to-end on real VCFs (PGP-UK); *3 coord Ensembl-GRCh38-verified; *3 EAS AF ~9.5% matches the thiopurine-toxicity spectrum. GeT-RM NUDT15 concordance = external wall (paper-supplement, like DPYD)",
        label_provenance="CPIC NUDT15 guideline (Relling 2019) allele-functionality + PharmVar NUDT15; deployment on PGP-UK PRJEB17529",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_decode_pgp_uk.py", incoming_data_gate="n/a",
        demotion_rule="*2 shares rs116855232 -> called *3 (SAME no-function phenotype, CPIC call unaffected); NO sentinel layer -> rarer non-core NUDT15 alleles called *1 (only *3/*2 change thiopurine dosing at v0)"),
    CellContract(
        cell_id="pgx:human:ugt1a1", track="pgx", route="dna-pgx", organism="human", target="ugt1a1",
        claim="UGT1A1 irinotecan-toxicity phenotype: CPIC activity-score over the SNP-callable *80 (rs887829, LD-tag for the *28 TA-repeat) + *6 (rs4148323) from a phased VCF",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="cpic_activity_score_ld_tag_proxy_star28_repeat_unassessed",
        validation_slice="v0 TAG-SNP tier: rs887829 (*80) EUR AF ~30% == the *28 frequency (confirms the LD tag); coords Ensembl-GRCh38-verified; decoded on real VCFs (PGP-UK). The direct *28 TA-repeat is a STRUCTURAL WALL (repeat-aware caller needed); GeT-RM UGT1A1 concordance = external wall",
        label_provenance="CPIC UGT1A1 guideline (Gammal 2016) + PharmVar UGT1A1; rs887829 as the validated *28 LD-tag; deployment on PGP-UK PRJEB17529",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_decode_pgp_uk.py", incoming_data_gate="n/a",
        demotion_rule="STRUCTURAL: *28 (promoter TA-repeat) is NOT directly called — rs887829 (*80) is an LD-tag PROXY (EUR r^2 ~0.9+, imperfect off-EUR); star28_ta_repeat_unassessed=True. *37/*36 repeat alleles + rarer non-core called *1"),
    CellContract(
        cell_id="pgx:human:vkorc1", track="pgx", route="dna-pgx", organism="human", target="vkorc1",
        claim="VKORC1 -1639G>A (rs9923231) warfarin-sensitivity genotype from a phased VCF",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="single_snp_genotype_to_sensitivity",
        validation_slice="direct genotype readout (minus-strand encoded); not a star/diplotype system",
        label_provenance="literature sensitivity assignment (no independent panel validation in-repo)",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="n/a", demotion_rule="n/a (deterministic single-SNP readout)"),
    CellContract(
        cell_id="pgx:human:cyp4f2", track="pgx", route="dna-pgx", organism="human", target="cyp4f2",
        claim="CYP4F2 *3 (rs2108622, V433M) warfarin dose-modifier genotype+function from a phased VCF (3rd warfarin gene with VKORC1+CYP2C9)",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="single_snp_genotype_to_function_readout",
        validation_slice="single-SNP rs2108622 genotype->CYP4F2 function readout (plus-strand genomic C>T == cDNA 433 Val>Met); AF-corroborated (*3 ~29% EUR / ~79% EAS); deployed on real VCFs (PGP-UK); trio-Mendelian consistency the only validation surface (no independent star truth)",
        label_provenance="CPIC warfarin guideline (Johnson 2017) CYP4F2*3 dose effect + dbSNP rs2108622; deployment on PGP-UK PRJEB17529",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_decode_pgp_uk.py", incoming_data_gate="n/a",
        demotion_rule="single-SNP *3 proxy (rs2108622 IS the *3 truth); a warfarin DOSE modifier, not a metabolizer phenotype; the dose direction is annotation only (NOT a clinical dose)"),
    CellContract(
        cell_id="pgx:human:abcg2", track="pgx", route="dna-pgx", organism="human", target="abcg2",
        claim="ABCG2 Q141K (rs2231142) rosuvastatin transporter-function genotype from a phased VCF (pairs with SLCO1B1 for statins)",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="single_snp_genotype_to_function_readout",
        validation_slice="single-SNP rs2231142 genotype->ABCG2 transporter function readout (plus-strand genomic G>T == cDNA 141 Gln>Lys); AF-corroborated (141K ~9% EUR / ~29% EAS); deployed on real VCFs (PGP-UK); trio-Mendelian consistency the only validation surface (no independent star truth)",
        label_provenance="CPIC rosuvastatin guideline (Cooper-DeHoff 2022) ABCG2 141K function + dbSNP rs2231142; deployment on PGP-UK PRJEB17529",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_decode_pgp_uk.py", incoming_data_gate="n/a",
        demotion_rule="single-SNP Q141K readout (rs2231142 IS the truth); ABCG2 transporter FUNCTION, not a metabolizer phenotype; rosuvastatin-specific (not all statins) — annotation only, NOT a clinical dose"),
    CellContract(
        cell_id="pgx:human:slco1b1", track="pgx", route="dna-pgx", organism="human", target="slco1b1",
        claim="SLCO1B1 c.521T>C (rs4149056, *5) transporter-function genotype -> simvastatin myopathy risk, from a phased VCF",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="single_snp_genotype_to_function",
        validation_slice="direct 521T>C genotype readout (plus-strand); NOT an independent star-diplotype concordance (single-SNP tautology); genotype+trio only",
        label_provenance="CPIC simvastatin function assignment (Cooper-DeHoff 2022); no independent panel validation in-repo (rs4149056 IS the truth for a 521 call)",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="n/a",
        demotion_rule="single-SNP proxy for *5/*15/*17; full SLCO1B1 star typing needs more variants (v0 scope-limit)"),
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

# --- Non-AMR trait cells (route dna-<trait>). Registered SEPARATELY from _TYPING_FINDER because each has an
# --- individually-earned tier, and collapsing them into the shared faithful-to-tool default would overstate
# --- two of them. These three shipped CLI-routable before this registration; the coverage guard caught it —
# --- which is the guard working as designed ("a new decoder cannot ship invisibly to the trust surface").
_TRAIT_CONTRACTS: list[CellContract] = [
    CellContract(
        cell_id="typing:arabidopsis:flowering", track="typing", route="dna-flowering",
        organism="Arabidopsis_thaliana", target="flowering",
        claim="flowering HABIT (winter-annual/late vs summer-annual/early) from FRI + FLC allele calls; "
              "a two-locus AND, NOT quantitative days-to-flower",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="fri_route_scored_in_distribution_structure_confounded",
        validation_slice=(
            "TWO runs, and the second closes the first's scope limit. (1) FRI ROUTE — Zhang & "
            "Jimenez-Gomez 2020 Table S3, N=854 phenotyped of 1,017 (scripts/flowering_tables3_score.py): "
            "pooled acc 0.733 vs 0.502 null, but the HONEST figure is the population-structure-weighted "
            "0.710 vs its own 0.676 null (+3.4pp; 7/9 ancestry groups beat their null, central_europe "
            "LOSES). Directional: FRI-LoF->early 93.9% (strong) vs FRI-functional->late 65.8% (weak) = "
            "necessary-not-sufficient. (2) FLC ROUTE — the distinctive two-locus claim, VALIDATED "
            "2026-07-17 (scripts/flowering_flc_route_test.py, wiki/flowering_flc_route_2026-07-17.md) by "
            "joining AraPheno phenotype 29 (measured FLC EXPRESSION, Atwell 2010) to S3 on n=106: ALL FOUR "
            "cells of the AND call their majority correctly (functional+strong 85% late; **functional+weak "
            "39% late = the Da(1)-12 class, a 46pp separation a FRI-only rule cannot see**; lof+strong 17% "
            "= the Lz-0 class, real but RARE at 1/6, which JUSTIFIES the MEDIUM cap; lof+weak 10%). FLC "
            "EARNS its place: net +5 calls fixed (14 rescued, 9 broken) on the 70 functional-FRI "
            "accessions; **within-ancestry two-locus 0.803 vs FRI-only 0.767 vs null 0.751 -> the FLC "
            "route roughly TRIPLES the within-ancestry advantage**"),
        label_provenance=(
            "FT16_mean (days to first flower, long days 16C) from 1001 Genomes via the paper's Table S3 "
            "(CC-BY 4.0); FRI status = the paper's own `deleterious_allele` call. IN-DISTRIBUTION: the cell's "
            "catalogue and this label both trace to the same literature — NOT an independent validation"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/flowering_tables3_score.py", incoming_data_gate="G2, G8",
        demotion_rule=(
            "the WITHIN-ANCESTRY gain is the real claim (+5.2pp two-locus / +1.6pp FRI-only), NOT the "
            "pooled one; if a structure-aware re-score drops it to <=0, demote to NOT_CENSUSED. **The FLC "
            "gain RIDES ON THE THRESHOLD** — measured across FLC-expression quantiles: q20 +0.028 / q30 "
            "+0.066 / q50 +0.047 / q60 +0.000 / **q70 -0.085** — so it holds only in the biologically "
            "plausible low-quantile range (Werner 2005: weak/null FLC alleles are RARE, which a median "
            "split cannot represent) and REVERSES if weak FLC is over-called. FLC EXPRESSION is a PROXY "
            "for allele status, not the same measurement. 16% of S3 lacks FT16 with NON-RANDOM dropout "
            "(9.8% deleterious among dropped vs 24% base rate) — a re-score on the full set may move it. "
            "NB the gate tags are BY ANALOGY: G2 is defined on source-study/submitter and G8 on Mash "
            "lineages, whereas the confounding grouping variable here is the STRUCTURE ancestry group — same "
            "shape (label confounded with a grouping variable; correcting for it collapses one group to a "
            "single class and shrinks the advantage), different variable"),
    ),
    CellContract(
        cell_id="finder:any:inverse", track="finder", route="dna-decode-inverse",
        organism="any", target="inverse",
        claim="proposes the edits at a target PERCENTILE of predicted molecular damage (Regime B) using the "
              "DMS-validated forward oracle as label-free ground truth -- a RANK, never a dose",
        evidence_tier=EvidenceTier.INDEPENDENT_MEASURED,
        claim_status="dms_measured_rank_inverse_regime_b_only",
        validation_slice=(
            "graded NON-circularly against MEASURED wet-lab ProteinGym DMS -- calibrate/select on disjoint "
            "POSITION splits, grade on the proposed variant's measured value, never the model's own "
            "re-score (the generating model grading its own proposals measures self-consistency). "
            "Beats an exact closed-form no-oracle null on **4/4 usable proteins across 4 kingdoms** "
            "(E. coli/human/yeast/Arabidopsis), ~2-5 percentile points at top-5. A 5th assay (CcdB) is "
            "EXCLUDED as censored (79.3% of variants tied at its ceiling -> percentile undefined). The "
            "magnitude round-trip separately PASSES on blaTEM (+53.0%, 6/6 paired splits) but is NOT "
            "deployable -- see demotion_rule. wiki/forward_inverse_{roundtrip,sweep,deployable}"
            "_2026-07-1{6,7}.md"),
        label_provenance=("ProteinGym DMS assays (free, published wet-lab per-variant fitness). The oracle "
                          "never sees a label; labels are used ONLY to grade the proposals"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/forward_inverse_deployable.py", incoming_data_gate="n/a",
        demotion_rule=(
            "FOUR rails, each measured. (1) IT RANKS, IT DOES NOT DOSE -- a magnitude claim needs a "
            "score->effect calibrator fit on the TARGET protein's own DMS (which would make the inverse "
            "unnecessary), and calibrators CANNOT transfer: the assays share no scale (CcdB's whole range "
            "[-9.00,-2.00] lies below TEM-1's minimum -3.56), so cross-protein magnitude is impossible by "
            "construction. The conformal interval is informative 0/6 splits on blaTEM -- it brackets while "
            "proving nothing, since coverage holds even for a useless model. (2) The LEARNED oracle earns "
            "its keep over plain BLOSUM62 on 3/4 hand-picked proteins -- but at SCALE (N=200 ProteinGym, "
            "wiki/proteingym_inverse_sweep_2026-07-17.md) the SHIPPED blosum62 default beats a random pick "
            "MATERIALLY on only 13.5% and is often WORSE than guessing, so the wheel-only default is NOT a "
            "reliable design tool -- the 4/4 headline was ESM (GPU). Utility also does NOT track forward "
            "rank (PTEN 0.5185 earns keep, RL40A 0.5190 does not), so a good Spearman does not license "
            "skipping the per-protein check. "
            "(3) REGIME B ONLY -- never clinical resistance, where this scorer class is BELOW CHANCE "
            "(0.454 vs the catalogue's 0.926). (4) top-1 is ~4x worse than best-of-5: the claim is "
            "'propose k, assay k, keep the best', not 'propose 1 and trust it'. Demote if a re-run stops "
            "beating the null on a majority of usable assays"),
    ),
    CellContract(
        cell_id="typing:human:pigment", track="typing", route="dna-pigment",
        organism="human", target="pigment",
        claim="IrisPlex eye-colour probability (blue/intermediate/brown) from 6 curated SNP genotypes",
        evidence_tier=EvidenceTier.FAITHFUL_TO_TOOL,
        claim_status="faithful_to_published_irisplex_model_not_scored",
        validation_slice=(
            "reproduces the published IrisPlex multinomial model (coefficients transcribed verbatim) and its "
            "reference anchors (rs12913832 GG -> blue 0.984; AA -> brown 0.846) via reference_integrity_ok(). "
            "NOT SCORED against measured phenotypes: the openSNP scoring run is still blocked on the host "
            "being down. Faithful-to-tool, not validated"),
        label_provenance="IrisPlex published model coefficients (Walsh et al.); no measured-phenotype cohort scored yet",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="n/a",
        demotion_rule=(
            "European-ancestry-derived model; predictive accuracy is known to degrade off that ancestry and "
            "intermediate is the weakest class. A measured-phenotype score (openSNP) would move this to a real "
            "tier — in either direction"),
    ),
    CellContract(
        cell_id="finder:any:forward", track="finder", route="dna-decode-forward",
        organism="any", target="forward",
        claim="molecular-effect RANK for a protein/CDS edit (Regime B: enzyme fitness/stability), with a "
              "conformal dosage interval — NEVER a clinical-resistance call (that routes to Regime A)",
        # DMS is a FREE, INDEPENDENT, per-variant WET-LAB measurement -- the molecular analogue of HIV
        # PhenoSense, and (as forward/README says) the one place this project's label wall does not bind.
        # The predictor never sees the label, so this is measured-independent, not faithful-to-tool.
        evidence_tier=EvidenceTier.INDEPENDENT_MEASURED,
        claim_status="dms_measured_rank_validated_regime_b_only",
        validation_slice=(
            "per-variant Spearman vs measured ProteinGym DMS fitness. CLI default `blosum62` (deterministic, "
            "wheel-only): TEM-1 0.3465 (n=4996) / PTEN 0.182 -- REAL but modest, and at SCALE (N=209 "
            "ProteinGym, wiki/forward_blosum_proteingym_2026-07-17.md) the shipped blosum62 default is "
            "|Spearman| MEDIAN 0.20 -- so TEM-1's 0.35 is top-13% (27/209 reach 0.30), NOT typical. "
            "Python-API `esm2` "
            "(ESM2-650M masked-marginal): TEM-1 **0.7315** / PTEN 0.518 / CcdB 0.5115; `alphamissense` PTEN "
            "0.539 (human-only). Genome-level nucleotide-edit path validated end-to-end on a real blaTEM CDS: "
            "**Spearman 0.7611** over 1,715 real single-nt-accessible variants "
            "(wiki/blatem_genome_demo_2026-07-14.json). Dosage head: conformal coverage calibrated 10/10 "
            "proteins, informative 7/10 (wiki/forward_dosage_sweep_2026-07-15.md)"),
        label_provenance=(
            "ProteinGym deep-mutational-scanning assays (free, published wet-lab per-variant fitness; "
            "BLAT_ECOLX_Stiffler_2015 + Firnberg_2014 + Deng_2012 + Jacquier_2013, PTEN_HUMAN_Mighell_2018, "
            "CCDB_ECOLI_Tripathi_2016, RL40A_YEAST, SR43C_ARATH). The predictor never sees the label"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/resistance_conservativeness_probe.py", incoming_data_gate="n/a",
        demotion_rule=(
            "SCOPE IS THE CLAIM, and it is narrow in three ways. (1) The validated quantity is a RANK "
            "correlation per protein -- 'ranks well' != 'pins the dose' (measured: CcdB-ESM2 ranks 0.49 yet "
            "does NOT narrow its magnitude interval), so a magnitude claim needs the dosage head's own "
            "informative flag, not the Spearman. (2) REGIME B ONLY: this must NEVER be read as a resistance "
            "predictor -- on antagonistically-selected resistance the same class of scorer is BELOW CHANCE "
            "(ESM2 0.454 vs the curated catalogue's 0.926; BLOSUM62 ranks real DRMs 4.0/19), which is why "
            "the router sends determinant hits to Regime A and organism-polygenic edits to ABSTAIN. (3) The "
            "shipped CLI default is blosum62 at 0.35/0.18, NOT the 0.73 headline -- the learned methods need "
            "a precomputed score table and stay in the Python API. Demote if a DMS re-score drops the rank "
            "materially, or if any path emits an organism-level or clinical call"),
    ),
]


# --- Mendelian (germline pathogenicity) cell (route dna-clinvar). Curated ClinVar catalog decoder. ---
_MENDELIAN_CONTRACTS: list[CellContract] = [
    CellContract(
        cell_id="mendelian:human:germline_pathogenicity", track="mendelian", route="dna-clinvar",
        organism="human", target="germline_pathogenicity",
        claim="curated ClinVar germline pathogenicity (P/LP + B/LB) over the ACMG SF v3.2 + carrier 86-gene panel, from a VCF",
        evidence_tier=EvidenceTier.FAITHFUL_TO_TOOL, claim_status="curated_clinvar_catalog_faithful_to_tool",
        validation_slice="deterministic ClinVar-catalog lookup; deployment demonstration on real PGP-UK individuals (N=5; 0 reportable pathogenic = expected ACMG-SF base rate, benign carrier load surfaced); faithful-to-ClinVar, no independent truth beyond the curated DB",
        label_provenance="ClinVar curated germline classifications (NCBI); P/LP + B/LB only; ACMG SF v3.2 (81) + 5 carrier genes",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="n/a",
        demotion_rule="VUS/conflicting excluded (deployable-claim tier only); bounded 86-gene panel -> out-of-panel = INDETERMINATE (absence != benign)"),
]


# --- HLA drug-hypersensitivity cells (route dna-hla). Tag-SNP LD-proxy carriage callers. ---
def _hla_contracts() -> list[CellContract]:
    from dna_decode.hla.catalog import CATALOG
    out: list[CellContract] = []
    for key, a in CATALOG.items():  # CATALOG now holds ONLY the validated cell(s) (b5701); failed tags demoted
        out.append(CellContract(
            cell_id=f"hla:human:{key}", track="hla", route="dna-hla", organism="human", target=key,
            claim=f"{a.allele} carriage (tag SNP {a.rsid}) -> {a.drug} {a.reaction} risk (CPIC)",
            evidence_tier=EvidenceTier.NEAR_INDEPENDENT,
            claim_status="tag_snp_ld_proxy_validated_vs_1000g_hla_truth",
            validation_slice=("sample-level concordance vs the free 1000G HLA truth (20140702_hla_diversity, "
                              "n=1103): sens 0.979 / spec 0.992 / PPV 0.855 — the deployed clinical abacavir "
                              "screen (rs2395029), independently measured"),
            label_provenance="1000G HLA types (20140702_hla_diversity) join rs2395029 tag genotypes; CPIC abacavir guideline",
            abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
            falsifier_ref="scripts/hla_concordance.py", incoming_data_gate="n/a",
            demotion_rule=("LD PROXY (not sequence-based typing) but VALIDATED vs real HLA truth (sens 0.979); "
                           "the sibling provisional tags (B*58:01 rs9263726 sens 0.61 weak; A*31:01 rs1061235 "
                           "not-paneled sens 0.0) FAILED validation and are demoted, NOT shipped")))
    return out


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
    """Every v0.1 cell contract (AMR projection + viral + PGx + HLA + Mendelian + typing/finder + traits)."""
    return (_amr_contracts() + _viral_contracts() + list(_PGX_CONTRACTS) + _hla_contracts()
            + list(_MENDELIAN_CONTRACTS) + _typing_finder_contracts() + list(_TRAIT_CONTRACTS))


def by_cell_id() -> dict[str, CellContract]:
    return {c.cell_id: c for c in cells()}


def amr_cells() -> list[CellContract]:
    return [c for c in cells() if c.track == "amr"]


def pgx_cells() -> list[CellContract]:
    return [c for c in cells() if c.track == "pgx"]


def mendelian_cells() -> list[CellContract]:
    return [c for c in cells() if c.track == "mendelian"]


def hla_cells() -> list[CellContract]:
    return [c for c in cells() if c.track == "hla"]


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
    from dna_decode.hla import HLA_ALLELES
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
        "dna-clinvar": {"germline_pathogenicity"},  # the Mendelian (ClinVar) single-decoder route
        "dna-hla": set(HLA_ALLELES),                 # HLA drug-hypersensitivity tag-SNP cells
        "traits": set(TRAITS) - {"amr", "pgx"},  # the typing/finder whole-tool traits
    }


def cli_routable_cell_ids() -> set[str]:
    return set(by_cell_id())
