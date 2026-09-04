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
    # which of the TEN rejection gates apply, or "n/a" (DECLARED, never executed here).
    # Was "8" until 2026-08-31: G9/G10 were added 2026-08-26 and 17 cells already declare them
    # (7 x "G9,G10", 7 x "G9", 3 x "G10"), so the count was stale against data in this same repo.
    # G1-G8 gate whether a usable LABEL exists; G9-G10 gate whether the decoder's own rule is
    # scoreable against a genotype at all. See wiki/negative_results_map_2026-06-13.md.
    incoming_data_gate: str
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


# HIV drugs that are NOT_CENSUSED for a DURABLE reason (no free label will ever appear), NOT pending work.
# delavirdine is a WITHDRAWN first-gen NNRTI absent from ALL Stanford PhenoSense datasets (only EFV/NVP/ETR/RPV/DOR
# present — verified 0/9 HIV datasets 2026-07-22). We KEEP it NOT_CENSUSED (the C1 anti-overclaim guardrail + its
# regression pin: a no-card-row drug must never claim to be measured; whether it should move to NO_FREE_SOURCE is a
# tiering decision left to the user) but ANNOTATE the reason so it is not misread as an open TODO.
HIV_NOT_CENSUSED_DURABLE_REASON = {
    "delavirdine": "withdrawn first-gen NNRTI absent from Stanford PhenoSense (0/9 datasets, 2026-07-22) — "
                   "no free fold-change label will appear; not pending work",
}


def _hcmv_card_drugs() -> dict[str, dict]:
    """{drug: card_row} from the PACKAGED HCMV report card (wheel-safe; force-included in pyproject).

    HCMV shipped 2026-07-23 with its own report card and five CLI-routable drugs, but no contract here --
    influenza NA is PROJECTED automatically from the frozen shipped surface, while HIV/SARS-CoV-2 are
    hand-declared in this block, and HCMV was neither. Data-driven from its card for the same reason HIV
    is: so the tier is read from the evidence rather than asserted.
    """
    from dna_decode.data import trust_surface
    card = trust_surface._load("hcmv_decoder_report_card.json")
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
            reason = HIV_NOT_CENSUSED_DURABLE_REASON.get(d)
            vslice = (f"CLI-routable; NOT in the HIV report card (uncensused) — {reason}" if reason
                      else "CLI-routable; NOT in the HIV report card (uncensused)")
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
    # HCMV (UL97/UL54/UL56) — the first herpesvirus cell. Tier read from its own card, never asserted.
    from dna_decode.data.hcmv_amr import all_supported_hcmv_drugs
    hcard = _hcmv_card_drugs()
    for d in sorted(all_supported_hcmv_drugs()):
        row = hcard.get(d)
        if row is None:
            out.append(CellContract(
                cell_id=f"viral:HCMV:{d}", track="viral", route="dna-amr", organism="HCMV", target=d,
                claim=f"HCMV antiviral-resistance call for {d}",
                evidence_tier=EvidenceTier.NOT_CENSUSED, claim_status="cli_routable_not_validated",
                validation_slice="CLI-routable; NOT in the HCMV report card (uncensused)",
                label_provenance="none (CLI-routable, not yet validated)",
                abstention_vocab=AbstentionVocab.NOT_CENSUSED, native_abstention="NOT_CENSUSED",
                falsifier_ref="none", incoming_data_gate="n/a",
                demotion_rule="re-tiers from the HCMV report card on validation"))
            continue
        genes = "/".join(row.get("genes") or [])
        out.append(CellContract(
            cell_id=f"viral:HCMV:{d}", track="viral", route="dna-amr", organism="HCMV", target=d,
            claim=f"HCMV {genes} target-site antiviral-resistance call for {d}",
            # IN_DISTRIBUTION on its own card -> KNOWLEDGE_BASELINE here. Its independence field records
            # the closure explicitly: HCMV phenotyping IS per-mutation recombinant marker-transfer and the
            # Chou compilations are its consensus, so no held-out per-isolate set disjoint from the
            # catalog exists. Structurally in-distribution, exactly like SARS-CoV-2 CoV-RDB -- NOT a
            # pending validation.
            evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="chou_recombinant_in_distribution",
            # NO scored metric exists for this cell -- the card carries catalog CENSUS counts
            # (n_resistance / n_benign) and no acc/sens/spec/n_scored. So `SCORED` would be an
            # overclaim; the slice says so plainly rather than letting a tier imply a number.
            validation_slice=(f"Chou recombinant fold-change on {genes}: catalog CENSUS of "
                              f"{row.get('n_resistance')} resistance + {row.get('n_benign')} "
                              "phenotyped-benign entries. NO acc/sens/spec exists for this cell"),
            label_provenance="Chou recombinant marker-transfer compilations (PMC3262590 / PMC5483911 / AAC 2018)",
            # NO_FREE_PHENOTYPE, not UNDERPOWERED: HCMV phenotyping is per-MUTATION marker transfer, so
            # there is no isolate-level phenotype source to be underpowered ON. The card's own
            # `independence` field records this as CLOSED for free data, not pending work.
            abstention_vocab=AbstentionVocab.NO_FREE_SOURCE, native_abstention="NO_FREE_PHENOTYPE",
            falsifier_ref="none", incoming_data_gate="G1",
            demotion_rule="a held-out clinical-isolate fold-change set would re-tier toward independent"))
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
        validation_slice="GeT-RM CYP2C8_getrm_ngs core-diplotype concordance 82/82 on real 1000G (UNCHANGED after the sentinel layer); 5 non-core samples now WITHHELD -> silent mis-calls 0/87 (the leak fully closed on this cohort)",
        label_provenance="GeT-RM consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) CYP2C8_getrm_ngs join 1000G",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="non-core *15/*16/*17/*18 now WITHHELD via sentinels (PharmVar-sourced, Ensembl-verified); function is substrate-dependent so NO PM/IM/NM is ever emitted; a non-core allele outside the 4-sentinel set would still mis-call *1 (v0.1)"),
    CellContract(
        cell_id="pgx:human:cyp3a5", track="pgx", route="dna-pgx", organism="human", target="cyp3a5",
        claim="CYP3A5 star-allele diplotype (*3/*6/*7) + CPIC expressor/non-expressor phenotype (tacrolimus) from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="calling_validated_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM CDC multi-lab consensus (CYP3A5_getrm_cons) 88/88 core-diplotype on 1000G-overlap (POWERED via the full GeT-RM Consolidated table, up from n=8); covers *1/*3/*6/*7 incl. *7 insertion",
        label_provenance="GeT-RM CDC Consolidated PGx/HLA table (363-sample; CYP3A4/CYP3A5 J Mol Diagn 2023) join 1000G",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="88 GeT-RM CYP3A5 samples now overlap 1000G (POWERED, was n=8; the underpowered flag is cleared); non-core *8/*9/*10/*11 sentinels POPULATED + Ensembl-verified + safe (no core-site collision, core 88/88 UNCHANGED) but 0 non-core carriers in this cohort so UNEXERCISED; non-core alleles outside the 4-sentinel set still mis-called *1 (v0.1)"),
    CellContract(
        cell_id="pgx:human:tpmt", track="pgx", route="dna-pgx", organism="human", target="tpmt",
        claim="TPMT COMPOUND star-allele diplotype (*3A=*3B+*3C) + CPIC thiopurine metabolizer phenotype from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="compound_calling_validated_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM CDC consolidated consensus 85/85 core-comparable on 1000G-overlap (truth *1/*3A/*3B/*3C), UNCHANGED after the sentinel layer; 10 non-core samples now WITHHELD (were silent *1) on the wider-region VCF; residual 3 silent = *6/*12/*40 NOT genotyped in the 1000G 30x panel (data limitation, not a code gap); compound *3A path exercised (6 *3A + 8 *3C samples)",
        label_provenance="GeT-RM CDC consolidated 363-sample PGx consensus (TPMT/NUDT15 J Mol Diagn 2022) join 1000G",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="first true compound-allele cell (>=2 SNPs in cis -> *3A); non-core alleles (*2/*8/*16/*40/*24/*32/*21/*12/*6/*33) now WITHHELD via the sentinel layer (PharmCAT-sourced, Ensembl-verified) rather than mis-called *1; a non-core allele NOT in the 10-sentinel set is still called *1 (v0.1 = extend the set)"),
    CellContract(
        cell_id="pgx:human:cyp2b6", track="pgx", route="dna-pgx", organism="human", target="cyp2b6",
        claim="CYP2B6 *6-proxy (516G>T signal) + CPIC efavirenz metabolizer phenotype from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="single_snp_proxy_calling_validated_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM CDC consolidated consensus 62/62 on clean *1/*6 truth on 1000G-overlap (UNCHANGED after the sentinel layer -- no false-withhold on core despite the *6-haplotype-shared 785 SNP); 18 non-core samples now WITHHELD (were silent mis-calls); SINGLE-SNP *6-proxy (516G>T) — cannot split *6/*9 (rs2279343/785A>G absent from 1000G 30x panel)",
        label_provenance="GeT-RM CDC consolidated 363-sample PGx consensus (CYP2B6) join 1000G",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="single-SNP *6-proxy: rs2279343 (785A>G) absent from the 1000G 30x panel so *6 can't be split from *9; non-core *2/*7/*18 now WITHHELD via distinctive-SNP sentinels (PharmCAT-sourced, Ensembl-verified; the *6-haplotype-shared 785 deliberately EXCLUDED to avoid false-withhold); absence-defined *4/*9 + rarer alleles outside the 3-sentinel set still mis-called (v0.1)"),
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
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="single_snp_ld_tag_validated_vs_getrm_repeat_unassessed",
        validation_slice="GeT-RM Consolidated single-SNP concordance: rs887829 (*80, *28/*37 LD-tag) 39/39 dosage-concordant + rs4148323 (*6) 39/39 on 1000G-overlap (ambiguous parenthetical/compound truth skipped); rs887829 EUR AF ~30% == the *28 frequency. The direct *28 TA-repeat is still a STRUCTURAL WALL (repeat-aware caller needed) -- the SNP validates the TAG, not the repeat length",
        label_provenance="CPIC UGT1A1 (Gammal 2016) + PharmVar; GeT-RM Consolidated table (363-sample) star-allele truth via scripts/pgx_single_snp_concordance.py; rs887829 as the *28/*37 LD-tag",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_single_snp_concordance.py", incoming_data_gate="n/a",
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
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="single_snp_validated_vs_getrm_star_truth",
        validation_slice="GeT-RM Consolidated single-SNP concordance: rs2108622 -> *3 dosage 54/54 (1.0) on 1000G-overlap (10 ambiguous skipped); plus-strand genomic C>T == cDNA 433 Val>Met; AF-corroborated (*3 ~29% EUR / ~79% EAS)",
        label_provenance="CPIC warfarin (Johnson 2017) CYP4F2*3 + GeT-RM Consolidated table star-allele truth via scripts/pgx_single_snp_concordance.py",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_single_snp_concordance.py", incoming_data_gate="n/a",
        demotion_rule="single-SNP *3 (rs2108622 IS the *3-defining variant -> 54/54 vs GeT-RM); a warfarin DOSE modifier, not a metabolizer phenotype; the dose direction is annotation only (NOT a clinical dose)"),
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
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT, claim_status="single_snp_validated_vs_getrm_star_truth",
        validation_slice="GeT-RM Consolidated single-SNP concordance: rs4149056 -> *5-family (521C: *5/*15/*17) dosage 87/87 (1.0) on 1000G-overlap (1 ambiguous skipped); plus-strand 521T>C readout",
        label_provenance="CPIC simvastatin (Cooper-DeHoff 2022) + GeT-RM Consolidated table star-allele truth via scripts/pgx_single_snp_concordance.py (rs4149056 is the 521C-defining variant of *5/*15/*17)",
        abstention_vocab=AbstentionVocab.SCORED, native_abstention="SCORED",
        falsifier_ref="scripts/pgx_single_snp_concordance.py", incoming_data_gate="n/a",
        demotion_rule="single-SNP 521C detector for *5/*15/*17 (87/87 vs GeT-RM); does NOT resolve WHICH of *5/*15/*17 (they share 521C); full star typing needs more variants (v0 scope-limit)"),
]

# --- Typing + determinant-finder whole-tool cells (route dna-<trait>). Faithful-to-tool curated-DB callers. ---
# (track, trait, organism-scope, one-line claim, native-abstention)
_TYPING_FINDER: list[tuple[str, str, str, str, str]] = [
    ("typing", "pathotype", "Escherichia_coli", "E. coli pathotype compatibility call + abstention (VirulenceFinder resolver)", "ABSTAIN"),
    # serotype MOVED to _TRAIT_CONTRACTS 2026-09-04: measured against wet-lab O:H labels, which found
    # and closed a live coverage-only selection bug. It has an individually-earned tier now.
    ("typing", "mlst", "bacteria", "multi-locus sequence type (PubMLST allele->profile->ST)", "ABSTAIN"),
    ("typing", "ktype", "Klebsiella", "Klebsiella K/O capsule type (Kaptive)", "ABSTAIN"),
    # salmserovar MOVED to _TRAIT_CONTRACTS 2026-09-04: it now has an individually-earned tier
    # (measured against a wet-lab label), and leaving it on the shared faithful-to-tool default would
    # both overstate its evidence class and HIDE that it underperforms the tool it wraps.
    # pneumoserotype MOVED to _TRAIT_CONTRACTS 2026-09-04: it was registered FAITHFUL_TO_TOOL while its
    # report card recorded INDEPENDENT Quellung validation (n=230). Under-claiming is as much a
    # trust-surface falsehood as over-claiming, so it gets its individually-earned tier.
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
            "proving nothing, since coverage holds even for a useless model. (2) FLOOR vs CEILING, both "
            "measured at scale (wiki/esm_at_scale_2026-07-17.md): the SHIPPED blosum62 default beats a random "
            "pick MATERIALLY on only 13.5% (N=200) and is often WORSE than guessing -- NOT a reliable design "
            "tool; the LEARNED method (ESM2-650M, GPU/precomputed-table) beats it on 72.9% (N=188), so the "
            "real capability is REAL and general but lives in the learned method, not the wheel-only "
            "default. Utility also does NOT track forward "
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
        claim="HIrisPlex-S visible-trait pigmentation: eye (blue/intermediate/brown, IrisPlex 6-SNP) + "
              "hair (blond/brown/red/black) + skin (very-pale..dark-black) probabilities from SNP genotypes",
        evidence_tier=EvidenceTier.FAITHFUL_TO_TOOL,
        claim_status="faithful_to_hirisplex_reference_tool_webtool_recovered_validated",
        validation_slice=(
            "EYE = the published IrisPlex model (coefficients in irisplex.py; reference anchors rs12913832 "
            "GG->blue / AA->brown via reference_integrity_ok()), POPULATION-VALIDATED on real 1000G (N=3202, "
            "Ensembl-pinned + strand-harmonized, 2026-07-29): known eye geography reproduced (EUR P(blue)=0.468, "
            "AFR/EAS/SAS brown ~1.0). HAIR (4-cat) + SKIN (5-cat) = the HIrisPlex-S deployed multinomial models "
            "RECOVERED from the erasmusmc webtool (2026-07-30) -- the papers publish the betas but NOT the "
            "intercepts (webtool-only), so the models were extracted by a designed-genotype-basis query + LS-fit "
            "and VALIDATED on 20 random held-out genotypes: max |ΔP| eye 6e-15 / hair 6e-16 / skin 9e-3, i.e. the "
            "offline models reproduce the deployed webtool to machine precision (eye/hair) / <1% (skin). "
            "Faithful-to-tool (the HIrisPlex-S webtool is the reference); population geography also confirmed on "
            "1000G. NOT per-individual concordance (openSNP, the free per-individual label source, deleted 2025-04-30)"),
        label_provenance=("EYE: IrisPlex published coefficients (Walsh 2011). HAIR/SKIN: HIrisPlex-S deployed "
                          "model recovered from hirisplex.erasmusmc.nl (user-authorized extraction; held-out-validated) "
                          "-- dna_decode/pigment/hirisplex_coefficients.json. Population geography vs 1000G; no per-individual cohort"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/pigment_1000g_hairskin_validate.py", incoming_data_gate="n/a",
        demotion_rule=(
            "European-ancestry-derived model; predictive accuracy is known to degrade off that ancestry and "
            "intermediate is the weakest class. The population geography holds on 1000G, but a PER-INDIVIDUAL "
            "measured score (needs a surviving openSNP mirror or a new consented genotype+phenotype cohort) "
            "would move this to a real independent tier — in either direction"),
    ),
    CellContract(
        cell_id="typing:dog:coatcolor", track="typing", route="dna-coatcolor",
        organism="Canis_lupus_familiaris", target="coatcolor",
        claim="dog coat colour (pigment type + eumelanin colour black/brown/blue/isabella + distribution "
              "solid/sable/agouti/tan-points) from the five classic OMIA loci E/K/A/B/D, resolved in fixed "
              "epistatic order — the first PHYSICAL/visible-trait animal cell",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="curated_epistatic_catalog_measured_black_only_substrate_limited",
        validation_slice=(
            "deterministic epistatic curated-catalog rule (Little 1957 / Schmutz & Berryere; OMIA causal loci "
            "MC1R/CBD103/ASIP/TYRP1/MLPH); reference_integrity_ok() pins breed genotypes -> colours incl. the "
            "E-locus EPISTASIS ANCHOR a naive rule mis-calls. MEASURED per-individual on the free Darwin's Ark "
            "cohort (Dryad doi:10.5061/dryad.83bk3jb4r; canFam4_gp-0.70_biallelic-SNV imputed, 3277 dogs x 29M "
            "SNVs; N=1930 owner coat colours; 2026-07-30): BLACK 160/161 = 0.994 (the eumelanin-default call is "
            "validated), blue/grey 11/31 = 0.355 (MLPH d1+d2 verified-present, partial). red/yellow + brown NOT "
            "SCORABLE — SUBSTRATE-LIMITED: the causal variants are indel/structural/low-freq (K/CBD103 delGGT, "
            "A/ASIP SINE, TYRP1 bs + MLPH d3 indels, MC1R e imputation-gap) and ABSENT from a biallelic-SNV panel; "
            "'red/yellow' is dominated by A^y (ASIP), unreachable without the SINE. Confirms the /probe: the cell "
            "is correct (black 0.994) but a SNP-only imputed substrate can't validate the full colour range. See "
            "wiki/dog_coat_darwins_ark_measured_2026-07-30.md"),
        label_provenance=("OMIA-curated causal loci (rule); Darwin's Ark/Dryad owner-reported coat colour "
                          "(measured); causal-variant coords OMIA canFam3.1 -> UCSC canFam3ToCanFam4 liftover -> .bim-verified"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/dog_coat_darwins_ark_validate.py", incoming_data_gate="G9,G10",
        demotion_rule=(
            "MEASURED but SUBSTRATE-LIMITED: black validated 0.994, other colours unscorable on the biallelic-SNV "
            "imputed panel (causal indels/SVs absent). A full-colour measured tier needs a substrate that "
            "genotypes CBD103/ASIP/TYRP1-bs directly (Embark/VGL panel or WGS), not imputed SNVs. Pattern loci "
            "(merle/spotting) ABSTAIN by design"),
    ),
    CellContract(
        cell_id="typing:dog:morphology", track="typing", route="dna-morphology",
        organism="Canis_lupus_familiaris", target="morphology",
        claim="dog body SIZE (relative rank toy/small..large/giant, additive polygenic score over "
              "IGF1/HMGA2/STC2/GHR) + EAR type (MSRB3 erect/drop) from pinned canFam4 causal SNP dosages — "
              "the quantitative/visible-trait sibling of the coat-colour cell",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="curated_pinned_catalog_measured_relative_height_and_ear",
        validation_slice=(
            "deterministic pinned + FUNCTIONALLY-VALIDATED catalog on the free Darwin's Ark cohort (Dryad "
            "doi:10.5061/dryad.83bk3jb4r; canFam4_gp-0.70_biallelic-SNV imputed, 3277 dogs x 29M SNVs; 2026-07-30): "
            "the 4-locus body-size polygenic score (IGF1/HMGA2/STC2/GHR, OMIA/lit canFam3.1 -> canFam4 liftover -> "
            ".bim-verified) tracks owner-reported height Q121 at r=+0.619 (R2=0.383, N=3276); the EAR lead MSRB3 "
            "chr10:8612500 tracks Q125 at r=+0.543 (N=2834), CLEANLY resolved from the HMGA2 body-size SNP (r=-0.13, "
            "the Morrill 2022 MSRB3-vs-HMGA2 confound). Unlike the coat indels/SVs, the body-size + ear causal SNPs "
            "ARE in-panel. RELATIVE size rank + ear axis, NOT calibrated absolute height (Q121 is a covariate-adjusted "
            "z-score); ear erect/drop NAMING is MSRB3-literature-anchored (Boyko 2010), not independently "
            "label-confirmed -> medium confidence. See wiki/dog_morphology_darwins_ark_validated_2026-07-30.md"),
        label_provenance=("OMIA/literature-curated causal loci + canFam4 gene windows (rule); Darwin's Ark/Dryad "
                          "owner-reported height Q121 + morphology Q125 (measured); coords lifted + .bim-verified"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/dog_morphology_darwins_ark_validate.py", incoming_data_gate="n/a",
        demotion_rule=(
            "MEASURED relative-signal (height polygenic r=0.619, ear r=0.543); RELATIVE rank not absolute inches. "
            "Coat length/curl (FGF5/KRT71), leg length (FGF4 SV), and the 4 covariate-adjusted rerun morph traits "
            "(Q124/127/128/245) ABSTAIN — no strong single-known-SNP mapping on this substrate (max |r|=0.21). A "
            "calibrated absolute-height tier or a label-confirmed ear polarity needs a raw-inches/codebook label"),
    ),
    CellContract(
        cell_id="typing:horse:horsecolor", track="typing", route="dna-horsecolor",
        organism="Equus_caballus", target="horsecolor",
        claim="horse coat colour (base chestnut/bay/black + cream dilution palomino/buckskin/cremello/perlino "
              "+ dun red-dun/grullo + progressive grey) from the five OMIA loci E/A/CR/D/G resolved in fixed "
              "epistatic order — a 2nd-organism visible-trait cell, the best-characterised animal coat system",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="deployed_rule_extended_validation_data_wall",
        validation_slice=(
            "deterministic epistatic curated-catalog rule; OMIA-sourced causal variants (MC1R p.Ser83Phe chestnut "
            "OMIA 001199-9796 / ASIP 11-bp-del black / SLC45A2 c.457G>A cream OMIA 001344-9796 / TBX3 dun "
            "Imsland 2016 / STX17 4.6-kb dup grey OMIA 001356-9796). The base E x A is the DEPLOYED VGL/Rieder-2001 "
            "rule (REUSES dna_decode.data.horse_coat.call_horse_base_colour); this cell EXTENDS it with cream/dun/"
            "grey. reference_integrity_ok() pins known genotypes -> colours incl. the TWO EPISTASIS ANCHORS a naive "
            "rule mis-calls: (1) e/e is CHESTNUT even when A/A bay (recessive epistasis); (2) a G/n horse GREYS out "
            "regardless of base (grey dominant + epistatic for the adult coat). KNOWLEDGE_BASELINE: the base rule's "
            "validator (scripts/horse_coat_validate.py) reports VALIDATION_DATA_WALL — no free INDEPENDENT-colour "
            "per-individual cohort (Dryad 3q111 is genotype-DERIVED=circular AND auth-gated; published contingencies "
            "Rieder/Synergy/Noma are PDF/paywalled), so no measured number (unlike the dog cells' Darwin's Ark)"),
        label_provenance=("OMIA-curated causal variants (rule); base E x A = deployed VGL/Rieder rule; no free "
                          "independent-colour per-individual cohort scored (documented data wall)"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/horse_coat_validate.py", incoming_data_gate="G10",
        demotion_rule=(
            "KNOWLEDGE_BASELINE curated catalog + validation-data wall. To reach a MEASURED tier, feed "
            "scripts/horse_coat_validate.py a TSV of mc1r,asip,INDEPENDENTLY-OBSERVED colour (not genotype-derived "
            "= circular); the free-substrate availability is the open risk (Darwin's-Ark-style browser-download "
            "wall). Cream/dun/grey extension has no validator yet. Dilution/pattern loci (champagne/silver/pearl/"
            "roan/tobiano/appaloosa) + sooty/flaxen shade ABSTAIN by design"),
    ),
    CellContract(
        cell_id="typing:cat:catcolor", track="typing", route="dna-catcolor",
        organism="Felis_catus", target="catcolor",
        claim="cat coat colour (base black/chocolate/cinnamon + dilute + X-linked orange red/cream + "
              "TORTOISESHELL/CALICO mosaic + tabby + colorpoint + white spotting/dominant white) from the OMIA "
              "loci W/O/A/B/D/C resolved in epistatic order — a 3rd-organism visible-trait cell, notable for the "
              "X-linked orange -> tortoiseshell X-inactivation mosaic",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="curated_epistatic_xlinked_catalog_no_free_validation_substrate",
        validation_slice=(
            "deterministic epistatic curated-catalog rule; OMIA-sourced causal variants incl. the 2025-identified "
            "X-linked ORANGE gene (ARHGAP36 5.1-kb intron-1 deletion, Toh/Kaelin Current Biology 2025), KIT "
            "dominant-white/spotting FERV1 (OMIA 000209/001737-9685, David 2014), ASIP agouti, TYRP1 brown, MLPH "
            "dilute, TYR albino-series (Siamese cs, OMIA 000202-9685, Lyons 2005). reference_integrity_ok() pins "
            "known genotypes -> colours incl. the THREE EPISTASIS ANCHORS a naive rule mis-calls: (1) W dominant-"
            "white masks ALL colour; (2) a female O/o is a TORTOISESHELL mosaic (X-inactivation), +white spotting "
            "= CALICO; (3) orange is EPISTATIC over brown (a b/b orange cat is red, not chocolate). The O locus is "
            "X-linked -> sex-dependent zygosity (1 allele male / 2 female). KNOWLEDGE_BASELINE — no free "
            "per-individual validation substrate"),
        label_provenance="OMIA-curated causal variants (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9,G10",
        demotion_rule=(
            "KNOWLEDGE_BASELINE curated catalog. To reach a MEASURED tier, score per-individual vs a public cat "
            "genotype+phenotype cohort (the open risk: a FREE such substrate may not exist — the Darwin's-Ark-"
            "style browser-download wall). Tabby sub-pattern (mackerel/classic/ticked), silver/inhibitor, "
            "wideband, karpati/roan + coat length ABSTAIN by design"),
    ),
    CellContract(
        cell_id="typing:chicken:plumage", track="typing", route="dna-plumage",
        organism="Gallus_gallus", target="plumage",
        claim="chicken plumage colour (eumelanin canvas extended-black/birchen/wheaten/partridge + Z-linked "
              "barring + Z-linked silver/gold + dominant white + blue/splash + lavender + recessive white) from "
              "the OMIA loci E/B/S/I/Bl/lav/c resolved in epistatic order — a 4th-organism (bird) visible-trait "
              "cell, notable for Z-LINKED loci where the FEMALE (ZW) is hemizygous (reversed from mammals)",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="curated_epistatic_zlinked_catalog_no_free_validation_substrate",
        validation_slice=(
            "deterministic epistatic curated-catalog rule; OMIA-sourced causal variants (MC1R E-locus series "
            "OMIA 000374-9031, E extended-black G274A/E92K; CDKN2A sex-linked barring OMIA 000102-9031, B1 V9D, "
            "Hellstrom 2010/Schwochow 2017; SLC45A2 silver OMIA 000370-9031, Y277C/L347M, Gunnarsson 2007; PMEL17 "
            "dominant white OMIA 000373-9031, 9-bp exon-10 insertion, Kerje 2004; Bl blue; MLPH lavender). "
            "reference_integrity_ok() pins the anchors a naive rule mis-calls: (1) EXTENSION is the canvas "
            "(barring/blue act on eumelanin, so barely show on a wheaten bird); (2) Z-LINKED barring/silver with "
            "REVERSED hemizygosity (the FEMALE is ZW-hemizygous -> 1 allele, the mirror of cat's X-linked orange); "
            "(3) dominant/recessive white mask eumelanin. KNOWLEDGE_BASELINE — no free per-individual validation substrate"),
        label_provenance="OMIA-curated causal variants (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9,G10",
        demotion_rule=(
            "KNOWLEDGE_BASELINE curated catalog. To reach a MEASURED tier, score per-individual vs a public "
            "chicken genotype+phenotype cohort (the open risk: a FREE such substrate may not exist — the Darwin's-"
            "Ark-style browser-download wall). Fine feather pattern (Columbian/mottling/pencilling/spangling), "
            "lacing, and comb/feather-structure genes ABSTAIN by design"),
    ),
    CellContract(
        cell_id="typing:rabbit:rabbitcolor", track="typing", route="dna-rabbitcolor",
        organism="Oryctolagus_cuniculus", target="rabbitcolor",
        claim="rabbit coat colour via the textbook A-E mammalian series (agouti/tan/self x black/chocolate x "
              "chinchilla/Himalayan/albino x dense/dilute x extension/steel/red) resolved epistatically",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA A-E series (A/ASIP, B/TYRP1, C/TYR, D/MLPH, E/MC1R) via the "
                          "shared mammalian-colour engine; reference_integrity_ok pins albino-masks / e-e-red-hides-"
                          "agouti / Ed-self-black. KNOWLEDGE_BASELINE — no free per-individual validation substrate"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9",
        demotion_rule=(
            "UNVALIDATABLE AS WRITTEN: none of its 5 loci record a causal variant, so no genotype file can be scored against this rule -- a cohort is NECESSARY and NOT SUFFICIENT. Curating the causal variants (OMIA/literature, per-locus, sourced) is the precondition; only then does a free cohort help. Spotting (En/Du) ABSTAIN"),
    ),
    CellContract(
        cell_id="typing:mouse:mousecolor", track="typing", route="dna-mousecolor",
        organism="Mus_musculus", target="mousecolor",
        claim="mouse coat colour via the foundational pigment loci (agouti x brown x albino x dilute x pink-eyed "
              "dilution x extension) resolved epistatically — the century-old model of mammalian colour genetics",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (A/ASIP, B/Tyrp1, C/Tyr, D/Myo5a, P/Oca2, E/Mc1r) via the "
                          "shared engine; reference_integrity_ok pins albino-masks / e-e-yellow. KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9",
        demotion_rule=(
            "UNVALIDATABLE AS WRITTEN: none of its 6 loci record a causal variant, so no genotype file can be scored against this rule -- a cohort is NECESSARY and NOT SUFFICIENT. Curating the causal variants (OMIA/literature, per-locus, sourced) is the precondition; only then does a free cohort help. Spotting (s/piebald) ABSTAIN"),
    ),
    CellContract(
        cell_id="typing:cattle:cattlecolor", track="typing", route="dna-cattlecolor",
        organism="Bos_taurus", target="cattlecolor",
        claim="cattle coat colour via MC1R Extension (dominant-black/wild/recessive-red) + incompletely-dominant "
              "PMEL/SILV dilution (Charolais/Highland -> dun/silver), resolved epistatically",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (MC1R ED>E+>e; PMEL Dc/Dh dosage dilution) via the shared "
                          "engine; reference_integrity_ok pins ED-dominant-black / e-e-red / PMEL-incomplete-dominant. "
                          "KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9",
        demotion_rule=(
            "UNVALIDATABLE AS WRITTEN: none of its 2 loci record a causal variant, so no genotype file can be scored against this rule -- a cohort is NECESSARY and NOT SUFFICIENT. Curating the causal variants (OMIA/literature, per-locus, sourced) is the precondition; only then does a free cohort help. Roan (KITLG)/spotting (MITF)/COPA-dom-red ABSTAIN"),
    ),
    CellContract(
        cell_id="typing:pig:pigcolor", track="typing", route="dna-pigcolor",
        organism="Sus_scrofa", target="pigcolor",
        claim="pig coat colour via KIT Dominant-White (epistatic, masks all) + MC1R Extension (dominant-black/wild/"
              "recessive-red), resolved epistatically",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (KIT dominant-white; MC1R OMIA 001199-9823 ED>E+>e) via the "
                          "shared engine; reference_integrity_ok pins KIT-masks / ED-black / e-e-red. KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9",
        demotion_rule=(
            "UNVALIDATABLE AS WRITTEN: none of its 2 loci record a causal variant, so no genotype file can be scored against this rule -- a cohort is NECESSARY and NOT SUFFICIENT. Curating the causal variants (OMIA/literature, per-locus, sourced) is the precondition; only then does a free cohort help. KIT belt/patch/roan sub-alleles ABSTAIN"),
    ),
    CellContract(
        cell_id="typing:sheep:sheepcolor", track="typing", route="dna-sheepcolor",
        organism="Ovis_aries", target="sheepcolor",
        claim="sheep coat colour via ASIP Agouti (A^Wt dominant white/tan from a 190kb duplication > a recessive "
              "black) + MC1R Extension (dominant-black overrides ASIP white), resolved epistatically",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (ASIP OMIA 000201-9940 dominant-white-duplication vs "
                          "recessive-black LOF; MC1R ED M73K/D121N) via the shared engine; reference_integrity_ok pins "
                          "ED-overrides-ASIP-white / dominant-white-tan / recessive-black. KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9",
        demotion_rule=(
            "UNVALIDATABLE AS WRITTEN: none of its 2 loci record a causal variant, so no genotype file can be scored against this rule -- a cohort is NECESSARY and NOT SUFFICIENT. Curating the causal variants (OMIA/literature, per-locus, sourced) is the precondition; only then does a free cohort help. Badgerface/spotting ASIP sub-alleles ABSTAIN"),
    ),
    CellContract(
        cell_id="typing:goat:goatcolor", track="typing", route="dna-goatcolor",
        organism="Capra_hircus", target="goatcolor",
        claim="goat coat colour via ASIP Agouti (A^Wt dominant white/tan, the CNV-driven many-pattern hub > a "
              "recessive nonagouti black) + TYRP1 brown, resolved epistatically",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (ASIP OMIA 000201-9925 dominant-white/tan CNV vs recessive-"
                          "black; TYRP1 brown) via the shared engine; goat MC1R association is incomplete in the "
                          "literature so ASIP is the modeled driver. KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9,G10",
        demotion_rule="KNOWLEDGE_BASELINE; MEASURED needs a free goat cohort. The ~11 ASIP pattern sub-alleles (badgerface/swiss/grey) ABSTAIN",
    ),
    CellContract(
        cell_id="typing:alpaca:alpacacolor", track="typing", route="dna-alpacacolor",
        organism="Vicugna_pacos", target="alpacacolor",
        claim="alpaca/llama fleece colour via MC1R Extension (E coloured; e/e = recessive WHITE regardless of ASIP — "
              "the camelid twist) + ASIP Agouti (A functional -> fawn/agouti > a loss-of-function -> black)",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (MC1R E/e camelid recessive-white; ASIP fawn-vs-black LOF) "
                          "via the shared engine; reference_integrity_ok pins ee-recessive-white (white regardless of "
                          "ASIP) / black-if-aa / fawn-if-A. KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9",
        demotion_rule=(
            "UNVALIDATABLE AS WRITTEN: none of its 2 loci record a causal variant, so no genotype file can be scored against this rule -- a cohort is NECESSARY and NOT SUFFICIENT. Curating the causal variants (OMIA/literature, per-locus, sourced) is the precondition; only then does a free cohort help. KIT grey/blue-eyed-white (open question) ABSTAINs"),
    ),
    CellContract(
        cell_id="typing:guineapig:guineapigcolor", track="typing", route="dna-guineapigcolor",
        organism="Cavia_porcellus", target="guineapigcolor",
        claim="guinea pig coat colour via the classic A/B/C/D/E series (ASIP agouti/non-agouti x TYRP1 brown x "
              "TYR x MLPH dilute x MC1R recessive-red), resolved epistatically",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (ASIP OMIA 000201-10141 c.181delTTCA non-agouti; MC1R "
                          "OMIA 001199-10141 e recessive-red Vidal 2018; TYRP1/TYR/MLPH classical) via the shared "
                          "engine. KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9,G10",
        demotion_rule="KNOWLEDGE_BASELINE; MEASURED needs a free guinea pig cohort. White spotting (OMIA 000214) ABSTAINs",
    ),
    CellContract(
        cell_id="typing:fox:foxcolor", track="typing", route="dna-foxcolor",
        organism="Vulpes_vulpes", target="foxcolor",
        claim="red/silver fox coat colour via the NON-epistatic silver-fox system — MC1R EA Alaska-Silver "
              "dominant-black + ASIP wild-red vs a Standard-Silver recessive-black (both routes reach dark)",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (MC1R EA Cys125Arg gain-of-function; ASIP OMIA 000201-9627 "
                          "Standard-Silver 166-nt exon-1 deletion; non-epistatic, Vage 1997) via the shared engine; "
                          "reference_integrity_ok pins Alaska + Standard silver both -> black via different loci, "
                          "wild red. KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9,G10",
        demotion_rule="KNOWLEDGE_BASELINE; MEASURED needs a free fox cohort. Platinum (KIT) + other farm morphs ABSTAIN",
    ),
    CellContract(
        cell_id="typing:donkey:donkeycolor", track="typing", route="dna-donkeycolor",
        organism="Equus_asinus", target="donkeycolor",
        claim="donkey coat colour via MC1R (e recessive-red) + ASIP (light-points/grey-dun vs solid-black-no-light-"
              "points) + TYR (Asinara white albinism), resolved epistatically",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (MC1R OMIA 001199-9793 e c.629T>C Abitbol 2014; ASIP OMIA "
                          "000201-9793 c.349T>C light-points Sun 2017; TYR Asinara albino c.604C>G) via the shared "
                          "engine. ASIP heterozygote variability -> one residual gene unidentified (documented). "
                          "KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="n/a",
        demotion_rule="KNOWLEDGE_BASELINE; MEASURED needs a free donkey cohort. KIT dominant-white/SLC45A2-cream/the residual ASIP-het gene ABSTAIN",
    ),
    CellContract(
        cell_id="typing:buffalo:buffalocolor", track="typing", route="dna-buffalocolor",
        organism="Bubalus_bubalis", target="buffalocolor",
        claim="water buffalo coat colour via ASIP A^W DOMINANT white (a 2809-bp LINE-1 insertion causing 10x ASIP "
              "overexpression) > a black; MC1R is monomorphic in buffalo so ASIP is the sole driver",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA locus (ASIP OMIA 000213-89462 LINE-1 dominant-white, Liang 2020 "
                          "-- convergent with cattle; MC1R monomorphic per Cruz 2020) via the shared engine. "
                          "KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal gene (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G10",
        demotion_rule="KNOWLEDGE_BASELINE; MEASURED needs a free buffalo cohort. KIT white-spotting (OMIA 001737) + the disputed ASIP-black SNP ABSTAIN",
    ),
    CellContract(
        cell_id="typing:camel:camelcolor", track="typing", route="dna-camelcolor",
        organism="Camelus_dromedarius", target="camelcolor",
        claim="dromedary camel coat colour via MC1R c.901C>T DOMINANT white (dominant-negative, heterozygote "
              "sufficient) + ASIP recessive black (exon-2 frameshift); wild = light brown",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (MC1R c.901C>T dominant-white, Almathen 2018/Alshanbari "
                          "2019; ASIP 23delT exon-2 frameshift recessive-black) via the shared engine. KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G10",
        demotion_rule="KNOWLEDGE_BASELINE; MEASURED needs a free camel cohort. KIT white-spotting + SLC45A2/TYR modifiers ABSTAIN",
    ),
    CellContract(
        cell_id="typing:mink:minkcolor", track="typing", route="dna-minkcolor",
        organism="Neogale_vison", target="minkcolor",
        claim="American mink coat colour via TYR albino/Himalayan + TYRP1 American-Palomino brown + MLPH "
              "Silverblue dilute; wild = dark",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA loci (TYR albino nonsense OMIA 000202-452646; TYRP1 Palomino "
                          "intron-2 insertion; MLPH Silverblue splice c.901+1G>A OMIA 000031-452646, Manakhov 2019) "
                          "via the shared engine. KNOWLEDGE_BASELINE"),
        label_provenance="OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9,G10",
        demotion_rule="KNOWLEDGE_BASELINE; MEASURED needs a free mink cohort. LYST Aleutian + MITF Hedlund-white + the 30+ other fur colours ABSTAIN",
    ),
    CellContract(
        cell_id="typing:roedeer:roedeercolor", track="typing", route="dna-roedeercolor",
        organism="Capreolus_capreolus", target="roedeercolor",
        claim="roe deer coat colour via ASIP c.33G>T p.Leu11Phe — A chestnut (phaeomelanin) > a recessive black",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE, claim_status="curated_epistatic_catalog_shared_engine",
        validation_slice=("deterministic curated OMIA locus (ASIP c.33G>T OMIA 000201-9858, Reissmann 2020: TT black "
                          "/ GG-GT chestnut) via the shared engine. KNOWLEDGE_BASELINE — a WILDLIFE cell"),
        label_provenance="OMIA-curated causal gene (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="n/a",
        demotion_rule="KNOWLEDGE_BASELINE; the single confirmed roe-deer colour variant is ASIP (chestnut vs black); other coat variation ABSTAINs",
    ),
    CellContract(
        cell_id="typing:pigeon:pigeoncolor", track="typing", route="dna-pigeoncolor",
        organism="Columba_livia", target="pigeoncolor",
        claim="pigeon plumage colour (base ash-red/blue/brown via the Z-linked B/TYRP1 series + recessive-red "
              "SOX10 + dilute SLC45A2 + wing pattern NDP T-check/checker/bar/barless) resolved epistatically — "
              "one of the best-characterised colour systems in any organism (Shapiro lab); a 2nd BIRD cell",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="curated_epistatic_zlinked_catalog_no_free_validation_substrate",
        validation_slice=(
            "deterministic epistatic curated-catalog rule; MOLECULARLY-CONFIRMED causal genes (TYRP1 B-locus "
            "ash-red/blue/brown, Domyan 2014 Curr Biol VAAST; SOX10 recessive-red; SLC45A2 dilute; NDP wing-pattern "
            "T-check>checker>bar>barless, Vickrey 2018 eLife). reference_integrity_ok() pins the anchors a naive "
            "rule mis-calls: (1) SOX10 e/e -> RED regardless of the TYRP1 base (epistatic); (2) Z-LINKED B/dilute "
            "with REVERSED hemizygosity (the FEMALE is ZW-hemizygous, same as chicken). KNOWLEDGE_BASELINE — no "
            "free per-individual validation substrate"),
        label_provenance="Shapiro-lab/OMIA-curated causal genes (rule only); no measured per-individual cohort scored",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="G9",
        demotion_rule=(
            "UNVALIDATABLE AS WRITTEN: none of its 4 loci record a causal variant, so no genotype file can be "
            "scored against this rule -- a cohort is NECESSARY and NOT SUFFICIENT. Curating the causal variants "
            "(OMIA/literature, per-locus, sourced) is the precondition; only then does a free cohort help (and the "
            "open risk remains that a FREE pigeon genotype+phenotype substrate may not exist). Modifiers (spread/"
            "grizzle/almond/indigo) + shade ABSTAIN by design"),
    ),
    CellContract(
        cell_id="typing:bacteriophage:phage", track="typing", route="dna-phage",
        organism="bacteriophage", target="phage",
        claim="bacteriophage host-RECEPTOR class (FhuA/BtuB/LPS_core/ECA/NfrA/LptD/...) from a phage genome "
              "(genome-homology transfer) or NCBI lineage (catalogue lookup) — the first non-AMR cell",
        # INDEPENDENT now: scored on a DIFFERENT LAB's measured receptors (LBNL/Arkin-Mutalik Phage Datasheets,
        # non-Bas isolates disjoint from the BASEL reference, K-12 BW25113 host). Covered classes only.
        evidence_tier=EvidenceTier.INDEPENDENT_MEASURED,
        claim_status="independent_measured_covered_classes_only",
        validation_slice=(
            "INDEPENDENT (LBNL/Arkin-Mutalik Phage Datasheets, github.com/mjohnson11/PhageDataSheets; measured by "
            "genome-wide genetic screens on K-12 BW25113; non-Bas isolates disjoint from the BASEL reference): on "
            "the classes the v0 catalogue COVERS (BtuB/LPS_core/LptD/ECA/NfrA) = 25/29 called = 0.862 (BtuB 22/26, "
            "LPS_core 3/3). In-distribution LOO on the 29 clade-conserved BASEL phages = 27/27. Overall independent "
            "= 25/86=0.291 because v0 does NOT model the RBP-variable classes (Tsx/OmpC/FhuA/OmpA/... 60+ phages) — "
            "out of v0 scope (the RBP-caller target), not a catalogue error. RBP-level caller (--rbp-fasta): "
            "within-LBNL LOO 0.975 BUT cross-lab (phageReceptor, Zhang 2020, independent) only 0.364 (4/11) — the "
            "within-study number does NOT generalize; T4-like phages mis-transfer to Tsx relatives (wiki/phage_rbp_crosslab_result_2026-07-24)"),
        label_provenance=(
            "BASEL catalogue (Maffei 2021 PLOS Biology 3001424) for the rule; INDEPENDENT test labels from LBNL "
            "Phage Datasheets (Moriniere et al.; measured on K-12 BW25113)"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/lbnl_independent_validate.py", incoming_data_gate="n/a",
        demotion_rule=(
            "the independent number is on COVERED classes only (0.862); expanding to the RBP-variable classes "
            "needs the RBP-level caller. A drop on a larger independent measured set would re-tier"),
    ),
    CellContract(
        cell_id="typing:klebsiella:kleb", track="typing", route="dna-kleb",
        organism="Klebsiella_pneumoniae", target="kleb",
        claim="Klebsiella phage depolymerase (enzymatic domain) -> host capsule KL-type, ranked top-K "
              "(cross-organism transfer of the E. coli phage-receptor paradigm) — FETCH-ONLY (no bundled data)",
        # in-distribution: prophage-LCA labels, clonality-corrected LOO; NOT independent wet-lab.
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="in_distribution_dpotropisearch_prophage_lca_labels",
        validation_slice=(
            "clonality-corrected leave-one-out (greedy-rep @0.90) on DpoTropiSearch depolymerase domains: "
            "top-1 ~0.45 / top-5 ~0.60 over 147-165 KL-types, lift +0.49 over a 0.10 prior null; the paradigm "
            "GENERALIZES cross-organism on modular depolymerase domains (harder problem than E. coli receptors, "
            "higher number). In-distribution (prophage-LCA labels), NOT independent wet-lab"),
        label_provenance=(
            "DpoTropiSearch (Concha-Eloko et al., Nat Commun 2025; Zenodo 10.5281/zenodo.14065540) — "
            "prophage-host-LCA-inferred KL-type; data NOT bundled (fetch-only; CC-BY record / repo "
            "Decapsulate Non-Commercial License v1.1 — user verifies their use)"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/klebsiella_depolymerase_crossorganism.py", incoming_data_gate="n/a",
        demotion_rule=(
            "an INDEPENDENT wet-lab depolymerase->KL-type test set (the 63 exp_validated set, currently "
            "blocked:external) would re-tier off in-distribution; a finer clonality collapse would move the number"),
    ),
    CellContract(
        cell_id="essentiality:any:essentiality", track="typing", route="dna-essentiality",
        organism="any", target="essentiality",
        claim="single-gene KO -> essential/non-essential via the deterministic conserved-core FUNCTION "
              "catalogue (translation/replication/transcription/envelope/division); label-independent, offline",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="conserved_core_validated_vs_gold_standard_in_distribution",
        validation_slice=(
            "E. coli AUROC 0.695 genome-wide vs the Goodall 2018 mBio TraDIS gold-standard (base rate 9.3%, "
            "sens 0.373 / spec 0.984 -- high-precision, conservative-recall); composition matches the known "
            "essentialome (208/4318, translation/envelope/replication-dominated). Cross-organism transfer to "
            "human (BAGEL CEG2/NEG) AUROC 0.580. The learned E3 complement lifts it (E. coli 0.795 / human 0.911)"),
        label_provenance=(
            "gold-standard essentiality: Goodall 2018 mBio Table S1 (E. coli TraDIS genome-wide, CC-BY) + "
            "BAGEL CEGv2/NEGv1 (human core-essential/non-essential reference, Hart lab). Free, independent screens"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/essentiality_e3_learned.py", incoming_data_gate="n/a",
        demotion_rule=(
            "the conserved-core is a deterministic PRIOR (high-precision, ~0.37 recall); the E3 learned complement "
            "is the accuracy tier. An independent per-organism essentiality screen would re-validate the transfer"),
    ),
    CellContract(
        cell_id="metabolic:escherichia_coli:metabolic", track="typing", route="dna-metabolic",
        organism="escherichia_coli", target="metabolic",
        claim="E. coli carbon-source utilization from gene/operon presence via the UPTAKE-GATED rule "
              "(utilizes iff catabolic enzymes present AND a transporter present AND transporter expressed "
              "under the O2 condition); the uptake-gate is what a naive has-the-genes rule misses; offline",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="curated_catalog_validated_vs_measured_k12_phenotypes",
        validation_slice=(
            "validated vs measured E. coli K-12 MG1655 phenotypes (EcoCyc / Neidhardt): lac+ ara+ mal+ xyl+ "
            "rha+ glc+, and the CITRATE anchor (Blount 2012 Nature LTEE) -- Cit- aerobic / Cit+ anaerobic, the "
            "case a naive has-the-genes rule mis-calls positive (K-12 carries the full TCA + citT yet the citT "
            "importer is anaerobic-only). Reads gene PRESENCE, not sequence integrity"),
        label_provenance=(
            "curated operon/transporter assignments from EcoCyc + Neidhardt (E. coli textbook physiology) + "
            "Blount et al. 2012 Nature for the LTEE Cit+ aerobic-citT anchor. Faithful-to-literature, not a new model"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="dna_decode/metabolic/carbon_catalog.py", incoming_data_gate="n/a",
        demotion_rule=(
            "v0 is presence-based E. coli carbon catabolism; a genome-mode that reads sequence integrity (not just "
            "gene presence) + cross-organism transfer + a measured Biolog/BV-BRC phenotype cohort would re-tier it"),
    ),
    CellContract(
        cell_id="motility:escherichia_coli:flagellar", track="typing", route="dna-motility",
        organism="escherichia_coli", target="motility",
        claim="flagellar SWIMMING motility from gene presence -- the first NON-metabolic trait catalog. "
              "MOTILE iff all 5 flagellar modules present (master flhDC / sigma-28 fliA / flagellin fliC-fljB / "
              "motor motAB / basal-body-export fliF-fliG-flhA-fliI); chemotaxis reported separately (a "
              "che-mutant still swims -- gating swimming on chemotaxis would be a biology error); offline",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="curated_catalog_validated_vs_literature_anchors",
        validation_slice=(
            "curated flagellar catalog vs literature-known anchors: E. coli K-12 MG1655 + Salmonella "
            "Typhimurium MOTILE (all modules) vs Shigella flexneri NON-motile (flagellar pseudogenes/deleted) "
            "+ any flhDC/fliC/motAB knockout non-motile. Presence-based DIRECTION (swim/no-swim), not speed"),
        label_provenance=(
            "curated flagellar-regulon assignments (class-1 flhDC master / class-2/3 cascade) from E. coli "
            "flagellar biology (Chevance & Hughes; EcoCyc). Faithful-to-literature, not a new model. No free "
            "genome-keyed swim-plate cohort is fetchable -> KNOWLEDGE_BASELINE, not a measured-cohort tier"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="dna_decode/motility/flagellar_catalog.py", incoming_data_gate="n/a",
        demotion_rule=(
            "v0 is presence-based E. coli flagellar swimming; a sequence-integrity genome-mode (to catch a "
            "present-but-IS-disrupted flhD, the K-12 non-motile case) + cross-organism transfer + a measured "
            "swim-plate/motility cohort would re-tier it. Twitching/gliding/swarming are separate traits (out)"),
    ),
    CellContract(
        cell_id="fba:escherichia_coli:growth_essentiality", track="typing", route="dna-fba",
        organism="escherichia_coli", target="fba",
        claim="a gene edit (single/double KO) -> a QUANTITATIVE cell-level trait (growth rate /h + "
              "essential/non-essential) via mechanistic flux-balance analysis over the iML1515 genome-scale "
              "model; GENERAL over all 1515 model genes (not a curated list); computes from stoichiometry + "
              "known biochemistry, so it sidesteps population-structure confounding by construction",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="mechanistic_model_validated_vs_keio_in_distribution",
        validation_slice=(
            "genome-wide single-gene-deletion essentiality (glucose M9 aerobic) vs the free Keio-collection "
            "mutant-fitness gold standard (Bernstein 2023 method, fitness<-2 = essential-on-glucose). Metrics on "
            "the assayable gene set + corroboration that FBA-essential genes have no viable Keio mutant. See "
            "wiki/fba_keio_validation_2026-08-03. METABOLIC traits only -- NOT virulence/regulation"),
        label_provenance=(
            "iML1515 genome-scale model (Monk 2017, Nat Biotechnol, free on BiGG) + Keio BW25113 RB-TnSeq mutant "
            "fitness (Baba 2006 / Wetmore-Price fitness browser, as vendored by Bernstein 2023). Free, published"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/fba_keio_validate.py", incoming_data_gate="n/a",
        demotion_rule=(
            "v0 is FBA on iML1515, glucose M9, in-distribution vs a knowledge baseline; a provenance-independent "
            "measured growth-rate/essentiality cohort (Biolog / a fresh TnSeq) + non-metabolic traits would re-tier it"),
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
            "|Spearman| MEDIAN 0.20 (TEM-1's 0.35 is top-13%, NOT typical); the LEARNED esm2 method is "
            "|Spearman| MEDIAN 0.49 at scale (N=194, wiki/esm_at_scale_2026-07-17.md) -- 2.5x the default. "
            "Python-API `esm2` "
            "(ESM2-650M masked-marginal): TEM-1 **0.7315** / PTEN 0.518 / CcdB 0.5115; `alphamissense` PTEN "
            "0.539 (human-only). Genome-level nucleotide-edit path validated end-to-end on a real blaTEM CDS: "
            "**Spearman 0.7611** over 1,715 real single-nt-accessible variants "
            "(wiki/blatem_genome_demo_2026-07-14.json). Dosage head: conformal coverage calibrated 10/10 "
            "proteins, informative 7/10 (wiki/forward_dosage_sweep_2026-07-15.md). "
            "LEAKAGE-FREE (the strongest evidence class this cell has; all the above is IN-DISTRIBUTION on "
            "the benchmark the methods were tuned against): on MaveDB assays whose gene is NOT in "
            "ProteinGym, esm2 median |Spearman| **0.478** over 2383 held-out assays (0.492 on 978 human; "
            "wiki/mavedb_full_esm2_2026-07-22.md), and alphamissense **0.502** over 57 held-out human "
            "assays (wiki/mavedb_am_holdout_2026-07-23.json). Held-out MODALITY comparison at scale (N=76 "
            "structurally-aligned, Kaggle T4, wiki/mavedb_holdout_hybrid_2026-07-23.md): ProSST (structure) "
            "median **0.596** > esm2 0.538, and the ESM2+ProSST `hybrid` (median 0.602) BEATS BOTH "
            "components PAIRED (70/76 vs esm2 +0.063; 52/76 vs ProSST +0.011, sign-test p=0.0009 -- "
            "SIGNIFICANT, confirmed by doubling N from 38; read the PAIRED delta, not the medians). "
            "So the deployable METHOD RANKING on held-out data is hybrid > prosst > esm2 > alphamissense "
            "> blosum62, while the CLI default stays blosum62 because it is the only wheel-only, "
            "no-model, no-structure option. ADDING EVOLUTION (the GEMME 3rd modality via the finalized "
            "Docker toolchain, TEM-1 0.719): on the held-out GEMME-covered subset (N=25) the 3-way "
            "ESM2+GEMME+ProSST beats the 2-way 21/25 (sign-p=0.0005) + beats GEMME-alone 22/25 -- adding "
            "evolution LIFTS the hybrid, and the three modalities (sequence/evolution/structure) are "
            "complementary (wiki/gemme_threeway_holdout_2026-07-23.md)"),
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
    CellContract(
        cell_id="typing:Salmonella:salmserovar", track="typing", route="dna-salmserovar",
        organism="Salmonella", target="salmserovar",
        claim="Salmonella enterica serovar from the antigenic formula (O:H1:H2) via blastn over the "
              "antigen allele DB + White-Kauffmann-Le Minor lookup",
        # The tier records EVIDENCE CLASS, not performance: this cell HAS now been measured against a
        # free, independent, wet-lab label. The headline carries the (poor) number.
        evidence_tier=EvidenceTier.INDEPENDENT_MEASURED,
        claim_status="measured_vs_wetlab_label_UNDERPERFORMS_the_tool_it_wraps",
        validation_slice=(
            "N=200 NCBI-PD isolates from reference public-health labs (CDC/PHE/FDA/USDA-FSIS/state "
            "health depts), 74 distinct serovars, 29 BioProjects, largest-source share 0.125 (CLEARS "
            "the project's own 0.60 diversity bar). Scored 2026-09-04 by "
            "scripts/salmserovar_validate.py with equivalence decided by "
            "dna_decode.salmserovar.equivalence (notation normalisation + the committed W-K-L formula "
            "table), applied IDENTICALLY to both callers. RESULT: ours 0.702 (99 hit / 42 miss / 59 "
            "no-call) vs the in-silico incumbent 0.925 (184/15/1) on the SAME isolates -- DELTA "
            "-0.222, and a 29.5% abstention rate. Diagnosed failure modes: phase-2 flagellin (H2) "
            "undetected in 33 of 59 no-calls (Salmonella is diphasic; '4:i:-' cannot resolve where the "
            "table wants '4:i:1,2'), O-antigen unresolved ('O?', e.g. Infantis/Rissen) or mis-grouped "
            "(Typhi called 1,3,19 not 9,12), plus a malformed DB antigen name ('22-gene2')"),
        label_provenance=(
            "NCBI-PD submitter `serovar` (traditional Kauffmann-White slide agglutination is the gold "
            "standard for this trait) restricted to reference labs; the in-silico comparator is PD's "
            "own `computed_types`. LABEL INDEPENDENCE IS CORROBORATED, not assumed: the incumbent "
            "scores 0.925 against these labels, matching published in-silico-vs-agglutination accuracy "
            "(~0.95) -- had the labels been copied from the tool it would score ~1.000"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/salmserovar_validate.py", incoming_data_gate="n/a",
        demotion_rule=(
            "ALREADY DEMOTED BY MEASUREMENT relative to its wrapper claim: it is WORSE than naive use "
            "of the tool it mimics (-0.222) and abstains on ~30%. Do NOT present it as a drop-in "
            "serovar caller. Residual circularity is bounded not eliminated (per-isolate agglutination "
            "provenance is unprovable), but both callers are scored on the SAME labels so the DELTA "
            "survives contamination even where the absolute levels are optimistic. Re-score after any "
            "H2/O-antigen DB fix; promote the claim only if the delta closes"),
    ),
    CellContract(
        cell_id="typing:Escherichia_coli:serotype", track="typing", route="dna-serotype",
        organism="Escherichia_coli", target="serotype",
        claim="E. coli O:H serotype from blastn over the SerotypeFinder allele DB",
        evidence_tier=EvidenceTier.INDEPENDENT_MEASURED,
        claim_status="measured_vs_wetlab_OH_labels_and_a_live_selection_bug_was_found_and_fixed",
        validation_slice=(
            "N=150 NCBI-PD isolates with a strict O:H label, 55 O-groups, 33 BioProjects, "
            "largest-source share 0.167 (clears the project's own 0.60 bar). Scored 2026-09-04 by "
            "scripts/serotype_oh_validate.py, O and H axes SEPARATELY (different denominators; never "
            "pooled). THE RUN FOUND A LIVE DEFECT: the caller selected alleles by COVERAGE ONLY -- the "
            "same bug the sibling Salmonella caller had already fixed with identity-primary selection, "
            "never propagated -- and fliC flagellin alleles cross-hybridize at near-full coverage, so "
            "the wrong H antigen won confidently. Confirmed by CODE INSPECTION before the fix, not "
            "inferred from numbers. Same 150 isolates, before -> after: H accuracy 0.770 -> 0.926 "
            "(misses 34 -> 11), O accuracy 0.931 -> 0.962, resolution UNCHANGED (the rule picks which "
            "allele, not whether). Misses were concentrated not diffuse (H21->H8 was 9 of 34 = 26%), "
            "which is what pointed at systematic allele confusion"),
        label_provenance=(
            "NCBI-PD submitter `serovar` parsed to a strict O:H shape; E. coli O:H serotyping is "
            "traditionally slide agglutination. WEAKER PROVENANCE THAN THE SALMONELLA CELL: PD does "
            "NOT populate `computed_types` for E. coli, so there is no in-silico incumbent -- neither "
            "as a comparator nor as the circularity probe the Salmonella run used"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="O?",
        falsifier_ref="scripts/serotype_oh_validate.py", incoming_data_gate="n/a",
        demotion_rule=(
            "The ABSOLUTE accuracy is weakly anchored (no incumbent); what is solid is the "
            "internally-controlled rule comparison. REPLICATED 2026-09-04 on 250 HELD-OUT isolates "
            "(seed 77, 31 discovery-overlapping accessions excluded, 47 BioProjects) with the "
            "prediction registered BEFORE the run: H accuracy 0.8171 -> 0.9228, gain +0.106, clearing "
            "the pre-registered +0.05 bar; no-call delta EXACTLY 0.0000, corroborating the mechanism "
            "(the rule picks WHICH allele wins, not WHETHER one does). QUOTE +0.106, NOT the "
            "discovery +0.155 -- the gain shrank by 0.050 on held-out data exactly as an unblinded "
            "choice predicts. The O axis is very slightly WORSE (-0.0086, 2 isolates); net is 26 H "
            "misses fixed vs 2 O misses introduced. Re-score if the allele DB changes. The SAME coverage-primary "
            "pattern is present in pneumoserotype and plasmid and was deliberately NOT fixed there "
            "(no validation cohort, and different biology: whole-locus cps / Inc replicon matching is "
            "not per-antigen flagellin) -- do not propagate the fix without measuring it"),
    ),
    CellContract(
        cell_id="typing:Streptococcus_pneumoniae:pneumoserotype", track="typing",
        route="dna-pneumo-serotype", organism="Streptococcus_pneumoniae", target="pneumoserotype",
        claim="pneumococcal capsular serotype from blastn of the assembly against cps reference loci",
        # Was registered FAITHFUL_TO_TOOL while its report card recorded independent Quellung
        # validation. The tier records EVIDENCE CLASS; the headline carries the resolution limit.
        evidence_tier=EvidenceTier.INDEPENDENT_MEASURED,
        claim_status="serogroup_validated_vs_wetlab_quellung_exact_serotype_is_a_v0_ceiling",
        validation_slice=(
            "GPS Poland cohort, wet-lab phenotypic serotype (Nat Commun 2025 Supplementary Data 1, "
            "`Phenotypic_serotype`), scored on GPS-deposited ENA assemblies -- label AND assembly both "
            "independent of this caller. 260 total = 230 SCORED + 25 assembly-unavailable + 5 no-call. "
            "On the 230: SEROGROUP concordance 0.939, EXACT serotype 0.661. Explicit-QUELLUNG-method "
            "subset n=42: serogroup 0.952 / exact 0.690 (consistent). The honest headline is the "
            "SEROGROUP number -- exact-serotype misses are systematically WITHIN-serogroup (9A/9V, "
            "6B/6E, 15B/15C), which is a single-best-reference v0 ceiling (the full GPS pipeline does "
            "allele-level within-serogroup resolution this v0 does not), NOT a bug. Selection-rule "
            "probe 2026-09-04 (scripts/pneumo_selection_rule_probe.py, 25 freshly-fetched assemblies): "
            "coverage-primary vs identity-primary flips 1 of 25 calls and that flip is wrong under "
            "BOTH orderings -- no evidence the E. coli fix transfers here, consistent with a different "
            "cause"),
        label_provenance=(
            "wet-lab phenotypic serotype (Quellung-explicit for n=42; method unspecified but still "
            "serology for the rest) -- NOT the in-silico Monocle field. The label is independent, but "
            "the cps DB and the serotype universe are a shared REFERENCE SYSTEM: reference-coupled, "
            "not circular"),
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="scripts/pneumo_gps_quellung_validate.py", incoming_data_gate="n/a",
        demotion_rule=(
            "Quote the SEROGROUP number (0.939) as this v0's resolution; exact-serotype 0.661 is the "
            "within-serogroup-limited LOWER bound and must not be presented as the cell's accuracy. "
            "~2.1% no-call rate is excluded from accuracy by construction (a utility fact, not an "
            "error). Promote only on allele-level within-serogroup typing (a v0.1). The cached GPS "
            "assemblies on D: are CORRUPT (HTTP 403 pages / truncated gzip) -- re-fetch with "
            "validation, never reuse them"),
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
    from dna_decode.data.routable_drugs import all_routable_amr_drugs
    # SHARED with the CLI's argparse choices (2026-09-01). This union used to be spelled out here as
    # well, and it drifted: HCMV's five drugs were added to the CLI and missed here, so they never
    # entered the routable set and the coverage test below could not notice they had no contracts.
    amr_drugs = all_routable_amr_drugs()
    per_target = {
        "dna-amr": {d.lower() for d in amr_drugs},
        "dna-pgx": set(PGX_GENES),
        "dna-clinvar": {"germline_pathogenicity"},  # the Mendelian (ClinVar) single-decoder route
        "dna-hla": set(HLA_ALLELES),                 # HLA drug-hypersensitivity tag-SNP cells
    }
    # "traits" = the WHOLE-TOOL (typing/finder) traits: every routable trait that does NOT already have a
    # per-target route key above. DERIVED from `per_target`, never hand-listed -- it used to read
    # `set(TRAITS) - {"amr", "pgx"}`, which silently went wrong the moment clinvar + hla became routable
    # traits (2026-08-23): both have their own per-target key here, yet both landed in "traits" and broke
    # the typing/finder coverage test. Adding a new per-target route now excludes its trait automatically.
    return {**per_target, "traits": set(TRAITS) - {k.removeprefix("dna-") for k in per_target}}


def cli_routable_cell_ids() -> set[str]:
    return set(by_cell_id())
