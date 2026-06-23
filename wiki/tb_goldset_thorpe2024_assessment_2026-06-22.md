# Thorpe et al. 2024 (Sci Rep 14:5201) as a TB gold set — assessment (2026-06-22)

**Verdict: NO-GO for a clean independent number from PUBLIC data — blocked on the per-isolate MEASURED pDST
(author contact). The genomes ARE genuinely independent; the blocker is the LABEL, not the data.**

## What was verified (EuropePMC full text PMC10908857 + Supplementary PDF MOESM1)
- **59 clinical isolates**, Thailand; raw reads on ENA (Illumina + ONT `ERR…` accessions, listed per-isolate
  in **Table S1**). 2024 deposit.
- **Genomes are independent of CRyPTIC:** leakage check (our `tb_goldset.assert_independent`) on the Table-S1
  ERR accessions = **0 overlap** with the CRyPTIC reuse table. Provenance-disjoint ✓.
- **Measured pDST EXISTS but only in aggregate:** *"Phenotypic drug susceptibility data was generated for
  rifampicin and isoniazid, and confirmed the genotypic resistance profiles."* Aggregate counts: RIF 7R/52S,
  INH 15R/44S (Table 1).
- **The per-isolate label in Table S1 is GENOTYPIC** (TB-Profiler-derived "Drug Resistance" category:
  Sensitive / RR-TB / MDR-TB / Pre-XDR-TB / Other). There is **no separate per-isolate measured-pDST column**.

## Why that's a NO-GO (the circularity gate G1)
Scoring our WHO-catalogue determinant rule against Thorpe's published per-isolate label = **a genome→
resistance rule (WHO catalogue) vs another genome→resistance rule (TB-Profiler catalogue)** = rule-vs-rule.
Agreement would be inflated BY CONSTRUCTION (overlapping catalogues) and would NOT measure performance vs
truth. That is exactly the project's #1 rejection gate (G1, circular label). The authors' statement that the
aggregate pDST "matched" the genotypic profile does not make the *published per-isolate* label a measurement.

## The unblock (one external action — yours)
Thorpe is the **closest viable candidate**: clean independent genomes + measured RIF/INH pDST was actually
performed. The ONLY missing piece is the per-isolate measured pDST table. Email the corresponding authors
for it; with it, this becomes a ready ~59-isolate independent gold set (genomes already leakage-clean).

### Drafted author request (you send)
> Subject: Request: per-isolate phenotypic RIF/INH DST for the 59 isolates in Sci Rep 14:5201
>
> Dear Dr Thorpe / Prof Clark (and colleagues),
>
> I'm validating a deterministic M. tuberculosis resistance caller and am looking for an independent,
> measured-phenotype benchmark. Your 2024 Scientific Reports paper (10.1038/s41598-024-55865-1) states that
> phenotypic DST for rifampicin and isoniazid was generated for the 59 isolates and confirmed the genotypic
> profiles. Supplementary Table S1 lists the ENA accessions and the genotypic resistance category, but I
> couldn't find the per-isolate *measured* pDST results (R/S per isolate for RIF and INH) tabulated.
>
> Would you be willing to share a per-isolate table linking each isolate's ENA accession to its measured
> RIF and INH pDST result (R/S, and method/critical concentration)? It would let me score my caller on a
> genuinely independent, phenotype-anchored set. Happy to acknowledge/cite as you prefer.
>
> Thank you, Farshad

## Once you have the per-isolate measured pDST
1. Build a candidate TSV: `strain_id, run_accession (ERR…), masked_vcf, regeno_vcf, rif_label, inh_label`
   (rif/inh = the MEASURED R/S from the authors).
2. `validate_tb_goldset_candidates.py` → `build_tb_goldset_manifest` (leakage check, already 0-overlap) →
   download the ERR FASTQs (to D:) → variant-call (TBProfiler / minimap2+bcftools) → masked VCFs →
   `score_tb_independent_goldset`. ~59 isolates (INH 15R/44S powered-ish; RIF 7R underpowered — report as such).

## Status vs the other candidates
Same wall as India `PRJNA1155695` (per-isolate measured label not public). Thorpe's advantage: **genomes
already leakage-clean + only 59 to fetch** (a true low-effort first cut) the moment the pDST table arrives.
TB Portals remains the no-author-contact powered option (sign the free DUA).
