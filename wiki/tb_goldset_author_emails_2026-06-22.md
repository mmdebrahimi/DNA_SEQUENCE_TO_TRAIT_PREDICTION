> **⛔ SUPERSEDED — DO NOT SEND (2026-06-28).** The independent TB number these emails would request
> ALREADY EXISTS, for free, no DUA, no author contact: RIF acc 0.937 / INH acc 0.914, N=2,845
> (`wiki/tb_independent_number_2026-06-23.md`), via the EBI AMR Portal. Sending these would request data
> the project already has. Kept for audit trail only. See `wiki/frontier_reassessment_2026-06-28.md`.

# Author-request emails for the TB independent gold set (2026-06-22)

What every one of these asks for (the ONE thing that's missing from the public record): a **per-isolate table
linking each isolate's ENA/SRA accession to its MEASURED (phenotypic) rifampicin and isoniazid DST result
(R/S)**. The genomes + genotypic calls are already public; only the measured per-isolate phenotype is not.

Reply-to / sign-off on all: **Farshad — mmdebrahimi@gmail.com**.

---

## Email A — Prof. Taane Clark (LSHTM) — covers BOTH Thorpe 2024 + Thailand 2025
**To:** taane.clark@lshtm.ac.uk
**Cc:** Angkana.cha@mahidol.ac.th
**Subject:** Request: per-isolate measured RIF/INH DST for two of your M. tuberculosis WGS datasets

Dear Prof. Clark (and Dr Chaiprasert),

I'm developing and independently validating an open, deterministic *Mycobacterium tuberculosis* drug-resistance
caller (a transparent rule derived from the WHO 2023 mutation catalogue). To test it honestly I need an
*independent, measured-phenotype* benchmark — i.e. genomes scored against laboratory DST results rather than
against another genome-based prediction, which would be circular.

Two of your datasets are ideal because the genomes are public and phenotypic DST was actually performed:

1. Thorpe et al. 2024, *Sci Rep* 14:5201 (Thailand, 59 isolates; ENA accessions in Table S1).
2. Thawong et al. 2025, *Sci Rep* (Thailand, ~2,005 isolates; ENA accessions in Table S1; phenotypic DST on
   n≈1,826 summarised in Table S3).

In both, the per-isolate *measured* RIF/INH results aren't tabulated alongside the accessions in the public
supplements. Would you be willing to share a simple per-isolate table linking each isolate's ENA/SRA accession
to its measured rifampicin and isoniazid DST result (R/S, and ideally the method / critical concentration)?
A CSV or spreadsheet is perfect — no need for anything formatted.

It would let me report a genuinely independent accuracy figure for the tool. I'm glad to acknowledge or cite
in whatever way you prefer, and to share the results with you. Thank you very much for considering it.

Best regards,
Farshad
mmdebrahimi@gmail.com

---

## Email B — Ethiopia childhood TB 2026 (Abebaw et al., BMC Infect Dis)
**To:** woldearegay.erku@aau.edu.et  (Dr Woldaregay Erku Abegaz, Dept. of Microbiology, Immunology &
Parasitology, College of Health Sciences, Addis Ababa University — corresponding author, verified)
**Subject:** Request: per-isolate measured RIF/INH DST for your childhood DR-TB WGS isolates (PRJNA1104194/1204469)

Dear Dr Erku Abegaz,

I'm developing and independently validating an open, deterministic *M. tuberculosis* drug-resistance caller
(a transparent rule from the WHO 2023 catalogue) and am looking for an independent, *measured-phenotype*
benchmark — genomes scored against laboratory DST, not against another genome-based prediction.

Your 2026 study ("First-line drug-resistant tuberculosis among children under 15 years in Ethiopia",
*BMC Infect Dis*) is an excellent fit: you performed MGIT 960 first-line phenotypic DST and deposited the
raw reads publicly (PRJNA1104194, PRJNA1204469). The public supplement gives genotypic mutations, but not the
per-isolate measured phenotype linked to the accessions.

Would you be willing to share a per-isolate table linking each isolate's sample/SRA accession to its measured
rifampicin and isoniazid DST result (R/S, with method/critical concentration if handy)? A CSV/spreadsheet is
ideal. I'd gladly acknowledge or cite as you prefer and share the validation results. Thank you for considering.

Best regards,
Farshad
mmdebrahimi@gmail.com

---

## Email C — India nationwide WGS (Tamilzhalagan et al., Front Microbiol 2024, PRJNA1155695)
**To:** umadevi.kr@icmr.gov.in  (Dr Uma Devi Ranganathan, National Institute for Research in Tuberculosis,
ICMR, Chennai — corresponding author, verified; paper PMC11750862, BioProject PRJNA1155695, 2,207 isolates)
**Subject:** Request: per-isolate measured RIF/INH DST linked to accessions (PRJNA1155695)

Dear Dr Ranganathan,

I'm developing and independently validating an open, deterministic *M. tuberculosis* drug-resistance caller
(a transparent rule from the WHO 2023 catalogue) and need an independent, *measured-phenotype* benchmark
(genomes vs laboratory DST, to avoid circularity with genome-based tools).

Your nationwide WGS study (BioProject PRJNA1155695) is attractive because it is large and recent and phenotypic
DST was performed; Supplementary Table S3 lists the per-isolate accessions, but I couldn't find the per-isolate
*measured* RIF/INH results tabulated against them.

Would you be willing to share a per-isolate table linking each isolate's SRA/BioSample accession to its measured
rifampicin and isoniazid DST result (R/S, plus method if available)? A CSV/spreadsheet is perfect. Happy to
acknowledge/cite as you prefer and to share results. Thank you very much.

Best regards,
Farshad
mmdebrahimi@gmail.com

---

## Email D (optional) — TB Portals clinical-data follow-up
**To:** _the AccessClinicalData@NIAID / TB Portals data-access contact you already submitted to_
**Subject:** Follow-up: clinical (Patient_Cases DST) access for independent resistance-caller validation

Dear TB Portals / AccessClinicalData team,

I recently requested access to the TB Portals clinical data. I'm validating an open, deterministic
*M. tuberculosis* resistance caller and specifically need the **Patient_Cases** phenotypic DST fields
(the culture-based `bactec_*` and `le_*` rifampicin/isoniazid results) joined to the `ncbi_biosample` /
`ncbi_sra` accessions, to score genomes against measured phenotype. Could you let me know the status of my
request or any additional step needed? Thank you.

Best regards,
Farshad
mmdebrahimi@gmail.com

---

## Email E — PLOS Glob. Public Health 2025 (Elton et al., UCL; PRJEB68143) — added 2026-06-27
**To:** linzy.elton@ucl.ac.uk  (Dr Linzy Elton, Centre for Clinical Microbiology, University College London —
corresponding author, verified)
**Cc:** t.mchugh@ucl.ac.uk  (Prof. Timothy D. McHugh, senior author — verify exact address before sending)
**Subject:** Request: per-isolate measured RIF/INH DST linked to accessions (PRJEB68143)

Dear Dr Elton,

I'm developing and independently validating an open, deterministic *M. tuberculosis* drug-resistance caller
(a transparent rule from the WHO 2023 catalogue) and need an independent, *measured-phenotype* benchmark —
genomes scored against laboratory DST rather than another genome-based prediction (which would be circular).

Your 2025 *PLOS Global Public Health* WGS drug-resistance/lineage pipeline study (ENA PRJEB68143) is a clean
fit: you performed phenotypic DST and deposited the reads, and the supplement appears to link isolates to both
sequencing accessions (S3) and phenotypic/genotypic results (S6). Would you be willing to share a simple
per-isolate table linking each isolate's ENA accession to its measured rifampicin and isoniazid DST result
(R/S, with method/critical concentration if handy)? A CSV/spreadsheet is perfect.

It would let me report a genuinely independent accuracy figure for the tool; even your ~17-isolate set is a
useful first-cut benchmark. I'm glad to acknowledge or cite as you prefer and to share results. Thank you.

Best regards,
Farshad
mmdebrahimi@gmail.com

---

### Notes
- **Send Email A first** — it's fully addressed and covers the two cleanest datasets in one go (and Prof. Clark
  is a central figure in public TB genomics, so a positive reply may also point to other measured-DST sets).
- **B and C are now finalised** with verified corresponding-author addresses (Dr Woldaregay Erku Abegaz, AAU,
  `woldearegay.erku@aau.edu.et`; Dr Uma Devi Ranganathan, NIRT/ICMR, `umadevi.kr@icmr.gov.in`). All four ready
  to send. Salutations adjustable if you prefer "Prof." — both titles weren't explicit in the articles.
- The moment any of them replies with a per-isolate measured-DST CSV, the built pipeline turns it into a scored,
  leakage-checked independent number in one command.
