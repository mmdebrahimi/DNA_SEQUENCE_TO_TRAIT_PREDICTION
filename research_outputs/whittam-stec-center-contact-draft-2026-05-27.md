# Whittam STEC Center Contact Draft (2026-05-27)

> Draft for user to send to stec@cvm.msu.edu (Manning lab, MSU). Covers Short-term Action 2-new from `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md`. Fires Q7-resolution.

---

**To:** stec@cvm.msu.edu
**Subject:** DECA reference collection — strain inventory + curation provenance request for an open-source pathotype CLI project

Dear Whittam STEC Center team,

I'm developing an open-source command-line tool for E. coli pathotype prediction from genome assemblies. The tool emits a multilabel virulence-cluster profile + derived pathotype call (EPEC / EHEC / ETEC / UPEC / EAEC / commensal-compatible + HYBRID + AMBIGUOUS + UNCLASSIFIED + STEC-non-LEE + tEPEC/aEPEC split) with explicit abstention rules, using pinned VirulenceFinder + ETECFinder gene calls as the underlying caller layer.

A central concern for the project is **substrate label independence** — the model's labels must come from sources whose pathotype assignments are independent of the same virulence-gene rules the tool itself applies, to avoid circular validation. The Whittam STEC Center's DEC and DECA collections are the strongest candidate "bedrock" substrate I've identified because the pathotype labels are anchored to clonal lineage (MLST/electrophoretic types) and curated literature, which predate the modern VirulenceFinder-style gene-rule pipelines.

I have three specific information requests:

1. **DECA collection per-pathotype strain counts.** The Hazen et al. 2012 paper (J Bacteriol 194:3026) references "the 15 most common diarrheagenic clones," but I have not been able to find a publicly accessible per-clone or per-pathotype breakdown. Could you share the per-pathotype counts (EPEC / EHEC / ETEC / EAEC / EIEC) and per-clone strain rosters?

2. **Curation methodology / provenance per strain.** For each strain in the DEC and DECA reference sets, how was the pathotype label originally assigned? Specifically, was each label derived from (a) clinical/epidemiological metadata at isolation, (b) MLST clonal-lineage assignment, (c) virulence-gene PCR or sequencing assays, (d) curated literature, or (e) some combination? If there is per-strain provenance metadata (e.g., a column in the strain registry indicating "label source"), I'd be very grateful to see the schema.

3. **MTA terms + access constraints.** I noticed the STEC Center has a $7/strain non-profit fee and a non-export-from-US restriction. Could you confirm what the MTA process looks like for a researcher working from Canada (Toronto) on an open-source software project, and whether the non-export restriction affects derivative data (e.g., assembled genome sequences + per-strain metadata used as training/test substrate, with no physical material crossing borders)?

For context: the project's source ledger is being maintained openly; I'm happy to share design notes, the working substrate-tier strategy, and acknowledge the Whittam STEC Center prominently in the resulting tool's documentation. The Phase 1 of the same project (E. coli antibiotic resistance prediction) closed in May 2026 with internal artifacts available for review on request; this is Phase 4 (pathotype prediction) of the same lineage of work.

Thank you for your time — and thanks for maintaining the DECA collection. It's a foundational resource for the field.

Best regards,
[your name]
[your affiliation, if applicable]
[your email]

---

## Notes for the user before sending

- **Customize the closing block** — fill in your name, affiliation (if any; "independent researcher" is fine), and reply-to email.
- **The mention of Phase 1 (AMR resistance) closeout** establishes credibility / prior work in the same domain. If you don't want to reference it, drop that sentence.
- **The "non-export restriction affects derivative data" question** is the load-bearing access question. If they say MTA allows derivative data, the project can proceed; if they say it requires physical strain transfer, the timeline + cost shift significantly.
- **Likely response time:** academic labs typically respond within 1-2 weeks. If no response in 14 days, follow up directly to Shannon Manning (the lab's PI; CV-tracked email contactable via MSU faculty page).
- **Optional cc:** Tracy Hazen (first author Hazen 2012 DECA paper) at the Institute for Genome Sciences, University of Maryland — she would have the original per-clone strain rosters if MSU's records are incomplete.
