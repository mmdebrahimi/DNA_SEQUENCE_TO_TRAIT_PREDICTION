# ASK #2 — the lab / MTA request (new measured labels = new capability)

**Status: DRAFT + verified target shortlist. Sending is the user's action — Soraya does not send.**

## The strategic point (from the path-B probe, 2026-07-16)

Free phenotype data is *abundant* but lives in the polygenic/ancestry-confounded regime and is already owned
by PGS Catalog — it does **not** clear the wall. The only acquisition worth relationship capital is
**mechanism-resolved, isolate-level, MEASURED phenotype paired with sequence**, for an organism/drug the
decoder doesn't yet cover. That is what this ask targets. (ARESdb — the *buyable* option — is verified
proprietary with no downloadable genomes; **money would not solve it**, so it is closed.)

## Verified target shortlist

These are groups whose paired WGS+measured-AST datasets this project has **already worked with or verified**
— so the ask is credible and specific, not cold. Each row's dataset is real and cited in a committed artifact.

| # | group / cohort | what they have | why them | status |
|---|---|---|---|---|
| 1 | **Oxford — Modernising Medical Microbiology** (Lipworth et al.; BioProject `PRJNA604975`) | ~2,875 *E. coli* bacteraemia WGS + **measured clinical MIC**, 8 drugs, UK | Their blood-culture WGS is already **openly deposited** — we've built the external-revalidation arm against it. The ask is for the **MIC table / the urine subset**, not the genomes. Warmest, most specific lead. | genomes public; MIC-table access = the ask |
| 2 | **CRyPTIC Consortium** (TB) | 12,287 *M. tuberculosis* + measured BMD MIC | We already use their public compendium; a relationship could unlock the **post-2023 / non-public** slice that would give a truly out-of-build-set TB number | public part already used |
| 3 | **Spain — PROBAC** (ENA `PRJEB62601`) | 224 *E. coli* + EUCAST BMD MIC, 16 drugs | Independent, non-US, measured MIC — a second external cohort | accession verified |
| 4 | **A local clinical micro lab** (your region) | routine AST + isolates | The only route to a **new organism/drug** cell; lowest data-sharing friction if there's a personal connection | requires your network |

> ⚠️ **Blanks I refuse to guess:** every **email address** and **current affiliation**. Take them from the
> paper's corresponding-author line or the group's website. **Do not send to an address I generated — I
> haven't generated any.**

> **Do the free check first:** for #1 and #3 the *genomes* are already public. Before asking anyone, confirm
> whether the MIC table is in the paper's supplement — you may need no email at all. (This is the same lesson
> that just downgraded the marine-Vibrio lead: `PRJNA1254427` turned out to have **0 public assemblies**,
> and `PRJNA645603` only **4** — so that path is dead on power, not on access.)

---

## Draft email (adapt per target)

**Subject:** Paired WGS + AST data — external validation of an open-source resistance decoder

Dear *[Dr X]*,

I maintain `dna-decode`, an open-source deterministic genotype-to-phenotype decoder (MIT-licensed, on PyPI).
It calls antibiotic resistance from a genome using curated determinant catalogues — no black-box model — and
it currently has 10 cells scored on provenance-disjoint NCBI Pathogen Detection cohorts, plus an independent
*M. tuberculosis* number and an HIV cell validated against Stanford's PhenoSense fold-change data.

I'm writing because the honest limit of the project is **labels, not methods**. The public genotype↔phenotype
space is saturated for the organisms I cover, and what I can't obtain freely is paired sequence + **measured**
AST for new organisms/drugs — which is exactly what *[your cohort / PRJNA…]* has.

**What I'd like to ask:** whether you'd be open to sharing *[the per-isolate MIC table for the deposited
genomes / a paired WGS+AST subset]*, under an MTA or data-sharing agreement if that's the right vehicle.

**What I'd do with it:** score the **already-frozen** decoder on it as a genuinely external cohort — no
retraining and no tuning, so the result is a real out-of-distribution test. I'd report per-cell
sensitivity/specificity with clonality correction and confidence intervals, **including the failures**, and
send the analysis back to you. I don't need patient identifiers or clinical metadata — only sequence, the AST
result, and the method used. Your data would not be redistributed, and nothing would be published without
your agreement.

I've attached a one-page summary of what the decoder does and how it's validated (including the results that
make it look bad — e.g. our TB sensitivity drops from 0.92 to 0.44 once we correct for clonal over-sampling,
and that corrected number is the one we headline).

If this isn't something you can share, no problem at all — and if there's a better-placed group or a route I
should be using instead, I'd be grateful for the pointer.

Best regards,
Farshad *[surname]*
*[email]*
`https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`

*[attach: `wiki/acquisition/decoder_credibility_sheet_2026-07-16.md`]*

---

**Why this framing:** it leads with what we *can't* do (labels are the limit) rather than a pitch; it asks for
a specific artifact; it names the exact analysis and the fact that we'd report failures; it explicitly
minimises what we're asking for (no PHI); and it gives an easy, face-saving out. The clonality anecdote is
load-bearing — it's the fastest way to signal to a microbiologist that we won't inflate their data.
