# ASK #1 — the flowering supplementary table (the smallest, highest-VOI request)

**Status: DRAFT. Sending is the user's action — Soraya does not send.**

## What this unblocks

One spreadsheet converts an external wall into a code wall: the flowering-habit cell's scoring harness is
**built and runs** (`scripts/flowering_arapheno_spotcheck.py`) — it already fetches free AraPheno DTF1
(934 accessions), derives the early/late split from the cohort median, joins by accession, and gates on a null
baseline. It is missing exactly one input: **per-accession FRI/FLC allele calls**. With that file, the cell
becomes the project's first **scored plant cell**, same day.

## Verified target (checked via Europe PMC API, 2026-07-16 — not from memory)

| field | verified value |
|---|---|
| paper | *Functional analysis of FRIGIDA using naturally occurring variation in Arabidopsis thaliana* |
| journal | **The Plant Journal**, 2020 · doi **10.1111/tpj.14716** · PMID **32022960** |
| authors | **Zhang L**, **Jiménez-Gómez JM** (only two) |
| corresponding | **José M. Jiménez-Gómez** (last author / PI) |
| affiliation *on the paper* | Dept. of Plant Breeding and Genetics, **Max Planck Institute for Plant Breeding Research** |
| open access | **No** (`isOpenAccess: N`) — independently confirms the paywall |
| what we want | **Table S1** — the per-accession FRI allele assignment across the 1,016 retained 1001-Genomes accessions |

> ⚠️ **Two blanks you must fill — I refuse to guess these:**
> 1. **Email address** — take it from the paper's corresponding-author line, or the current institute
>    directory. **PIs move**: the MPIPZ affiliation is what the 2020 paper says; verify where he is now.
> 2. **His current institution** in the greeting, if it has changed.
>
> **Free alternatives if he doesn't reply** (same data, different papers): Shindo et al. 2005 *Plant Physiol*
> (176 accessions: FRI haplotype + flowering time + FLC level) or Werner et al. 2005 *Genetics* (145
> accessions genotyped for both common FRI lesions). Also worth one try first: **a university librarian or
> any colleague with institutional access** can pull Table S1 in about two minutes — no email needed.

---

## Draft email

**Subject:** Request: Table S1 (per-accession FRI alleles) from your 2020 Plant Journal FRIGIDA paper

Dear Dr Jiménez-Gómez,

I'm building an open-source, deterministic genotype-to-phenotype decoder (`dna-decode`, MIT-licensed) as an
independent project. It applies curated causal-locus catalogues to call traits from a genome — currently
antibiotic resistance across bacteria/viruses/fungi, plus human pharmacogenomics — and I've recently added a
first plant cell: an *Arabidopsis* flowering-habit caller built on the FRI/FLC mechanism (functional FRI AND
strong FLC → winter-annual; loss of either → summer-annual, with the FRI-independent cases explicitly
flagged as outside its scope).

The cell reproduces the canonical literature anchors, but I want to do the honest thing and **score it against
measured data** rather than leave it as a restatement of published assignments. I have the free AraPheno DTF1
flowering-time values for ~934 accessions; what I'm missing is the per-accession FRI allele assignment — i.e.
**Table S1 from your 2020 Plant Journal paper** (doi 10.1111/tpj.14716). I don't have institutional access to
the Wiley supporting information.

Would you be willing to share that table? I'd use it purely to evaluate the caller's accuracy on the
accessions where measured flowering time exists, and I'd of course cite the paper as the source of the allele
assignments. If it's easier, even a subset (accession name + FRI allele/functional status) would be enough. I'm
happy to send back whatever the evaluation shows — including if it shows the two-locus rule performs poorly,
which is a result I'd publish either way.

Thank you for considering it, and for the allele resource itself — the 103-allele characterisation is exactly
the kind of curated ground truth this approach depends on.

Best regards,
Farshad *[surname]*
*[email]*
`https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`

---

**Why this framing:** it's short, it asks for one specific thing, it explains the use, it offers reciprocity,
and it signals we'll report a negative result honestly — which is what makes a researcher trust a stranger
with their data. It does **not** overclaim: the cell is described as reproducing anchors, not as validated.
