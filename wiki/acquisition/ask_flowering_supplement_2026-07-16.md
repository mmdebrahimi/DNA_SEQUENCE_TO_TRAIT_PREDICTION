# ASK #1 — the flowering supplementary table (the smallest, highest-VOI request)

> ## ⚠️ CORRECTION 2026-07-16 — THIS EMAIL IS PROBABLY UNNECESSARY. TRY THE BROWSER FIRST.
>
> **The paper is OPEN ACCESS.** An Unpaywall lookup on the DOI returns **`is_oa: True`, `oa_status: hybrid`**,
> with a publisher OA PDF (`.../doi/pdfdirect/10.1111/tpj.14716`) and an MPG repository copy
> (`http://hdl.handle.net/21.11116/0000-0007-03AB-5`).
>
> **The earlier "HTTP 402 Payment Required → money gate" conclusion in this file's first draft was WRONG.**
> Two compounding errors: (a) the 402 came from the `/doi/full/` *HTML* route, not the OA PDF route; and
> (b) `pdfdirect` returns **403 = BOT-BLOCKING**, which is not a payment wall — a human in a browser gets it
> free. Europe PMC's `isOpenAccess: N` was also misleading: it tracks *PMC deposit*, not OA status.
>
> ## ⚠️ CORRECTION #2 (same day) — THE TARGET IS **TABLE S3**, NOT TABLE S1.
>
> I fetched the **publisher-version PDF free** from the MPG repository REST API (no browser, no paywall:
> handle → `item_3252572` → `component/file_3252619/content`; saved at
> `data/arabidopsis/zhang2020/zhang_tpj14716_publisher.pdf`). Its own Supporting-Information caption list
> settles what each file actually is:
>
> | file | actual contents |
> |---|---|
> | Table S1 | **Variants found in the FRI locus** — the 171 mutations. **NOT the per-accession mapping.** |
> | Table S2 | Variant comparison with previous works |
> | **Table S3** | **"Allele id and flowering time for each accession"** ← **THE FILE WE NEED** |
> | Data S1/S2 | FRI allele DNA sequences / ClustalW alignment |
>
> My earlier "Table S1 = 1,016 accessions × FRI allele" came from a **search-engine summary** and was wrong.
> **Bonus:** S3 carries **allele id AND flowering time**, so scoring may not need AraPheno at all — the join
> key and the label ship in one file.
>
> **Do this (2 min, free, no favour owed) — a BROWSER is genuinely required:**
> 1. Open `https://onlinelibrary.wiley.com/doi/10.1111/tpj.14716` → **Supporting Information** →
>    download **Table S3** (and S1 too if it's one bundled file — no harm).
> 2. Save it anywhere; hand Soraya the path.
>
> **Why a human must do this step:** Wiley's Cloudflare returns **403 to every scripted request** — article
> HTML, `pdfdirect`, and `downloadSupplement` alike — even with a full browser User-Agent. That is
> **bot-blocking, not a paywall** (the paper is OA). The MPG repository holds **only** the article PDF (one
> component, no supplements), so there is no non-browser route to S3. This is a real capability boundary,
> not an unexplored option.
>
> **Only if the SI itself turns out to be genuinely unavailable** does the email below become the move.
>
> **Reusable rule:** an HTTP 402/403 from `curl`/WebFetch is **not** evidence of a paywall — publishers
> bot-block. **Check Unpaywall by DOI (`is_oa`) before ever calling something a money gate.**

**Status: DRAFT — DEMOTED TO FALLBACK by the correction above. Sending is the user's action; Soraya does not send.**

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
