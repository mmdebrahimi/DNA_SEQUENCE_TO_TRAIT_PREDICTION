# Zhang & Jiménez-Gómez 2020 — FRIGIDA natural variation

Source of the per-accession FRI allele calls needed to **score** the flowering-habit cell
(`dna_decode/organism_rules/arabidopsis_flowering.py`). Until the table below lands, that cell is
**faithful-to-catalogue but NOT SCORED**.

| | |
|---|---|
| paper | *Functional analysis of FRIGIDA using naturally occurring variation in Arabidopsis thaliana* |
| journal | The Plant Journal (2020) · doi **10.1111/tpj.14716** · PMID 32022960 |
| open access | **YES** — Unpaywall `is_oa: true`, `oa_status: hybrid` |

## The article PDF — free, scripted, no browser (this works)

Wiley's Cloudflare 403s every scripted request, but the Max Planck repository holds the
**publisher version** and its REST API serves it plainly:

```bash
# handle 21.11116/0000-0007-03AB-5 -> pubman item_3252572 -> component file_3252619
curl -L -o zhang_tpj14716_publisher.pdf \
  https://pure.mpg.de/rest/items/item_3252572/component/file_3252619/content
```

596,981 bytes. **Gitignored** — regenerable in one call, and we don't redistribute a publisher binary.

## The supplementary tables — browser required (this is the open item)

The PuRe record has **exactly one component** (the PDF above). **No supplements are deposited there**, and
Wiley 403s `downloadSupplement` like everything else. So a human with a browser is genuinely required:

> `https://onlinelibrary.wiley.com/doi/10.1111/tpj.14716` → **Supporting Information** → download **Table S3**.

**403 here is bot-blocking, not a paywall** — the paper is OA and a human gets the file free. See
[[feedback_buildable_from_literature_but_scorable_only_behind_paywall]].

## Which table (corrected 2026-07-16 — read this before fetching)

Straight from the PDF's own Supporting-Information caption list — **not** from a search summary:

| file | actual contents | want? |
|---|---|---|
| Table S1 | Variants found in the FRI locus (the 171 mutations) | no |
| Table S2 | Variant comparison with previous works | no |
| **Table S3** | **"Allele id and flowering time for each accession"** | **YES** |
| Data S1 / S2 | FRI allele DNA sequences / ClustalW alignment | no |

An earlier draft of the acquisition ask named **Table S1** as the per-accession mapping. That came from a
search-engine summary and was **wrong**; S1 is the variant list. The PDF settles it.

**Bonus:** S3 carries **allele id AND flowering time** in one file — both the join key and the label. So
scoring may not need the free AraPheno DTF1 values at all (`scripts/flowering_arapheno_spotcheck.py`
currently fetches those and joins by accession).

## Once Table S3 is in hand

Point the spot-check at it. Note the null-baseline gate is load-bearing: the AraPheno-median split already
produced `DEGENERATE_NO_LABEL_VARIATION` once — the cell went 4/5 while a constant-`early` predictor went
5/5. A score is only meaningful on a set where both classes are actually observed.
