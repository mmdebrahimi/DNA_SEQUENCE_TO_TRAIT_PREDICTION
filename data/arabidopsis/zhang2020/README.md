# Zhang & Jiménez-Gómez 2020 — FRIGIDA natural variation (the flowering cell's substrate)

**Table S3 is here, and the flowering cell is now SCORED on it** (N=854 phenotyped of 1,017) —
see `wiki/flowering_tables3_score_2026-07-16.md` and `scripts/flowering_tables3_score.py`.

| | |
|---|---|
| paper | *Functional analysis of FRIGIDA using naturally occurring variation in Arabidopsis thaliana* |
| journal | The Plant Journal 103:154-165 (2020) · doi **10.1111/tpj.14716** · PMID 32022960 |
| open access | **YES** — Unpaywall `is_oa: true`, `oa_status: hybrid` |
| **licence** | **CC-BY 4.0** (Crossref: `creativecommons.org/licenses/by/4.0/`, content-version `vor`, delay 0) |

**The licence is why the supplements are committed here.** CC-BY permits redistribution with attribution,
and these files are *browser-only* (below) — so committing them is what makes the scoring reproducible for
anyone else. Cite the paper (see `CITATION` in the scorer) with any reuse.

## What's committed

| file | contents | used for |
|---|---|---|
| **`tpj14716-sup-0012-TableS3.tsv`** | **allele id + FRI status + FT16 per accession** (1,017 rows) | **the scored substrate** |
| `tpj14716-sup-0010-TableS1.xlsx` | the 171 FRI variants × 1,016 accessions genotype matrix (`prot_change`, REF, ALT, `funct_consequ`) | **v0.1 genome-mode input** (not yet used) |
| `tpj14716-sup-0011-TableS2.xlsx` | variant comparison with previous works | — |
| `tpj14716-sup-0003-figs3.pdf` | % functional/non-functional FRI per STRUCTURE group, Fisher's test | evidence for the population-structure confound the scorer corrects for |
| `tpj14716-sup-0013-DataS1.docx` | DNA sequences of all cloned FRI alleles | v0.1 reference material |
| `tpj14716-sup-0014-DataS2.docx` | ClustalW alignment of cloned alleles | v0.1 reference material |

*Not committed:* the article PDF (gitignored — it's one API call away, below).

## Table S3 columns

`accession_id`, `cloned`, **`deleterious_allele`** (TRUE/FALSE — the FRI functional call), `allele_group`
(a001…a103), `name`, `CS_number`, `count` (**mislabeled — it holds the country code**), `latitude`,
`longitude`, **`group`** (STRUCTURE population group), **`FT16_mean`** (days to first flower, long days 16 °C).

**Two traps, both live:**
1. **Missing FT16 is the string `NA`, not blank** — 163 of 1,017 (16%). A naive `float()` crashes; a naive
   truthiness check *passes it through*. The dropout is **not random**: only 9.8% of dropped accessions carry
   a deleterious FRI vs a 24% base rate, so functional-FRI (late-candidate) accessions are preferentially
   unphenotyped. `phenotype_attrition()` reports this.
2. **`group` is `NA` for 130 accessions** — all of which also lack FT16, so they vanish from scoring anyway.

`load_table_s3()` verifies the file against the paper's own stated counts (1,017 rows / 245 deleterious /
103 allele groups) and **refuses to score** on a mismatch.

## The article PDF — free, scripted, no browser

Wiley's Cloudflare 403s every scripted request, but the Max Planck repository serves the publisher version:

```bash
# handle 21.11116/0000-0007-03AB-5 -> pubman item_3252572 -> component file_3252619
curl -L -o zhang_tpj14716_publisher.pdf \
  https://pure.mpg.de/rest/items/item_3252572/component/file_3252619/content
```

596,981 bytes. Gitignored — regenerable in one call.

## The supplements — browser-only (a real capability boundary)

PuRe holds **exactly one component** (the article PDF); no supplements are deposited. Wiley 403s
`downloadSupplement` like everything else — **bot-blocking, not a paywall** (the paper is OA/CC-BY). A human
with a browser gets them free from the article's *Supporting Information* section; that's how these arrived.
See [[feedback_buildable_from_literature_but_scorable_only_behind_paywall]].

**Which table is which** (from the PDF's own caption list — an earlier draft named Table S1 as the
per-accession mapping, sourced from a search-engine summary, and was **wrong**):

- **Table S1** = *variants* in the FRI locus (the 171 mutations) — not per-accession phenotype
- **Table S3** = *"Allele id and flowering time for each accession"* — the scored substrate

## Result summary

Pooled accuracy **0.733** vs a 0.502 constant-predictor null — but the honest figure is the
**population-structure-weighted 0.710 vs its own 0.676 null (+3.4 pp)**, because FRI genotype correlates
with ancestry. The rule is strongly predictive in one direction only (FRI loss → early, 93.9%) and weak in
the other (FRI intact → late, 65.8%): functional FRI is **necessary but not sufficient**, exactly as the
cell's two-locus rule says. Full detail + all caveats: `wiki/flowering_tables3_score_2026-07-16.md`.
