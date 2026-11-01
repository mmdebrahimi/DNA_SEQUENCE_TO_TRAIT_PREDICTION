# PromoterAtlas is identified, the circularity flag resolves, and my own "worse lead" is withdrawn

<!-- doc-status: LIVE as of 2026-10-31. Supersedes the UNRESOLVED half of
     wiki/promoteratlas_circularity_2026-10-31.md (written earlier the same day), which recorded
     "no source found uses the name PromoterAtlas" and "the arm names a model whose identity and
     provenance are both unpinned". Both are now false. That memo's Promoter Calculator half
     (CLEARED vs GSE108535) stands unchanged. -->

## CORRECTIONS FIRST — two claims of mine from earlier today are retracted

**1. "No source found uses the name PromoterAtlas."** False. It is a published model with public
MIT-licensed weights, and the query that found it was one word away from the one I ran. I searched
the model name paired with the repo's internal head label (`urtecho2018`); the name alone, paired
with "pretrained weights", returns the paper, the journal version, the GitHub repository and the
HuggingFace training corpus in the first five hits. **This is the recorded pattern of reaching for
an external blocker before exhausting a readable source** — the same error the salmserovar O-axis
arc logged three times in one session. The `ecolipromoterdb.com` TLS failure was real and was
irrelevant.

**2. The "worse lead" is WITHDRAWN as NOT SUPPORTED.** I proposed that the arm might be circular on
the **measured** tile/fragment rung via the Urtecho **genome-wide** study (321,123 sheared fragments
= GSE144621), which would have made the arm **inadmissible as a ceiling comparator** rather than
merely low-merit. The weights manifest is complete and lists five prediction heads — `urtecho2018`,
`lafleur2022`, `hossain2020`, `yu2021` and `kosuri2013` (translation) — and **no genome-wide or
tiling Urtecho head exists.** The base corpus is genomic sequence, not MPRA. The 321,123 match was a
numerical coincidence with the same lab's **later** study, which this model does not use.

I recorded it as *"a LEAD, not a finding"* at the time, which is the only reason withdrawing it
costs nothing. Had it been asserted, the arm would have been wrongly foreclosed.

## What PromoterAtlas is

Coppens et al., *"PromoterAtlas: decoding regulatory sequences across Gammaproteobacteria using a
transformer model"* — bioRxiv `2025.07.11.664312` → *Nature Communications* `s41467-026-72837-3`
(PMC13377108, PubMed 42140912).

| | |
|---|---|
| weights | `github.com/LucasCoppens/PromoterAtlas`, `trained_model_weights/`, **MIT** |
| base model | `promoteratlas-base.pt` — 1.8M-param custom **DNATransformer** (conv filters + rotary attention + residual feed-forward) |
| base corpus | **~9M regulatory sequences from >3,000 gammaproteobacterial species** |
| corpus artifact | HuggingFace `LCoppens/PromoterAtlas-data`, named in the repo README as *"the data used to train the base model"* |
| heads | `trspred-{urtecho2018, lafleur2022, hossain2020, yu2021}`, `tslpred-kosuri2013`, plus `annotation` |

The arm declares `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` — **the base model's
frozen embeddings**, not any prediction head. That distinction carries the whole result.

## The flag resolves in TWO parts, and collapsing them inverts the answer

Two committed artifacts carried the same open flag: *"whether the Promoter Calculator's 346
parameters or PromoterAtlas's urtecho2018 head were fit on GSE108535."*

### Part 1 — the flag's own question is answered YES

A head `promoteratlas-trspred-urtecho2018.pt` exists, described verbatim as *"Transcription
prediction based on Urtecho et al. 2018"*. **"Urtecho et al. 2018" is the σ70 variant library** —
bioRxiv 207332 posted 2018-01-28, *Biochemistry* online December 2018, 58(11):1539–1551, 10,898
designed promoters combining 8 −35 × 8 −10 × 3 UP × 8 spacer × 8 background elements = **GSE108535**.

So **that head is circular on the sigma70 rung and must never be used as a comparator there.** The
flag was right to exist and was pointing at something real.

### Part 2 — the arm does not use that head, and the base model is clean of GSE108535

The base corpus on HuggingFace is **3,371 genome FASTA files** keyed by RefSeq assembly accession
(`fasta/GCF_*.fasta`, matching the paper's ">3,000 species"), plus exactly two non-FASTA files: one
packed `sequence_dataset.h5` (~2.3 GB) and one `segmentation_data.csv`. **There is no expression
file, no MPRA file, no count matrix and no GSE accession of any kind.** The dataset card states the
provenance in its own words: *"The sequences were extracted from GenBank files of gammaproteobacterial
genomes."* Pretraining is self-supervised on genomic sequence.

**Verdict: `cleared_vs_GSE108535_for_the_BASE_model`.** The arm as declared carries no label
contamination from the sigma70 library. Spending it is now a **merit** call, no longer an
admissibility one.

## The residual that is real — sequence exposure, which is NOT what the flag asked

The base corpus is gammaproteobacterial genome assemblies. The peak-tile and fragment substrates
(GSE144621) are tiles and fragments **of the E. coli genome**. If an E. coli assembly sits among
those 3,371, the frozen encoder has seen those **sequences** during pretraining — never their
**labels**.

That is the ordinary situation for any pretrained DNA LM used as a frozen feature extractor and is a
strictly weaker concern than circularity. But it does mean **the arm may not be described as
evaluating on unseen sequence**, and that caveat would ship with its number.

**Measured:** the K-12 MG1655 reference assembly **`GCF_000005845` is ABSENT** from the corpus,
checked directly against the file list, as are `GCF_000008865` (O157:H7 Sakai) and `GCF_000019645`.
That is consistent with the repo's own demo annotating `NC_000913.1` — K-12 MG1655 — as a *held-out*
genome.

**Left OPEN deliberately, with the route recorded.** Whether some *other* E. coli assembly is among
the 3,371 needs an accession→organism resolution of the full list. The working route is the NCBI
Datasets v2 **GET** endpoint
(`/genome/accession/{acc,...}/dataset_report?returned_content=ASSEMBLY_INFO`); the **POST** form
returns 404, and the GET form returns **400 at 150 accessions per URL** while working at 2, so this
needs small chunks and roughly 340 requests. **I chose not to spend them:** this is a caveat on a
number that does not exist, for an arm that is `planned` and that the standing merit-negative below
argues against ever running. If the arm is ever funded, settle this first — it is cheap in absolute
terms, just not worth it now.

## A defect in my own measurement, recorded because it nearly shipped

The first accession→organism resolution used the wrong endpoint shape and resolved **0 of 3,371**
accessions — then printed, in its own output, **`Escherichia assemblies in the base corpus: 0`**.

That line reads exactly like a finding. It is a broken lookup wearing the costume of a result, and
it would have supported a *stronger* clean bill of health than the evidence allows. This is the
project's own recorded pattern (`feedback_a_test_that_errors_looks_like_a_finding`; the serotype
scorer's `exit 3` when both axes resolve zero) firing in a script I wrote one command after reading
that rule.

The rewritten script **refuses any presence/absence claim below a 90% resolution rate**. It was then
re-run, hit the same wall by a different cause (HTTP 400 on long URLs), resolved 0% again, and
**printed the refusal instead of the absence** — so the guard is not just written, it is
demonstrated.

## What this does NOT change

**The standing INFERRED merit-negative, untouched.** PromoterAtlas would be a **third
sequence-learned channel** beside the `conv` already in the stack — the failure mode `gc`
demonstrated at 0.4217 against the hybrid's 0.4254. On this rung every sequence-statistics arm
clusters **0.2691–0.3234** while the one **thermodynamic** arm sits at **0.3968**, so the next
channel should carry a different **PRIOR**, not a different neural net. The orthogonality
measurement that would settle it needs per-tile vectors the arm has not produced, so this stays an
inference.

**The arm stays `planned`, not killed.** Under R6 an arm dies only by its own executed evidence, and
an admissibility clearance is not executed evidence in either direction. `ceiling_verdict` stays
`blocked:untried-arm` with two planned arms (`promoteratlas`, `substrate-change-eukaryotic`).

## Honest limits

- The head→dataset identification rests on the repo's own one-line description plus the year, and on
  there being exactly one Urtecho σ70 library published as a 2018 preprint. I did **not** read the
  paper's methods section — bioRxiv returned HTTP 429 twice and the journal version sits behind an
  auth redirect. A methods read could still reveal a head trained on something else.
- "No expression file in the corpus" is read from the **file manifest** plus the dataset card's own
  provenance sentence, not from the contents of `sequence_dataset.h5`. A packed h5 holding labels
  would not show up as a filename.
- Clearing admissibility says nothing about whether the arm's number would be any good. The merit
  inference above still argues it would not.
