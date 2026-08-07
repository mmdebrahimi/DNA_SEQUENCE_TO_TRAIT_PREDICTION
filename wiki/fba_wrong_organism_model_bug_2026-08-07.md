# FBA cross-organism registry ships TWO WRONG-ORGANISM models (2026-08-07)

**Severity: HIGH — user-facing, in a PUBLISHED package (v0.11.0 live on PyPI). v0.12.0 is tagged but
NOT yet published; do NOT publish it as-is.**

## The defect

`dna_decode/fba/model.py::_BIGG_MODELS` maps two organism aliases to models of a DIFFERENT organism.
Verified against the BiGG API (`http://bigg.ucsd.edu/api/v2/models/<id>`) 2026-08-07:

| `--organism` alias | loads model | BiGG's actual organism | verdict |
|---|---|---|---|
| `escherichia_coli` / `ecoli` | `iML1515` | *Escherichia coli* K-12 MG1655 | CORRECT |
| `saccharomyces_cerevisiae` / `yeast` | `iMM904` | *Saccharomyces cerevisiae* S288C | CORRECT |
| `staphylococcus_aureus` / `saureus` | `iYS1720` | **Salmonella pan-reactome** | **WRONG ORGANISM** |
| `pseudomonas_aeruginosa` / `paeruginosa` | `iJN1463` | **Pseudomonas putida KT2440** | **WRONG ORGANISM** |

Corroborating evidence (not just the API label): the cached `iYS1720.xml.gz` has **1262 of 1707 gene
products prefixed `STM`** — the *Salmonella* Typhimurium LT2 locus-tag prefix. The v0.11.0-era note in
`essentiality_labels.py` ("iYS1720 gene ids are STM#### (Salmonella-style locus tags)... needs a
gene-name crosswalk") OBSERVED this smell but misdiagnosed it as an ID-convention quirk rather than
as "this is the wrong organism's model".

## What it invalidates

- The v0.11.0 / v0.12.0 **"cross-organism engine generalizes"** claim: it generalizes to
  *E. coli*, yeast, **Salmonella**, and ***P. putida*** — not to S. aureus or P. aeruginosa.
- Ledger row 631's `LABEL_WALLED` framing for S. aureus + P. aeruginosa. The binding constraint was
  never the LABEL — it was that the MODEL was the wrong organism. A gold standard for either organism
  could never have joined, no matter how good the label.
- Claim surface to correct: `README.md`, `CHANGELOG.md`, `dna_decode/cli.py`,
  `dna_decode/data/shipped_decoder_surface.py`, `dna_decode/fba/cli.py`,
  `dna_decode/fba/essentiality_labels.py`, `dna_decode/fba/model.py`, plus ledger rows 629/631.

E. coli (`iML1515`, Keio-validated, MCC 0.652) and yeast (`iMM904` vs SGD, MCC 0.252) results are
**UNAFFECTED** — both models are correctly assigned.

## The fix (available, verified)

- **S. aureus** — BiGG HAS real models: `iYS854` (*S. aureus* USA300_TCH1516, 866 genes, gene ids
  `USA300HOU_####` / `USA300HOU_RS#####`) and `iSB619` (*S. aureus* N315, 619 genes, ids `SA####`).
  `iYS854` is the better pick: USA300 is the strain background of the **NTML / Nebraska transposon
  library**, so the essentiality label and the model are the same strain lineage. Caveat: NTML is
  USA300 **JE2/FPR3757** (`SAUSA300_####`) while iYS854 is USA300 **TCH1516** (`USA300HOU_####`) —
  near-1:1 orthologs but still a crosswalk, NOT a free join.
- **P. aeruginosa** — **BiGG has NO P. aeruginosa model at all** (queried the full model list;
  zero hits on "aeruginosa"). Options: (a) drop the alias and fail loudly, (b) source a P. aeruginosa
  GEM from outside BiGG (e.g. iPae1146 / PA14 reconstructions), (c) re-label the alias as
  *P. putida* (honest, since `iJN1463` IS a valid P. putida model).

## Bonus: the P. aeruginosa gold standard IS fetchable (but keyed to PA14)

The PLOS Comput Biol 2026 supplementary **downloaded cleanly** from this host today (the bot-block hit
on 2026-08-03 was transient):
`https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1013945.s011&type=supplementary`
→ XLSX with sheets **`GOLD_84`** (84 genes) and **`GOLD_115`** (115 genes), columns
`Locus tag | Old Locus tag | Gene name | Product description | Function PseudoCAP | ...`.

**CORRECTION to `wiki/label_wall_data_sources_2026-08-03.md`:** that memo claimed GOLD_84 uses PAO1
`PA####` tags giving "a DIRECT join to the shipped iJN1463 model, no crosswalk". **Both halves are
wrong.** The tags are **PA14** (`PA14_RS00005` / old `PA14_00010`), and `iJN1463` is *P. putida*.
There is no clean join; there is not even a correct model yet.

## Reusable lesson

A gene-ID prefix that doesn't match the organism you think you loaded is a **wrong-organism** signal,
not an ID-convention quirk. Verify a fetched model's organism against the source registry's own
metadata BEFORE building a validation claim on it — the same class as
`feedback_verify_cited_accessions_resolve_before_build`.
