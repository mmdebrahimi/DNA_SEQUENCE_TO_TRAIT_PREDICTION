# What this project actually is — read from the record, not from memory

Written after mis-stating the project's own state twice in one session. This is distilled from the four
session transcripts (2026-07-09 → 2026-08-29), the full git history, and the live registries — **not** from
recollection.

## Scale, measured

| | |
|---|---:|
| commits | **1,533** |
| span | 2026-05-11 → 2026-08-29 (**16 weeks**) |
| conversation sessions on record | 4 (202 MB of transcript) |
| user turns | 1,121 |
| distinct steering directives | **469** |
| routine `--advance` invocations | 99 |
| tests | 3,953 passing |

## What actually ships — and this is the part I did not have

**`dna-decode` is a PUBLISHED tool at v0.13.1**, not a research repo.

- **48 console entry points**, **46 CLI traits**
- **110 registered cells** in `cell_registry`

Release arc: v0.9.0 productization → v0.10.0 animal fleet + microbial/viral → v0.11.0 FBA metabolic cell →
v0.12.0 FBA cross-organism + motility → v0.13.0 reachability (clinvar/hla routable, `dna-forward` takes real
DNA) → **v0.13.1 published**.

### Evidence distribution across all 110 cells

| tier | cells |
|---|---:|
| `INDEPENDENT_MEASURED` | **28** |
| `NEAR_INDEPENDENT` | 24 |
| `KNOWLEDGE_BASELINE` | 33 |
| `FAITHFUL_TO_TOOL` | 13 |
| `NO_FREE_SOURCE` | 11 |
| `NOT_CENSUSED` | 1 |

**The correction that matters most:** the validation report card's *"27 cells / 10 SCORED"* — which I
treated all session as the project's validated surface — is the **AMR provenance-disjoint arm only**. The
tool's actual evidence surface is 110 cells with **28 independently measured**. I was reporting a subset as
if it were the whole.

## Where the 1,533 commits went

| theme | commits | share |
|---|---:|---:|
| AMR decoder (bacterial rules/cells) | 202 | 13% |
| Validation surface / report card | 189 | 12% |
| Infra / CLI / tests / packaging | 185 | 12% |
| Planning artifacts (plans, retros, ledger) | 100 | 7% |
| FBA / metabolic | 89 | 6% |
| Forward / variant-effect / inverse | 83 | 5% |
| TB | 80 | 5% |
| Staleness auditor / docs hygiene | 79 | 5% |
| Foundation models / embeddings | 70 | 5% |
| HIV / viral | 51 | 3% |
| Pathotype / virulence / genome-map | 51 | 3% |
| Human clinical (pgx / clinvar / hla) | 37 | 2% |
| Colour / pigment trait cells | 33 | 2% |
| Eukaryote (Arabidopsis / BXD / yeast) | 23 | 2% |
| Fungal | 15 | 1% |

## Arcs I had no awareness of this session

Found only by reading the record:

- **Phage** (`dna-phage`) — including a cross-lab RBP result recorded as **DATA-BLOCKED**, and an overnight
  diagnostic that **REFUTED** its own 0.364 caveat.
- **Multimodal** — a Family A/B falsifier that **killed Family B and promoted Family A**.
- **A whole molecular-typing suite** — `plasmid`, `serotype`, `mlst`, `ktype`, `salmserovar`,
  `pneumo-serotype`, `resfinder`, `pointfinder`, `disinfinder`.
- **Human clinical cells** — `pgx`, `clinvar`, `hla`.
- **`motility`**, **`concordance`**, **`coloc`**, **`profile`** — first non-metabolic trait catalog.
- **A 19-cell animal colour/plumage fleet**, later frozen.

## The genotype→phenotype picture, corrected

Both this session's strategic errors were the same shape: a scoped result compressed into an unscoped
label.

| regime | status |
|---|---|
| natural population + zero-shot embedding | **closed negative** (0-for-5, de-confounded) |
| constructed variation → molecular phenotype | **works** — TEM-1 genome-edit path, **Spearman 0.761** vs measured ampicillin fitness |
| constructed variation → organism phenotype (yeast cross) | **works** — **12/12 traits, r 0.46–0.80** |
| constructed variation → organism phenotype, per-condition | **works** — FBA iML1515 **MCC 0.70–0.74** |
| constructed variation → organism phenotype, **condition-switch** | **OPEN**, ~null, bottleneck measured |

**The discriminating variable is population design, not organism complexity.**

## Why this document exists

Across 16 weeks and 1,533 commits, no single context window holds this project. That is a structural fact,
not a failure of attention — and it means **every strategic claim must be re-derived from the artifacts**,
because a confident summary from memory will be wrong in ways that sound right.

Two concrete instances, both this session:

1. Called organism-level g→p a *"closed negative"* — while a clean **12/12, r 0.46–0.80** positive sat in
   the repo. (Third occurrence of that specific compression.)
2. Proposed the FBA/Keio benchmark as an *unexplored* direction — while it is the **deepest line in the
   repo**, ~25 artifacts with pre-registrations, an adversarial review, and a retraction.

Both groundings were on disk the whole time.

**Reproduce this document:** parse `~/.claude/projects/C--Users-Farshad-PythonProjects-dna-decode/*.jsonl`
for user turns; `git log --pretty=format:'%ad|%s'` for what landed;
`dna_decode.data.cell_registry.cells()` for the live evidence surface;
`wiki/decoder_validation_report_card.json` for the AMR arm specifically.
