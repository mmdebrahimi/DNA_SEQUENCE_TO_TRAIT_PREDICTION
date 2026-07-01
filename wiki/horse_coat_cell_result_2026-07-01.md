# Horse coat-colour cell — rule built; validation hits the anti-circular data wall (2026-07-01)

Greenlit build of the label-first scan's winner (`wiki/label_first_scan_2026-07-01.md`): the horse base
coat-colour decoder — the first off-pathogen cell targeting a **non-human, observed (not self-report)**
label. Outcome: the RULE cell is built + tested; the non-circular real-data SCORE is blocked by a data wall.

## What shipped (durable)
- **`dna_decode/data/horse_coat.py`** — the deployed Mendelian two-locus rule (Rieder 2001 / UC Davis VGL):
  MC1R E/e epistatic to ASIP A/a → `e/e → chestnut`, `E_ A_ → bay`, `E_ a/a → black`. Molecular basis pinned
  (chestnut = MC1R c.901C>T; black = ASIP 11-bp exon-2 deletion).
- **`scripts/horse_coat_validate.py`** — validation harness that scores a user-provided **observed-colour**
  TSV (mc1r/asip/observed_colour) and otherwise reports `VALIDATION_DATA_WALL` (never a circular/synthetic
  number). Surfaces discordant cells (rule-breakers) — the only informative signal for a deployed rule.
- 7 tests (`tests/test_horse_coat.py` 5 + `tests/test_horse_coat_validate.py` 2). Frozen AMR surface
  byte-unchanged (leak guard 9/9).

## The wall: no headlessly-fetchable NON-CIRCULAR per-individual dataset

The anti-circular gate (colour must be *independently observed*, not assigned *from* genotype) eliminated
every reachable dataset:

| Source | Why it fails |
|---|---|
| Dryad `10.5061/dryad.3q111` (N=215) | colour was *"determined based on the genotypes"* → **CIRCULAR**; also **auth-gated** (401 headless, even via the v2 API) |
| Rieder 2001 (120 horses, observed) / 709-horse Synergy / Noma CIE-Lab | genuinely observed + non-circular, but **PDF/paywalled supplementary tables** — not headlessly fetchable as structured per-individual data |
| Figshare / Zenodo / GitHub | no matching open observed-colour dataset found |
| Arabidopsis AraPheno (open API) | the open traits are **quantitative** (flowering/trichome), not a clean binary Mendelian rule |

## Verdict + what this proves

This **re-confirms the label wall at fine grain**: even the *best* off-pathogen gate-candidate — the cleanest
Mendelian rule there is — lacks a clean OPEN, NON-CIRCULAR, per-individual observed-phenotype dataset. The
one clean open CSV is circular by construction; the non-circular data is PDF-locked. This is fork #2's
premise, now demonstrated rather than assumed.

**Honest classification:** the rule cell is DEPLOYED-RULE INTEGRATION (VGL sells these tests) — like
IrisPlex/ABO — so even a successful score would re-confirm a known tool, not a novel finding. The build
honours the greenlight; the score is an EXTERNAL wall.

## Unblock (if pursued later)
Supply an independently-observed dataset: `uv run python scripts/horse_coat_validate.py --data <TSV:
mc1r,asip,observed_colour>`. Realistic sources: transcribe the Rieder 2001 / 709-horse observed contingency
from the papers (manual, small), or a future open observed-colour repository. NOT the Dryad file (circular).

## Recovery attempts EXHAUSTED (2026-07-01, second pass) — the wall is demonstrated across 5 avenues

Per the user's "find a way", a second pass actively hunted a non-circular observed dataset to score. All 5
reachable avenues fail the clean-non-circular-per-individual bar:

| Avenue | Fails because |
|---|---|
| Dryad `10.5061/dryad.3q111` per-individual CSV | colour genotype-DERIVED (circular) + auth-gated (401) |
| 709-horse Synergy open contingency (PMC6447268, Table 2) | **marginal, not joint** (MC1R + ASIP counts reported SEPARATELY per colour → per-horse pairing unrecoverable) **+ non-functional 2-SNP diplotype coding** |
| Rieder 2001 / Noma (joint FUNCTIONAL genotype × observed colour) | PDF/paywalled — not headlessly fetchable |
| Figshare / Zenodo / GitHub | no matching open observed-colour dataset |
| Arabidopsis AraPheno (open API) | traits quantitative — no clean binary Mendelian rule |

**The 709-study finding is decisive + informative:** its colour IS independently observed (photos +
registration, non-circular), but its genotype notation (`EeEe`/`AaAa` = heterozygous at TWO linked SNPs) is
NOT the functional extension/agouti allele — proven because **all 31 observed BLACK horses are `AaAa`**
(would be functional `a/a`) and **all 123 CHESTNUT are `EeEe`** (would be functional `e/e`): het-at-both ≠
homozygous-recessive. So it cannot be mapped to the C901T / 11-bp-deletion functional rule. Even the best
OPEN, NON-CIRCULAR, observed contingency is unusable for the functional rule.

## Recommendation
**Bank (fork #1) — the wall is now demonstrated, not assumed.** Across 5 distinct open-data avenues, the
cleanest off-pathogen Mendelian rule has NO clean, non-circular, per-individual, functional-genotype ×
observed-colour open dataset. This is the LABELS wall at its finest grain. The deterministic decoder's honest
off-pathogen ceiling: eye/ABO pilots (self-report) + the built, tested, but data-walled horse rule. The ONLY
remaining unblock is a USER-supplied joint functional-genotype × observed-colour TSV (e.g. hand-transcribed
from Rieder 2001 Table 2 via the PDF) → `horse_coat_validate.py --data`. Not worth autonomous pursuit.
