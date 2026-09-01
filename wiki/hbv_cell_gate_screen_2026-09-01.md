# HBV screened as a 6th viral cell — NO-BUILD, on G1

**Verdict: do not build it.** The engine fits, the biology fits, the catalog is curatable — and there is
**no free measured-phenotype label source** to validate it against. Every free HBV resource is either an
*interpretation system* (a rule) or a *prevalence annotation*, both of which fail the circular-label
gate. An HBV cell would ship unvalidatable from day one.

Screened 2026-09-01 after the idea surfaced from a hepatitis-B explainer video. **The video is not
evidence for any of this** — it prompted the question; the screen answered it.

---

## Why it looked promising

HBV replicates through **reverse transcription** — its polymerase carries an RT domain. Nucleos(t)ide-
analogue resistance therefore sits at catalogued RT positions, which is **structurally identical** to the
HIV RT cell that is this project's strongest: same target-site shape, and the shipped
`HIVTargetClass` engine plus the codon-mapper would carry over almost unchanged. It would be the 5th
viral kingdom after HIV-1, SARS-CoV-2, influenza NA and HCMV.

Confirmed absent first: `grep` for hepatitis/HBV across `dna_decode/` and `scripts/` returns nothing,
and the registry's viral track is HIV-1 + SARS-CoV-2 + HCMV only.

## The screen: G1 (circular label)

| free resource | what it actually provides | verdict |
|---|---|---|
| geno2pheno[HBV], HIV-GRADE HBV tool | rules-based genotype→interpretation | **G1 FAIL** — a rule, not a measurement |
| Stanford HBVseq / HBVrtDB | mutation **prevalence** by genotype and treatment | **G1 FAIL** — prevalence ≠ phenotype; also dormant (founded 2010, last update 2012) |
| HBVdb (Lyon) | sequence knowledge base + annotation/genotyping | not a phenotype source |
| EASL / AASLD guidance tables | expert consensus signature mutations | a catalog SOURCE, not a validation label |

**The field states the cause directly:** phenotypic drug-resistance testing for HBV is especially
labor-intensive because there is no simple cell-culture system, so the genotype–phenotype correlation
dataset that underpins a Stanford-style scoring algorithm is far thinner for HBV than for HIV.

That is the whole difference. The HIV cell is this project's only `INDEPENDENT_MEASURED` viral arm
**solely** because PhenoSense fold-change exists as a free, independent, isolate-level *wet-lab*
measurement, and validating against it — never against HIVDB's own Sierra interpretation — is what made
it non-circular. HBV has no such artifact.

## Why this is a NO-BUILD and not a "build it as a knowledge baseline"

I predicted HBV would resemble **HCMV**: a curated measured-fold-change compilation (Chou's recombinant
marker-transfer tables) giving an honest `KNOWLEDGE_BASELINE` cell. **That prediction is wrong**, and the
difference matters: HCMV at least has *measured* per-mutation fold-changes behind its catalog. HBV's free
tier is rules and prevalence, so an HBV cell would be one tier lower again — curated from guidance tables
with nothing measured behind it.

That is exactly the pattern the **colour-cell freeze** exists to stop: 19 cells built to
`KNOWLEDGE_BASELINE` before anyone asked whether they *could* be validated, 7 of which turned out
unvalidatable as written. The lesson recorded there is to screen G9/G10 — and G1 — **before** building.
This is that screen doing its job on the first candidate since.

## What would change the verdict

A free, isolate-level, **measured** susceptibility set — in vitro EC50/fold-change per RT variant, or a
clinical virologic-response cohort with per-isolate genotypes — released openly. The published
per-mutation in vitro work exists but is scattered across primary literature rather than compiled into a
fetchable table the way Chou's HCMV compilations or Stanford's HIV datasets are. **Compiling one is a
curation project, not a decoder project**, and it carries the fabrication hazard the colour-cell memo
names: entries must be sourced per-mutation or not written.

## Honest limits of this screen

Bounded, not exhaustive. Two search framings were safeguard-blocked and rephrased, so the sweep is
narrower than intended; a paid/registration-gated phenotype resource would not have surfaced. I did not
attempt to enumerate primary-literature in vitro studies. The verdict is therefore *"no free compiled
measured-phenotype source was found"*, not *"none exists anywhere"*.

## Sources

- [Hepatitis B Virus Drug Resistance Tools: One Sequence, Two Predictions](https://karger.com/int/article/57/3-4/232/178641/Hepatitis-B-Virus-Drug-Resistance-Tools-One) — the commercial-vs-free interpretation-system landscape
- [Stanford HBVseq / HBVrtDB](https://hivdb.stanford.edu/HBV/HBVseq/development/HBVseq.html) — prevalence-annotation scope; dormant since 2012
- [HBVdb: a knowledge database for Hepatitis B Virus](https://academic.oup.com/nar/article/41/D1/D566/1051781)
- [geno2pheno[ngs-freq]](https://academic.oup.com/nar/article/46/W1/W271/4990638)
- [Rationale and Uses of a Public HIV Drug-Resistance Database](https://pmc.ncbi.nlm.nih.gov/articles/PMC2614864/) — what HIV's ~600-study correlation base rests on, and why HBV lacks the equivalent
