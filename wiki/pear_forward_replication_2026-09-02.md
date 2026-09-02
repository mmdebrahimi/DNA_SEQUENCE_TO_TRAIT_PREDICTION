# The forward cell transfers to a second β-lactamase — at less than half strength

**External replication of the one learned regime that works.** The forward cell was validated on TEM-1 +
ampicillin (genome-edit path, ESM2, Spearman **0.761** over 1,715 variants). PEAR is the same shape on a
**different β-lactamase (CTX-M-14)** with **different drugs (cefotaxime, ceftazidime)**, measured by an
independent lab. Nothing was re-fit: `predict_genome_edit` was called exactly as shipped.

| method | cefotaxime (CTX) | ceftazidime (CAZ) | n missense |
|---|---|---|---|
| BLOSUM62 | 0.198 | 0.020 | 1,513 |
| **ESM2-650M** | **0.352** | 0.078 | 1,513 |
| ProSST-2048 (structure) | **−0.040** | 0.025 | 1,513 |
| ESM2 + ProSST rank-average hybrid | 0.204 | 0.083 | 1,513 |
| *TEM-1 + ampicillin, ESM2 (the cell's own validation)* | *0.761* | — | *1,715* |

**Two things are true at once, and both matter.** The learned model earns its keep — ESM2 nearly doubles
BLOSUM62 (0.352 vs 0.198), exactly as the regime map predicts for constructed→molecular. And **0.761 does
not transfer**: on a second enzyme in the same family, the same shipped path reaches less than half of it.

**And a third: the modality hybrid does NOT transfer either.** The published finding is that a naive
rank-average of orthogonal modalities beats ESM2-650M on 84–90% of ProteinGym proteins. Here it **loses**,
0.204 vs 0.352 — because the premise fails: **ProSST alone is at chance (−0.040)**, so rank-averaging it
in halves the signal ESM2 had. Structure is not an orthogonal modality on this protein; it is noise.

That was worth measuring rather than assuming. The prior expectation was that the hybrid would *raise*
0.352 and might change this memo's conclusion. It did the opposite.

---

## The pipeline validated itself before any correlation was read

A coordinate error here produces plausible garbage rather than an exception, so the run asserts something
the data must satisfy regardless of whether the predictor works:

| consequence | n | median CTX | median CAZ |
|---|---|---|---|
| silent | 538 | **+0.0009** | +0.0052 |
| missense | 1,513 | −0.0457 | +0.0180 |
| nonsense | 63 | **−0.1286** | −0.0221 |

Silent ≈ exactly neutral, nonsense strongly deleterious, missense in between — on both drugs. No
frame-shifted or offset mapping produces that ordering by accident.

**Coordinates were established by measurement, not assumption.** PEAR's `C648T` notation is 1-based on the
authors' own 795-nt reference: all **2,114/2,114** variants satisfy `ref[pos-1] == wt`, while the `+81`
convention used by their own Figure-2 axis matches only **452/2,114** (chance). That reference is the
**mature** protein — the gene's 81-nt signal peptide is trimmed — so residue numbers from this run are
**not Ambler numbering** and must be converted before being compared to the literature.

## Ceftazidime's near-zero is structural, not a shrug

CAZ ρ=0.078 looks like failure and is better read as a mismatch of question. CTX-M-14 is a
**cefotaximase**: it hydrolyses cefotaxime efficiently and ceftazidime poorly. The measured distribution
says exactly that — CAZ has a **tighter core** than CTX (IQR 0.080 vs 0.180) with a **much longer positive
tail** (max +3.68 vs +1.28).

So most CAZ variants have little to discriminate, and the variance that exists is concentrated in rare
**gain**-of-function events — which is the known route by which CTX-M enzymes extend their spectrum. **A
damage predictor predicts loss of function; it cannot predict gain.** That is a property of the question,
not of ESM2, and it is a real limit on what the forward cell can be pointed at.

## What the ceiling is — unknown, and I checked

Before reading 0.352 as "weak", the obvious question is how much signal is achievable at all. I tested it
by correlating the two published tables (Figure 2B per-nucleotide vs Figure 3A per-variant) as a proxy for
measurement reproducibility.

**Spearman = 1.0000 exactly, over all 2,114 variants.** That is not reproducibility — it means the two
tables are the *same numbers under a monotone transform*, a different representation rather than a
replicate. **No noise-ceiling estimate is available from the published artifacts**, so 0.352 stands
uncorrected against 0.761 and the gap cannot be attributed to assay noise on this evidence.

## What this changes

- **It bounds a shipped claim.** "The forward cell reaches Spearman 0.761 on genome edits" is now known to
  be **protein-specific**, not a general property of the path. The honest statement of the cell's range is
  0.35–0.76 across two β-lactamases, and the low end came from an independent lab.
- **It does not overturn the regime.** Constructed variation → molecular endpoint still works: the
  correlation is positive, significant in size, and the learned model beats the substitution-matrix
  baseline. That is the regime's claim, and it holds.
- **It is genuinely independent.** Different enzyme, different drugs, different lab, no re-fitting, and a
  substrate that cleared all ten rejection gates first.

## Why ProSST is at chance — the structure is exonerated, the toolchain is not

The obvious suspect is a bad structure. AlphaFold has **no model for any CTX-M-14 UniProt entry** (all 9
that contain the mature sequence 404 — AlphaFold DB covers reference proteomes and these are plasmid-borne
entries from clinical isolates), so the structure was **folded with ESMFold v1 from the mature sequence**.
That is also what removed the numbering risk entirely: the structure is of exactly the sequence being
scored, so the offset is 0 by construction rather than a cross-database alignment to verify.

**The fold is high-confidence: mean pLDDT 0.954, with 95.5% of residues above 0.90 and 99.6% above 0.70**
(only the N-terminal residue is low, 0.53). A bad structure does not explain a chance result.

**What is NOT excluded:** `torch_cluster` failed to install on the Kaggle image (`OSError`), so
`torch_geometric` fell back for graph construction inside ProSST's quantizer. The quantizer path itself
was validated locally on GRB2 (self-quantized == ProteinGym's pre-quantized tokens, 217/217), but **not
with that fallback and not on this protein**. So the honest statement is *ProSST scored at chance here,
and a quantizer-fallback artefact is the leading un-excluded alternative to "structure adds nothing on
CTX-M-14"*. Resolving it needs a working `torch_cluster` build, or ProteinGym-style pre-quantized tokens.

## Honest limits

- **The hybrid was RUN and it lowered the number** — 0.352 is not a floor, as this memo first assumed.
  GEMME (the evolution modality) was still not run; it needs an MSA pipeline.
- **ProSST's chance result carries the quantizer-fallback caveat above.**
- **No noise ceiling**, as above.
- **Aggregated effect sizes**, ~2,100 per drug — not the ~23,000 raw barcoded strains.
- **Comparability caveat:** the TEM-1 number came from a different DMS with its own normalization. Both
  are Spearman on single-nucleotide-accessible missense sets, which is the closest available match, but
  they are not the same assay.

## Reproduce

```bash
uv run python scripts/pear_forward_replication.py --method blosum62
HF_HOME=D:/hf_cache uv run --with torch --with transformers \
  python scripts/pear_forward_replication.py --method esm2 \
  --esm-table D:/dna_decode_cache/esm/esm2_t33_650M_UR50D__ctxm14_mature.json
```

A trap the runner now refuses rather than reports: JSON round-trips the ESM table's position keys as
**strings**, `predict_effect` raises `KeyError` per variant, and the script records those as *skips* — so
an uncoerced table silently scores only the 538 synonymous variants and prints `missense n=0` instead of
failing. It now exits 3 with the skip reason.
