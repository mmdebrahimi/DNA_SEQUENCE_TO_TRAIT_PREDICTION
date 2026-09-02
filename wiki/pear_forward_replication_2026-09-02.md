# The forward cell transfers to a second β-lactamase — at less than half strength

**External replication of the one learned regime that works.** The forward cell was validated on TEM-1 +
ampicillin (genome-edit path, ESM2, Spearman **0.761** over 1,715 variants). PEAR is the same shape on a
**different β-lactamase (CTX-M-14)** with **different drugs (cefotaxime, ceftazidime)**, measured by an
independent lab. Nothing was re-fit: `predict_genome_edit` was called exactly as shipped.

| method | cefotaxime (CTX) | ceftazidime (CAZ) | n missense |
|---|---|---|---|
| BLOSUM62 | 0.198 | 0.020 | 1,513 |
| **ESM2-650M** | **0.352** | 0.078 | 1,513 |
| *TEM-1 + ampicillin, ESM2 (the cell's own validation)* | *0.761* | — | *1,715* |

**Two things are true at once, and both matter.** The learned model earns its keep — ESM2 nearly doubles
BLOSUM62 (0.352 vs 0.198), exactly as the regime map predicts for constructed→molecular. And **0.761 does
not transfer**: on a second enzyme in the same family, the same shipped path reaches less than half of it.

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

## Honest limits

- **ESM2-650M only.** ProSST/GEMME/hybrid were not run; the modality-hybrid work found rank-averaging
  orthogonal modalities beats ESM2 alone on 84–90% of proteins, so 0.352 is a floor for this substrate,
  not the best achievable.
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
