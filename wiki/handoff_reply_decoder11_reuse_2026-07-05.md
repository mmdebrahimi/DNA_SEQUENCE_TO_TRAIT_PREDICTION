# Handoff-back → Soraya V3-5 fact-check session: "decoder 11" resolved + §5 tightening inputs

**Reply to:** `D:\PythonProjects\ideas_research\research_outputs\handoff-fact-check-g0xl0oijkb4-2026-07-05.md` §5.
**Verified against the LIVE dna_decode repo** (`C:\Users\Farshad\PythonProjects\dna_decode`, HEAD `63a0492`,
main == origin) — the exact gap your handoff flagged ("repo grep timed out; couldn't pin decoder 11").
Every claim below re-derived from committed code this session. **Additive, not a re-do** of your video
fact-check — this only supplies the repo-side facts you couldn't reach from D:.

## 1. "decoder 11" resolved — it is NOT a session number

**"decoder 11" = the dna_decode project, named by its umbrella ledger `project_state/dna-decode-2026-05-11.md`.**
The **"11" is the DATE suffix** (`2026-05-11`), not a session/decoder index. There is no literal
"decoder 11" / "session 11" / "decoder_11" artifact anywhere in the repo (the only repo-wide string hit is an
unrelated `scipy v0.11.0` URL in site-packages).

**Real §5 targets (cite these):**
- Umbrella ledger: `project_state/dna-decode-2026-05-11.md`
- Human PGx arm (see §3): `dna_decode/pgx/` + `wiki/pgx_report_card.md`
- Deployed-surface roll-up: `wiki/decoder_validation_report_card.md`

## 2. Your verdict HOLDS — "NO direct reuse" is correct

The domain-mismatch verdict is right and gets *stronger* when grounded in the real repo. VERVE-102
(human PCSK9 base-EDIT → LDL drop, a therapeutic gene-therapy trial readout) shares no substrate with
dna_decode: no free measured genotype↔phenotype label of that kind, and it crosses the universal
"NOT a clinical tool" rail (**15 occurrences across `dna_decode/pgx/*.py` alone**, verified).

## 3. ONE correction that tightens §5 (your §5.1 premise is now outdated)

Your §5.1 says the single-nucleotide-variant→phenotype analogy is **microbial only**
(`pointfinder`/`codon_map`, `gyrA S83L → cipro R`). Those files exist as you cited
(`dna_decode/pointfinder/`, `dna_decode/typing/codon_map.py` — verified). **But dna_decode now has a HUMAN
SNV→phenotype arm** — 7 shipped PGx cells, several literally human single-variant → phenotype decoders:
- `dna_decode/pgx/vkorc1.py` — `VKORC1 rs9923231 → warfarin sensitivity`
- `dna_decode/pgx/slco1b1.py` — `SLCO1B1 rs4149056 → statin-myopathy risk`
- `cyp2c19 / cyp2c9 / cyp2c8 / cyp3a5 / cyp2b6 / tpmt` catalogs — star-allele → CPIC phenotype
  (`tpmt_catalog.py` + `compound_caller.py` = the first true compound-allele resolver, *3A=*3B+*3C)

So the human mirror you imagined **already exists in-repo** — a more accurate illustration than the
microbial one for a human-cardiovascular video.

**Why no-reuse STILL holds (sharper reasons):**
1. **Calling ≠ prediction.** Every cell (microbial + human PGx) is *curated-catalog CALLING* — a KNOWN
   variant → a KNOWN phenotype label. VERVE-102 needs predicting the *effect of a NOVEL base-edit* on LDL —
   learned therapeutic edit-effect prediction, a different target class.
2. **The learned-predictor arm is a CLOSED NEGATIVE by design** (embeddings 0-for-4 across the kingdom
   boundary; the shipped product is deterministic decoders + the `genome_map/` label-free honesty report).
   dna_decode deliberately does not do learned human variant-effect prediction.
3. **"NOT a clinical tool"** is enforced on every PGx cell (VERVE-102 is a clinical-trial readout).

**Net §5 edit:** "microbial-only SNV analogy" → "microbial **+ human-PGx** SNV analogy, but all
*catalog-calling not learned-prediction*, validated on free measured labels, non-clinical → still no reuse."

## 4. Provenance / honesty
Read-only cross-session verification (R4 value-add lane). No dna_decode source changed; this note is the only
artifact. All file paths + counts re-derived from HEAD `63a0492` this session (2026-07-05). The frozen AMR
surface + all decoder cells are untouched.
