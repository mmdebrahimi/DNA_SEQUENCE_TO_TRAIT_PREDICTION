# Dog "world model" masked-genome reconstruction — F1 engine + first real result (2026-07-31)

> **⚠️ DEPRECATED HEADLINE (2026-07-31, superseded by F1′).** The `−0.15` delta below is **not a clean
> negative** and must not be cited. An adversarial review found the comparison biased against NT on TWO
> independent axes: (1) the order-5 Markov baseline was fit on the SAME slice it scored, so each masked
> target's own count leaked into the table predicting it (transductive leakage → inflated Markov accuracy);
> (2) NT was scored by its single argmax 6-mer, discarding the per-base marginal distribution (true-token
> prob ~0 does NOT prove the true BASES have low marginal prob). Both bias against NT; the sign may flip.
> The engine + degeneracy guard remain valid. See the F1′ clean re-run
> (`wiki/dog_masked_reconstruct_clean_2026-07-31.md`): disjoint/LOO Markov + per-base marginal NLL.


**Question (user):** take one data-rich animal (dog) and see if our "world model" can mask parts of the
genome and guess them. **Framing (ratified):** dog / Framing A (a DNA masked-LM reconstructs masked
sequence — the model's native objective) / NT-v2-100M. **Why it matters strategically:** masked
reconstruction is *self-supervised* — the genome is its own label — so it **dodges the label wall** that
closed every prior track (AMR/embedding: only dog ever cleared the free-independent-label bar).

## Result — NT does NOT beat a cheap Markov baseline (honest, verified)

Substrate: **canFam4 chr1** `NC_051804.1:20,000,000–20,006,000` (6,001 bp, 100% ACGT). Smoke window
600 bp; **40 masked 6-mer tokens (240 bases)**. NT-v2-100M masked-LM (real weights), CPU.

| metric | value |
|---|---|
| **HEADLINE delta (NT − Markov-5)** | **−0.150** (NT LOSES) |
| NT per-base accuracy | 0.3375 |
| NT token (exact-6mer) accuracy | **0.0000** (0/40) |
| NT mean true-token prob | 0.0007 (~3× uniform) |
| order-5 Markov per-base accuracy (teacher-forced, fit on 6,001 bp) | 0.4875 |
| null (uniform) | 0.2500 |

Artifact: `wiki/dog_masked_reconstruct_smoke_2026-07-31.json`. **Raw accuracy is never the claim** — the
delta vs the cheap baseline is (a repetitive genome is reconstructed well by Markov with no
"understanding").

## The result is VERIFIED, not a pipeline artifact (the load-bearing check)

The `token_accuracy = 0.0` + `true_prob = 0.0007` signature is exactly what a **`trust_remote_code`
version-drift degenerate load** looks like (near-random logits). So before believing "NT loses", the
degeneracy guard `scripts/nt_mlm_sanity_probe.py` ran controls a working DNA LM must pass:

| control | true | NT pred | true-token prob | verdict |
|---|---|---|---|---|
| poly-A (`AAAA…`) | AAAAAA | **AAAAAA** | **0.9728** | HIT |
| poly-AT (`ATAT…`) | ATATAT | **ATATAT** | **0.9659** | HIT |
| real canFam4 token | TCTACC | GGATAT | 0.0000 | (near-chance) |

**VERDICT: OK — NT is NOT degenerate.** It nails periodic sequence near 1.0, so the harness alignment is
correct and the near-chance behavior on *real* sequence is a genuine capability ceiling. My initial
suspicion triggered the right check; the check *cleared* the −0.15 rather than invalidating it.

## Interpretation (honest, scoped)

- **NT-v2-100M works, but at 6-mer-token resolution it is near-chance on real non-repetitive canFam4
  sequence and loses to order-5 Markov.** Predicting a specific 6-mer (1 of 4,104) from flanking 6-mers
  is genuinely hard; the model's low argmax confidence (0.03) reflects real uncertainty, not brokenness.
- **Consistent with the project's standing finding** that these DNA foundation models track composition /
  structure more than specific sequence — now shown in the **pure self-supervised reconstruction** regime,
  with no phenotype label anywhere. It is NOT another embedding-vs-phenotype negative; it is a *native-task*
  measurement.
- **Named measurement caveat (fair-comparison, F2 to fix):** NT emits ONE argmax 6-mer token (all 6 bases
  from a single token choice), while Markov predicts each base independently from the TRUE left context
  (teacher-forced — the strongest cheap baseline). This token-vs-base asymmetry partly explains the gap; a
  fairer NT per-base metric marginalizes the token logits to per-base distributions (F2 hardening).
- **Scope:** this is a SMOKE (40 tokens, one 600 bp window). A region-stratified full-chromosome sweep
  (coding vs intergenic vs conserved; the coat/size causal loci) is F2/F3.

## Environment finding (repo knowledge)

**NT-v2-100M is UNLOADABLE under the repo's `transformers==5.14.1`** — the vendored `modeling_esm.py`
(trust_remote_code) imports `find_pruneable_heads_and_indices` (removed) and reads `config.is_decoder`
(gone) → cascading failures. Shimming the import only surfaced the next break. The real-surface run
required an **isolated `transformers==4.30.2` env** (`uv run --isolated --with …`), reusing the cached
weights at `D:/hf_cache`. Do NOT downgrade the repo env (the `forward`/ESM cells need 5.x). This matches
the standing `trust_remote_code library-version drift` lesson — the degeneracy guard is the mandated
defense, and it fired correctly here.

## F1 bar — MET

- ✅ harness runs end-to-end on real canFam4 via NT-v2-100M → per-position reconstruction table (real surface)
- ✅ headline = NT − Markov delta (raw accuracy explicitly not the claim); null + token-acc + true-prob reported
- ✅ offline tests green with the mock model (`tests/test_masked_reconstruct.py`, 7/7)
- ✅ does not assert NT wins — it in fact loses, honestly reported
- ✅ degeneracy guard proves the real number is trustworthy (poly-A 0.97), per the R3 real-surface + version-drift discipline
