# Hybrid Learned + Deterministic Decoder — multiphase plan (2026-07-03)

**Goal (user, verbatim intent):** build the best version of the "hybrid" idea — combine masked-prediction /
representation-learning (CLIP/MLM/world-model style) with the validated deterministic decoder so the whole is
stronger than either alone.

**Framing (load-bearing — this plan is EVIDENCE-GROUNDED, not speculative):** the naive forms are already
tested and bounded. This plan is built on three established results, so it targets only the forms that can
still win:

1. **Pretrained-FM embeddings → phenotype = 0-for-5 de-confounded** (learn population/lineage structure).
   Genotype-level self-supervised representation is confounded by LD. → the genotype world-model is deferred,
   pure-SSL form rejected.
2. **Masked protein-LM (ESM) zero-shot resistance signal is PARTIAL + POSITION-MECHANISM-DEPENDENT**
   (`wiki/esm_hiv_resistance_matched_test_2026-07-03.md`): AUC 0.82 (NRTI/3TC, active-site) & 0.72
   (INSTI/RAL, core) where resistance rides structural constraint; AUC 0.24 (NNRTI/EFV, tolerant pocket)
   where it rides drug-pocket evasion. So a fitness masked-LM is a **mechanism-gated partial** resistance
   signal — NOT a blanket win, NOT a blanket blind. (The earlier Mpro blanket-negative was an overclaim from
   an underpowered probe; corrected + superseded.)
3. **Masked-LM / conservation WINS on molecular FUNCTION** (DMS ~0.5) — already captured
   (conservation/AlphaMissense arc). Not re-litigated.

**North-star fit:** the deterministic decoder is the deployed product; every learned component here must BEAT
the deterministic baseline on INDEPENDENT labels to earn inclusion — else it stays out. The frozen AMR
surface is the floor and is never regressed.

---

## The three value-add hypotheses (what a learned layer could add to the deterministic decoder)

| # | Deterministic-decoder weakness | Learned component | Evidence it could work |
|---|---|---|---|
| **V1** | Catalog INCOMPLETENESS — misses novel variants at KNOWN positions (SARS Mpro FN; HIV mutant-level gaps) | a learned variant-effect SCORER (supervised on fold) that scores uncatalogued variants | fair test: zero-shot already 0.72–0.82 at conserved-position resistance; supervised should extend to pocket classes |
| **V2** | ABSTAIN on uncallable/missing input (TB regeno callability, ABO `--` no-calls, low-cov genome-map) | a masked-genotype IMPUTATION model to fill uncallable determinant positions | LD-learning is a FEATURE for imputation (the one place the confound is desirable) |
| **V3** | Position-based OVER-CALL (HIV v0 over-calls benign pocket polymorphisms → low spec) | a learned R-vs-benign discriminator at known positions | the v0.1 mutant-specific catalogs did this manually; a model could generalize |

---

## Phases (dependency DAG + per-phase falsifier + pre-committed verdict branches)

### Phase 1 — Fair zero-shot signal test  ✅ DONE (2026-07-03)
- **Claim:** does zero-shot masked-LM representation carry ANY resistance signal, fairly measured?
- **Result:** YES but PARTIAL + mechanism-gated (best powered AUC 0.821; NNRTI anti-predictive 0.244).
  `wiki/esm_hiv_resistance_matched_test_2026-07-03.md`. **Gate for Phase 2: PASSED** (signal exists to build on).

### Phase 2 — Supervised resistance head (THE decisive learned-scoring test)  ✅ DONE 2026-07-03 → PARTIAL
- **Result:** `PARTIAL` (`wiki/esm_supervised_head_result_2026-07-03.md`). Leave-one-position-out head BEATS
  zero-shot ESM (pooled AUC 0.691 vs 0.587) but does NOT beat the deterministic catalog (balacc 0.590 vs
  0.783). On KNOWN positions, curated knowledge > learning-from-embeddings; the head does not recover NNRTI
  pocket resistance to a deployable level. → **V1 = novel-variant FALLBACK only, not a catalog replacement.**
  Phase 4 integrates it fail-closed on the catalog-SILENT subset only; Phase 5 stays deferred.

- **Terminal claim:** a supervised head on ESM embeddings (logreg / shallow MLP, OR LoRA-fine-tuned ESM)
  trained on PhenoSense fold BEATS both (a) the deterministic catalog and (b) zero-shot ESM, on HELD-OUT
  variants — INCLUDING recovering the NNRTI pocket signal zero-shot misses.
- **Substrate:** HIV single+multi-mutant PhenoSense (rich R+S, all 4 classes); features = ESM per-residue
  embeddings at the variant position(s) (± window).
- **Method:** leave-one-POSITION-out and leave-one-DRUG-out CV (NOT random split — must prove generalization
  to unseen resistance sites, not memorization). Baselines scored on the SAME folds: deterministic catalog +
  zero-shot ESM `-LLR`.
- **Falsifier (pre-committed):**
  - **WIN** (supervised beats catalog by ≥3pp balacc AND beats zero-shot ESM on NNRTI) → V1 greenlit; proceed
    to Phase 4 integration.
  - **PARTIAL** (beats zero-shot but not the curated catalog) → learned layer is a novel-variant FALLBACK only
    (score variants the catalog is silent on), not a replacement.
  - **FAIL** (does not beat the catalog on held-out positions) → learned scoring does not add deployable value;
    V1 closed; the deterministic catalog remains the sole scorer.
- **Cost:** ESM embeddings already computable (torch+esm installed); training is a shallow head (CPU-fine).
- **Independent-label note:** in-distribution first (HIVDB); an external-cohort fold confirmation mirrors the
  decoder's own provdisjoint/external arc IF Phase 2 wins.

### Phase 3 — Imputation hybrid (input completeness; INDEPENDENT of Phase 2)  [class d/e]
- **Terminal claim:** a masked-genotype imputation model reduces the deterministic decoder's ABSTAIN rate on
  uncallable determinant positions WITHOUT lowering accuracy vs abstaining.
- **Substrate:** TB regeno-callability gaps (ABSTAIN vs S-by-absence) + ABO `--` no-calls + low-coverage
  genome-map determinant windows. Imputation model = a standard genotype imputation approach (Beagle-class
  or a small masked-genotype net trained on the cohort's own variant matrix — LD IS the desired signal here).
- **Falsifier (pre-committed):** on held-out MASKED-known positions, imputation accuracy ≥ a stated bar AND
  imputing-then-calling ≥ abstaining on downstream R/S accuracy. If imputation accuracy is below the bar → the
  decoder keeps abstaining (honest ABSTAIN > wrong impute); V2 closed.
- **Why LD is OK here:** V2 uses LD to fill MISSING GENOTYPE, not to predict phenotype — the confound that
  kills phenotype-prediction is exactly the signal that makes imputation work.

### Phase 4 — Mechanism-aware hybrid decoder + end-to-end validation  [class d; gated on P2 and/or P3]
- **Terminal claim:** a hybrid = deterministic catalog (primary) + learned scorer (ONLY where Phase 2 proved
  it beats baseline, mechanism-gated — e.g. conserved-position novel variants) + imputation (input) ≥ the pure
  deterministic decoder on INDEPENDENT cohorts, and NEVER worse (the frozen surface is the floor).
- **Method:** the hybrid is a ROUTING layer, not a replacement: deterministic call stands where the catalog is
  confident; the learned scorer only fires on catalog-silent variants AND only in mechanism classes where
  Phase 2 showed a win (NOT NNRTI-pocket). Validated on provdisjoint + external cohorts via the existing
  harness.
- **Falsifier (pre-committed):** hybrid balacc ≥ deterministic balacc on every scored cell (non-inferiority),
  with ≥1 cell strictly improved. Any cell where hybrid < deterministic → that route is disabled (fail-closed
  to deterministic). If no cell improves → the hybrid adds no value; ship deterministic-only (honest negative).

### Phase 5 — Genotype world-model  [DEFERRED / conditional; low priority]
- Built ONLY if Phases 2–4 show learned representations add deployable value AND a specific gap remains that a
  GENOTYPE-level (not protein-level) model uniquely fills. The pure self-supervised form stays rejected
  (LD-confound). The only viable shape is a supervised/phenotype-conditioned genotype model — evaluate its
  necessity after Phase 4, not before.

**Dependency DAG:** P1(done) → P2 ; P3 (parallel, independent) ; (P2 win ∪ P3 win) → P4 ; P4 → P5(conditional).

---

## Cross-cutting invariants (non-negotiable)
- **Frozen AMR surface byte-unchanged** (`amr_rules.py` + `calibrated_amr_rules.json`; leak guard 9/9) — every
  learned component is ADDITIVE + fail-closed to the deterministic call.
- **Beat-the-baseline gate** — no learned component ships unless it beats the deterministic catalog on
  INDEPENDENT labels (the validate-wrapper-vs-underlying-tool rail). In-distribution wins are necessary, not
  sufficient.
- **Held-out-by-position/drug CV** — never random-split (would memorize resistance sites).
- **Honest verdict** — every phase pre-commits WIN/PARTIAL/FAIL branches before results; a FAIL ships the
  deterministic-only decoder, which is already the validated product.

## Success criterion (whole plan)
Either (a) ≥1 learned component (V1 supervised scorer OR V2 imputation) ships as an additive, fail-closed,
baseline-beating hybrid layer with an independent-label win; OR (b) a documented negative that the learned
components do not beat the deterministic decoder on independent labels — a real finding that closes the
hybrid question with evidence, leaving the deterministic product as the terminal.

## State of the plan (2026-07-03)
Phase 1 ✅ (partial zero-shot signal, mechanism-gated) → Phase 2 ✅ **PARTIAL** (head beats zero-shot, loses to
the catalog on known positions). **The evidence-grounded conclusion so far:** a learned protein-level scorer
is a BOUNDED FALLBACK (catalog-silent variants), NOT a catalog replacement — the deterministic decoder stays
primary. The grand "learned beats deterministic" form is closed on independent-label held-out tests.

**Remaining options (user's call — the hybrid's honest ceiling is now known):**
- **Phase 4 (fallback-only)** — evaluate the head on the catalog-SILENT subset specifically; wire it as a
  fail-closed fallback ONLY if it beats the K/N null there. Modest, bounded value.
- **Phase 3 (imputation)** — the OTHER, independent value-add (reduce ABSTAIN rate); untouched by the Phase-2
  result and arguably the higher-VOI remaining move (LD-learning as a feature, not competing with the catalog).
- **Bank** — the hybrid question is substantively answered: learned scoring does not beat curated knowledge
  on known positions; the deterministic decoder is the terminal product. Phase 3 imputation is the only
  remaining branch that could add deployable value.
