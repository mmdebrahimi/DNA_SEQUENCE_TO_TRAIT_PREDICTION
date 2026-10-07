# The GLM: project-family decomposition + requirements flow-down

**Date:** 2026-10-06
**Mission:** build the generative genome-edit engine — *name a trait, get the genome edits*.
**Design basis:** `wiki/glm_architecture_design_2026-10-06.md` (written first, deliberately).
**Standing constraint:** no catalog comparison. See orientation error 0 in `CLAUDE.md`.

---

## 0. The design decision that shapes every family: retrieve the CDS, GENERATE the regulation

Four alternative architectures were steelmanned before fixing the decomposition. Three were rejected for
stated reasons; the fourth changed the design.

| Alt | Architecture | Verdict |
|---|---|---|
| A | End-to-end trait→edit net trained on (edit, measured phenotype) pairs | **Rejected — correctly label-blocked.** This is the one thing the label wall genuinely stops. |
| B | Frontier-LLM agent reasoning over genome edits with tools | **Rejected as the engine** (no validation path, hallucinates biology) — but it IS the right tool for the trait-resolver family, where the task is symbolic reasoning over literature. |
| C | RL (GRPO/PPO) fine-tuning a DNA LM directly against the oracle | **Deferred as OPTIMIZATION, not v0.** FBA is seconds per solve and ESM2 scoring is slow, so neither survives an RL inner loop at scale. Best-of-N reranking first — which is exactly the user's *"then we will think of optimization"*. |
| D | Discrete diffusion over DNA instead of autoregressive | **Noted.** Plausible and active; binding constraint is weight availability, not merit. |
| **E** | **Retrieval-augmented genome design** | **ADOPTED — it changes the generative requirement.** |

**Alt E, stated plainly.** For *"make this host produce protein X"*, some organism already makes X. So the
engine should **retrieve** the natural cassette and **adapt** it to the host, not invent an enzyme. De-novo
CDS generation is the rare hard case, not the common one.

Consequences, and they are all favourable:

1. **The generative surface becomes SHORT sequence** — promoters, RBS, 5'UTRs, codon-adapted CDS variants.
   Short, high-variance, and the part that actually governs expression.
2. **It fits the available compute.** A 4 GB CC-5.0 laptop cannot host a 7B genome model; generating and
   scoring ~100 bp regulatory elements is tractable locally, and the free T4 covers the rest.
3. **Retrieval is already in-repo**: BLAST+ 2.17.0 is installed natively, `genome_map/` annotates,
   `dna-identify` types an organism, and `fba/gapfill.py` already does donor-model retrieval at the
   *reaction* level — this extends the same instinct to the *sequence* level.
4. **It is honestly scoreable.** A retrieved-and-adapted cassette has a real wild-type reference to compare
   against, which a de-novo invention does not.

---

## 1. Project-family decomposition

| # | Family (slug) | Scope (one line) | Why a distinct family |
|---|---|---|---|
| 1 | `glm-edit-representation` | ONE `GenomeEdit` type expressing all four edit classes as actual DNA changes on a host genome: apply, invert, round-trip, and emit the constraint-layer `ctx` | Keystone data contract. Pure, no network/GPU/labels. Wrong here poisons every downstream number silently — this repo's most-repeated failure shape |
| 2 | `glm-sequence-generator` | A generative model over short DNA (promoter/RBS/UTR/codon variants) that runs on available compute, with a self-supervised quality floor vs a Markov baseline | External weights + licence + VRAM questions; GPU-bound; own falsifier that needs no labels |
| 3 | `glm-expression-oracle` | Predict the expression level a regulatory sequence yields in a host, from a measured dataset if one is free, else a biophysical model | A *scoring* problem with its own data-availability question; separable from generation and gates the honesty of 4 and 5 |
| 4 | `glm-construct-assembler` | reaction/gene/protein target → a buildable cassette (promoter + RBS + codon-adapted CDS + integration locus), via retrieval-then-adapt | The integration layer that converts the repo's existing organism-level design into an actual answer to "give me the edits" |
| 5 | `glm-design-loop` | The propose→score→select loop: oracle committee, trust region, disagreement reporting, and the oracle-hacking guards | Where the novel contribution AND the central risk live. Its falsifier is a held-out oracle, not a metric |
| 6 | `glm-trait-resolver` | trait → scoreable target, for traits that are not a metabolite or a named protein | PROPOSED, not accepted. This is the climb to animals; its blocker is conceptual, and families 1-5 are its prerequisites |

Families 1-5 are proposed as ACCEPTED. Family 6 is explicitly parked so it does not consume attention
before its prerequisites exist.

---

## 2. Requirements flow-down + dependency graph

| Family | Binding requirement | Depends on | Critical-path position | Cheapest next test |
|---|---|---|---|---|
| `glm-edit-representation` | must express knockout, substitution, insertion AND expression-tuning without a special case, and survive apply→invert on a real genome | — | **Gate A (keystone)** | round-trip all four classes on the committed MG1655 reference; assert the constraint layer REFUSES a frame-breaking edit |
| `glm-sequence-generator` | downloadable weights + a licence permitting this use + inference inside 4 GB local or 16 GB T4 | — | **Gate A (parallel)** | generate N promoters; measure constraint-pass rate + codon/GC realism **against a 3-mer Markov baseline** — if Markov matches it, the model adds nothing |
| `glm-expression-oracle` | a free measured regulatory-sequence→expression dataset, or an accepted biophysical fallback | — | **Gate A (parallel, highest info)** | resolve whether such a dataset is public and downloadable **before** building anything on it |
| `glm-construct-assembler` | a retrieved cassette must pass constraints AND the expression oracle AND reproduce a known case | 1, 2, 3 | **Gate B** | emit a construct for the measured `FFSD` sucrose case (`gapfill` took growth 0.000 → 1.7798 /h) and check it is constraint-clean |
| `glm-design-loop` | must beat a random-proposal baseline on the committee **and** survive a held-out oracle that never scored during search | 1, 2, 3, 4 | **Gate C** | random-vs-generated proposal A/B under the committee, with the held-out oracle as the audit |
| `glm-trait-resolver` | a trait must map to something an existing oracle can score | 1-5 | **Phase 2 (parked)** | — |

**Critical path:** `Gate A {1 ∥ 2 ∥ 3} → Gate B {4} → Gate C {5} → Phase 2 {6}`

**Cross-cutting (not a family):** the oracle-hacking guard discipline — hard constraints first, committee
over a single scalar, trust region, held-out oracle, report disagreement — applies inside 4 and 5 rather
than being its own workstream.

---

## 3. VOI-ranked first moves

Ranked by information retired per unit cost:

1. **Resolve the expression-dataset question (family 3's unknown).** Highest info: a *yes* gives the whole
   stack a real supervised target and a genuine held-out number; a *no* changes family 3's shape and
   weakens the honest claim of 4 and 5. Cost: a search. **In flight.**
2. **Resolve which generator is actually runnable (family 2's unknown).** A *no* across all candidates
   would force the design onto retrieval + biophysics only — still useful, but a different claim.
   Cost: a search plus one load test. **In flight.**
3. **Build the keystone (`GenomeEdit`).** Depends on neither answer, unblocks everything, and is pure
   enough to be fully tested offline. **This is the first thing to BUILD.**
4. Everything else is downstream of those three.

Moves 1 and 2 are already running as parallel searches, so move 3 is where execution starts.

---

## 4. Terminal condition

**Mission terminal (mechanical part):** families 1-5 each meet their recorded MVP criteria, and the
`glm-design-loop` family additionally produces a run where the **held-out oracle agrees with the committee's
selection** above the random-proposal baseline — i.e. the engine proposes edits that survive a scorer it
never optimized against.

**Not auto-closable.** Whether the engine is *good enough to act on* is the user's call, and no wet-lab
claim is implied by any of the above: a proposal is a hypothesis for the bench, the same stamp
`fba/design.py` already puts on every record.

**Explicitly excluded from the terminal condition:** any comparison against the deterministic catalog.
