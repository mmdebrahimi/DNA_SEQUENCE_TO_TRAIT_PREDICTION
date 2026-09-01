# Five shipped decoders had no evidence contract — and the guard that exists to catch that was blind

**HCMV's 5 CLI-routable drugs were absent from `cell_registry` from 2026-07-23 until today.** The
coverage test that exists precisely to stop a decoder shipping without a contract was correct the whole
time; its **input set** was wrong.

Registry **110 → 115 cells**; viral track **29 → 34**. The AMR report card is byte-unchanged.

---

## The gap

`dna-amr --drug ganciclovir` (plus `valganciclovir`, `cidofovir`, `foscarnet`, `letermovir`) routed and
produced calls, while `cell_registry.cells()` — the live evidence surface — contained nothing for HCMV.
It has its own report card (`wiki/hcmv_decoder_report_card.json`, 5 cells) and a v0 memo, so it was
*validated*; it just wasn't *contracted*, and the registry is what the trust surface reads.

**Why it slipped, structurally.** There are three ways a cell reaches the registry, and HCMV took none:

| kingdom | route into the registry |
|---|---|
| bacterial / fungal / antimalarial / **influenza NA** | PROJECTED verbatim from the frozen `shipped_decoder_surface` |
| **HIV-1 / SARS-CoV-2** | hand-declared in `_viral_contracts()` |
| **HCMV** | **neither** |

## The real defect: a drift-proof guard with a hand-maintained input

`cell_registry.cli_routable_manifest()` built its AMR drug set from a **hand-enumerated union of six
catalog imports**. `all_supported_hcmv_drugs()` was not among them — while the CLI's own argparse
`choices` at `amr/cli.py:437` built the same union and *did* include it.

So HCMV drugs never entered the "CLI-routable" set, and
`test_every_cli_amr_drug_has_a_contract` compared a set that already excluded them.

The comment two lines below that union congratulates itself for deriving `traits` from `per_target`,
"**never hand-listed**" — while the drug union directly above it was hand-listed. **Fifth instance of
this bug class in this repo.**

**Fix:** `dna_decode/data/routable_drugs.py::all_routable_amr_drugs()` is now the single definition, and
both the CLI parser and the registry import it. A seventh catalog cannot be added to one and missed by
the other.

**Confirmation the diagnosis was right:** fixing only the union — before adding any contract — made the
pre-existing guard fail and name all five drugs. The test never needed changing.

## The contracts, derived not asserted

Data-driven from the packaged HCMV card (the same pattern HIV uses), so the tier is read from evidence:

- **`KNOWLEDGE_BASELINE`** — the card self-reports `IN_DISTRIBUTION`, which is exactly this tier's
  definition ("literature/catalogue assignment, in-distribution").
- **`NO_FREE_PHENOTYPE` abstention**, not `SCORED` and not `UNDERPOWERED`. The card carries **no
  acc/sens/spec/n_scored at all** — it is a catalog *census* (`n_resistance` / `n_benign`), so `SCORED`
  ("a real validated number exists") would be an overclaim. And `UNDERPOWERED` is wrong in kind: HCMV
  phenotyping is per-**mutation** recombinant marker transfer, so there is no isolate-level source to be
  underpowered *on*. The card's own `independence` field records this as **CLOSED for free data**, not
  pending work.
- The `validation_slice` says plainly that no acc/sens/spec exists, so a tier can never imply a number.

I first wrote `native_abstention="IN_DISTRIBUTION"` — an invented term — and the vocabulary guard
rejected it. Second guard to fire correctly in this run.

## Scope

Touches **no frozen file**. HCMV is hand-declared on the `viral` track exactly like HIV/SARS-CoV-2, so
`shipped_decoder_surface.py` is untouched and the **2026-08-31 gentamicin v2 lock stays valid**
(`verify_lock` re-checked). The AMR report card is unchanged: 27 cells, zero field diffs.

## Honest limits

This makes five decoders **visible**, not better validated — their tier is a knowledge baseline and
their independence is closed for free data. Whether an HCMV cell should ship at all with no scored
metric is a scope question this does not answer; it answers only "is it declared?".

Two of the three registry entry points remain hand-maintained (`_viral_contracts` is a literal block),
so a *sixth* viral kingdom could still be added to the CLI and missed here — the new coverage guard
would now catch it, which is the point, but the declaration itself is still manual.

## Artifacts

`dna_decode/data/routable_drugs.py` · `_hcmv_card_drugs()` + the HCMV loop in `cell_registry.py` ·
4 tests in `tests/test_routable_drugs_single_source.py` (including resolving the CLI's choices through
the **real parser**, and a per-catalog non-vacuity check).
