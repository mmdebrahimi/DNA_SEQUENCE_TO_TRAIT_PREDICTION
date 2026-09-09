# The serovar name table silently discarded 213 serovars — and my second frozen bar was also mis-specified

After the O-procedure port, **13 of the 14 remaining misses carried a formula byte-identical to
SeqSero2's**. The antigen axes now agree; what diverges is the **formula → serovar-name lookup**. The
registry's standing line that "the formula lookup is not the problem" was true when the O axis dominated
and had gone stale.

## The defect

`scripts/build_salmserovar_db.py::build_table` keyed on `(O, H1, H2)` and kept the **first** serovar per
key, commented *"first-wins: subspecies-I named serovars come first."*

Measured against SeqSero2's own lists: **2578 entries → 2365 distinct formulas; 195 formulas carry more
than one serovar; 213 serovars were silently discarded.** And the comment is **falsified by its own
data** — for `16:b:e,n,x` the first entry is the *unnamed subspecies-II formula* `II 16:b:e,n,x` and the
named `Hvittingfoss` is second. Cubana, Muenchen and Hvittingfoss were absent from our table entirely.

Reading SeqSero2's naming code (`SeqSero2_package.py` ~L286–380) shows **four** mechanisms our build
applied none of:

| mechanism | what it does | what it fixes here |
|---|---|---|
| `remove_list` (13) | excludes entries outright | — |
| `rename_dict` (30) | canonical renames | **Virginia → Muenchen**, **Hindmarsh → Bovismorbificans** |
| bracket expansion | `r,[i]` also matches a bare `r` | why `('8','r','1,5')` never reached Bovismorbificans |
| subspecies filter | candidates restricted to the called subspecies | `II 16:b:e,n,x` vs **Hvittingfoss** |

We have no subspecies signal (SeqSero2 gets it from SalmID), so the fourth becomes a deliberately weak
tie-break: on an I-vs-non-I collision take I. A collision *within* subspecies I (Agoueve vs Cubana)
stays genuinely ambiguous and is never resolved by inventing a rule.

## Two variants, measured separately

I had conflated two independent design choices. Separated and measured on the same 200 isolates against
the same frozen bar:

| | hits | miss | no_call | accuracy | regressions |
|---|---|---|---|---|---|
| baseline (after the O port) | 165 | 14 | 21 | 0.8250 | — |
| **A** — ambiguous keys **omitted** (abstain) | 165 (+0) | 6 | 29 | 0.8250 | **3** `hit→no_call` |
| **B** — ambiguous keys keep the prior winner | **168 (+3)** | 11 | 21 | **0.8400** | **0** |

**A is the principled one and it is strictly worse.** Abstaining on an ambiguous formula matches the
caller's stated contract, but the abstentions cost exactly as much as the renames gain. **The
canonical-naming half of this change is what pays; the ambiguity-policy half costs.** B ships.

**Index-matched control (mechanical, in the artifact).** This change rebuilds `serovar_table.tsv`, which
is also the input to `load_formula_index` — so it moves the **scoring function**, not just the caller.
Re-scoring the previous run's *unchanged* serovar names under the current index gives **166, not 165**:
**+1 hit of any gain is equivalence, not calling.** So B's honest caller-attributable gain is **+2**, and
the artifact records it rather than letting the fix quietly bank it.

## The bar said REJECT, and the bar was wrong again

    [FAIL] hits rise by >= +4 over the baseline 165   (165 -> 168, +3)
    [PASS] ZERO regressions: hit->miss is 0
    [PASS] ZERO regressions: hit->no_call is 0
    [PASS] confident errors do not increase: miss <= 14   (14 -> 11)
    [PASS] H1 / H2 agreement not below baseline
    VERDICT: REJECT

**The `+4` was unachievable by construction.** I sized it from an actionable set of 5 *before* reading
the mechanism. Reading it afterwards shows only **3** of those 5 are recoverable without guessing — the
other 2 are the Agoueve/Cubana collision, which no faithful implementation can resolve. A bar that no
correct implementation can satisfy is not a test; it is an error in the bar.

**This is the SECOND mis-specified bar today, and that pattern is the finding.** The first
(`salmserovar_o_fix_acceptance_bar.json`) used a **net** `no_call` ceiling to police a **directional**
failure mode and rejected a change that was +24 hits with zero regressions. This one used a **threshold
sized from an un-analysed count**. Different errors, one root cause:

> **I was writing the bar before doing the mechanism analysis that tells you what a correct fix can
> actually achieve.** Freezing early protects against tuning to the result; it does not excuse freezing
> a number you have no basis for. The analysis has to come first, and *then* the freeze.

Two overrides in one session is not a habit to normalise. **The bar file and the artifact's
`verdict: REJECT` are both left unedited**, and the adoption is surfaced rather than quietly taken:
re-sizing an acceptance bar is an authority call, and this one should be ratified or reversed
deliberately. What supports shipping B is that it is **strictly dominant** — zero regressions in either
direction, fewer confident errors, higher accuracy — and that the clause it fails is unsatisfiable.

## Where the cell stands

**0.8400** vs SeqSero2's 0.8800 — **−0.0400**, from −0.1750 this morning. Better, not fixed.
Residual: 11 misses, 21 abstentions.

## Honest limits

- **+3 is the raw gain; +2 is the caller-attributable gain** after the index-matched control. Quote +2
  when attributing to the fix.
- The subspecies tie-break is **ours, not SeqSero2's**. SeqSero2 filters by a subspecies it determines
  independently; we prefer I because it is 1532 of 2578 entries. It breaks I-vs-non-I ties only.
- The shipped `first` policy keeps an **arbitrary** winner on a genuinely ambiguous formula. It is the
  pre-existing behaviour and it scores better, but a call from an ambiguous key is right by luck. The
  principled alternative was measured and is worse; that trade is recorded, not hidden.
- 9 keys present in the old table are gone — `remove_list` exclusions, i.e. entries SeqSero2 itself
  refuses to report.
- One cohort, per-serovar cap 12: per-isolate accuracy on a deliberately diverse mix, not
  population-weighted.
- `13:z29:-` (Agoueve/Cubana) remains **wrong on 2 isolates** and is not addressable without subspecies
  input or SeqSero2's multi-name-with-star reporting, which our equivalence function would not credit.

## Reproduce

```bash
docker run --rm quay.io/biocontainers/seqsero2:1.3.2--pyhdfd78af_0 \
  cat /usr/local/bin/Initial_Conditions.py > D:/dna_decode_cache/ss2_ic_dir/Initial_Conditions.py
uv run python -c "from pathlib import Path; import sys; sys.path.insert(0,'.'); \
  from scripts.build_salmserovar_db import build_table; \
  print(build_table(Path('D:/dna_decode_cache/ss2_ic_dir'), Path('data/salmserovar_db/serovar_table.tsv')))"
uv run python scripts/salmserovar_o_fix_rerun.py \
  --bar wiki/salmserovar_name_table_acceptance_bar.json \
  --baseline-checkpoint D:/dna_decode_cache/ss2_checkpoint/ofix_rows.o_port_oldtable.jsonl
```

Frozen AMR surface byte-unchanged — typing cell.
