# SeqSero2's O decision procedure, read from source — and what it corrects

> **PARTLY SUPERSEDED THE SAME DAY — read this header before trusting the threshold section below.**
> This memo is kept VERBATIM as the pre-implementation record (the acceptance bar was frozen against
> it), so nothing below is edited. But its central threshold-translation argument was then MEASURED,
> and it **fails in one place**: see `wiki/salmserovar_o_fix_result_2026-09-09.md`.
>
> The coverage↔k-mer-score mapping holds for the branch **gates** (`wbaV > 70`, `wzy > 40` — different
> antigens, where presence-strength is the question) and **breaks for the `tyr` O9-vs-O2 refinement**,
> where the two references are ~99% identical to each other so both align full length under BLAST.
> Implementing that comparison on coverage fired `O-2` on 16 isolates for **13 miss / 3 no-call / zero
> hits**. The refinement is now REFUSED rather than approximated. The "Honest limits" section below
> already said the mapping was argued and not measured and should be validated first — it has been, and
> that is the result.

The review's recommended first step was to **read** SeqSero2's O logic rather than infer rules from
header names, since inferring-from-name had already produced three wrong causal claims here. This is
that reading, from `SeqSero2_package.py` 1.3.2 (`call_O_and_H_type`, lines 1070–1131; scoring in
`get_kmer_dict`, lines 1053–1069), extracted from the pinned biocontainer.

**It corrects two claims — one from the review, one of mine.**

## The scores are coverage-of-reference, not alignment scores

```python
score = (len(lib_dict[h] & input_Ks) / len(lib_dict[h])) * 100
```

The score is **the percentage of the reference entry's own k-mers found in the input**. That is a
*coverage-of-reference* quantity. Our BLAST caller computes `cov = (alignment_length / qlen) * 100`
with the allele as query — **the same semantic axis**.

Admission to `O_dict` requires `score > 25`.

| SeqSero2 condition | our equivalent |
|---|---|
| in `O_dict` | `percent_coverage > 25` |
| `wbaV > 70` | `percent_coverage > 70` |
| `wzy > 40` | `percent_coverage > 40` |

**So the branch thresholds ARE translatable** — the open question "what if k-mer scores can't be
translated" largely dissolves. Caveat: k-mer *set overlap* is not contiguous alignment, so a scattered
match could score differently. Analogous, not identical — but a principled mapping, not a blind
approximation. Note our current coverage cut of 40 is **stricter** than SeqSero2's admission bar of 25.

## The algorithm: three branches, in order

**Branch A — the O9 family, gated on `wbaV > 70`:**
```
if wbaV__1002 > 70 or wbaV-from-II-...__1002 > 70:
    if wzy__1191 > 40 or wzy_partial__216 > 40:  -> O-9,46
    elif O-9,46,27_partial_wzy__1019 present:    -> O-9,46,27
    else:                                        -> O-9        # DEFAULT within the branch
         if tyr-O-2 score > tyr-O-9 score:       -> O-2
```

**Branch B — the 3,10 / 1,3,19 family:**
```
elif O-3,10_wzx__1539 present and O-9,46_wzy__1191 present:
    if O-3,10_not_in_1,3,19__1519 present:       -> O-3,10
    else:                                        -> O-1,3,19
```

**Branch C — generic argmax, with two guards:**
```
else:
    winner = argmax(O_dict by score), antigen = key.split("_")[0]
    guard 1: if winner is a wbaV entry and NEITHER wzy nor wzy_partial present -> O-9
    guard 2: if antigen == "O-1,3,19": redo the argmax EXCLUDING 'O-1,3,19_not_in_3,10__130'
```

## Correction 1 — to the review: the plan is NOT unimplementable

The review's critical issue said the O9 branch cannot be implemented because our build drops
`tyr-O-9__609`. Reading the source shows **`tyr` only splits O9 from O2**, both *inside* Branch A. The
determinant that reaches `O-9` is **`wbaV`**, and plain `O-9` is the branch's **default**. So O9 is
reachable without any O9-specific allele and without `tyr`.

Checked by length against our DB — **every entry the algorithm names is already present**:

| SeqSero2 key | ours |
|---|---|
| `O-9,46_wzy__1191` | `O__9,46__66` |
| `O-9,46_wbaV__1002` (×2 variants) | `O__9,46__128`, `O__9,46__239` |
| `O-9,46_wzy_partial__216` | `O__9,46__362` |
| `O-9,46,27_partial_wzy__1019` | `O__9,46,27__63` |
| `O-3,10_wzx__1539` | `O__3,10__92` |
| `O-3,10_not_in_1,3,19__1519` | `O__3,10__269` |
| `O-1,3,19_not_in_3,10__130` | `O__1,3,19__126` |

**Branches A and B are implementable with the sequences we already ship.** The only genuine addition is
`tyr-O-9`/`tyr-O-2`, needed *solely* for the O9-vs-O2 refinement.

## Correction 2 — to my own Round-2 claim

I asserted "the determinant needed to call O9 is not in our DB at all" and treated that as a scope
increase. That was **overstated**: `wbaV` is the determinant, we have both variants, and `tyr` is only
the O2 refinement. I made that claim from the *docstring's* exclusion list rather than from the
algorithm — the same inferring-without-reading pattern the reading step exists to stop, committed while
arguing that the reading step is load-bearing.

## What is actually blocking implementation

Not sequences. **Role labels.** Our schema collapses all four `9,46` entries to the antigen string
`9,46`, so the caller cannot tell `wbaV` from `wzy` from `wzy_partial` — and every branch condition is
phrased in terms of *which entry* hit. The build (`build_salmserovar_db.py:71`) writes only
`{axis}__{antigen}__{index}`.

So the minimum change is: **carry the role**, then implement the three branches against coverage
thresholds. No biological sequence is authored.

## Both `not_in` markers, precisely

They are used in **opposite** directions, which is why a symmetric rule would have been wrong:

- **`O-3,10_not_in_1,3,19`** — used **positively**: its presence *selects* O-3,10 (Branch B).
- **`O-1,3,19_not_in_3,10`** — **never** produces a positive call. It is not consulted in Branch B at
  all, and in Branch C guard 2 explicitly re-runs the argmax with that exact key excluded.

`O-1,3,19` is therefore reached only by the **absence** of the reciprocal marker inside Branch B. Our
caller lets the 130 bp marker win outright, which is precisely the 13 over-calls.

## Honest limits

- Read for **k-mer mode (`-m k`)**, which is what we run (`-m k -t 4`). The allele/microassembly modes
  have separate paths (`decide_O_type_and_get_special_genes`, `predict_O_and_H_types`) and are **not**
  covered here.
- The coverage↔k-mer-score mapping is an **argued equivalence on the same axis**, not a measured one.
  It should be validated empirically before the numbers are trusted.
- Branch C's argmax compares scores **across different reference lengths**; whether our coverage
  behaves the same way under that comparison is untested.
- This reading does not, by itself, improve any number. Nothing is implemented.

## Reproduce

```bash
docker run --rm quay.io/biocontainers/seqsero2:1.3.2--pyhdfd78af_0 \
  cat /usr/local/bin/SeqSero2_package.py > ss2_source.py
sed -n '1053,1131p' ss2_source.py
```
