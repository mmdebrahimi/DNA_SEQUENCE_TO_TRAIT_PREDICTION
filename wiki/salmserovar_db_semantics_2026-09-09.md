# The O-axis root cause is a lost header qualifier — and my "external wall" was wrong

One run ago I diagnosed the serovar O-axis gap as **two reference-DB content defects** and concluded the
fix needed *sourced* reference sequences, i.e. an **external wall**. That conclusion is **wrong**, and
this corrects it. The sequences were never missing. The **semantics** were.

## What the DB actually is

Our `salmonella_antigens.fasta` is **SeqSero2's antigen database, re-headered**. Measured, not assumed:

> **All 62 of our O alleles are byte-identical to a SeqSero2 entry.**

SeqSero2 names its entries with a **semantic qualifier**; our conversion maps everything to
`O__<antigen>__<index>`, which has no field for one. Six SeqSero2 O entries carry such a qualifier:

| SeqSero2 header | what the qualifier means |
|---|---|
| `O-1,3,19_not_in_3,10__130` | a **differential marker** — present in 1,3,19, absent in 3,10 |
| `O-3,10_not_in_1,3,19__1519` | the reciprocal differential marker |
| `O-9,46_wbaV__1002` | **`wbaV`** — the gene that distinguishes O9 from O9,46 |
| `O-9,46_wbaV-from-II-9,12:z29:1,5-SRR1346254__1002` | the same gene, second source |
| `O-9,46,27_partial_wzy__1019` | explicitly **partial** |
| `O-9,46_wzy_partial__216` | explicitly **partial** |

The two alleles I flagged as anomalous are exactly two of these, confirmed by byte-identical sequence:

- our `O__1,3,19__126` (130 bp) **is** `O-1,3,19_not_in_3,10__130`
- our `O__9,46__362` (216 bp) **is** `O-9,46_wzy_partial__216`

They are not defective sequences. They are **correctly-labelled special-purpose entries whose labels our
build discarded.**

## One root cause, both symptoms

Our caller treats every entry as a standalone positive antigen allele, so:

1. **`not_in_3,10` becomes a positive `1,3,19` call whenever it hits.** That is exactly the measured 60%
   hit-rate at 28× its antigen's true prevalence, and the 13 `1,3,19` over-calls. The marker is doing its
   job; the caller is reading it wrong.
2. **The `wbaV` entries — the O9-vs-O9,46 discriminator — flatten into two more ordinary `9,46` alleles.**
   With no way to express "wbaV present ⇒ O9", plain O9 becomes **unemittable**, which is the 15 `9,46→9`
   over-calls.

So the "two DB content defects" are **one defect**: the conversion dropped the qualifiers.

## The correction to my own claim

I said this needed sourced sequences and was therefore external. **It is not external — it is
code/data-build closable.** Every sequence required is already present and byte-correct. What is missing
is (a) preserving the qualifier in our header schema and (b) honoring it in the caller.

**This is the third time this session I called a wall too early** — after `ktype` ("cohort-blocked", it
wasn't) and my own coverage threshold ("the likely cause", it wasn't). The pattern is consistent enough
to name: I reach for *external* before exhausting *readable*.

## What the fix is, and why it is not attempted here

Honoring these qualifiers means implementing **differential O typing**, not a patch:

- a `not_in_X` marker may only be used to *discriminate between* two candidate antigens, never to assert
  one on its own;
- a `partial` entry should not win a selection against a full-length allele on coverage terms;
- `wbaV` presence/absence must map to the O9 / O9,46 distinction.

That changes a shipped caller's semantics and must be measured against a **pre-registered** bar, the way
the coverage-40 change was. Attempting it in the tail of a long run — with no bar registered and no
re-measurement budget — is exactly the under-measured change this project keeps catching. **Sized and
specified here; deliberately not implemented.**

**Expected ceiling if implemented:** the two symptoms account for ~28 of the 33 O disagreements, and the
O axis is 100% of the caller's −0.175 gap against SeqSero2. That is an estimate from disagreement counts,
**not a measured post-fix result**.

## Honest limits

- Byte-identity is measured over the **O axis** (62/62). I did not check H1/H2, so "our DB is SeqSero2's"
  is established for O and merely likely elsewhere.
- The qualifier *meanings* are read from SeqSero2's header names (`not_in_3,10`, `partial`, `wbaV`). That
  is a strong reading of self-documenting names, and it matches the measured behaviour, but I did not
  read SeqSero2's source to confirm exactly how it consumes each one. **Before implementing, read that
  code** — inferring the rule from a name is how the previous three wrong causes happened.
- Whether our DB was built from SeqSero2 deliberately (with the qualifiers judged unnecessary) or the
  loss was accidental is **not established** by this evidence.
- No fix is applied, so nothing here is a measured improvement.

## Reproduce

```bash
docker run --rm quay.io/biocontainers/seqsero2:1.3.2--pyhdfd78af_0 \
  sh -c "grep '^>O-' /usr/local/lib/python*/site-packages/seqsero2_db/H_and_O_and_specific_genes.fasta"
uv run python scripts/salmserovar_db_audit.py
```

Frozen AMR surface byte-unchanged — typing cell.
