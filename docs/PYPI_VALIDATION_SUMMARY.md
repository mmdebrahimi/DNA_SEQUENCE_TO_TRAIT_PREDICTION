# PyPI Validation Summary

## Principle

`dna-decode` should not present all surfaces as equally validated.

Every output should make it clear:

- what the tool called
- why it called it
- what the trust tier is
- what the main blind spots are

## Trust-tier shape

Recommended public badges:

- `INDEPENDENT_MEASURED`
- `INDEPENDENT_WETLAB`
- `IN_DISTRIBUTION`
- `FAITHFUL_TO_TOOL`
- `ABSTAINS_BY_DESIGN`

## Compact validation summary

| Surface family | Preferred trust tier | What that means |
|---|---|---|
| bacterial AMR | `INDEPENDENT_MEASURED` | evaluated against measured AST on provenance-disjoint cohorts |
| M. tuberculosis | `INDEPENDENT_MEASURED` | deterministic rule tested against measured AST |
| HIV-1 | `INDEPENDENT_WETLAB` | evaluated against isolate-level wet-lab fold-change data |
| SARS-CoV-2 / influenza | often `IN_DISTRIBUTION` | useful, but some cells are underpowered or non-independent |
| fungal AMR | scope-limited | do not overstate as a universal cross-lineage binary surface |
| typing / pathotype | often `FAITHFUL_TO_TOOL` | faithful to a reference DB/tool, not independent phenotype validation |
| human PGx | mixed | genotype calling may be independently checkable while phenotype mapping stays faithful-to-guideline |
| abstaining surfaces | `ABSTAINS_BY_DESIGN` | refusal is the honest output when the surface is out of scope |

## Good public wording

Prefer:

- "reports a call plus its trust badge and provenance"
- "validation tier varies by decoder surface"
- "typing and some PGx surfaces may be faithful-to-tool rather than independently phenotype-validated"
- "some surfaces abstain by design rather than overclaim"

Avoid:

- "validated everywhere"
- "works across all organisms"
- "better than existing tools" without a cell-specific benchmark

## Deeper links

The short PyPI page should link to deeper validation/report-card material rather than embedding the full research history.

Useful link buckets:

- cross-kingdom validation summary
- bacterial AMR report card
- TB report card
- HIV report card
- fungal scope-limit memo

## Release hygiene

Before publishing a PyPI README update, check:

- no stale version strings
- every mentioned trust tier has a real linked artifact behind it
- no in-distribution or faithful-to-tool surface is accidentally described as independent measured validation
