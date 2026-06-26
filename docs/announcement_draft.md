# Announcement draft (for YOU to post — not auto-posted)

These are prepared drafts. Posting to LinkedIn / Reddit / Biostars, emailing labs, submitting to Bioconda /
JOSS, and minting a Zenodo DOI are **your** actions (your accounts + identity + a release decision) — they
are intentionally not automated. Edit freely before posting.

## Positioning (use this, not hype)

> Research-use, interpretable microbial genomics CLI for supported AMR/typing surfaces.

**Strongest hook:** `dna-decode` doesn't just output a call — it outputs **the call, the biological
evidence, the blind spots, and the validation tier.**

## Short announcement

> I released **dna-decode**, a research-use Python CLI for interpretable microbial genome→trait calls.
> It takes genome assemblies, AMRFinder outputs, or observed mutations and reports:
> - AMR R/S calls where supported (bacteria, M. tuberculosis, C. auris, HIV-1, SARS-CoV-2, influenza)
> - the exact genes/mutations driving each call
> - pathotype / typing outputs for supported organisms
> - a validation tier, caveats, and provenance — and **abstention** when the mechanism is out of scope
>
> ```bash
> pip install dna-decode
> dna-decode amr --drug efavirenz --observed RT:K103N
> dna-decode list
> ```
>
> Output includes the call, the determinant evidence, and a trust badge. **Not a clinical diagnostic tool.**
> Looking for feedback from microbial-genomics / AMR users on usability, missing organisms, validation
> surfaces, and docs.

## Channel checklist (your call, in rough priority order)

- [ ] GitHub: tag a release + add the repo **topics** (bioinformatics, genomics, amr, antimicrobial-resistance, microbial-genomics, python, cli, resfinder, pointfinder, mlst, pathotyping, genotype-phenotype).
- [ ] PyPI: the metadata is prepared in `pyproject.toml`; a re-publish is a version bump + your token (parked — say the word).
- [ ] Bioconda: `grayskull pypi dna-decode` to scaffold a recipe (needs adequate tests + a stable sdist + checksum) — high value for bioinformatics adoption.
- [ ] Zenodo: enable the GitHub↔Zenodo integration → next release mints a citable DOI (feeds `CITATION.cff`).
- [ ] Reddit r/bioinformatics / Biostars: a tool-announcement + "feedback wanted", no overselling.
- [ ] Direct outreach: 10–20 AMR/genomics labs/tool maintainers.
- [ ] (later) JOSS: once docs + tests + statement-of-need + examples + citation metadata are mature.

## What NOT to say (regulatory / credibility)

Avoid: "AI DNA decoder" · "predicts phenotype from DNA" · "breakthrough" · "better than AMRFinder"
(over-broad) · "clinical AMR prediction". These are over-broad or carry medical/regulatory risk.
