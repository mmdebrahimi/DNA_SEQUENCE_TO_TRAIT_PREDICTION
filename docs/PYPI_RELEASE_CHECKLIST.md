# PyPI Release Checklist

## Goal

Publish a short, accurate PyPI surface without regressing packaging truth.

## 1. Metadata sync

Before build:

- bump package version everywhere
- remove stale version strings from public-facing docs
- point `[project].readme` at `README.pypi.md`
- confirm `description` is a one-line package summary, not a project-history block
- confirm `project.urls` includes:
  - Homepage
  - Repository
  - Documentation
  - Issues
  - Changelog
  - Validation
  - FeatureMatrix

## 2. README discipline

`README.pypi.md` should:

- explain the package in under 30 seconds
- put install first
- show at least one command that really works
- keep the research-use / not-clinical warning near the top
- avoid long internal history
- use absolute GitHub links, not repo-relative links

## 3. Feature / dependency clarity

Before release, confirm the PyPI docs answer:

- which commands work after base install
- which need BLAST+
- which need DB downloads
- which need Docker
- which need precomputed AMRFinder or other tool outputs

## 4. Packaging truth

Do not ship a PyPI page that promises commands the wheel does not expose.

Check:

- console scripts in wheel match advertised commands
- package data needed for trust/report cards is included
- optional dependency groups match the real command surface
- `[project.urls]` points at real public docs/issues paths

## 5. Build and validation

Run:

```bash
python -m build
twine check dist/*
```

Then do an installed-wheel smoke check in a clean env:

```bash
pip install dist/dna_decode-<version>-py3-none-any.whl
dna-decode list
```

And at least one lightweight command that should work without hidden source-tree assumptions.

## 6. Publish boundary

The release branch should be based on fresh `origin/main`.

Do not publish from a stale long-lived fork when:

- PyPI has already advanced
- `origin/main` has moved
- release-facing docs were drafted elsewhere and not yet reconciled
