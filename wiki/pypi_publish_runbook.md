# PyPI / TestPyPI publish runbook (2026-06-24)

The package is **publish-ready + metadata-valid** (`twine check` PASSED; the wheel ships the deterministic
decoder + the 4 trust cards + all console entry points incl. the new `dna-ktype`). What remains is a
**USER-authority + credential** step — the agent cannot upload (no token, and publishing is outward/
semi-permanent). This runbook is the exact path; run it yourself with your token.

## Before you upload — 3 user decisions
1. **Name availability.** Check `dna-decode` / `dna_decode` is free on the index. On TestPyPI:
   `https://test.pypi.org/project/dna-decode/` (404 = available). If taken, rename in `pyproject.toml`
   (`[project] name = "..."`).
2. **Version is one-shot + permanent.** Once `0.5.0` is uploaded to (Test)PyPI you can NEVER re-upload that
   exact version (you can only bump to `0.5.1`). So upload a version you're happy to freeze.
3. **It's public + under your name.** The scientific claims (the trust badges) go out attributed to you.
   "Not a clinical tool" matters. TestPyPI is the sandbox; do it first.

## Step 1 — TestPyPI dry-run (sandbox, recommended first)
```bash
# (a) get a TestPyPI token: register at https://test.pypi.org/ -> Account -> API tokens -> create
# (b) the artifacts are already built + checked (dist/dna_decode-0.5.0*); rebuild if you changed anything:
uv build
uv run --with twine python -m twine check dist/dna_decode-0.5.0*

# (c) upload ONLY the 0.5.0 artifacts to TestPyPI (the sandbox):
uv run --with twine python -m twine upload --repository testpypi dist/dna_decode-0.5.0*
#   username: __token__
#   password: <your TestPyPI API token, pypi-...>

# (d) verify a clean install FROM TestPyPI works (deps come from real PyPI):
uv venv /tmp/tpv && uv pip install --python /tmp/tpv/Scripts/python.exe \
  --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ dna-decode
/tmp/tpv/Scripts/dna-amr.exe --drug efavirenz --observed RT:K103N --sample-id t   # expect the badge
```

## Step 2 — real PyPI (only after TestPyPI is clean; a SEPARATE explicit decision)
```bash
# get a real PyPI token at https://pypi.org/ ; then:
uv run --with twine python -m twine upload dist/dna_decode-0.5.0*
#   username: __token__ ; password: <your PyPI API token>
```

## What ships (verified)
- 0.3 MB wheel: the deterministic decoder (9 console commands incl. `dna-ktype`) + the 4 validation report
  cards (`dna_decode/report_cards/`). NO genomes / caches / private data (verified).
- Genome-mode AMR + the typing decoders need Docker / BLAST+ / external DBs pip can't install; the
  observed-substitution + cached-run paths work out of the box (see `QUICKSTART.md`).

## Honest status
Agent did: build wheel+sdist, `twine check` (PASSED), verify ktype + cards ship, write this runbook.
Agent did NOT (and will not without your say-so + token): upload to any index. Publishing is yours.
