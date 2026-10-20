# This host: BLAST+ is absent, and two torch processes crash the suite

**Date:** 2026-10-20 · Both findings are measurements on `C:/Users/Laptop Portable`, the host CLAUDE.md
describes. Held here rather than in CLAUDE.md because that file sits against a word ceiling
(`tests/test_claude_md_growth_ceiling.py`), and a pointer is its prescribed compression mechanism.

## 1. BLAST+ 2.17.0 is NOT installed here — CLAUDE.md asserted it for months

CLAUDE.md's `vf_diff` bullet (2026-06-04) said *"NCBI BLAST+ 2.17.0 (blastn + makeblastdb) is installed
natively at `C:/Users/Farshad/ncbi-blast/bin` on this host (Q1 ratified)"*. Measured through the DEPLOYED
resolver, not by inference — `dna_decode/pathotype/vf_runner.py::find_blastn()`:

| resolver branch | result |
|---|---|
| `$BLASTN_BIN` | unset |
| `shutil.which("blastn")` | not on PATH |
| `~/ncbi-blast/bin/blastn.exe` (= `C:/Users/Laptop Portable/...`) | **absent** |
| `C:/Users/Farshad/ncbi-blast/bin/blastn.exe` | **absent** (profile exists, 167 entries, no `ncbi-blast`) |

`find_blastn()` returns **None**. The claim was true of the host it was written on; the binary is
gitignored-class so it never synced, and the `Farshad` vs `Laptop Portable` profile mismatch is why it
reads as plausible here. Independently found by a peer session and re-verified here before acting.

**Consequence:** every blastn-dependent path takes its offline-safe branch — `vf_diff` / the virulence
overlay, the fungal ERG11 caller, HIV RT/PR/IN/CA genome-FASTA mode, SARS-CoV-2 Mpro genome mode, and the
typing suite (serotype / ktype / salmserovar / pneumoserotype / plasmid / resfinder / pointfinder /
disinfinder).

**NOT claimed: that any committed typing validation is invalid.** Those ran where blastn resolved. What is
measured is only that they cannot be RE-RUN here until BLAST+ is installed — `winget install --id
NCBI.BLAST`, or the NCBI win64 installer into `C:/Users/Laptop Portable/ncbi-blast/bin`, which satisfies
the resolver's 3rd branch with no env var. Installing it is a USER decision, not taken.

**I propagated this claim before checking it** (the 2026-10-19 clonality memo cited the install as the
reason a native AMRFinder route was worth testing) — the exact error CLAUDE.md's own header warns about.

## 2. Do NOT run the full suite while a torch model job is loaded

`uv run pytest tests/ -q` died mid-run with **`Windows fatal exception: access violation`** at
`tests/test_genome_jepa_prototype.py` while `scripts/glm_generator_falsifier.py` held
GENERator-1.2B in a second process. pytest died outright — **two Windows Error Reporting records at
22:36:47/48**, no summary line, and the harness recorded no completion, so it is indistinguishable from
a test failure at a glance.

Two hypotheses were measured WRONG first:

- **Not memory.** Pagefile PEAK 466 MB against a 2432 MB allocation, 8.6 GB RAM free.
- **Not the test file.** Isolated, nothing else running: **6 passed, 3 runs for 3, zero faults.**

**What is established vs inferred.** Established by measurement: not memory, and not the test file.
The remaining attribution — two torch processes at once — rests on **ONE crash event** plus three clean
isolated runs. That is consistent with concurrency and is NOT proof of it; `n=1` on the failing side
cannot separate "concurrency" from "intermittent native fault." Treat *do not run them together* as a
sound precaution, not a demonstrated mechanism. Installed torch is **2.7.1+cu118**, a CUDA build pinned
deliberately (`[tool.uv.sources] torch = { index = "pytorch-cu118" }`) for a **GTX 860M (Maxwell, CC 5.0)**
that modern torch no longer supports — a plausible common root for native instability here, also unproven.

**A SECOND, DIFFERENT native crash, found 2026-10-20 while testing the generator alone.** Loading the GLM
segfaults (exit 139, zero stdout/stderr) with **nothing else running**, so it is NOT the crash above.
`faulthandler` localises it exactly: `glm/generate.py:233 -> load() -> AutoTokenizer.from_pretrained ->
hf_hub_download -> httpx/httpcore` — a **network fetch during tokenizer load**. The model file is intact
(`model.safetensors` 4,648,274,384 B = 4.33 GiB) and E: has 911 GB free, so it is not corruption or space.
This is the same Hub re-fetch path `setup_tokenizer` exists to work around: the remote modeling code
re-fetches its own tokenizer with no `cache_dir` and no `local_files_only`. **`HF_HUB_OFFLINE=1` bypasses
it** (the load then proceeds past the tokenizer instead of faulting instantly). Use that for any local GLM
run on this host.

**Why it only started biting:** the file dates from 2026-07-11 and is guarded by
`pytest.importorskip("torch")`, so for three months it SKIPPED. Installing torch during a session converts
the skip into real execution — which is also the 6625 → 6647 pass jump.

**Rule: serialise torch work**, or accept that a concurrent suite run is uninterpretable.

## 3. Clean baseline, measured with nothing concurrent

**16 failed / 6655 passed / 185 skipped** in 394 s. Zero failures in changed files. The 16 are absent-data
and absent-extra cases, consistent with the standing accounting in CLAUDE.md's reproducibility-freeze
bullet.
