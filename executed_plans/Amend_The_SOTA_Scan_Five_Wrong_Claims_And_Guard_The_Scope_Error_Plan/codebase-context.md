**The two artifacts under amendment.** `wiki/ai_genomics_sota_scan_2026-10-22.md` (prose, ~230 lines)
and its `.json` sidecar. Verified line anchors in the `.md`: the F1 heading at **line 115** and the
"independently reproduce them" sentence at **line 68** carry C-1; the D2 block begins at **line 40**;
the component miscount is at **line 33**; the "unasked" assertion at **line 111**; the "retires"
clause at **line 178**. Lines **58** and **92** also contain "reproduces" but refer to the
fold-range check (16,227 vs 16,226) — that usage is **correct and must not be edited**.

**JSON fields under amendment** (verified by parse): `decisions_needed[].id` = `["D1","D2"]`;
`findings[0].headline` ends "…reproduce this project's encoder-gate negative";
`findings[0].why_they_agree`; `findings[1].consequence` carries the "retires" claim;
`decisions_needed[0].components_computable` / `components_structurally_absent` carry the 4/3 split
with the two composites absent; `measured_here.open_gap_the_field_left` carries C-3.

**Pinned-model ground truth** (`grep -roE "GenerTeam/[A-Za-z0-9._-]+"`): five sites, all
`GenerTeam/GENERator-v2-prokaryote-1.2b-base` — `dna_decode/glm/generate.py` (×2),
`scripts/kaggle/build_glm_sampling_sweep_kernel.py`, `scripts/kaggle/glm_sampling_sweep_kernel.py`,
`tests/test_glm_generate.py`. Local cache: one dir,
`D:/hf_cache/hub/models--GenerTeam--GENERator-v2-prokaryote-1.2b-base`, snapshot
`8b2f768b0d293953518ff91d34600f9322ef1f94`.

**Component ground truth** (`scripts/promcalc_decompose.py:67-68`): `DG_FIELDS = ("dG_total",
"dG_10", "dG_35", "dG_ext10", "dG_spacer", "dG_disc", "dG_ITR", "dG_UP", "dG_bind")` — nine.

**Wiki replay — two HIGH-salience precedents that bear directly.**
- `decisions-log.md:193` — an independent-cohort write-up that *"overclaimed 'single-cohort artifact /
  floor not crossable'"* and was **corrected after a /brainstorm audit**. Same shape as this plan:
  the review caught an over-claim in a shipped artifact, and the fix was to narrow the claim to what
  the denominator supports rather than to re-run anything.
- `decisions-log.md:132` — *"When a /brainstorm surfaces 'no critical issues' + flags continued
  review-cycling as the next risk, the discipline-honoring move is to ship the artifact, not to
  invoke the next planning skill."* **Checked and does not block this plan:** the review returned
  three Critical findings, not "no critical issues". But it is a live caution against expansion, and
  it is why Step 3 is the only item beyond the review's literal recommendation (see Risk Flags).

**Test infrastructure.** `pytest` via `uv run pytest tests/ -q`. Prose-guard precedent exists and is
the reuse target (below). JSON-parseability precedent: `tests/test_fba_solver_audit.py:213`
(`test_every_fba_wiki_artifact_is_parseable_json`), plus hand-edited-JSON invalidity is a recorded
CLAUDE.md gotcha — a wiki JSON was invalid for a day because nothing read it.

**External surfaces used, and how each was grounded.**
- HuggingFace model API — `GET https://huggingface.co/api/models/<id>`, read-only; `sha` and
  `lastModified` fields observed directly on both candidate ids. Grounded by invocation, not docs.
- `arXiv 2604.06549` / `2506.10271v3` / `2608.05329` — read via their arXiv HTML renderings. Their
  metric definitions (LL-derived) are the load-bearing facts and were read from the papers, not
  inferred. `[grounded]`
- NAR `gkaf863` — read in two windows covering chars 0–200,000 of 209,274. **Tail unread; this is
  exactly C-3.** `[inferred]`

### Reusable-Code Survey

- **`tests/test_essentiality_overclaim.py` — the direct reuse target, adopted wholesale as a
  pattern.** It is already a prose-**relationship** guard for an over-claim of this class, with the
  four parts this plan needs: a `SURFACES` tuple of files to scan; a `PRE_FIX_STRINGS` tuple kept
  *only* to prove non-vacuity (its own comment: *"If the guard passed on these it would be testing
  nothing"*); `_segments_mentioning()` splitting on blank lines so the scanned unit is the paragraph
  a human actually reads; and `_violations()` checking a required marker's presence in the same
  segment. Its closing test also pins the *premise* against the artifact rather than asserting it
  from memory. Step 3 mirrors all four.
- **`tests/test_build_validation_report_card.py:405`
  (`test_prospective_prose_names_the_ACTIVE_lock_and_never_a_retired_one`)** — second instance of the
  same pattern, confirming it is a repo convention rather than a one-off, and the source of the
  "assert the right name appears AND no other candidate does" idiom.
- **`tests/test_fba_solver_audit.py:213`** — reuse for Step 4's JSON-parseability assertion rather
  than writing a new parse check.
- Searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/utils`, `src/lib`, `src/common`,
  `src/helpers`, `utils`, `lib`, `common`, `helpers` (none exist), `tests/` for prose-guard and
  JSON-parse precedent (two and one hits respectively).
