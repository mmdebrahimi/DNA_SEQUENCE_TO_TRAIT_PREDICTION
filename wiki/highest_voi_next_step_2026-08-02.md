# Highest-VOI next step — strategic analysis (2026-08-02)

**Ask:** "spend an entire run thinking about the best / highest-VOI next step." This is that analysis:
state read, parallel-session ownership checked (R4), forks ranked adversarially.

## Headline recommendation

**Ship the tool. Cut and publish the prepared v0.9.0 release (or a fresh cut folding the post-0.7.23
cells).** Confidence: HIGH. It is the most direct realization of the project's own north star — *"AI DNA
decoder TOOL, not papers"* — and the largest block of unrealized value in the repo.

## The evidence that reframes the whole question

The decisive fact isn't about any single cell — it's the **built-vs-delivered gap**:

| fact | value |
|---|---|
| Current built version | **v0.9.0** (version bumped in `d5d2c18 "release: v0.9.0"`) |
| Last **public tag** | **v0.7.0** (2026-07-11) |
| Commits built-but-undelivered | **458** |
| Unreleased minor versions | 0.8.0 / 0.8.1 / 0.9.0 (none ever tagged/published) |
| Test suite | **3203 tests** collected clean |
| CHANGELOG `[0.9.0]` entry | **already written** (2026-07-23) |
| Publish mechanics | `~/.pypirc` present; `uv publish` free (token-based, no money) |
| CLI cells shipped in the build | **20** (`dna-amr`, `dna-forward`, `dna-inverse`, `dna-pathotype`, `dna-coatcolor`/horse/cat/chicken/pigeon, `dna-morphology`, `dna-essentiality`, `dna-metabolic`, `dna-plasmid`, `dna-serotype`, `dna-phage`, `dna-kleb`, `dna-flowering`, `dna-pigment`, pgx …) |

**v0.9.0 is a fully-prepared release that stalled at the final publish gate.** Version bumped, changelog
written, 458 commits of validated+tested capability — sitting one `git tag` + one `uv publish` from users.
Everything else on the board is *adding cells to a tool users can't get past v0.7.0*. Delivering closes
that gap in a single action.

## VOI ranking of the real forks

| # | Move | VOI | Reversible / free / mine | Verdict |
|---|---|---|---|---|
| **1** | **Ship v0.9.0** (prepared cut, or roll a fresh 0.9.1/0.10.0 folding the coat-colour fleet + confound-free arm) | **HIGH** | free; publish = Care-pause (user ratifies); mine to PREP | **Recommended #1** |
| 2 | Acquire a FREE independent clinical/wet-lab label (clears the binding LABELS wall — the one lever that unblocks learned-decode expansion) | HIGH (future) | USER-authority; partly in-progress (sibling pgx/GeTRM) | User fork; draft anchor exists |
| 3 | genome-map deferred tiers (KEGG/Pfam/eggNOG homology; non-E.coli VF DBs) | MED | reversible/free/mine | Real product work; heavy infra |
| 4 | TB independent gold-set (Thorpe 2024 ~59 ENA isolates) | LOW | external wall (network+Docker), likely underpowered | Not a bounded batch |
| 5 | Prospective-lock accrual re-sweep | LOW now | free but data hasn't accrued (eligible=0) | Waiting game |

## #1 in detail — what's mine vs the gate

**My executor role (all reversible, free, non-colliding):**
1. Confirm the full 3203-test suite is green (running now) — R3 real-surface gate; don't claim
   "release-ready" without it.
2. Confirm HEAD is a clean release point; decide 0.9.0-as-prepared vs a fresh cut that promotes the
   `[Unreleased]` changelog section (folds in the 14-organism coat-colour fleet, the confound-free
   cross-kingdom arm, pgx certification).
3. Update CHANGELOG + bump version if rolling a fresh cut.

**The gated final steps (USER ratifies):**
- `git tag v0.X.0 <commit> && git push origin` — reversible-outward (a tag can be deleted); a release tag is
  still a semantic authority act → surfaced.
- `uv publish` — **genuinely-irreversible-outward** (a PyPI version cannot be unpublished) → **Care-PAUSE.**
  I do not cross this without explicit user go.

**Timing / coordination (R4):** a sibling session has in-flight pgx *validation wiki docs* (uncommitted:
`pgx_report_card`, `getrm_concordance`, `certification_capstone`). These are docs, **not package code** —
they don't ship in the wheel, so they do NOT block a release. The pgx *cell* itself is already committed and
in the build. So a clean-state release is safe today; the only real choice is **0.9.0-as-prepared vs a
fresh cut** — a value/timing judgment for the user.

## Pre-mortem (why #1 survives)

- *"Shipping is motion, not signal."* → No — 458 commits of validated capability undelivered is the
  opposite of motion; it realizes built value. The "decode novel signal" frontier is a closed negative
  (organism-level DL) or owned/blocked.
- *"It'll break for users."* → Gated on a full green suite + frozen surfaces byte-locked + prospective-lock
  intact.
- *"Shipping mid-pgx-work causes a half-cell."* → pgx code is committed; only validation *docs* are
  in-flight and don't ship in the wheel.

## Honest plateau statement

The **research** frontier that is mine + reversible + free is at a genuine plateau: AMR banked (freeze
2026-06-13), forward molecular cell at a genuine terminal (0.478 over 2383 held-out), organism-level DL a
closed negative, and this session's confound-free cross-kingdom arm just reached a clean terminal
(architecture-vs-power settled). The forward/ProSST + pgx-certification frontiers are the **sibling
session's** owned lane. When every high-VOI *research* move is banked, closed, owned-elsewhere, or
externally gated, the highest-VOI move flips from "decode more" to **"deliver what's decoded"** — the
release. That is not a retreat; it is the north star.

## What I need from you (the one decision)

Ratify the release path: **(A)** publish 0.9.0 as-prepared, **(B)** roll a fresh cut (0.9.1/0.10.0) folding
the post-0.7.23 cells first [my lean], or **(C)** hold. I'll do all reversible prep + verification and stop
at the `uv publish` Care-pause for your explicit go.
