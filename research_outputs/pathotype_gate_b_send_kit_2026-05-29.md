# Pathotype Gate B — Consolidated SEND KIT — 2026-05-29

> **Single source of truth for executing Gate B outreach.** Merges this laptop's named 5-contact shortlist (2026-05-27) with the workhorse's post-Gate-A bundle (2026-05-29). Supersedes the pre-Gate-A framing in `research_outputs/gate_b_cold_email_packet_2026-05-27.md` for send purposes.
>
> **Why this exists:** two machines produced overlapping Gate B materials. This file resolves them into one send-ready packet so you don't juggle near-duplicates.

## Status going in

- **Gate A: PASS** (workhorse, 2026-05-28) — real VirulenceFinder runtime (`python -m virulencefinder` + cloned VF DB + local blastn), 5-strain frozen panel all called correctly:
  - `ehec_edl933 → EHEC_COMPATIBLE`, `etec_h10407 → ETEC_COMPATIBLE`, `eaec_042 → EAEC_COMPATIBLE`, `upec_cft073 → UPEC_COMPATIBLE`, `commensal_mg1655 → COMMENSAL_LOW_MARKER_BURDEN`
- **Gate B: NOT STARTED** — this kit is the means to run it.
- **PASS rule: ≥2 credible yes-to-pilot within 60 days.** (`reports/pathotype_gate_b/pathotype_gate_b_response_rubric_2026-05-29.md`)

## Binding action is YOURS

Claude cannot send cold emails. This kit prepares the send + scores replies. You: verify 5 names (~5-10 min total), send, log replies.

---

## The 5 contacts (named; verify before send)

Send order (highest-leverage first): **4 → 3 → 5 → 2 → 1**

| # | Profile | Named contact (verify) | Stable URL | Verification query | Why |
|---|---|---|---|---|---|
| 4 | Bioinformatics maintainer | **Torsten Seemann** (ABRicate/mlst/shovill, Melbourne) | https://github.com/tseemann | `"Torsten Seemann" abricate contact` | Maintainer of the gene-presence frontend this tool would adopt-from; strongest possible "fills a gap" signal |
| 3 | Clinical-micro PI | **Tracy H. Hazen** (DECA author, IGS, U. Maryland) | https://www.medschool.umaryland.edu/profiles/hazen-tracy/ | `"Tracy Hazen" site:umaryland.edu` | Authored DECA; dominant independent-label source in Horesh 2021; knows the substrate-circularity problem |
| 5 | Strategic pivot path | **PHAC-NML ECTyper team** (Catherine Yoshida / current lead) | https://github.com/phac-nml/ecoli_serotyping | `site:github.com phac-nml ecoli_serotyping CODEOWNERS` | Most-likely acquirer if Gate B fails — audit layer is a natural ECTyper PR |
| 2 | nf-core pipeline | **nf-core/funcscan or /bacass maintainer** (J. Fellows Yates et al.) | https://nf-co.re/funcscan | `site:nf-co.re funcscan maintainers` | nf-core acceptance = realistic v0 distribution path; adoption + channel in one reply |
| 1 | PH surveillance | **PulseNet CDC sequencing lead** (Heather Carleton historic / current) | https://www.cdc.gov/pulsenet/contact-us/index.html | `site:cdc.gov pulsenet bioinformatics lead` | Daily WGS on E. coli outbreaks; auditable layer fits reproducibility discipline; slowest reply |

**Conflict rule:** Profile 3 = Hazen (UMD), NOT Shannon Manning (MSU) — Manning is reserved for the separate Whittam STEC substrate-access ask (`research_outputs/whittam-stec-center-contact-draft-2026-05-27.md`). Do not double-send.

---

## Email body (post-Gate-A; copy-paste, personalize salutation)

**Subject:** `Quick pilot question: local audit companion for E. coli pathotype interpretation`

```text
Hi [Name],

I'm working on a local-first CLI for E. coli pathotype interpretation. The narrow goal is not to replace existing callers broadly, but to add a more honest audit/QA layer for awkward cases — returning explicit HYBRID / AMBIGUOUS / UNCLASSIFIED states plus provenance, instead of forcing overconfident labels.

We've already passed a small local sanity gate on a frozen 5-strain panel (EHEC, ETEC, EAEC, UPEC-compatible, commensal). We have not started the heavy cohort build yet because we're gating that on whether real users would actually pilot this.

If this were available as a local CLI that fits beside an existing workflow, would you consider piloting it within 60 days?

Simple answer is enough:
- yes
- yes, with conditions
- no

If "with conditions," the most useful detail would be what exact use case or pain point would make it worth testing.

Thanks,
[Your name]
[independent researcher]
```

**Follow-up (once, after 5 business days, no reply):**

```text
Hi [Name],

Following up once on the note below. We're using this outreach as a go / no-go gate before investing in a larger benchmark build.

Even a one-line answer is useful:
- yes
- yes, with conditions
- no

Thanks,
[Your name]
```

**If asked "how is this different?":**
> The narrow difference is the audit layer: explicit ambiguous/hybrid/unclassified states, provenance, and a workflow aimed at honest interpretation in edge cases rather than forcing a clean single label every time.

**If asked "what stage is it in?":**
> Early pilot stage. We passed a small local sanity gate, but we have deliberately not started the larger build yet because we want to know whether real users would actually pilot it.

---

## One-page concept attachment

Send as a separate `.md`/PDF: use `reports/pathotype_gate_b/pathotype_gate_b_one_pager_2026-05-29.md` verbatim (post-Gate-A version; lists the 11-class surface + wedge + ideal pilot shape).

---

## Reply scoring (the real Gate B layer)

Score each reply with `reports/pathotype_gate_b/pathotype_gate_b_response_rubric_2026-05-29.md`. Log in `reports/pathotype_gate_b/pathotype_gate_b_response_log_template_2026-05-29.csv`.

| Outcome | Next move |
|---|---|
| **2+ credible yes** | Gate B PASS → authorize Horesh substrate build |
| 1 credible yes | HOLD → more outreach or refine wedge |
| 0 credible yes | FAIL → stop / pivot to ECTyper PR before heavy build |

`credible yes` = relevant archetype + explicit pilot/eval willingness + near-term + not mere politeness. Vague "interesting, keep me posted" does NOT count.

---

## Hard "do not" (until ≥2 yes land)

- do not rerun Gate A
- do not start the Horesh heavy build
- do not reopen AMR work
- do not broaden Phase 4 implementation before replies land

## Downtime work while awaiting replies (bundle-sanctioned, non-overlapping)

- commensal-substrate follow-up
- pathotype audit-policy / opacity proposal

(Note: `research_outputs/horesh2021-file-f1-label-provenance-audit-2026-05-29.md` is substrate groundwork — fine to refine, but do NOT escalate into the heavy build before the gate clears.)

## Provenance

- Contacts: `research_outputs/gate-b-target-contact-shortlist-2026-05-27.md` (this laptop, 2026-05-28)
- One-pager + drafts + rubric + CSV + digital proxy + Gate A result: `reports/pathotype_gate_b/*` (workhorse bundle, 2026-05-29)
- Superseded for send: `research_outputs/gate_b_cold_email_packet_2026-05-27.md` (pre-Gate-A framing; kept for history)
