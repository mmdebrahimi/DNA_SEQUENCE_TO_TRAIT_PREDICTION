# Gate B Cold-Email Packet — pre-Gate-A concept validation (2026-05-27)

> **Calibration:** this packet is **concept-validation-tier**, NOT a polished v0 one-pager. Gate A (workhorse-driven manual VF + ETECFinder sanity check on 5 strains) has not fired yet. If Gate A surfaces install / DB schema / decision-table-uncomputable friction, the user-facing promise will change — so this packet stays at the "would you consider using this kind of tool" level, not "here's what we're shipping." Re-issue polished v0 one-pager AFTER Gate A passes.
>
> **Purpose:** validate adoption-signal demand for a deterministic auditable pathotype-call layer on top of pinned VirulenceFinder + ETECFinder. Per /idea-validation-council Gate B threshold: PASS = ≥2/5 contacts say yes; FAIL = ≤1/5 say yes (→ KILL or PIVOT to ECTyper PR per T2 trade-off in project ledger).
>
> **Discipline:** keep contact list to 5. More invitations dilute the signal + increase the "no reply" rate.

---

## Suggested 5 target-contact profiles

Pick contacts who span both demand sides + tool-maintainer sides; the 5 slots are not rigid (substitute equivalents you actually know):

| # | Profile | Why this contact | Example targets |
|---|---|---|---|
| 1 | **Public-health surveillance lab** (state or national level) | Daily E. coli pathotype calls; existing ECTyper / VirulenceFinder workflow; institutional pull on tool adoption | PulseNet (CDC) sequencing lead; ECDC E. coli surveillance contact; Public Health England (UKHSA) E. coli reference lab; state-level food-safety lab (e.g., Texas DSHS, Minnesota MDH) |
| 2 | **Food-safety pathogen-detection contact** | Audit-grade provenance + DB-version-pin discipline matters for traceback investigations + recall decisions | FDA GenomeTrakr team; USDA FSIS sequencing contact; nf-core/funcscan or nf-core/bacass maintainer |
| 3 | **Clinical-microbiology research PI** doing E. coli pathotype work | Direct user; publishes pathotype papers; will say if the audit-packet abstraction would replace their existing manual cluster-collapsing | Tracy Hazen (DECA author, U. Maryland); Mark Achtman / EnteroBase team; Shannon Manning (Whittam STEC Center successor) — **NOTE: same contact as Whittam draft email; consolidate** |
| 4 | **Bioinformatics pipeline maintainer** building E. coli typing pipelines | Will say if this slots into Snakemake/Nextflow / Galaxy / IRIDA workflows; tool-distribution gate | nf-core E. coli pipeline maintainer; Torsten Seemann (ABRicate maintainer); PHAC-NML ECTyper team (could be the natural acquirer if Gate B fails — see T2 pivot path) |
| 5 | **Wildcard / cross-organism analog** | If the audit-packet pattern resonates outside E. coli, the scaffold has wider scope than the current project frames | Kleborate maintainer (Klebsiella; sibling ontology); SeqSphere/Ridom team (commercial pathogen typing); WHO sequencing-for-surveillance contact |

**Do NOT cc all 5 on a single thread** — that creates social-proof gaming. Send 5 individual emails.

**Conflict avoidance:** Profile 3 may overlap with the Whittam STEC Center draft email (`research_outputs/whittam-stec-center-contact-draft-2026-05-27.md`). If you've already sent Whittam the substrate-access ask, do NOT also send Gate B to Shannon Manning — pick a different Profile 3 contact (Tracy Hazen at Maryland is the cleanest alternative).

---

## Email body (~150 words, copy-paste ready)

> **Subject:** brief question — would a deterministic auditable E. coli pathotype CLI be useful to you?

> Hi [Name],
>
> I'm exploring whether to build an open-source CLI that turns CGE VirulenceFinder + ETECFinder gene calls into a deterministic, auditable, abstention-aware E. coli pathotype label — with explicit `HYBRID` / `AMBIGUOUS` / `UNCLASSIFIED` output classes, full DB-version-pin provenance per call, and a side-by-side diff against the underlying caller output.
>
> The pitch is not a new gene caller — it's a thin, transparent rule layer on top of an existing one, so the call is reproducible and the abstention is explicit when the evidence is mixed.
>
> Before I commit to the build, I'm asking 5 people who would be plausible users two questions:
>
> 1. Would you install + try this in a real workflow within ~60 days of it shipping?
> 2. What existing tool or workflow (if any) would it replace for you?
>
> If the second question is uncomfortable to answer in detail, even a one-line gut reaction is useful.
>
> Thanks for your time.
>
> Best,
> [your name]
> [your affiliation, if any; "independent researcher" is fine]
> [your email]

---

## One-page concept attachment (copy as a separate `.md` or PDF)

> # E. coli Pathotype Decoder — Concept (Pre-Gate-A Validation)
>
> ## Pitch
>
> Open-source CLI: FASTA in → multilabel virulence-cluster profile + derived pathotype call + audit packet, on top of pinned CGE VirulenceFinder 3.2.0 + ETECFinder gene calls. Deterministic + auditable + abstention-aware. **Not a new gene caller. Not an ML classifier.** Just an honest rule + reporting layer.
>
> ## What the output looks like
>
> Default surface (honest 11-class):
> `EHEC_COMPATIBLE`, `STEC_NON_LEE`, `tEPEC_COMPATIBLE`, `aEPEC_COMPATIBLE`, `ETEC_COMPATIBLE`, `EAEC_COMPATIBLE`, `UPEC_COMPATIBLE`, `HYBRID`, `AMBIGUOUS`, `UNCLASSIFIED`, `COMMENSAL_LOW_MARKER_BURDEN`
>
> `--legacy-6class` flag collapses to {EPEC, EHEC, ETEC, UPEC, EAEC, commensal} for downstream pipelines that expect the legacy surface.
>
> Each call carries a JSON audit packet: per-cluster %ID + contig position + accession + source DB + DB version + DB checksum + parameters. The `caller_is_independent_baseline: false` honesty flag is set explicitly — the tool does not pretend to be an independent baseline against VirulenceFinder.
>
> ## What's different vs existing options
>
> | vs | What this adds |
> |---|---|
> | VirulenceFinder web/CLI | a deterministic resolver layer on top + explicit abstention class + DB checksum in provenance |
> | ECTyper v2.0 | explicit `HYBRID` / `AMBIGUOUS` / `UNCLASSIFIED` / `_COMPATIBLE` honesty qualifier; UPEC/ExPEC coverage; per-call DB checksum |
> | ABRicate + hand-rolled rules | ships the rules as code + tests + version pin so you don't re-implement them each time |
>
> ## What it does NOT add
>
> - new gene calls (uses pinned VirulenceFinder + ETECFinder)
> - a trained classifier (deferred to v0.1+, conditional on independent-label substrate)
> - a web UI (CLI only)
>
> ## Pre-build gates (this is why I'm asking before building)
>
> 1. **Gate A** — manual sanity check (5 strains): tool installs + decision table is computable from raw caller outputs.
> 2. **Gate B (this email)** — adoption-signal: ≥2 of 5 people I'm asking would install + try it in a real workflow within 60 days.
>
> Both gates have to pass before I commit 3 months to the v0 build. If they fail, the work pivots to a contribution against an existing tool (ECTyper PR most likely) or stops entirely.
>
> ## Status
>
> - Substrate availability validated: Horesh 2021 curated 10,146-genome collection has 2,077 records with publication-extracted (non-circular) pathotype labels — enough for ExPEC/EPEC/ETEC class floors. STEC/EHEC/EAEC/COMMENSAL supplements pending direct contact with Whittam STEC Center.
> - Architecture decided: deterministic multilabel cluster resolver + abstention. No classifier in v0.
> - Output schema drafted (23-marker catalog, 11-class decision table, abstention rules, provenance JSON).
> - Not yet written a line of v0 code.
>
> ## What I'd most want to hear back
>
> - "Yes, would use it" → confirms demand; thank you for the signal
> - "No, would not use it" → please share what you'd use instead so I can route to the actual unmet need
> - "Would use it but only if X" → most valuable answer; tells me what to scope-cut or scope-add

---

## Tracking sheet template (for the user)

| # | Contact | Profile | Sent date | Replied date | Q1: would install/try? | Q2: what would it replace? | Notes |
|---|---|---|---|---|---|---|---|
| 1 | | | | | | | |
| 2 | | | | | | | |
| 3 | | | | | | | |
| 4 | | | | | | | |
| 5 | | | | | | | |

Gate B verdict (after replies in, or 7 days after last send, whichever first):
- **PASS:** ≥2 of 5 said "yes" to Q1 → proceed to v0 build conditional on Gate A also PASSing.
- **FAIL:** ≤1 of 5 said "yes" → KILL the v0 build OR pivot to ECTyper PR per project ledger T2 trade-off.

---

## Sending tips

- **Send to one contact at a time, not bcc-blast.** Personalize the salutation. Cold-email reply rates are typically 5-15%; quality matters more than volume at N=5.
- **Window:** allow 5-7 business days for replies before declaring Gate B done. Send on Tue-Thu; avoid Mon (deluge) and Fri (lost).
- **Follow-up once max:** if no reply after 5 business days, send a single short bump ("any thoughts on the below?"). Do not chase further.
- **If a contact wants a call instead of an email reply:** great signal; even 1 call counts as a strong "yes" for Q1 + much richer answer to Q2.
- **Wildcard contact (profile 5) may be a low-yield slot** — keep expectations low; if no reply, that's fine and doesn't move the verdict against you.

## What this packet is NOT

- Not a polished v0 product sheet — explicitly labeled pre-Gate-A concept validation.
- Not a sales pitch — the framing leads with "before I commit to the build."
- Not a request for funding / collaboration — just adoption-signal validation.
- Not for cc / mass-list — N=5 individual asks.
- Not the Whittam contact — that's a separate email at `research_outputs/whittam-stec-center-contact-draft-2026-05-27.md` requesting DECA collection access (different ask, different framing).
