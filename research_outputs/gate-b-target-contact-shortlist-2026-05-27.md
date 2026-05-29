# Gate B Target-Contact Shortlist (2026-05-27)

> 5-contact shortlist for the user to fire the Gate B cold-email packet at (`research_outputs/gate_b_cold_email_packet_2026-05-27.md`). Each row: role + organization + stable public contact URL + suggested-name + verification-search-query. **Names are training-knowledge-anchored and need 1-line user verification before send — the search query is the verification path.**
>
> Discovery-machine constraint: rate-limited at WebSearch, so this shortlist is built from durable public org URLs + role anchors. User has 5-10 min of verification work per row to confirm current named contact + email.

## Why these 5 (profile rationale)

Per the Gate B packet's contact-profile guidance:
- Profile 1: Public-health surveillance lab (PulseNet CDC / ECDC / UKHSA)
- Profile 2: Food-safety pathogen-detection contact (FDA GenomeTrakr / USDA FSIS / nf-core maintainer)
- Profile 3: Clinical-microbiology research PI doing E. coli pathotype work
- Profile 4: Bioinformatics pipeline maintainer building E. coli typing
- Profile 5: Wildcard / cross-organism analog

The shortlist below assigns 1 contact per profile. Whittam STEC Center (Shannon Manning) is intentionally NOT used for Profile 3 — that conflicts with the separate `whittam-stec-center-contact-draft-2026-05-27.md` substrate-access ask. Tracy Hazen (DECA paper first author at Maryland) is the cleanest Profile 3 alternative.

## Shortlist

| # | Profile | Role + organization | Stable public contact URL | Suggested named contact (verify) | Verification search query | 1-line rationale |
|---|---|---|---|---|---|---|
| 1 | PH surveillance lab | PulseNet sequencing lead — CDC Enteric Diseases Laboratory Branch | https://www.cdc.gov/pulsenet/contact-us/index.html | Heather Carleton (historic role) OR current PulseNet bioinformatics lead | `site:cdc.gov pulsenet bioinformatics lead site:cdc.gov OR linkedin "PulseNet" "sequencing"` | PulseNet runs daily WGS on E. coli outbreak isolates; auditable pathotype-call layer fits their reproducibility discipline; ECTyper is already in their toolchain |
| 2 | Food-safety + bioinformatics pipeline | nf-core/funcscan or nf-core/bacass maintainer | https://nf-co.re/funcscan + https://nf-co.re/bacass (each pipeline page lists named maintainers) | Jasmin Frangoulidis OR James Fellows Yates (nf-core/funcscan maintainers per current pipeline page) | `site:nf-co.re funcscan maintainers` and `site:github.com nf-core funcscan README maintainers` | nf-core acceptance is the realistic v0 distribution path; if a funcscan maintainer says "I'd add this as a module", that's adoption signal AND distribution path in one reply |
| 3 | Clinical-micro research PI | Tracy H. Hazen (DECA paper first author) — Institute for Genome Sciences, Univ. Maryland School of Medicine | https://www.medschool.umaryland.edu/profiles/hazen-tracy/ (faculty profile; public email link) | Tracy Hazen directly | `"Tracy Hazen" site:umaryland.edu` and `"Tracy Hazen" "Diarrheagenic Escherichia coli Collection"` | Authored the DECA 2012 paper; familiar with the substrate-circularity problem; her grant-funded follow-ups (Hazen 2013, Hazen 2016) are the dominant independent-label sources in Horesh 2021 (n=139 records); strong "would-use-or-route-to-actual-user" signal |
| 4 | Bioinformatics pipeline maintainer | Torsten Seemann — ABRicate / mlst / shovill maintainer, Melbourne | https://github.com/tseemann (public GitHub profile with contact info) + http://thegenomefactory.blogspot.com/ | Torsten Seemann directly | `"Torsten Seemann" abricate contact` and `"Torsten Seemann" Melbourne` | ABRicate is the most-installed gene-presence frontend in bacterial bioinformatics; if Torsten Seemann says "this fills a gap ABRicate users hand-roll", that's strongest possible signal from the maintainer of the user-base the tool would adopt-from |
| 5 | Wildcard | PHAC-NML ECTyper team — Public Health Agency of Canada, National Microbiology Laboratory | https://github.com/phac-nml/ecoli_serotyping (GitHub repo; lists current maintainers + issues for direct contact) | Catherine Yoshida OR current PHAC-NML lead per repo's CODEOWNERS / README | `site:github.com phac-nml ecoli_serotyping CODEOWNERS` and `"ECTyper" phac-nml lead` | Most likely *acquirer* if Gate B fails — the audit-packet/abstention layer is a natural ECTyper PR; getting an early signal from the actual maintainer team is high-leverage. If they say "we'd merge this PR", the project might pivot from standalone-tool to ECTyper-extension entirely |

## Conflict-avoidance reminder (from Gate B packet §3 line 1)

> "Profile 3 may overlap with the Whittam STEC Center draft email. If you've already sent Whittam the substrate-access ask, do NOT also send Gate B to Shannon Manning — pick a different Profile 3 contact (Tracy Hazen at Maryland is the cleanest alternative)."

This shortlist uses Tracy Hazen (UMD), not Shannon Manning (MSU), for Profile 3 — exactly per the cleanup rule.

## Send-order recommendation (not mandatory)

Suggested order: 4 (Seemann) → 3 (Hazen) → 5 (PHAC-NML) → 2 (nf-core) → 1 (PulseNet). Rationale: contact 4 is the highest-leverage / highest-likely-to-reply (open-source maintainer with public-facing presence); 3 is the next-highest (academic; UMD email + active on professional networks); 5 is the strategic pivot signal (do-they-want-this-as-a-PR question); 2 is the distribution-channel question; 1 is the institutional public-health question (slowest reply, most-bureaucratic). Maximizes early signal; lets you stop after 2-3 high-confidence "no/yes" answers.

## Verification time per row

~5-10 min total. Each search query above is single-domain-scoped (CDC / nf-core / UMD / GitHub / PHAC-NML); a quick visit to each row's URL + email confirmation is the verification path. No multi-step research needed.

## Use with the Gate B email body

Per `research_outputs/gate_b_cold_email_packet_2026-05-27.md` email body — paste verbatim, personalize the salutation per row, attach the 1-page concept attachment from the packet. Track replies in the packet's tracking sheet template. Gate B verdict thresholds (≥2/5 yes = PASS) unchanged.

## What this shortlist is NOT

- Not auto-verified — names + role assignments are training-knowledge-anchored; the search queries are the verification path
- Not a replacement for sending the emails — discovery-machine cannot send; user-driven send is the binding action
- Not authoritative on current emails — institutional emails (CDC, ECDC, UKHSA, UMD, MSU, PHAC) rotate as people change roles. Each row's public URL is the durable anchor; named individuals need a 30-second verification step
- Not Whittam — that's a separate substrate-access ask with its own draft at `research_outputs/whittam-stec-center-contact-draft-2026-05-27.md`
