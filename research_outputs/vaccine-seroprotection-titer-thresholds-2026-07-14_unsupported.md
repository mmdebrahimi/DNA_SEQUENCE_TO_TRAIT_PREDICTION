# Vaccine seroprotection titer thresholds — unsupported claims (V1 invocation)

> Slug: vaccine-seroprotection-titer-thresholds-2026-07-14. Captured 2026-07-14.

## Rejected rows

| Row content | Rejection reason | Suggested follow-up |
|---|---|---|
| Measles protection against severe disease = ≥200 mIU/mL | Low confidence: the directly-fetched JID systematic review returned VALUE NOT PRESENT for a separate severe-disease threshold; the ≥200 figure is commonly cited but not confirmed by the primary this run | Fetch Chen et al. (original measles challenge/outbreak data) or a WHO measles position paper for the ≥200 mIU/mL severe-disease claim verbatim |
| Poliovirus neutralizing seroprotection = ≥1:8 | Low confidence: standard convention but NOT fetched this run | Fetch a WHO polio position paper / peer-reviewed correlate review for the ≥1:8 neutralizing titer verbatim |

## Summary

- Total rejected: 2 · Low confidence (unconfirmed/not-fetched): 2 · Missing locator: 0 · Mapping-floor: 0 · Banned phrase: 0

## Documented absences of evidence

- No established serologic correlate of protection exists for **pertussis, mumps, HPV, RSV** (efficacy-based, no accepted antibody threshold) — omitted rather than fabricated.
- Rabies (0.5 IU/mL, RFFIT) carried at LOW in the supported memo: WHO-attributed but the confirming search hit an AUP content filter this run → direct WHO-position-paper fetch recommended before use.
- CDC HepB MMWR (rr6210a1.htm) returned HTTP 403 on direct fetch → the 10 mIU/mL HepB row rests on a search-summary of that page + broad corroboration (carried medium, not rejected).
