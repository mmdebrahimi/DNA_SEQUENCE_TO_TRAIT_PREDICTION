<!-- memo-schema: 0.4 -->
# Unsupported / verify-needed rows — ΔΔG→fold-change calibration (2026-07-09)

Rows held back from the supported memo. Not fabricated — real claims, but with a locator/provenance defect the intake floor flags.

| Claim | Value | Why held back | How to promote |
|---|---|---|---|
| FoldX stability ΔΔG on the S461 benchmark | PCC 0.30, RMSE 1.91 kcal/mol | URL↔claim MISMATCH: the S461 FoldX figure surfaced in a benchmark *discussion*, but the URL attached at capture (a DDGemb/Bioinformatics article) is not verified as the S461-FoldX primary source. Claim is plausible + consistent with FoldX's known 0.3–0.5 PCC range, but the exact locator is unconfirmed. | WebFetch the actual S461/S669 benchmark paper (e.g. the ThermoMutDB / S461 source) and confirm the FoldX row verbatim. |
| DDGemb multi-point (PTmul-NR) stability ΔΔG | PCC 0.59, RMSE 2.16, MAE 1.59 kcal/mol | Real DDGemb result but full-text not WebFetch-verified (captured from search summary); the URL may point to the DDGemb article but the PTmul-NR table cell was not directly read. | WebFetch the DDGemb Bioinformatics 2025 article; confirm the PTmul-NR row. Relevant because resistance is often multi-mutation and accuracy degrades there (RMSE 2.16 → ~30× fold error). |

**Note:** the FEP+/Abl (88%), kinase MUE (~1.0), and FoldX/Rosetta homology-sensitivity rows are in the SUPPORTED memo at *medium* (real primaries, abstract-via-search provenance downgrade), not here — they have coherent locators and quotes; only full-text WebFetch was blocked by paywall/redirect.
