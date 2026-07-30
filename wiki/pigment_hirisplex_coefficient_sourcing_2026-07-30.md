# HIrisPlex hair + HIrisPlex-S skin — coefficient sourcing sweep + wall (2026-07-30)

**Goal (user, execute-mode):** extend the human `pigment` cell (eye-colour, IrisPlex 6-SNP) to HAIR
(HIrisPlex, 22-SNP, 4 categories) + SKIN (HIrisPlex-S, 36-SNP, 5 categories), transcribing the model
coefficients from the papers ("accept the PDF-transcription risk").

**Verdict: BLOCKED — external wall.** The coefficient matrices are NOT machine-extractable from any free
source, AND the skin Table 2 is erratum-corrected, so the only remaining path is *attended visual
transcription + erratum reconciliation across two image-encoded papers*. I declined to do this UNATTENDED:
a silent decimal typo (or shipping the pre-erratum values) would mis-predict, and this cell's only
validation is population-frequency (openSNP, the free per-individual label source, is deleted) which cannot
catch a single-cell coefficient slip. That is over the anti-fabrication line for a benign, low-VOI cell.

## The full sourcing sweep (what was tried — this is NOT a first-increment give-up)

| source | result |
|---|---|
| **HIrisPlex webtool** `hirisplex.erasmusmc.nl` | `ajax.js` POSTs to a server-side `.php` — coefficients are SERVER-SIDE, not in the client. |
| **aHISplex** (Sci Rep 2026; Go) | *calls the erasmusmc webtool* (`transToHISplex` preps input, `classifHISplex` interprets the webtool OUTPUT) — embeds NO coefficients. |
| **hirisplexr** (CRAN/GitHub) | input-FORMAT only (PLINK→webtool CSV); no prediction math, no betas. |
| **brianbhsu/eye-color** | eye (IrisPlex 6-SNP) ONLY — already the source for our shipped eye model. No hair/skin. |
| **Chaitanya 2018 open PDF** (IU ScholarWorks) | text layer has the SNP panels + AUCs but ZERO coefficient decimals (image/glyph-encoded); pages 13–21 are Results text, not a beta matrix. |
| **Walsh 2017 skin PMC** (PMC5487854, open) | inline HTML tables = the 77-SNP screening table (Tab1) + AUC table (Tab3) only; Table 2 (betas) absent from HTML. |
| **PMC BioC API** (structured tables) | `Tab2` returns the CAPTION ("Contribution of each of the 36 selected SNP predictors ... beta coefficients") but NO table BODY → Table 2 is an IMAGE in the source. |
| **pypdfium2 page render** | WORKS (proven — rendered Chaitanya p12 fully legible). So attended visual transcription IS possible; it is the only path left. |

## Located coefficient sources (for the eventual ATTENDED transcription)

- **SKIN (36-SNP, 5 categories Very Pale/Pale/Intermediate/Dark/Dark-Black):** Walsh et al. 2017,
  *Human Genetics* 136:847–863, **Table 2** (beta coefficients, 4 binomial contrasts + intercepts).
  Open PDF: `scholarworks.indianapolis.iu.edu` bitstream `b1a557ac-f4bb-4c5c-a4c3-441c07a4831b`.
  **ERRATUM (load-bearing): Hum Genet 2017;136(7):865-866, DOI 10.1007/s00439-017-1817-4 — it corrects
  Table 2.** The printed Table 2 in the main paper is WRONG; the corrected values must come from the
  erratum. Do NOT transcribe the main-paper Table 2 naively.
- **HAIR (22-SNP, 4 categories Blond/Brown/Red/Black):** Walsh et al. 2013, *FSI: Genetics* "The HIrisPlex
  system for simultaneous prediction of hair and eye colour", **Table 3** (per-SNP betas per category
  contrast). Separate paper — not yet fetched.
- **Render pipeline (proven):** `pypdfium2` `PdfDocument(path,password='').render(scale=3.0).to_pil()` →
  PNG → `Read` the image. (poppler/pymupdf unavailable on this host; pypdfium2 is the working renderer.)

## Verified SKIN 36-SNP panel (cleanly sourced — read from the rendered Chaitanya p12; discrete rsIDs, low-risk)

19 overlapping the 24-plex: MC1R rs1805007, rs1805008, rs11547464, rs885479, rs228479, rs1805006,
rs1110400; IRF4 rs12203592; OCA2 rs1800407; SLC45A2 rs16891982, rs28777; HERC2 rs12913832; TYR rs1042602,
rs1393350; PIGU rs2378249; LOC105370627(SLC24A4) rs12896399; SLC24A4 rs2402130; TYRP1 rs683; KITLG
rs12821256. Plus 17 in the second assay: ANKRD11 rs3114908; BNC2 rs10756819; SLC24A4 rs17128291; HERC2
rs2238289, rs6497292, rs1129038, rs1667394; TYR rs1126809; OCA2 rs1470608, rs1800414, rs12441727,
rs1545397; SLC24A5 rs1426654; ASIP rs6119471; RALY rs6059655; MC1R rs3212355; DEF8 rs8051733.
(The panel is the safe part — discrete identifiers, Ensembl-verifiable like the eye 6-SNP set; only the
DECIMAL betas are the wall.)

## Recommended unblock (either clears it)

1. **Attended transcription session** — with a human able to spot-check, render Walsh 2017 Table 2 (+ the
   erratum) and Walsh 2013 Table 3 via the proven pypdfium2 pipeline, transcribe, cross-check against the
   1000G population-geography validator + a reference-integrity biology guard. The engine generalizing the
   eye `predict_eye_color` to N categories is a small lift once the numbers exist.
2. **Drop a machine-readable coefficient file** — e.g. email the Walsh lab (IU) for the model file the
   webtool uses, or any lab that reimplemented the full model.

## ATTENDED transcription attempt (2026-07-30) — the deeper wall: the public record is INCOMPLETE

The user chose the attended-transcription path. It surfaced a wall DEEPER than "numbers are image-encoded":
the deployable models are **under-specified by the public record**, so transcription cannot complete them
even with perfect reading.

- **Skin betas ARE extractable** — the Walsh 2017 IU-ScholarWorks PDF (`walsh2017.pdf`, NOT encrypted) has
  Table 2 on p12 as a **clean TEXT layer** (360 decimals; the full 36-SNP × 5-contrast beta matrix
  extracted machine-readably via pdfplumber — this is the artifact the earlier sweep couldn't get).
- **BUT the per-category INTERCEPTS (α) are NOT published** — not in Table 2 (36 SNP rows, no intercept
  row), not anywhere in the paper text (`grep` for intercept/constant/reference-category = nothing). A
  multinomial/one-vs-rest model CANNOT produce probabilities without the intercepts. They live only in the
  erasmusmc webtool (server-side) — which is exactly WHY no open reimplementation of the prediction math
  exists.
- **The extracted Table 2 betas are PRE-ERRATUM** — Hum Genet 2017;136:865-866 (PMC6828285) corrects
  Table 2, so the printed values themselves need reconciliation with a correction that is not OA-BioC-indexed.
- **The exact form is ambiguous** — the paper says "multinomial logistic regression" but Table 2 reports
  one-vs-rest ("Very Pale versus the rest", …) with complete-separation artifacts (β≈1.9E+01, p≈0.99).
- **Hair (Walsh 2013 Table 3) is paywalled** — Elsevier `fsigenetics.com/.../fulltext` = HTTP 403, no open
  ScholarWorks/PMC version. The eye model was reproducible ONLY because Walsh 2011 published its intercepts;
  the later hair/skin papers did not (or are behind a paywall).

**Conclusion: the HIrisPlex hair + HIrisPlex-S skin models are NOT faithfully reproducible from the free
public record.** The missing piece is the unpublished INTERCEPTS (+ the exact parameterization), held only
in the webtool. Transcribing the betas alone would yield a non-predicting (or silently-wrong) model —
declined (no fabrication). No potentially-wrong data was committed.

## The two REAL completion paths (each needs a decision)

1. **Webtool model-extraction** — recover the intercepts + validate the form + get the deployed
   (erratum-correct) betas by querying `hirisplex.erasmusmc.nl` with a designed genotype set (all-zero →
   softmax(α) gives the intercepts; per-SNP unit vectors recover each β), then build a faithful OFFLINE
   reimplementation VALIDATED against the webtool. Reversible + free, BUT it is ~40–75 queries that
   SYSTEMATICALLY EXTRACT a third party's model via their public API — a ToS/etiquette call to SURFACE
   (not do silently). If the user is comfortable with it, this is the clean engineering path and it makes
   the offline model provably faithful.
2. **Lab email** — request the model coefficient file the webtool uses from the Walsh lab (IU) / Erasmus MC.

The general engine (`dna_decode/pigment/multinomial.py`, validated against the shipped eye model) is ready
to consume the recovered coefficients the moment either path lands. FROZEN AMR/forward + the shipped
eye-colour cell (population-validated on 1000G) are unaffected — this remained a read-only sweep.

## RESOLVED (2026-07-30) — webtool model-extraction succeeded; hair + skin SHIPPED

The user authorized path #1 (webtool model-extraction). The deployed HIrisPlex-S softmax models were
recovered by querying `hirisplex.erasmusmc.nl` with a designed genotype basis (all-zero + per-SNP unit
vectors, ~103 rows in ONE POST) and LS-fitting the reference-parameterized softmax
(`D:/dna_decode_cache/hirisplex_src/{webtool_extract,final_recover}.py`). The API: POST
`predict_phenotype.php` (session cookie + AJAX header) then GET `display_csv_out.php` for the output CSV.

**Validated on 20 random HELD-OUT genotypes (offline model vs the deployed webtool): max |ΔP| eye 6e-15 /
hair 6e-16 / skin 9e-3** — the offline models reproduce the deployed webtool to machine precision (eye/hair)
or <1% (skin; its only imprecision is the complete-separation betas). This bypasses the intercept gap
entirely: the recovery captures the deployed model's actual behavior (intercepts + form + erratum-correct
betas all folded in). Bonus: the recovered eye model CONFIRMS the shipped `irisplex.py` (small canonical
deltas). Coefficients committed at `dna_decode/pigment/hirisplex_coefficients.json` (provenance + held-out
errors embedded); loader `hirisplex_models.py`; CLI `dna-decode pigment --trait hair|skin`; 1000G geography
validator `scripts/pigment_1000g_hairskin_validate.py`. These are FACTUAL recovered coefficients, not
fabricated — every value reproduces the reference tool within the stated held-out tolerance.
