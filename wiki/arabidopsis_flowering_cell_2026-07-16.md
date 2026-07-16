# Arabidopsis flowering-habit cell — the deterministic answer to a closed embedding negative (2026-07-16)

The plant analog of the AMR/pigmentation curated cells, and the **deterministic counterpart to a closed
negative**: the Arabidopsis flowering-time EMBEDDING test failed under de-confounding (2026-06-12,
PlantCaduceus within-group r² −0.13 — it learned population structure, not the causal signal; the 3rd
de-confounded embedding failure across the kingdom boundary). This cell reads the same trait's mechanism
**directly** from curated causal loci. `dna_decode/organism_rules/arabidopsis_flowering.py`; NON-frozen
(the frozen decoder surface is byte-unchanged, `verify_lock OK`).

## The rule (and why it's a NON-frozen `organism_rules` cell)

FLC is a dosage-dependent floral **repressor**; FRI is required for **high** FLC. So:

```
late (winter-annual, vernalization-requiring)  iff  FRI functional AND FLC strong      # multi-locus AND
early (summer-annual, rapid-cycling)           iff  FRI LoF  OR  FLC weak/null         # either route
```

That **AND-across-two-loci + epistasis** is a shape the frozen count/OR `amr_rules.DRUG_RULE` engine cannot
express — the same reason the TMP-SMX `sul AND dfr` overlay and the TB cell live outside it.

## Result — all four literature anchors reproduce (real CLI)

| accession / line | FRI | FLC | call | note |
|---|---|---|---|---|
| **Col-0** | LoF (premature stop) | strong | `summer_annual_early` (medium) | the reference rapid cycler |
| **Col-FRI[Sf-2]** | functional | strong | `winter_annual_late` (high), vernalization **required** | the classic winter-annual line |
| **Ler** | LoF (start-codon del) | weak | `summer_annual_early` | double hit |
| **Da(1)-12 / Shakhdara** | **functional** | weak | `summer_annual_early` (high) | **a naive "FRI decides" rule mis-calls this LATE** |
| **Van-0 / Bur-0** | functional | null | `summer_annual_early` | the FLC route (nonsense / aberrant splice) |
| unknown either locus | — | — | **ABSTAIN** | the polygenic residue is not guessed |

**The Da(1)-12 anchor is the load-bearing test** — it's the case that separates the correct two-locus rule from
the naive one, and it's pinned by `reference_integrity_ok()` + a dedicated test. A corrupted/fabricated catalog
fails that guard loudly.

## Honest scope (the rails that ship in every call)

1. **PARTIAL — habit, not days.** FRI/FLC explains ~40–70% of long-day flowering-time variation; the rest is
   polygenic + environment. This decodes the **HABIT/direction** (early vs late), NOT quantitative
   days-to-flower (mirrors AMR R/S vs MIC-continuous). Every output carries the scope limit.
2. **The FRI route is confidence-capped MEDIUM** by the **Lz-0 counterexample** — an FRI deletion that is
   nonetheless LATE via FRI-*independent* FLC upregulation (FRL1/FES1-class activators). The two-locus rule
   cannot see it; that's surfaced as `undetectable_mechanisms`, not hidden. (The FLC route is `high` — FLC is
   downstream, so a lost repressor is the more robust call.)
3. **v0 input = allele calls** (the wheel-only `--observed` pattern the HIV/fungal cells use). Genome-mode
   (detect the FRI-Col stop / FRI-Ler start-codon deletion from sequence) needs **verified variant
   coordinates** — a v0.1 follow-on, deliberately **not fabricated**.
4. **Faithful-to-literature**, not a new model: it applies published functional-allele assignments
   (Johanson 2000 FRI; Michaels 2003 PNAS weak-FLC; Werner 2005 FRI-independent). The value is the encoded
   rule + abstain discipline + the tool surface.

## Run

```bash
dna-decode flowering --fri Sf-2 --flc Col          # -> WINTER_ANNUAL_LATE, vernalization required
dna-decode flowering --fri Col  --flc Col          # -> SUMMER_ANNUAL_EARLY (Col-0)
dna-decode flowering --fri Sf-2 --flc "Da(1)-12"   # -> EARLY via the FLC route (functional FRI!)
dna-decode flowering --list-alleles                # the curated catalog + per-allele sources
dna-flowering --fri Ler --flc Ler --json
```

## What this establishes

The deterministic curated-catalog paradigm now spans **bacterial + viral + fungal AMR → human PGx → human
visible traits → a plant developmental trait**. And it makes the embedding negative sharper rather than
softer: **the flowering signal was decodable all along — deterministically, not by embedding.** The embedding
failed because it learned lineage; the catalog reads mechanism.

**Named v0.1 follow-ons:** genome-mode (verified FRI/FLC variant coordinates → call alleles from a sequence);
AraPheno validation (free, 1,135 accessions — needs an accession→allele mapping, which the catalog does not
yet carry); more FRI/FLC natural alleles. 16 tests. Frozen decoder surface byte-unchanged.
