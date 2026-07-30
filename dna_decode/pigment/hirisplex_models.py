"""HIrisPlex-S hair + skin (+ eye) multinomial models, loaded from coefficients RECOVERED from the
deployed erasmusmc webtool (user-authorized model extraction, 2026-07-30).

WHY recovered, not transcribed: the published papers give the beta matrices but NOT the per-category
INTERCEPTS (webtool-server-side only) — so the models are under-specified by the public record and cannot
be reproduced by transcription (see wiki/pigment_hirisplex_coefficient_sourcing_2026-07-30.md). Instead the
deployed softmax models were recovered by querying the webtool with a designed genotype basis (all-zero +
per-SNP unit vectors) and LS-fitting the reference-parameterized softmax, then VALIDATED on 20 random
held-out genotypes: max |ΔP| eye 6e-15 / hair 6e-16 / skin 9e-3 — i.e. the offline models reproduce the
deployed webtool to machine precision (eye/hair) or <1% (skin, whose complete-separation betas are the only
imprecision). These are FACTUAL model coefficients (not fabricated); provenance + validation ride in the JSON.

The eye model here is the webtool's canonical IrisPlex; the shipped `irisplex.py` eye cell is unchanged
(this recovery CONFIRMS it — small canonical deltas). Consumers of visible-trait prediction use these via
`dna_decode.pigment.multinomial.predict`. Benign visible-trait genetics, NOT a forensic tool.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from dna_decode.pigment.multinomial import PigmentModel, PigmentSNP

_COEF_PATH = Path(__file__).parent / "hirisplex_coefficients.json"
TRAITS = ("eye_colour", "hair_colour", "skin_colour")


@lru_cache(maxsize=1)
def load_hirisplex_models() -> dict[str, PigmentModel]:
    """{trait -> PigmentModel} for eye/hair/skin, built from the recovered coefficient table."""
    d = json.loads(_COEF_PATH.read_text(encoding="utf-8"))
    out: dict[str, PigmentModel] = {}
    for trait in TRAITS:
        m = d[trait]
        ref = m["reference"]
        ordered = (ref,) + tuple(c for c in m["categories"] if c != ref)  # PigmentModel: categories[0]=ref
        counted = m["counted_alleles"]
        betas = m["betas"]
        snps = tuple(
            PigmentSNP(rsid=rsid, counted_allele=allele,
                       betas={c: betas[c][rsid] for c in betas if rsid in betas[c]})
            for rsid, allele in counted.items()
        )
        out[trait] = PigmentModel(
            trait=trait, categories=ordered,
            intercepts={c: m["intercepts"][c] for c in m["intercepts"]},
            snps=snps, source="recovered from HIrisPlex-S webtool 2026-07-30 (held-out-validated)")
    return out


def provenance() -> dict:
    return json.loads(_COEF_PATH.read_text(encoding="utf-8"))["_provenance"]
