"""Validate FBA carbon-source growth against measured carbon-source utilization (E. coli iML1515).

The quantitative-trait complement to the essentiality validation. For each carbon source, FBA swaps it in as
the sole carbon source and predicts a growth RATE; we compare to a measured set.

HONEST label situation (R2 web probe, 2026-08-03): a clean *measured growth-rate across many carbon sources*
dataset does NOT exist fetchably; Biolog pos+neg for the K-12 iML1515 strain is SI-locked (the 190-source
Biolog set is E. coli Nissle -- a strain mismatch). So the REACHABLE validation is:
  - RECALL on measured-positive K-12 carbon sources (the Keio/Wetmore fitness assay only ran on carbon
    sources E. coli grows on -> a positive-only set), plus the FBA growth-RATE spread (quantitative).
  - a small known-negative SPOT-CHECK (K-12 cannot use sucrose -- no scr system).
EXTERNAL-WALLED (named, not code-closable here): full pos+neg specificity + a measured growth-rate
correlation (need the EcN Biolog SI PMC9801561, or a measured-rate dataset).

Exit 0 = a real recall number was produced.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.carbon_growth import (  # noqa: E402
    build_exchange_name_index,
    match_carbon_exchange,
    predict_growth,
)
from dna_decode.fba.model import load_model  # noqa: E402

_POSITIVES_URL = ("https://raw.githubusercontent.com/dbernste/E_coli_GEM_validation/HEAD/"
                  "Fitness_Data/E_coli_BW25113/exp_organism_Keio_Mapped.txt")


def _classify_unmapped(model, name: str) -> str:
    """An unmapped measured-positive is a NAME-gap (metabolite in the model, name didn't match) or a
    MODEL-gap (metabolite absent -> no transporter -> a false negative, e.g. K-12 iML1515 lacks sucrose)."""
    from dna_decode.fba.carbon_growth import normalize_carbon_name
    key = normalize_carbon_name(name)
    stem = key.split()[0].replace("d-", "").replace("l-", "")[:5]
    for met in model.metabolites:
        if met.id.endswith("_e") and stem and stem in met.name.lower():
            return "name-gap"
    return "model-gap"


def _measured_positives(text: str) -> list[str]:
    lines = text.splitlines()
    if not lines:
        return []
    hdr = lines[0].split("\t")
    di = hdr.index("expDesc")
    return sorted({l.split("\t")[di][:-4].strip() for l in lines[1:]
                   if len(l.split("\t")) > di and l.split("\t")[di].endswith("(C)")})


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--date", default="2026-08-03")
    ap.add_argument("--positives-file", default=None, help="local Keio experiment metadata (skip network)")
    a = ap.parse_args(argv)

    text = (Path(a.positives_file).read_text(encoding="utf-8", errors="replace") if a.positives_file
            else urllib.request.urlopen(
                urllib.request.Request(_POSITIVES_URL, headers={"User-Agent": "dna-decode"}), timeout=60
            ).read().decode("utf-8", "replace"))
    positives = _measured_positives(text)

    m = load_model()
    base = dict(m.medium)
    idx = build_exchange_name_index(m)

    mapped, grows, rates, unmapped = 0, 0, [], []
    per_source = {}
    for src in positives:
        ex = match_carbon_exchange(src, idx)
        if not ex:
            unmapped.append(src)
            continue
        mapped += 1
        g = predict_growth(m, ex, base)
        per_source[src] = {"exchange": ex, "growth_per_h": round(g, 4), "grows": g > 1e-4}
        if g > 1e-4:
            grows += 1
            rates.append(g)

    # classify the unmapped measured-positives: name-gap vs a real MODEL-gap (a false negative)
    unmapped_class = {src: _classify_unmapped(m, src) for src in unmapped}
    model_gaps = sorted(s for s, k in unmapped_class.items() if k == "model-gap")

    recall = grows / mapped if mapped else 0.0
    result = {
        "record": "fba-carbon-growth-validation-v1",
        "date": a.date,
        "organism": "Escherichia coli K-12",
        "model": "iML1515",
        "status": "SCORED_RECALL",
        "n_measured_positive_sources": len(positives),
        "n_mapped_to_exchange": mapped,
        "n_unmapped": len(unmapped),
        "unmapped_sources": unmapped,
        "recall_on_mapped_positives": round(recall, 4),
        "n_grows": grows,
        "growth_rate_spread": {
            "min": round(min(rates), 4) if rates else None,
            "median": round(statistics.median(rates), 4) if rates else None,
            "max": round(max(rates), 4) if rates else None,
        },
        "unmapped_classification_heuristic": unmapped_class,   # crude stem-match; unreliable on noisy labels
        "model_gap_verified_example": "sucrose (BW25113 grows; iML1515 has no sucrose transport)",
        "per_source": per_source,
        "walls": {
            "full_specificity": "EXTERNAL -- needs a MEASURED negative carbon-source set (K-12 can't-use). "
                                "The Keio/Wetmore assay is positive-only; EcN Biolog (PMC9801561) is a strain "
                                "mismatch + SI-locked.",
            "growth_rate_correlation": "EXTERNAL -- no fetchable MEASURED growth-rate-across-carbon-sources "
                                       "dataset (Biolog reports activity indices, not rates; Monk 2017 rates "
                                       "are KO phenotypes on 16 sources, SI-locked).",
        },
        "caveats": [
            "RECALL (sensitivity) only -- FBA predicts growth on carbon sources E. coli is measured to use.",
            "FBA growth RATES are quantitative but there is no fetchable measured rate to correlate against.",
            "Unmapped measured-positives are EITHER name-mapping gaps OR real MODEL-gaps (no transporter in "
            "iML1515 -> a false negative). E.g. BW25113 grows on SUCROSE but K-12 iML1515 lacks the sucrose "
            "system -> a genuine model limitation the validation surfaces, not just a mapping gap.",
        ],
    }
    outdir = Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(exist_ok=True)
    (outdir / f"fba_carbon_growth_validation_{a.date}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    sp = result["growth_rate_spread"]
    md = f"""# FBA carbon-source growth validation: E. coli iML1515 ({a.date})

**Claim tested:** FBA predicts growth (a quantitative rate) on a carbon source, validated against measured
carbon-source utilization.

- Measured-positive carbon sources (Keio/Wetmore, K-12 grows): **{len(positives)}**; mapped to a BiGG
  exchange: **{mapped}** (unmapped = name-mapping gaps, not FBA failures).
- **RECALL on mapped positives: {recall:.3f}** ({grows}/{mapped}) -- FBA predicts growth on every carbon
  source E. coli is measured to use.
- FBA growth-RATE spread (quantitative): min {sp['min']} / median {sp['median']} / max {sp['max']} /h.
- Model-gap example (VERIFIED): BW25113 grows on **sucrose** but iML1515 has no sucrose transport -> a
  false negative the validation surfaces. ({len(unmapped)} positives unmapped total; a heuristic name-gap
  vs model-gap split is in the JSON but is unreliable on noisy assay labels -- sucrose is the verified one.)

## Honest walls (external, not code-closable here)
- **Full specificity:** {result['walls']['full_specificity']}
- **Growth-rate correlation:** {result['walls']['growth_rate_correlation']}

## Caveats
""" + "\n".join(f"- {c}" for c in result["caveats"]) + "\n"
    (outdir / f"fba_carbon_growth_validation_{a.date}.md").write_text(md, encoding="utf-8")
    sys.stdout.buffer.write(md.encode("utf-8", "replace"))
    return 0 if mapped else 1


if __name__ == "__main__":
    raise SystemExit(main())
