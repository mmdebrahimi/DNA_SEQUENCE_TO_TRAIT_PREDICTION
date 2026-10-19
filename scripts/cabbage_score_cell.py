"""Score the FROZEN decoder on a CABBAGE cell, using the portal's PRE-COMPUTED AMRFinderPlus.

This is the first scored cell on the reopened public-label AMR track (user authority decision
2026-10-09). It is possible at this scale only because the EMBL-EBI AMR Portal ships AMRFinderPlus
output per assembly at ~1 KB, where running it locally costs ~95 s/strain through Docker.

THREE RAILS, each from a measured failure in this repo:

1. HEADER DIALECT. The portal writes `Element_symbol`; `eval/amr_rules` reads `Element symbol`. Fed raw,
   two of three drugs score a confident SUSCEPTIBLE -- including `gyrA_T86I`, the determinant the
   project's own Campylobacter cell rests on. Every file goes through
   `dna_decode.data.amr_portal.prepare`, which normalizes the header and REFUSES an unrecognised schema.

2. LEAKAGE, FAIL-CLOSED. 1,011 of the repo's 1,509 manifest accessions appear in CABBAGE. Every
   accession in `cohort_manifest.prior_accessions` is EXCLUDED (version-insensitively), and an
   INCOMPLETE manifest aborts rather than scoring -- an incomplete manifest cannot prove disjointness.

3. LABELS VIA THE REPO'S OWN TIERS. Phenotypes come from `external_mic_labels.tier_for_isolate` on the
   operator-aware MIC, not from the categorical column. Only STRICT tiers are scored; everything else
   (BORDERLINE / AMBIGUOUS / CENSORED_EXCLUDED) is EXCLUDED and counted.

Run:  uv run python scripts/cabbage_score_cell.py --smoke
      uv run python scripts/cabbage_score_cell.py --max-isolates 1500
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import random
import sys
import time
import urllib.error
import urllib.request
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

PHENO_PARQUET = "D:/dna_decode_cache/amr_portal/phenotype.parquet"
CACHE = Path("D:/dna_decode_cache/amr_portal/amrfinder")
DELAY_S = 0.12


def fetch_amrfinder(accession: str) -> Path | None:
    """Cached fetch of one pre-computed AMRFinderPlus file. None on any failure (recorded, not silent).

    Integrity: the gzip is FULLY decompressed before the cache entry is kept -- a truncated gzip
    decompresses its opening blocks perfectly and raises only at the end, so a head read would keep
    exactly the files that later fail.
    """
    from dna_decode.data.amr_portal import genome_url

    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / f"{accession}.tsv"
    if dest.exists():
        return dest
    try:
        url = genome_url(accession)
    except ValueError:
        return None
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "dna_decode-research/0.1"})
        with urllib.request.urlopen(req, timeout=60) as fh:
            raw = fh.read()
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError):
        return None
    except Exception:                                 # noqa: BLE001
        return None
    try:
        text = gzip.decompress(raw).decode("utf-8", errors="replace")
    except Exception:                                 # noqa: BLE001 - truncated/corrupt
        return None
    tmp = dest.with_suffix(".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(dest)
    return dest


def build_labels(species: str, drug: str):
    """{accession: 'R'|'S'} from the repo's STRICT MIC tiers, plus the exclusion census."""
    import pandas as pd

    from dna_decode.data.external_mic_labels import (EXCLUDED_TIERS, STRICT_TIERS,
                                                     tier_for_isolate)
    df = pd.read_parquet(PHENO_PARQUET)
    s = df[(df.species == species) & (df.antibiotic_name.str.lower() == drug.lower())]
    s = s[s.assembly_ID.notna() & s.measurement.notna()].copy()
    s["tok"] = (s.measurement_sign.fillna("=").replace({"==": "="}).astype(str)
                + s.measurement.astype(str))
    gm = s.groupby("assembly_ID")["tok"].apply(list)
    gc = s.groupby("assembly_ID")["resistance_phenotype"].apply(
        lambda x: set(v for v in x.dropna()))
    srcs = s.drop_duplicates("assembly_ID").set_index("assembly_ID")["database"].to_dict()
    yrs = s.drop_duplicates("assembly_ID").set_index("assembly_ID")["collection_year"].to_dict()
    ctry = s.drop_duplicates("assembly_ID").set_index("assembly_ID")["country"].to_dict()

    labels, tiers = {}, collections.Counter()
    for acc in gm.index:
        t = tier_for_isolate(gm[acc], gc.get(acc, set()), drug)
        tiers[t] += 1
        lab = STRICT_TIERS.get(t)
        if lab:
            labels[acc] = lab
    return labels, {"tier_census": dict(tiers),
                    "n_with_asm_and_mic": int(len(gm)),
                    "n_strict_labelled": len(labels),
                    "n_excluded_tiers": int(sum(n for t, n in tiers.items()
                                                if t in EXCLUDED_TIERS))}, srcs, yrs, ctry


def main(argv=None) -> int:
    from dna_decode.data.amr_portal import PortalHeaderIncompatible, prepare
    from dna_decode.eval.amr_rules import call_resistance
    from dna_decode.eval.clonality import wilson_ci
    from dna_decode.eval.cohort_manifest import build_manifest, prior_accessions

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", default="Campylobacter jejuni")
    ap.add_argument("--drug", default="ciprofloxacin")
    ap.add_argument("--organism", default="Campylobacter", help="AMRFinder -O / registry organism")
    ap.add_argument("--max-isolates", type=int, default=1500,
                    help="0 = no cap. A cap is applied as a SEEDED RANDOM sample, never an "
                         "accession-ordered prefix: the prefix is measurably easier (Campylobacter "
                         "cipro first-1,444 sens 0.9971 vs full-4,265 0.9794).")
    ap.add_argument("--sample-seed", type=int, default=20261010,
                    help="seed for the random subset when --max-isolates caps the cohort")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--allow-incomplete-manifest", action="store_true")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    if a.smoke:
        a.max_isolates = 40

    # --- leakage gate, FAIL CLOSED -------------------------------------------------------------
    man = build_manifest()
    if getattr(man, "incomplete", False) and not a.allow_incomplete_manifest:
        print("INCOMPLETE_MANIFEST -- refusing to score; an incomplete manifest cannot prove "
              "disjointness. Re-run with --allow-incomplete-manifest to stamp a degraded artifact.")
        return 2
    prior = {x.split(".")[0] for x in prior_accessions(man, exclude_cohort="__none__")}
    print(f"leakage gate: {len(prior):,} prior accessions (manifest_incomplete="
          f"{getattr(man, 'incomplete', None)})", flush=True)

    print(f"building labels: {a.species} x {a.drug} ...", flush=True)
    labels, census, srcs, yrs, ctry = build_labels(a.species, a.drug)
    print(f"  {census}", flush=True)

    kept = {acc: lab for acc, lab in labels.items() if acc.split(".")[0] not in prior}
    n_leak = len(labels) - len(kept)
    print(f"  excluded for leakage: {n_leak}  -> scoreable: {len(kept):,}", flush=True)

    # DETERMINISTIC RANDOM subset, not an accession-ordered prefix. Measured, not assumed: on
    # Campylobacter x ciprofloxacin the first 1,444 accessions scored sens 0.9971 while the full
    # 4,265 scored 0.9794, so the prefix is EASIER and a capped prefix reports an optimistic number.
    # Accession order tracks deposition batch, which tracks study, so it is nowhere near random.
    # A seeded shuffle keeps the run reproducible while making the cap unbiased.
    order = sorted(kept)
    if a.max_isolates and len(order) > a.max_isolates:
        rng = random.Random(a.sample_seed)
        rng.shuffle(order)
        order = sorted(order[: a.max_isolates])
        print(f"  capped to {a.max_isolates:,} of {len(kept):,} by SEEDED RANDOM sample "
              f"(seed={a.sample_seed}); an accession-ordered prefix is measurably easier",
              flush=True)
    print(f"scoring {len(order):,} isolates (fetch + normalize + call_resistance) ...", flush=True)

    conf = collections.Counter()
    fetch_fail, dialect_fail, call_fail = [], [], []
    scored_src, scored_year, scored_ctry = collections.Counter(), [], collections.Counter()
    t0 = time.time()
    for k, acc in enumerate(order, 1):
        cached = (CACHE / f"{acc}.tsv").exists()
        p = fetch_amrfinder(acc)
        if not cached:
            time.sleep(DELAY_S)
        if p is None:
            fetch_fail.append(acc)
            continue
        try:
            norm = prepare(p, CACHE / "norm" / f"{acc}.tsv")
        except PortalHeaderIncompatible:
            dialect_fail.append(acc)
            continue
        try:
            res = call_resistance(norm, a.drug, organism=a.organism)
        except Exception as e:                        # noqa: BLE001
            call_fail.append(f"{acc}: {type(e).__name__}")
            continue
        pred = res.get("prediction")
        truth = kept[acc]
        conf[(truth, pred)] += 1
        scored_src[str(srcs.get(acc, "?"))] += 1
        if yrs.get(acc) == yrs.get(acc):
            scored_year.append(yrs.get(acc))
        scored_ctry[str(ctry.get(acc, "?"))] += 1
        if k % 200 == 0:
            el = time.time() - t0
            print(f"  {k}/{len(order)}  scored={sum(conf.values())}  "
                  f"fetch_fail={len(fetch_fail)}  {el:.0f}s", flush=True)

    tp = conf[("R", "R")]; fn = conf[("R", "S")]
    tn = conf[("S", "S")]; fp = conf[("S", "R")]
    n = tp + fn + tn + fp
    sens = tp / (tp + fn) if (tp + fn) else None
    spec = tn / (tn + fp) if (tn + fp) else None
    acc_ = (tp + tn) / n if n else None
    top_src = scored_src.most_common(1)[0] if scored_src else ("?", 0)
    top_ctry = scored_ctry.most_common(1)[0] if scored_ctry else ("?", 0)
    db_share = top_src[1] / n if n else None
    ct_share = top_ctry[1] / n if n else None
    # The repo's own diversity bar, applied to the repo's own new cohort.
    DIVERSITY_SHARE_BAR = 0.60          # scripts/source_diverse_validate.py
    DIVERSITY_MIN_SOURCES = 5
    gate_axes = {
        "database_label": {"n_distinct": len(scored_src), "largest_share": db_share},
        "country": {"n_distinct": len(scored_ctry), "largest_share": ct_share},
    }
    for ax in gate_axes.values():
        ax["passes"] = bool(ax["largest_share"] is not None
                            and ax["n_distinct"] >= DIVERSITY_MIN_SOURCES
                            and ax["largest_share"] <= DIVERSITY_SHARE_BAR)

    out = {
        "schema": "cabbage-scored-cell-v1", "generated": date.today().isoformat(),
        "cell": {"species": a.species, "drug": a.drug, "amrfinder_organism": a.organism},
        "smoke": bool(a.smoke),
        "label_source": ("CABBAGE / EMBL-EBI AMR Portal release 2026-07, scored through this repo's "
                         "external_mic_labels.tier_for_isolate on operator-aware MIC (STRICT tiers only)"),
        "genotype_source": ("portal PRE-COMPUTED AMRFinderPlus, header-normalized via "
                            "dna_decode.data.amr_portal.prepare"),
        "rule": "FROZEN eval/amr_rules.call_resistance (byte-unchanged)",
        "label_census": census,
        "leakage": {"n_prior_accessions": len(prior), "n_excluded_for_leakage": n_leak,
                    "manifest_incomplete": bool(getattr(man, "incomplete", False))},
        "attempted": len(order),
        "subset": {"n_scoreable": len(kept), "max_isolates": a.max_isolates,
                   "selection": ("seeded_random_sample" if a.max_isolates
                                 and len(kept) > a.max_isolates else "all_scoreable"),
                   "sample_seed": a.sample_seed,
                   "why_not_prefix": ("an accession-ordered prefix is EASIER -- measured on "
                                      "Campylobacter x ciprofloxacin, first 1,444 sens 0.9971 vs "
                                      "full 4,265 sens 0.9794 -- because accession order tracks "
                                      "deposition batch and therefore study")},
        "confusion": {"tp": tp, "fn": fn, "tn": tn, "fp": fp, "n_scored": n},
        "metrics": {
            "accuracy": round(acc_, 4) if acc_ is not None else None,
            "sensitivity": round(sens, 4) if sens is not None else None,
            "specificity": round(spec, 4) if spec is not None else None,
            "sens_wilson_ci": [round(x, 4) for x in wilson_ci(tp, tp + fn)] if (tp + fn) else None,
            "spec_wilson_ci": [round(x, 4) for x in wilson_ci(tn, tn + fp)] if (tn + fp) else None,
        },
        "failures": {"n_fetch_failed": len(fetch_fail), "n_dialect_refused": len(dialect_fail),
                     "n_call_failed": len(call_fail),
                     "fetch_failed_sample": fetch_fail[:10], "call_failed_sample": call_fail[:5]},
        "source_concentration_of_scored_set": {
            "bar": {"max_largest_share": DIVERSITY_SHARE_BAR,
                    "min_distinct_sources": DIVERSITY_MIN_SOURCES,
                    "provenance": "scripts/source_diverse_validate.py (this repo's own bar)"},
            "n_distinct_database_labels": len(scored_src),
            "largest_database_share": round(db_share, 4) if db_share is not None else None,
            "largest_database_label": top_src[0],
            "n_distinct_countries": len(scored_ctry),
            "largest_country_share": round(ct_share, 4) if ct_share is not None else None,
            "largest_country": top_ctry[0],
            "collection_year_range": ([int(min(scored_year)), int(max(scored_year))]
                                      if scored_year else None),
            "n_distinct_collection_years": len(set(scored_year)),
            "axes": {k: {kk: (round(vv, 4) if isinstance(vv, float) else vv)
                         for kk, vv in v.items()} for k, v in gate_axes.items()},
            "verdict": ("source_diverse" if all(ax["passes"] for ax in gate_axes.values())
                        else "source_concentrated"),
            "reading": ("A cohort that fails this bar cannot reveal a determinant-family blind spot it "
                        "does not contain -- the mechanism that made the gentamicin rmt gap "
                        "structurally invisible. It does NOT make the metrics wrong; it bounds what "
                        "they could ever have detected."),
        },
        "honest_limits": [
            "NOT CLONALITY-CORRECTED. These are isolate-level counts; Campylobacter surveillance "
            "collections are clonally structured, so read the metrics as isolate-level. The lineage "
            "layer (G8) is the open gate on this substrate and is NOT closed by this run.",
            "The `database` column is an overlapping AGGREGATOR label, not a BioProject, so the "
            "source-concentration figures here are NOT this repo's per-BioProject diversity bar.",
            "Campylobacter cipro resistance is overwhelmingly one well-catalogued determinant "
            "(gyrA_T86I), so a high score is MECHANISTICALLY EXPECTED and is not evidence of broad "
            "catalog coverage.",
            "A fetch failure means the portal's genomes/ tree lacks that assembly; those isolates are "
            "EXCLUDED and counted, never scored as wrong. Tree coverage is uncensused (see "
            "constraint `cabbage-genomes-tree-coverage-uncensused`).",
            "Leakage exclusion is accession-string based (version-insensitive). It does not resolve "
            "BioSample identity, so an isolate deposited under two different assemblies could survive.",
            "GEOGRAPHICALLY CONCENTRATED: see source_concentration_of_scored_set.verdict. Do NOT quote "
            "this cell as source-diverse. What it DOES add is INDEPENDENCE: the provenance-disjoint "
            "Campylobacter cipro cell is PRJNA560409 (Univ. Melbourne, N=40) and the prospective cell "
            "is PRJNA563067 (PulseNet Canada, N=426, 100% one BioProject). Each fails the diversity "
            "bar alone; that they are three MUTUALLY INDEPENDENT national collections is the part "
            "worth keeping.",
        ],
    }
    # POWERING GATE. Mirrors MIN_PER_CLASS in independent_cohort_validate / the external-cohort arm
    # rather than inventing a bar. The n==0 and class==0 branches REFUSE rather than write: an
    # artifact of null metrics is indistinguishable from an artifact of a genuinely empty cohort, and
    # the empty reading is the reassuring one. A clean-looking verdict computed over silently-zero
    # input is exactly what this guards.
    MIN_PER_CLASS = 10
    n_R, n_S = tp + fn, tn + fp
    if n == 0:
        print(f"REFUSING to write an artifact: ZERO isolates scored of {len(order):,} attempted "
              f"({len(fetch_fail)} fetch-failed, {len(dialect_fail)} dialect-refused, "
              f"{len(call_fail)} call-failed). Nothing was measured.")
        return 3
    if n_R == 0 or n_S == 0:
        print(f"REFUSING to write an artifact: one class is EMPTY (R={n_R}, S={n_S}). A sensitivity "
              "or specificity over an absent class is not a measurement.")
        return 3
    out["powering"] = {
        "min_per_class": MIN_PER_CLASS, "n_R": n_R, "n_S": n_S,
        "degraded": bool(n_R < MIN_PER_CLASS or n_S < MIN_PER_CLASS),
        "note": ("below MIN_PER_CLASS the corresponding metric is reported but should not be quoted; "
                 "the bar is imported from this repo's external-cohort arm, not asserted here"),
    }
    if out["powering"]["degraded"]:
        print(f"WARNING degraded powering: R={n_R} S={n_S} against MIN_PER_CLASS={MIN_PER_CLASS}")

    if a.out is None:
        tag = f"{a.species.split()[0].lower()}_{a.drug.lower()}"
        a.out = REPO / "wiki" / f"cabbage_scored_cell_{tag}_{date.today().isoformat()}.json"
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print()
    print(f"SCORED {n:,} isolates  tp={tp} fn={fn} tn={tn} fp={fp}")
    print(f"  accuracy={acc_}  sens={sens}  spec={spec}")
    print(f"  fetch_failed={len(fetch_fail)}  dialect_refused={len(dialect_fail)}  "
          f"call_failed={len(call_fail)}")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
