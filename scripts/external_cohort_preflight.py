"""Gate-0 preflight for external-cohort re-validation — the wave-0 go/no-go.

Before any genome download or scoring, this answers three questions for a candidate
external project accession (e.g. Oxford PRJNA604975 / Spain PROBAC PRJEB62601):

  (a) ASSEMBLY-AVAILABILITY — how many cohort BioSamples have a downloadable GCA/GCF
      (FREE pilot) vs reads-only (ASSEMBLY-REQUIRED, excluded from the free pilot).
  (b) BIOSAMPLE-LEVEL LEAKAGE — resolve the decoder's TUNING accessions (every
      data/raw/*/selected.tsv + data/processed/*.parquet, via cohort_manifest) to
      BioSamples and intersect with the cohort BioSamples. This closes the
      accession-string leakage gate's blind spot (a GCA and an ERR run that are the
      same physical isolate). FAIL-CLOSED if ANY overlap, OR >5% of tuning accessions
      are unresolved-to-BioSample, OR Entrez/ENA disagree on one (disagreement counts
      as unresolved).
  (c) MIC-OPENNESS — a manual human-confirmed flag (the one fact code can't settle):
      is the per-isolate MIC table actually downloadable, or MTA-gated.

Emits `wiki/external_preflight_<cohort>_<date>.json` with an overall PASS/FAIL verdict.
The pure verdict functions are unit-tested offline; only `preflight()` touches network.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.eval.biosample_resolver import BioSampleResolver
from dna_decode.eval.cohort_manifest import build_manifest, prior_accessions

UNRESOLVED_THRESHOLD = 0.05  # >5% of tuning accessions unresolved -> FAIL-CLOSED


# --------------------------------------------------------------------------- #
# Pure verdict logic (no network) — unit-tested directly.
# --------------------------------------------------------------------------- #
def classify_availability(bs_to_assemblies: dict[str, list[str]]) -> dict:
    """Partition cohort BioSamples into FREE (>=1 assembly) vs ASSEMBLY-REQUIRED."""
    free = sorted(bs for bs, gcas in bs_to_assemblies.items() if gcas)
    assembly_required = sorted(bs for bs, gcas in bs_to_assemblies.items() if not gcas)
    return {
        "n_biosamples": len(bs_to_assemblies),
        "n_free": len(free),
        "n_assembly_required": len(assembly_required),
        "free": free,
        "assembly_required": assembly_required,
    }


def leakage_verdict(tuning_acc_to_bs: dict[str, str | None], cohort_biosamples: set[str],
                    unresolved_threshold: float = UNRESOLVED_THRESHOLD) -> dict:
    """Decide the BioSample leakage verdict from precomputed resolutions.

    `tuning_acc_to_bs` maps each tuning accession -> its resolved BioSample, or None
    (unresolved: not found OR Entrez/ENA disagreement). Overlap of ANY resolved
    tuning BioSample with the cohort set, OR an unresolved fraction above the
    threshold, FAILS the gate (cannot prove disjointness).
    """
    n_total = len(tuning_acc_to_bs)
    resolved = {a: bs for a, bs in tuning_acc_to_bs.items() if bs}
    unresolved = [a for a, bs in tuning_acc_to_bs.items() if not bs]
    overlap = sorted({bs for bs in resolved.values() if bs in cohort_biosamples})
    unresolved_frac = (len(unresolved) / n_total) if n_total else 0.0
    fail_overlap = bool(overlap)
    fail_unresolved = unresolved_frac > unresolved_threshold
    return {
        "n_tuning": n_total,
        "n_resolved": len(resolved),
        "n_unresolved": len(unresolved),
        "unresolved_fraction": round(unresolved_frac, 4),
        "unresolved_threshold": unresolved_threshold,
        "overlap_biosamples": overlap,
        "fail_overlap": fail_overlap,
        "fail_unresolved": fail_unresolved,
        "passed": not (fail_overlap or fail_unresolved),
    }


def overall_verdict(leakage: dict, availability: dict, mic_open: bool | None) -> dict:
    """Combine the three signals into a PASS/FAIL with reasons."""
    reasons: list[str] = []
    if leakage["fail_overlap"]:
        reasons.append(f"BioSample overlap with tuning: {leakage['overlap_biosamples']}")
    if leakage["fail_unresolved"]:
        reasons.append(f"unresolved tuning fraction {leakage['unresolved_fraction']} "
                       f"> {leakage['unresolved_threshold']}")
    if mic_open is False:
        reasons.append("MIC table is not openly downloadable (MTA-gated)")
    if mic_open is None:
        reasons.append("MIC-openness not confirmed (pass --mic-open / --mic-gated)")
    if availability["n_free"] == 0:
        reasons.append("no FREE (assembly-resolvable) BioSamples — free pilot N is zero")
    verdict = "PASS" if not reasons else "FAIL"
    return {"verdict": verdict, "reasons": reasons}


# --------------------------------------------------------------------------- #
# Orchestrator (network) — composes the resolver + manifest + verdicts.
# --------------------------------------------------------------------------- #
def preflight(project: str, cohort_name: str, *, mic_open: bool | None,
              resolver: BioSampleResolver | None = None,
              wiki_dir: str | Path = "wiki", write: bool = True) -> dict:
    resolver = resolver or BioSampleResolver()

    # cohort BioSamples (ENA read_run for the project)
    runs = resolver.runs_for_project(project)
    cohort_biosamples = sorted({bs for _run, bs in runs})

    # (a) assembly availability per cohort BioSample
    bs_to_assemblies = {bs: resolver.biosample_to_assemblies(bs) for bs in cohort_biosamples}
    availability = classify_availability(bs_to_assemblies)

    # (b) BioSample-level leakage vs the decoder's tuning accessions
    manifest = build_manifest()
    tuning_accs = sorted(prior_accessions(manifest, exclude_cohort=cohort_name))
    tuning_acc_to_bs = {a: resolver.assembly_to_biosample(a) for a in tuning_accs}
    leakage = leakage_verdict(tuning_acc_to_bs, set(cohort_biosamples))

    resolver.save_cache()

    verdict = overall_verdict(leakage, availability, mic_open)
    artifact = {
        "_schema": "external-preflight-v1",
        "date": _date.today().isoformat(),
        "project": project,
        "cohort_name": cohort_name,
        "verdict": verdict["verdict"],
        "reasons": verdict["reasons"],
        "n_cohort_biosamples": len(cohort_biosamples),
        "assembly_availability": availability,
        "leakage": leakage,
        "mic_open": mic_open,
        "manifest_complete": not manifest.incomplete,
        "manifest_n_cohorts": len(manifest.cohorts),
    }
    if write:
        out = Path(wiki_dir) / f"external_preflight_{cohort_name}_{_date.today().isoformat()}.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
        artifact["_artifact_path"] = str(out)
    return artifact


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project", required=True, help="ENA/NCBI project accession (e.g. PRJNA604975)")
    ap.add_argument("--cohort-name", required=True, help="cohort slug (used in output filenames)")
    mic = ap.add_mutually_exclusive_group()
    mic.add_argument("--mic-open", dest="mic_open", action="store_true", default=None,
                     help="per-isolate MIC table confirmed openly downloadable")
    mic.add_argument("--mic-gated", dest="mic_open", action="store_false",
                     help="per-isolate MIC table is MTA-gated")
    a = ap.parse_args()
    art = preflight(a.project, a.cohort_name, mic_open=a.mic_open)
    print(f"\nPREFLIGHT {art['verdict']} — {art['cohort_name']} ({a.project})")
    print(f"  cohort BioSamples: {art['n_cohort_biosamples']} "
          f"(FREE {art['assembly_availability']['n_free']} / "
          f"ASSEMBLY-REQUIRED {art['assembly_availability']['n_assembly_required']})")
    print(f"  leakage: overlap={art['leakage']['overlap_biosamples']} "
          f"unresolved={art['leakage']['unresolved_fraction']}")
    for r in art["reasons"]:
        print(f"  - {r}")
    return 0 if art["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
