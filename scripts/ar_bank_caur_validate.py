"""AR Isolate Bank C. auris fluconazole scorer — the G1-validated fungal cell vs CDC measured MIC.

Independent re-validation of the NON-FROZEN fungal ERG11 cell (`dna_decode/data/fungal_amr` +
`scripts/fungal_erg11_caller`) on the CDC AR Isolate Bank's C. auris isolates. The fungal cell was
G1-validated on a de-confounded 24-genome cohort (clades I Y132F + III F126L); this scores it on a
SEPARATE, provenance-disjoint (0 overlap vs the G1 cohort) CDC C. auris set.

Pipeline per genome (NO AMRFinder — fungi have none): download the GCA assembly ->
`fungal_erg11_caller.observed_substitutions(genome, ERG11_cds_ref)` (BLAST) ->
`fungal_amr.call_from_observed_substitutions("fluconazole", observed)` -> R/S.

Labels: fluconazole MIC (from the AR Bank page; INT is blank because C. auris breakpoints are tentative)
-> `fungal_amr.mic_to_phenotype("fluconazole", mic)` (CDC tentative >=32 -> R). Phenotype side is
independent of the ERG11 genotype prediction, so NOT circular. Only fluconazole has a configured
breakpoint (voriconazole has none -> unlabelable). CURATED_NONFROZEN; frozen surface untouched.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data import fungal_amr as F
from dna_decode.data import refseq
from scripts.amr_portal_score_independent import wilson_ci
from scripts.external_cohort_revalidate import _read_selected
from scripts.fungal_erg11_caller import observed_substitutions
from scripts.independent_cohort_validate import _conf

ERG11_REF = "data/fungal_ref/Cauris_ERG11_cds.fna"
SPEC_FLOOR = 0.85
ORGANISM = "Candida auris"


def predict_fluconazole(gca: str, cache: Path) -> tuple[str, list[str]]:
    """Download the genome, BLAST ERG11, call fluconazole R/S. Returns (prediction, determinants)."""
    # Some BioSample resolutions return a versionless accession (e.g. GCA_003014415); NCBI Datasets
    # needs a version -> try the given accession, then .1/.2 as a fallback.
    candidates = [gca] if "." in gca else [gca + ".1", gca + ".2", gca]
    last_err = "no_candidate"
    for acc in candidates:
        try:
            refseq.download_genome(acc, cache)
            gca = acc
            break
        except Exception as e:  # noqa: BLE001
            last_err = type(e).__name__
    else:
        return "INDETERMINATE", [f"download_failed:{last_err}"]
    fasta = refseq.fasta_path(gca, cache)
    if not Path(fasta).exists():
        return "INDETERMINATE", ["no_genome_fna"]
    obs = observed_substitutions(str(fasta), ERG11_REF)
    if obs is None:
        return "INDETERMINATE", ["blast_unavailable_or_no_hit"]
    call = F.call_from_observed_substitutions("fluconazole", obs)
    return call.prediction, list(call.determinants)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--labels-dir", default="data/raw/ar_bank_caur_extval_fluconazole")
    ap.add_argument("--cohort-manifest", default="wiki/cohort_manifest_external_arbank_caur_2026-07-20.json")
    ap.add_argument("--cache", default="data/raw/ar_bank_caur/genomes")
    ap.add_argument("--run-id", default=f"caur_{_date.today().isoformat()}")
    ap.add_argument("--min-per-class", type=int, default=5)
    a = ap.parse_args()

    labels = _read_selected(Path(a.labels_dir) / "selected_strict.tsv")
    from dna_decode.data.external_cohort_genomes import resolve_cohort_genomes
    from dna_decode.eval.biosample_resolver import BioSampleResolver
    resolver = BioSampleResolver()
    genomes = resolve_cohort_genomes(set(labels), resolver)
    resolver.save_cache()
    free = genomes["free"]
    cache = Path(a.cache)

    records, n_excluded = [], 0
    for bs, rs in sorted(labels.items()):
        gca = free.get(bs)
        if gca is None:
            n_excluded += 1
            continue
        pred, dets = predict_fluconazole(gca, cache)
        records.append({"biosample": bs, "gca": gca, "prediction": pred, "label": rs,
                        "y": 1 if rs == "R" else 0, "determinants": dets})
        print(f"  {bs} {gca}: label={rs} pred={pred} {dets}")

    Path(a.labels_dir).mkdir(parents=True, exist_ok=True)
    scored = [(r["prediction"], r["y"]) for r in records if str(r["prediction"]).upper() in ("R", "S")]
    conf = _conf(scored)
    n_R = conf["tp"] + conf["fn"]
    n_S = conf["tn"] + conf["fp"]
    n_indet = sum(1 for r in records if str(r["prediction"]).upper() not in ("R", "S"))
    powered = n_R >= a.min_per_class and n_S >= a.min_per_class
    endorsed = bool(powered and conf["spec"] is not None and conf["spec"] >= SPEC_FLOOR)
    artifact = {
        "_schema": "ar-bank-caur-validation-v1", "date": _date.today().isoformat(), "run_id": a.run_id,
        "organism": ORGANISM, "drug": "fluconazole", "gene": "ERG11",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "not_in_shipped_surface": True,
        "label_source": "ar_bank_MIC_cdc_tentative_breakpoint (fluconazole >=32 -> R)",
        "independence_tier": ("CDC AR Isolate Bank measured MIC; provenance-disjoint (0 overlap vs the "
                              "fungal G1 tuning cohort); genotype = BLAST ERG11 target-site; phenotype = "
                              "CDC MIC. NOT methodology-independent (rule is the curated ERG11 catalog)."),
        "binary": conf, "n_R": n_R, "n_S": n_S, "n_indeterminate": n_indet, "n_excluded_no_assembly": n_excluded,
        "sens_wilson95": wilson_ci(conf["tp"], n_R), "spec_wilson95": wilson_ci(conf["tn"], n_S),
        "powered": powered, "spec_floor": SPEC_FLOOR,
        "headline": "SCORED_ENDORSED" if endorsed else ("UNDERPOWERED" if not powered else "SCORED_NOT_ENDORSED"),
        "frozen_surface_changed": False,
    }
    (Path(a.labels_dir) / "predictions.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    out = Path(f"wiki/ar_bank_caur_validation_fluconazole_{a.run_id}_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nRESULT C. auris fluconazole: n={len(scored)} ({n_R}R/{n_S}S) acc={conf['acc']} "
          f"sens={conf['sens']} spec={conf['spec']} indet={n_indet} -> {artifact['headline']}")
    print(f"  artifact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
