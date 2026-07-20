"""Power the AR-Bank C. auris fluconazole cell by combining the 8 assembled predictions with the 5
SRA-read-mapped fluconazole-S consensuses, then recomputing sens/spec + the powering verdict.

The assembled arm (`ar_bank_caur_validate`) scored 5R/3S (0 errors) but S landed at 3 < the 5-per-class
powering floor. The 5 remaining fluconazole-S isolates are SRA-reads-only (no downloadable assembly);
`assemble_sra_cohort --method map` produced a targeted ERG11 consensus FASTA per isolate. Here we run the
SAME fungal cell (`fungal_erg11_caller.observed_substitutions` -> `fungal_amr.call_from_observed_substitutions`)
on each consensus, combine with the assembled records, and recompute. Expected: 5 WT-ERG11 -> S, giving
5R + 8S = POWERED at >=5-per-class. Frozen surface untouched (NON-FROZEN fungal cell).

  uv run python -m scripts.ar_bank_caur_powered_finalize
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data import fungal_amr as F
from scripts.amr_portal_score_independent import wilson_ci
from scripts.fungal_erg11_caller import observed_substitutions
from scripts.independent_cohort_validate import _conf

ERG11_REF = "data/fungal_ref/Cauris_ERG11_cds.fna"
SPEC_FLOOR = 0.85
ORGANISM = "Candida auris"


def _read_assembled_tsv(path: Path) -> list[dict]:
    """Parse the assemble_sra_cohort out-tsv (header row -> isolate_id + genome_fasta column)."""
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    idx = {h: i for i, h in enumerate(header)}
    fasta_col = "genome_fasta" if "genome_fasta" in idx else header[-1]
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        rows.append({"isolate_id": parts[idx["isolate_id"]], "genome_fasta": parts[idx[fasta_col]]})
    return rows


def predict_from_consensus(consensus_fna: str) -> tuple[str, list[str], str]:
    obs = observed_substitutions(consensus_fna, ERG11_REF)
    if obs is None:
        return "INDETERMINATE", ["blast_unavailable_or_no_hit"], "NA"
    call = F.call_from_observed_substitutions("fluconazole", obs)
    return call.prediction, list(call.determinants), call.confidence


def _confidence_from_determinants(determinants: list[str]) -> str:
    """Recompute the confidence tier for an already-scored record from its stored ERG11 determinants."""
    obs: dict[str, set[str]] = {}
    for d in determinants:
        if ":" in d:
            gene, sub = d.split(":", 1)
            obs.setdefault(gene, set()).add(sub)
    return F.call_from_observed_substitutions("fluconazole", obs).confidence


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assembled-predictions",
                    default="data/raw/ar_bank_caur_extval_fluconazole/predictions.json")
    ap.add_argument("--sra-tsv", default="data/raw/ar_bank_caur/caur_S_sra.assembled.tsv")
    ap.add_argument("--sra-cohort", default="data/raw/ar_bank_caur/caur_S_sra_cohort.tsv",
                    help="cohort TSV mapping isolate_id -> fluconazole_label (all S)")
    ap.add_argument("--run-id", default=f"caur_powered_{_date.today().isoformat()}")
    ap.add_argument("--min-per-class", type=int, default=5)
    a = ap.parse_args()

    assembled = json.loads(Path(a.assembled_predictions).read_text(encoding="utf-8"))
    assembled_bs = {r["biosample"] for r in assembled}

    # SRA cohort labels (isolate_id == BioSample; all fluconazole-S here)
    labels: dict[str, str] = {}
    cohort_lines = [ln for ln in Path(a.sra_cohort).read_text(encoding="utf-8").splitlines() if ln.strip()]
    chdr = {h: i for i, h in enumerate(cohort_lines[0].split("\t"))}
    for ln in cohort_lines[1:]:
        p = ln.split("\t")
        labels[p[chdr["isolate_id"]]] = p[chdr.get("fluconazole_label", len(p) - 1)].strip().upper() or "S"

    sra_rows = _read_assembled_tsv(Path(a.sra_tsv))
    sra_records = []
    for row in sra_rows:
        iso = row["isolate_id"]
        if iso in assembled_bs:          # guard: never double-count an isolate already scored assembled
            print(f"  SKIP {iso}: already in assembled predictions")
            continue
        fna = row["genome_fasta"]
        if not fna or not Path(fna).exists():
            print(f"  {iso}: consensus FASTA missing ({fna!r}) -> INDETERMINATE")
            pred, dets, conf = "INDETERMINATE", ["no_consensus_fna"], "NA"
        else:
            pred, dets, conf = predict_from_consensus(fna)
        rs = labels.get(iso, "S")
        sra_records.append({"biosample": iso, "gca": f"SRA:{fna}", "prediction": pred, "label": rs,
                            "y": 1 if rs == "R" else 0, "determinants": dets, "confidence": conf,
                            "source": "sra_erg11_map"})
        print(f"  {iso} (SRA ERG11 map): label={rs} pred={pred} conf={conf} {dets}")

    combined = [{**r, "source": r.get("source", "assembly"),
                 "confidence": r.get("confidence") or _confidence_from_determinants(r.get("determinants", []))}
                for r in assembled] + sra_records
    # Lineage-marker disclosure: R calls driven ONLY by the clade-background haplotype (non-discriminative).
    low_conf = [r for r in combined if r.get("confidence") == "LOW_LINEAGE_ONLY"]
    lineage_fp = [r for r in low_conf if r["label"] == "S"]
    lineage_tp = [r for r in low_conf if r["label"] == "R"]

    # HIGH-confidence (mechanism-attributable) subset: exclude the clade-background-only LOW calls.
    hi = [r for r in combined if r.get("confidence") != "LOW_LINEAGE_ONLY"
          and str(r["prediction"]).upper() in ("R", "S")]
    hi_scored = [(r["prediction"], r["y"]) for r in hi]
    hi_conf = _conf(hi_scored)
    hi_R = hi_conf["tp"] + hi_conf["fn"]
    hi_S = hi_conf["tn"] + hi_conf["fp"]
    scored = [(r["prediction"], r["y"]) for r in combined if str(r["prediction"]).upper() in ("R", "S")]
    conf = _conf(scored)
    n_R = conf["tp"] + conf["fn"]
    n_S = conf["tn"] + conf["fp"]
    n_indet = sum(1 for r in combined if str(r["prediction"]).upper() not in ("R", "S"))
    powered = n_R >= a.min_per_class and n_S >= a.min_per_class
    endorsed = bool(powered and conf["spec"] is not None and conf["spec"] >= SPEC_FLOOR)

    artifact = {
        "_schema": "ar-bank-caur-powered-validation-v1", "date": _date.today().isoformat(), "run_id": a.run_id,
        "organism": ORGANISM, "drug": "fluconazole", "gene": "ERG11",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "not_in_shipped_surface": True,
        "label_source": "ar_bank_MIC_cdc_tentative_breakpoint (fluconazole >=32 -> R)",
        "genotype_source": ("BLAST ERG11 target-site; 8 isolates from downloadable assembly + "
                            "5 fluconazole-S from SRA-read-mapped ERG11 consensus (assemble_sra_cohort --method map)"),
        "independence_tier": ("CDC AR Isolate Bank measured MIC; provenance-disjoint (0 overlap vs the fungal "
                              "G1 tuning cohort); NOT methodology-independent (rule is the curated ERG11 catalog)."),
        "binary": conf, "n_R": n_R, "n_S": n_S, "n_indeterminate": n_indet,
        "n_assembled": len(assembled), "n_sra_mapped": len(sra_records),
        "sens_wilson95": wilson_ci(conf["tp"], n_R), "spec_wilson95": wilson_ci(conf["tn"], n_S),
        "powered": powered, "spec_floor": SPEC_FLOOR,
        "headline": "SCORED_ENDORSED" if endorsed else ("UNDERPOWERED" if not powered else "SCORED_NOT_ENDORSED"),
        "lineage_marker_disclosure": {
            "clade_iv_haplotype": ["ERG11:K177R", "ERG11:N335S", "ERG11:E343D"],
            "n_low_confidence_R": len(low_conf),
            "low_confidence_R_that_are_FP": [r["biosample"] for r in lineage_fp],
            "low_confidence_R_that_are_TP": [r["biosample"] for r in lineage_tp],
            "finding": ("clade IV ERG11 haplotype is non-discriminative (identical genotype in >=1 R + >=1 S "
                        "isolate) -> +TP and +FP wash -> lineage marker, not causal. R calls preserved "
                        "(sensitivity) but flagged LOW_LINEAGE_ONLY. Fungal analogue of the QRDR-vs-lineage "
                        "confound."),
        },
        "high_confidence_subset": {
            "note": ("mechanism-attributable calls only (Y132F/K143R/F126L causal markers); the clade IV "
                     "haplotype-only isolates abstain to LOW_LINEAGE_ONLY and are excluded here."),
            "binary": hi_conf, "n_R": hi_R, "n_S": hi_S,
            "sens_wilson95": wilson_ci(hi_conf["tp"], hi_R), "spec_wilson95": wilson_ci(hi_conf["tn"], hi_S),
        },
        "frozen_surface_changed": False,
    }
    outdir = Path("data/raw/ar_bank_caur_extval_fluconazole")
    (outdir / "predictions_powered.json").write_text(json.dumps(combined, indent=2), encoding="utf-8")
    out = Path(f"wiki/ar_bank_caur_powered_validation_fluconazole_{a.run_id}_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nRESULT C. auris fluconazole (POWERED): n={len(scored)} ({n_R}R/{n_S}S) acc={conf['acc']} "
          f"sens={conf['sens']} spec={conf['spec']} indet={n_indet} -> {artifact['headline']}")
    print(f"  HIGH-confidence (mechanism-attributable) subset: n={hi_conf['n_scored']} ({hi_R}R/{hi_S}S) "
          f"acc={hi_conf['acc']} sens={hi_conf['sens']} spec={hi_conf['spec']} "
          f"(excluded {len(low_conf)} clade-IV-haplotype-only: {len(lineage_tp)}R-TP + {len(lineage_fp)}S-FP)")
    print(f"  artifact: {out}")
    return 0 if powered else 1


if __name__ == "__main__":
    raise SystemExit(main())
