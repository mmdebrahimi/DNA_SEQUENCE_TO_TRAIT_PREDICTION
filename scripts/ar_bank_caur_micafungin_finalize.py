"""Finalize the AR-Bank C. auris micafungin (FKS1 echinocandin) cell: combine the assembled micafungin-S
genomes with the SRA-read-mapped micafungin-R FKS1 consensuses, score, and DISCLOSE uncatalogued FKS1
variants (the ERG11-lineage-marker lesson applied to the echinocandin target).

Why the disclosure matters: the first R consensus (SAMN38094230, MIC>8) carries FKS1 W691L -- a clean,
full-length, non-canonical variant (S639/F635/R1354 all WT) that the catalog does NOT list. Silently
scoring it S would hide that the target-site scan FOUND a variant it doesn't recognize. We therefore record
the RAW FKS1 substitutions (catalogued + uncatalogued) for BOTH R and S isolates and flag
`uncatalogued_fks1_variants` -- a variant seen ONLY in R is a candidate blind-spot determinant; one seen in
BOTH R and S is a lineage/background marker (non-causal), exactly the ERG11 clade IV pattern. We do NOT add
any uncatalogued variant to the catalog (1-isolate, no independent causal evidence = the over-call the
ERG11 lesson warns against).

  uv run python -m scripts.ar_bank_caur_micafungin_finalize
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

FKS1_REF = "data/fungal_ref/Cauris_FKS1_cds.fna"
SPEC_FLOOR = 0.85
ORGANISM = "Candida auris"
CATALOGUED = F.resistance_mutations_for("micafungin")["FKS1"]   # {S639F, S639P, S639Y, F635del, R1354H}


def _read_assembled_tsv(path: Path) -> list[dict]:
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    idx = {h: i for i, h in enumerate(header)}
    fasta_col = "genome_fasta" if "genome_fasta" in idx else header[-1]
    return [{"isolate_id": p[idx["isolate_id"]], "genome_fasta": p[idx[fasta_col]]}
            for p in (ln.split("\t") for ln in lines[1:])]


def _resolve_genome(cache: Path, gca: str) -> Path | None:
    """Genome path with versionless-accession fallback (predictions.json may store GCA_x w/o .1/.2)."""
    for acc in ([gca] if "." in gca else [gca + ".1", gca + ".2", gca]):
        p = cache / acc / "genome.fna"
        if p.exists():
            return p
    return None


def score_one(fasta: str) -> dict:
    """Raw FKS1 obs + micafungin call + catalogued/uncatalogued split. INDETERMINATE if BLAST fails."""
    obs = observed_substitutions(fasta, FKS1_REF, gene="FKS1")
    if obs is None:
        return {"prediction": "INDETERMINATE", "fks1_subs": [], "determinants": [], "uncatalogued": []}
    subs = sorted(obs.get("FKS1", set()))
    call = F.call_from_observed_substitutions("micafungin", obs)
    uncat = [s for s in subs if s not in CATALOGUED]
    return {"prediction": call.prediction, "fks1_subs": subs,
            "determinants": list(call.determinants), "uncatalogued": uncat}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--labels-dir", default="data/raw/ar_bank_enterococcus_faecium_extval_micafungin")
    ap.add_argument("--mica-labels", default="data/raw/ar_bank_caur_extval_micafungin/selected_strict.tsv")
    ap.add_argument("--assembled-preds", default="data/raw/ar_bank_caur_extval_micafungin/predictions.json")
    ap.add_argument("--genome-cache", default="data/raw/ar_bank_caur/genomes")
    ap.add_argument("--sra-tsv", default="data/raw/ar_bank_caur/caur_R_mica_sra.assembled.tsv")
    ap.add_argument("--sra-assemblies", default="D:/dna_decode_cache/caur_fks1_sra/assemblies",
                    help="dir of per-isolate FKS1 consensus FASTAs (fallback if --sra-tsv absent)")
    ap.add_argument("--run-id", default=f"caur_mica_{_date.today().isoformat()}")
    ap.add_argument("--min-per-class", type=int, default=5)
    a = ap.parse_args()

    from scripts.external_cohort_revalidate import _read_selected
    labels = _read_selected(Path(a.mica_labels))          # {biosample: R/S}

    records = []
    # S side: re-derive raw FKS1 from the already-downloaded assembled genomes (predictions.json has GCAs)
    assembled = json.loads(Path(a.assembled_preds).read_text(encoding="utf-8")) \
        if Path(a.assembled_preds).exists() else []
    for r in assembled:
        bs, gca = r["biosample"], r["gca"]
        fasta = _resolve_genome(Path(a.genome_cache), gca)
        sc = score_one(str(fasta)) if fasta else \
            {"prediction": "INDETERMINATE", "fks1_subs": [], "determinants": [], "uncatalogued": ["no_genome"]}
        records.append({"biosample": bs, "source": "assembly", "label": labels.get(bs, r.get("label")),
                        **sc})
    seen = {r["biosample"] for r in records}
    # R side: SRA FKS1 consensuses (from the out-tsv if present, else the assemblies dir)
    sra_rows = _read_assembled_tsv(Path(a.sra_tsv)) if Path(a.sra_tsv).exists() else \
        [{"isolate_id": p.stem, "genome_fasta": str(p)}
         for p in sorted(Path(a.sra_assemblies).glob("*.fna"))]
    for row in sra_rows:
        bs = row["isolate_id"]
        if bs in seen:
            continue
        fasta = row["genome_fasta"]
        sc = score_one(fasta) if Path(fasta).exists() else \
            {"prediction": "INDETERMINATE", "fks1_subs": [], "determinants": [], "uncatalogued": ["no_consensus"]}
        records.append({"biosample": bs, "source": "sra_fks1_map", "label": labels.get(bs, "R"), **sc})

    for r in records:
        r["y"] = 1 if r["label"] == "R" else 0
        print(f"  {r['biosample']} ({r['source']}): label={r['label']} pred={r['prediction']} "
              f"fks1={r['fks1_subs']} uncat={r['uncatalogued']}")

    scored = [(r["prediction"], r["y"]) for r in records if str(r["prediction"]).upper() in ("R", "S")]
    conf = _conf(scored)
    n_R = conf["tp"] + conf["fn"]
    n_S = conf["tn"] + conf["fp"]
    powered = n_R >= a.min_per_class and n_S >= a.min_per_class
    endorsed = bool(powered and conf["spec"] is not None and conf["spec"] >= SPEC_FLOOR)

    # Uncatalogued-variant disclosure: which uncatalogued FKS1 variants, and are they R-only or R+S?
    uncat_by_variant: dict[str, dict] = {}
    for r in records:
        for v in r["uncatalogued"]:
            if v in ("no_genome", "no_consensus"):
                continue
            d = uncat_by_variant.setdefault(v, {"in_R": [], "in_S": []})
            (d["in_R"] if r["label"] == "R" else d["in_S"]).append(r["biosample"])
    def _interpret(d):
        if not d["in_S"]:                       # R-only uncatalogued variant
            return ("CANDIDATE_BLIND_SPOT (R-only, uncatalogued -> possible determinant OR R-lineage "
                    "marker; needs independent evidence, NOT auto-added to catalog)")
        if not d["in_R"]:                        # S-only uncatalogued variant
            return "BENIGN_POLYMORPHISM (present only in susceptible isolate(s) -> not a resistance determinant)"
        return "LINEAGE_MARKER (present in BOTH R and S -> non-discriminative, not causal)"
    disclosure = {v: {"in_R": d["in_R"], "in_S": d["in_S"], "interpretation": _interpret(d)}
                  for v, d in uncat_by_variant.items()}

    artifact = {
        "_schema": "ar-bank-caur-micafungin-validation-v1", "date": _date.today().isoformat(),
        "run_id": a.run_id, "organism": ORGANISM, "drug": "micafungin", "gene": "FKS1",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "not_in_shipped_surface": True,
        "label_source": "ar_bank_MIC_cdc_tentative_breakpoint (micafungin >=4 -> R)",
        "genotype_source": ("BLAST FKS1 (annotated GSC1) target-site; assembled S + SRA-read-mapped R "
                            "FKS1 consensus (assemble_sra_cohort --method map)"),
        "binary": conf, "n_R": n_R, "n_S": n_S,
        "n_indeterminate": sum(1 for r in records if str(r["prediction"]).upper() not in ("R", "S")),
        "sens_wilson95": wilson_ci(conf["tp"], n_R), "spec_wilson95": wilson_ci(conf["tn"], n_S),
        "powered": powered, "spec_floor": SPEC_FLOOR,
        "uncatalogued_fks1_variants": disclosure,
        "catalog_note": ("uncatalogued FKS1 variants are DISCLOSED, never auto-added to the catalog "
                         "(1-isolate, no independent causal evidence = the ERG11 clade IV over-call trap)."),
        "headline": "SCORED_ENDORSED" if endorsed else ("UNDERPOWERED" if not powered else "SCORED_NOT_ENDORSED"),
        "frozen_surface_changed": False,
    }
    outdir = Path("data/raw/ar_bank_caur_extval_micafungin")
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "predictions_fks1_final.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    out = Path(f"wiki/ar_bank_caur_micafungin_validation_{a.run_id}.json")
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nRESULT C. auris micafungin (FKS1): n={len(scored)} ({n_R}R/{n_S}S) acc={conf['acc']} "
          f"sens={conf['sens']} spec={conf['spec']} -> {artifact['headline']}")
    if disclosure:
        print("  uncatalogued FKS1 variants:", json.dumps(disclosure))
    print(f"  artifact: {out}")
    return 0 if powered else 1


if __name__ == "__main__":
    raise SystemExit(main())
