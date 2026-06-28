"""Score the FROZEN bacterial AMR rule on the AMR Portal's provenance-disjoint, measured-phenotype isolates
— the first genuinely-INDEPENDENT numbers for the deployed cells, no Docker / no genome fetch.

The AMR Portal ships a per-isolate GENOTYPE table in AMRFinderPlus format (`amr_element_symbol`, `class`,
`subclass`, `element_subtype` AMR/POINT) — the EXACT fields the frozen `amr_rules.call_resistance` consumes
from an AMRFinder `main.tsv`. So for each provenance-disjoint isolate (BioSample/ERS/GCA not in CRyPTIC or
our cohorts) that has a MEASURED `resistance_phenotype`, we reconstruct a faithful main.tsv from its
determinants and call the FROZEN rule unchanged → confusion vs the measured label → sens/spec + Wilson CI.

HONESTY RAILS:
 - INDEPENDENT vs the frozen cells' tuning cohorts: leakage is accession-disjoint (BioSample/ERS/GCA). This
   is an ACCESSION-STRING upper bound on independence; BioSample cross-archive resolution (biosample_resolver)
   would only tighten it. The phenotype is wet-lab MIC/disk (non-circular).
 - The genotype is the AMR PORTAL's own AMRFinder run (a different operator/possibly-different version than
   our pinned image) — that makes it MORE independent (different pipeline), but the AMRFinder version is a
   named caveat (point-mutation calling can differ across versions).
 - The frozen rule is applied UNCHANGED via `call_resistance` (organism=None → the validated DRUG_RULE path,
   exactly as the frozen E. coli cells were validated). FROZEN surface byte-unchanged.

Pure logic (`wilson_ci` / `genotype_to_main_tsv` / `confusion`) unit-tested; the parquet load + frozen-rule
call is the live part.
"""
from __future__ import annotations

import math
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.eval.amr_rules import call_resistance  # noqa: E402  (FROZEN rule, called unchanged)
from scripts.amr_portal_feasibility import (  # noqa: E402  (reuse the verified leakage + SIR helpers)
    DEFAULT_CRYPTIC, DEFAULT_PARQUET, binarize_sir, iso_leaked, _load_leak_set,
)

GENO_PARQUET = Path("D:/dna_decode_cache/data files donwload/amr_portal/genotype.parquet")
# Per-organism rule routing: AMR-Portal organism name -> the `organism=` passed to call_resistance.
#   None             -> the validated DRUG_RULE default (E. coli; Shigella shares the E. coli rules).
#   "Klebsiella"/"Salmonella" -> the OPT-IN calibrated registry (cipro: Kleb qrdr_point+oqxAB-exclusion,
#                    Salmonella broad) — IN-SAMPLE-calibrated (N~30) configs whose own provenance asks for an
#                    INDEPENDENT cohort before becoming default; the AMR Portal disjoint set IS that cohort.
RULE_ORGANISM = {
    "Escherichia coli": None, "Shigella sonnei": None, "Shigella flexneri": None,
    "Klebsiella pneumoniae": "Klebsiella", "Salmonella enterica": "Salmonella",
    # Campylobacter: routes to the endorsed calibrated cipro rule (qrdr_point, LOO balacc 1.0); added
    # 2026-06-28 (it has a deployed rule + is SCORED on the NCBI-PD card, and AMR Portal has it powered).
    "Campylobacter jejuni": "Campylobacter", "Campylobacter coli": "Campylobacter",
}
_DRUGS = ("ciprofloxacin", "ceftriaxone", "tetracycline", "gentamicin", "meropenem")
# Scope contract (2026-06-28): the 5 deployed-default/calibrated organisms get ALL 5 drug rules;
# Campylobacter is endorsed for ciprofloxacin ONLY (its sole calibrated registry rule). The generic
# DRUG_RULE default for cef/mero/tet/gent is E. coli-derived + NOT validated for Campylobacter, so it is
# NOT sprayed across it (the intrinsic-gene over-call guardrail). TB uses a different rule (WHO catalogue
# on VCF) -> scored separately, not here.
_FULL_DRUG_ORGS = ("Escherichia coli", "Shigella sonnei", "Shigella flexneri",
                   "Klebsiella pneumoniae", "Salmonella enterica")
_CIPRO_ONLY_ORGS = ("Campylobacter jejuni", "Campylobacter coli")
CELLS = ([(org, drug) for org in _FULL_DRUG_ORGS for drug in _DRUGS]
         + [(org, "ciprofloxacin") for org in _CIPRO_ONLY_ORGS])
_MAIN_TSV_HEADER = "Element symbol\tMethod\tClass\tSubclass\tElement name\t% Identity to reference"


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float] | None:
    """Wilson score interval for a binomial proportion k/n. None if n==0."""
    if n == 0:
        return None
    phat = k / n
    denom = 1 + z * z / n
    centre = (phat + z * z / (2 * n)) / denom
    half = (z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def genotype_to_main_tsv(dets: list[dict]) -> str:
    """AMR Portal genotype rows for one isolate -> AMRFinder main.tsv text (the columns the frozen parser
    reads). element_subtype POINT -> Method 'POINTX' (the parser's `'POINT' in Method` gate); else 'EXACTX'."""
    lines = [_MAIN_TSV_HEADER]
    for d in dets:
        method = "POINTX" if (d.get("element_subtype") or "").upper() == "POINT" else "EXACTX"
        sym = (d.get("amr_element_symbol") or d.get("gene_symbol") or "").strip()
        name = (d.get("evidence_description") or d.get("gene_symbol") or "").strip()
        lines.append("\t".join([sym, method, (d.get("class") or ""), (d.get("subclass") or ""), name, ""]))
    return "\n".join(lines) + "\n"


def confusion(pred: str, label: str) -> str | None:
    """('R'/'S' pred, 'R'/'S' label) -> 'TP'/'FP'/'TN'/'FN'; None if pred not R/S."""
    if pred not in ("R", "S"):
        return None
    if pred == "R":
        return "TP" if label == "R" else "FP"
    return "TN" if label == "S" else "FN"


def _load_genotype(parquet: Path, organisms: set[str]) -> dict[str, list[dict]]:
    import pyarrow.parquet as pq
    cols = ["BioSample_ID", "organism", "amr_element_symbol", "gene_symbol", "class", "subclass",
            "element_subtype", "evidence_description"]
    by_bs: dict[str, list[dict]] = defaultdict(list)
    pf = pq.ParquetFile(str(parquet))
    for batch in pf.iter_batches(batch_size=200_000, columns=cols):
        d = batch.to_pydict()
        for i in range(batch.num_rows):
            if d["organism"][i] not in organisms and not any(o in (d["organism"][i] or "") for o in organisms):
                continue
            bs = (d["BioSample_ID"][i] or "").strip()
            if bs:
                by_bs[bs].append({k: d[k][i] for k in cols})
    return by_bs


def _load_phenotype_cells(parquet: Path, cells, leak: set[str]):
    """-> {(org,drug): [(biosample, leaked_bool, SIR)]} for disjoint+leaked isolates with a measured label."""
    import pyarrow.parquet as pq
    from scripts.amr_portal_feasibility import DRUG_ALIASES
    want = {(o, d) for o, d in cells}
    cols = ["organism", "antibiotic_name", "BioSample_ID", "SRA_accession", "assembly_ID", "resistance_phenotype"]
    seen: set = set()
    out = defaultdict(list)
    pf = pq.ParquetFile(str(parquet))
    for batch in pf.iter_batches(batch_size=200_000, columns=cols):
        dd = batch.to_pydict()
        for org, ab, bs, sra, asm, rp in zip(dd["organism"], dd["antibiotic_name"], dd["BioSample_ID"],
                                              dd["SRA_accession"], dd["assembly_ID"], dd["resistance_phenotype"]):
            drug = DRUG_ALIASES.get((ab or "").strip().lower())
            if drug is None or (org, drug) not in want:
                continue
            sir = binarize_sir(rp)
            if sir is None:
                continue
            bs = (bs or "").strip()
            key = (org, drug, bs or sra or asm)
            if key in seen:
                continue
            seen.add(key)
            leaked = iso_leaked({bs, (sra or "").strip(), (asm or "").strip()}, leak)
            out[(org, drug)].append((bs, leaked, sir))
    return out


def score(cells=CELLS, pheno_parquet: Path = DEFAULT_PARQUET, geno_parquet: Path = GENO_PARQUET,
          cryptic: Path = DEFAULT_CRYPTIC, min_per_class: int = 10) -> dict:
    leak = _load_leak_set(cryptic)
    organisms = {o for o, _ in cells}
    geno = _load_genotype(geno_parquet, organisms)
    pheno = _load_phenotype_cells(pheno_parquet, cells, leak)
    results = {}
    tmp = Path(tempfile.gettempdir()) / "_amrportal_main.tsv"
    for (org, drug), isolates in pheno.items():
        conf = {"TP": 0, "FP": 0, "TN": 0, "FN": 0}
        n_no_geno = n_indet = 0
        for bs, leaked, sir in isolates:
            if leaked:
                continue                              # provenance-disjoint only
            dets = geno.get(bs)
            if dets is None:
                # no genotype row for this BioSample = AMRFinder found NO AMR determinants -> empty main.tsv
                dets = []
            tmp.write_text(genotype_to_main_tsv(dets), encoding="utf-8")
            # FROZEN rule; organism routes E. coli/Shigella -> DRUG_RULE default, Kleb/Salmonella -> calibrated.
            pred = call_resistance(tmp, drug, organism=RULE_ORGANISM.get(org))["prediction"]
            if pred not in ("R", "S"):
                n_indet += 1
                continue
            conf[confusion(pred, sir)] += 1
        tp, fp, tn, fn = conf["TP"], conf["FP"], conf["TN"], conf["FN"]
        n_R, n_S = tp + fn, tn + fp
        sens = tp / n_R if n_R else None
        spec = tn / n_S if n_S else None
        acc = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) else None
        results[f"{org}|{drug}"] = {
            "organism": org, "drug": drug, "confusion": conf, "n_R": n_R, "n_S": n_S,
            "sens": sens, "spec": spec, "accuracy": acc, "n_indeterminate": n_indet,
            "sens_wilson95": wilson_ci(tp, n_R), "spec_wilson95": wilson_ci(tn, n_S),
            "powered": n_R >= min_per_class and n_S >= min_per_class,
            "status": "PROVENANCE_DISJOINT_INDEPENDENT_ACCESSION_LEVEL",
        }
    return results


def main(argv=None) -> int:
    import argparse, json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pheno", type=Path, default=DEFAULT_PARQUET)
    ap.add_argument("--geno", type=Path, default=GENO_PARQUET)
    ap.add_argument("--cryptic", type=Path, default=DEFAULT_CRYPTIC)
    a = ap.parse_args(argv)
    for p in (a.pheno, a.geno):
        if not p.exists():
            print(f"ERROR: AMR Portal parquet not found: {p}", file=sys.stderr)
            return 2
    res = score(pheno_parquet=a.pheno, geno_parquet=a.geno, cryptic=a.cryptic)
    print(f"{'cell':<40} {'nR':>5} {'nS':>5} {'sens':>6} {'spec':>6} {'acc':>6} powered")
    for k, r in res.items():
        s = lambda x: f"{x:.3f}" if isinstance(x, float) else "  -  "
        print(f"{k:<40} {r['n_R']:>5} {r['n_S']:>5} {s(r['sens']):>6} {s(r['spec']):>6} {s(r['accuracy']):>6} "
              f"{'YES' if r['powered'] else 'no'}")
        for m in ("sens", "spec"):
            ci = r[f"{m}_wilson95"]
            if ci:
                print(f"    {m} 95% CI: [{ci[0]:.3f}, {ci[1]:.3f}]")
    out = REPO / "wiki" / "amr_portal_independent_scores.json"
    out.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
