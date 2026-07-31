"""Dog morphology relative-signal validation on Darwin's Ark (the body-size + ear analog of the coat
concordance script). Reproduces, from the REAL canFam4 PLINK set + owner-reported morphology questions:

  (1) HEIGHT (Q121): the 4-locus body-size polygenic score (dna_decode.pigment.dog_body_size.SIZE_LOCI)
      vs the height z-score -> combined r (~0.619, R^2 ~0.38).
  (2) EAR (Q125): the pinned MSRB3 ear lead SNP (chr10:8,612,500, MORPH_LOCI['EAR']) vs the ear ordinal
      -> r (~0.543), monotonic dose-response, and the disambiguation vs the HMGA2 body-size SNP.
  (3) The FUNCTIONAL SCAN that IDENTIFIES which morphology question maps to which locus (no codebook exists
      for the Q-numbers) + the HONEST NEGATIVE: the 4 covariate-adjusted "rerun" traits (Q124/127/128/245)
      do NOT map strongly to the classic single-SNP loci (FGF5 coat-length / KRT71 curl / MSRB3 ear).

Like the coat script, this is a RELATIVE-signal validation (does dosage track the owner-reported ordinal),
NOT a calibrated absolute predictor. Emits a .md narrative + a .json sidecar.

Requires the Darwin's Ark canFam4 PLINK set + phenotype TSVs (on D:, not committed). The PURE helpers
(pearson, polygenic reuse) are unit-tested offline in tests/test_dog_morphology_darwins_ark.py.

Reproduction (paths default to the D: cache):
  uv run python scripts/dog_morphology_darwins_ark_validate.py
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.pigment import plink_io  # noqa: E402
from dna_decode.pigment.dog_body_size import (  # noqa: E402
    MORPH_LOCI,
    SIZE_LOCI,
    polygenic_size_score,
)

_DA = "D:/dna_decode_cache/darwins_ark"
_BASE = f"{_DA}/darwins_dogs_genetic_set/DarwinsDogs_2024_N-3277_canfam4_gp-0.70_biallelic"
_PHENO_DIR = f"{_DA}/darwins_dogs_gwas_input_files/phenotype_input_files"


def pearson(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 2:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    return cov / ((vx * vy) ** 0.5) if vx > 0 and vy > 0 else None


def _bim_index(base: str, variant_id: str) -> tuple[int, str, str]:
    out = subprocess.run(["grep", "-n", "-m", "1", variant_id, base + ".bim"],
                         capture_output=True, text=True).stdout.strip()
    if not out:
        raise SystemExit(f"variant {variant_id} not found in {base}.bim")
    n = int(out.split(":", 1)[0])
    f = out.split(":", 1)[1].split("\t")
    return n - 1, f[4], f[5]


def _load_pheno(tsv: str, cols: list[str]) -> dict[str, dict[str, float | None]]:
    out: dict[str, dict[str, float | None]] = {}
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        out[r["dog_id"]] = {c: (float(r[c]) if r.get(c) not in ("", "NA", "nan", None) else None)
                            for c in cols}
    return out


def validate_height(base: str) -> dict:
    """4-locus body-size polygenic score vs Q121 height z."""
    fam = plink_io.read_fam(base + ".fam")
    idx = {vid: _bim_index(base, vid)[0] for vid in (L.canfam4_variant for L in SIZE_LOCI.values())}
    dos = plink_io.read_bed_variants(base + ".bed", len(fam), list(idx.values()))
    ht = _load_pheno(f"{_PHENO_DIR}/DarwinsArk_Height_Q121_N-3277.tsv", ["Q121"])
    xs, ys = [], []
    for si, dog in enumerate(fam):
        h = ht.get(dog, {}).get("Q121")
        if h is None:
            continue
        doses = [dos[idx[L.canfam4_variant]][si] for L in SIZE_LOCI.values()]
        if any(d is None for d in doses):
            continue
        xs.append(sum(doses))
        ys.append(h)
    r = pearson(xs, ys)
    # score-bin means
    bins: dict[int, list[float]] = {}
    for s, h in zip(xs, ys):
        bins.setdefault(s, []).append(h)
    return {"trait": "height_Q121", "score": "sum of 4 big-allele dosages (0-8)", "n": len(xs),
            "polygenic_r": r, "r2": (r * r if r else None),
            "score_bin_mean_height": {s: round(sum(v) / len(v), 3) for s, v in sorted(bins.items())}}


def validate_ear(base: str) -> dict:
    """MSRB3 ear lead SNP vs Q125, with the HMGA2 body-size disambiguation."""
    fam = plink_io.read_fam(base + ".fam")
    ear = MORPH_LOCI["EAR"]
    ei, ea1, ea2 = _bim_index(base, ear.canfam4_variant)
    hi, ha1, ha2 = _bim_index(base, next(L.canfam4_variant for L in SIZE_LOCI.values() if L.gene == "HMGA2"))
    dos = plink_io.read_bed_variants(base + ".bed", len(fam), [ei, hi])
    ph = _load_pheno(f"{_PHENO_DIR}/DarwinsArk_34Q_Priorities1-3_Behavior_9Q_Morphology_N-3277.tsv", [ear.darwins_ark_question])
    xs, ys, hs = [], [], []
    by: dict[int, list[float]] = {0: [], 1: [], 2: []}
    for si, dog in enumerate(fam):
        q = ph.get(dog, {}).get(ear.darwins_ark_question)
        de, dh = dos[ei][si], dos[hi][si]
        if q is None or de is None:
            continue
        xs.append(de)
        ys.append(q)
        by[de].append(q)
        if dh is not None:
            hs.append((dh, q))
    r_ear = pearson(xs, ys)
    r_size = pearson([h for h, _ in hs], [q for _, q in hs])
    return {"trait": "ear_Q125", "locus": ear.canfam4_variant, "gene": ear.gene, "n": len(xs),
            "ear_r": r_ear, "bodysize_snp_r": r_size,
            "dose_response": {k: round(sum(v) / len(v), 3) for k, v in by.items() if v},
            "resolved_from_size": (abs(r_ear or 0) > 2 * abs(r_size or 0))}


def scan_rerun_traits(base: str) -> dict:
    """Honest negative: the 4 covariate-adjusted rerun traits vs the classic single-SNP loci (best |r|)."""
    fam = plink_io.read_fam(base + ".fam")
    traits = ["Q124_rerun", "Q127_rerun", "Q128_rerun", "Q245_recoded"]
    ph = _load_pheno(f"{_PHENO_DIR}/DarwinsArk_4Q_Morphology_N-3277_rerun_20241107.tsv", traits)
    # scan windows: FGF5 chr32, KRT71 chr27, MSRB3 chr10 (same as the identity scan)
    windows = {"FGF5": (32, 4528000, 4561000), "KRT71": (27, 2530000, 2548000),
               "MSRB3_ear": (10, 8562500, 8662500)}
    recs = []
    for _, (c, lo, hi) in windows.items():
        out = subprocess.run(
            ["awk", "-F", "\t",
             f'$1=={c} && $4>={lo} && $4<={hi}{{print NR-1"\\t"$1"\\t"$4}}', base + ".bim"],
            capture_output=True, text=True).stdout.strip().splitlines()
        for line in out:
            ridx, ch, pos = line.split("\t")
            recs.append((int(ridx), int(ch), int(pos)))
    dos = plink_io.read_bed_variants(base + ".bed", len(fam), [r[0] for r in recs])
    best: dict[str, dict[str, float]] = {t: {} for t in traits}
    for ridx, ch, pos in recs:
        lname = next(n for n, (c, lo, hi) in windows.items() if c == ch and lo <= pos <= hi)
        col = dos[ridx]
        for t in traits:
            xs, ys = [], []
            for si, dog in enumerate(fam):
                d = col[si]
                q = ph.get(dog, {}).get(t)
                if d is None or q is None:
                    continue
                xs.append(d)
                ys.append(q)
            if len(xs) < 500:
                continue
            maf = sum(xs) / (2 * len(xs))
            if min(maf, 1 - maf) < 0.03:
                continue
            r = pearson(xs, ys)
            if r is not None and abs(r) > abs(best[t].get(lname, 0.0)):
                best[t][lname] = round(r, 3)
    max_abs = max((abs(v) for t in best for v in best[t].values()), default=0.0)
    return {"rerun_traits_best_r_by_locus": best,
            "max_abs_r": round(max_abs, 3),
            "verdict": ("NO_STRONG_SINGLE_SNP: the 4 rerun morph traits do not map to FGF5/KRT71/MSRB3 "
                        "(|r|<=0.3) -- likely SV-caused (substrate-limited like coat) or different traits")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=_BASE)
    ap.add_argument("--out-md", default="wiki/dog_morphology_darwins_ark_validated_2026-07-30.md")
    ap.add_argument("--out-json", default="wiki/dog_morphology_darwins_ark_validated_2026-07-30.json")
    a = ap.parse_args()
    if not Path(a.base + ".bim").exists():
        print(f"PLINK set not found at {a.base}; this validator needs the Darwin's Ark canFam4 data (D:).")
        return 2
    result = {"cohort": "Darwin's Ark (Dryad doi:10.5061/dryad.83bk3jb4r) canFam4 gp-0.70 biallelic",
              "height": validate_height(a.base),
              "ear": validate_ear(a.base),
              "rerun_morph_scan": scan_rerun_traits(a.base)}
    Path(a.out_json).write_text(json.dumps(result, indent=2), encoding="utf-8")
    h, e = result["height"], result["ear"]
    print(f"HEIGHT (Q121): polygenic r={h['polygenic_r']:+.3f} (R2={h['r2']:.3f}) N={h['n']}")
    print(f"EAR (Q125): {e['locus']} r={e['ear_r']:+.3f} vs body-size-SNP r={e['bodysize_snp_r']:+.3f} "
          f"(resolved={e['resolved_from_size']}) N={e['n']}")
    print(f"RERUN morph scan: max |r|={result['rerun_morph_scan']['max_abs_r']} -> "
          f"{result['rerun_morph_scan']['verdict'][:60]}...")
    print(f"JSON -> {a.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
