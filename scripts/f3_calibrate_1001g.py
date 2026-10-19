"""F3: run the de-confounder calibration for real, with 1001-Genomes as the positive control.

This is the measurement `wiki/f3_deconfound_calibration_2026-10-09.md` could not make. That run
returned `blocked:positive_control_not_fold_compatible` because `dgrp` -- the only cached natural
population -- has no grouping that is both fold-compatible (largest cluster <= 1/k) and not mostly
singletons. 1001-Genomes does (measured: p5, 108 clusters, largest 9.5%), so the control pair is now
complete.

WHAT IS REUSED, DELIBERATELY: the FROZEN verdict machinery (`compute_verdict`, `sign_test`) is imported
from `f3_deconfound_calibration`, not reimplemented -- the verdict states were named before any number
existed and must not drift. The negative control is the SAME bloom arm, through the same `run_arm`.

THE CROSS-CHECK IS PART OF THE DELIVERABLE, not an afterthought. The 2026-10-09 grouping came from the
1001G project's OWN kinship matrix (IBS at MAC>=5), converted to a distance. This run ALSO derives the
grouping from this repo's own `panel_structure.ibs_distance_matrix` over the raw genotypes and reports
both. A disagreement is a finding to diagnose, not a number to choose between.

Run (long; ~tens of minutes at the defaults):
    uv run python scripts/f3_calibrate_1001g.py --smoke          # fast correctness check first
    uv run python scripts/f3_calibrate_1001g.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

GENO_H5 = "D:/dna_decode_cache/arabidopsis_1001g/1001_SNP_MATRIX/imputed_snps_binary.hdf5"
KINSHIP_H5 = "D:/dna_decode_cache/arabidopsis_1001g/1001_SNP_MATRIX/kinship_ibs_mac5.hdf5"
PHENO = {
    "FT10": "data/arabidopsis/FT10_pheno_261.csv",
    "FT16": "data/arabidopsis/FT16_pheno_262.csv",
    "FLC": "data/arabidopsis/flowering_1001g/pheno_29_FLC.csv",
}

MAF_FLOOR = 0.05          # asserted, standard for a GP ridge
N_BLOCKS = 60             # chunk-aligned blocks of 1000 SNPs to read before filtering
N_MARKERS = 5000          # markers kept after MAF filter + subsample
N_SPLITS = 5
N_PERM = 100


def load_1001g_genotypes(*, n_blocks: int = N_BLOCKS, n_markers: int = N_MARKERS,
                         maf_floor: float = MAF_FLOOR, seed: int = 0):
    """(X, ids) with X as (individuals, markers), MAF-filtered and subsampled.

    The HDF5 is (10.7M SNPs, 1135 accessions) chunked (1000, 1135), so a random row gather would
    read a whole chunk per row. Reading CHUNK-ALIGNED BLOCKS spread across the genome costs ~n_blocks
    chunk reads instead of ~n_markers, and spreading them keeps the sample off any single region.
    """
    import h5py

    rng = np.random.default_rng(seed)
    with h5py.File(GENO_H5, "r") as f:
        snps = f["snps"]
        ids = [s.decode() for s in f["accessions"][:]]
        n_snp = snps.shape[0]
        starts = np.linspace(0, n_snp - 1000, n_blocks).astype(int) // 1000 * 1000
        blocks = [np.asarray(snps[s:s + 1000, :], dtype=np.int8) for s in starts]
    M = np.vstack(blocks)                       # (n_blocks*1000, 1135)
    X = M.T.astype(np.float64)                  # -> (1135, markers)

    # MAF filter on the observed coding (imputed binary: 0/1)
    p = X.mean(axis=0)
    maf = np.minimum(p, 1.0 - p)
    keep = maf >= maf_floor
    X = X[:, keep]
    if X.shape[1] > n_markers:
        pick = rng.choice(X.shape[1], size=n_markers, replace=False)
        X = X[:, np.sort(pick)]
    return X, ids, {"n_snps_read": int(M.shape[0]), "n_after_maf": int(keep.sum()),
                    "n_markers_used": int(X.shape[1]), "maf_floor": maf_floor,
                    "n_blocks": n_blocks}


def load_1001g_traits(ids, max_traits: int = 3, *, min_distinct: int = 20) -> dict:
    """Traits aligned to the genotype accessions, preferring the ranked AraPheno cache.

    `min_distinct` keeps the target CONTINUOUS: a ridge on a binary/ordinal score is a different
    problem, and 20 distinct levels is the same floor the repo's `assay_degeneracy` screen uses.
    Falls back to the three committed CSVs when the cache is absent.
    """
    import pandas as pd

    idset = {str(i) for i in ids}
    pos = {str(i): k for k, i in enumerate(ids)}

    def _vec(pairs):
        y = np.full(len(ids), np.nan)
        hit = 0
        for acc, val in pairs:
            if acc is None or val is None:
                continue
            try:
                k = str(int(acc)); v = float(val)
            except (TypeError, ValueError):
                continue
            if k in idset:
                y[pos[k]] = v
                hit += 1
        return y, hit

    out = {}
    rank = sorted((REPO / "wiki").glob("arapheno_trait_overlap_*.json"))
    cache = Path("D:/dna_decode_cache/arapheno")
    if rank and cache.exists():
        ranked = json.loads(rank[-1].read_text(encoding="utf-8")).get("usable", [])
        # ONE TRAIT PER STUDY. Measured: 163 usable traits span only 8 studies, and 121 of them come
        # from a single rainfall-manipulation factorial -- so a naive top-N by overlap would hand the
        # sign test six correlated replicates of one experiment and call them six observations. Within
        # each study the highest-overlap trait wins (`ranked` is already overlap-sorted).
        seen_studies: set[str] = set()
        for row in ranked:
            if len(out) >= max_traits:
                break
            study = str(row.get("study", ""))[:60]
            if study in seen_studies:
                continue
            f = cache / f"values_{row['phenotype_id']}.json"
            if not f.exists():
                continue
            try:
                vals = json.loads(f.read_text(encoding="utf-8"))
            except Exception:                         # noqa: BLE001
                continue
            y, hit = _vec((r.get("accession_id"), r.get("phenotype_value")) for r in vals)
            fin = y[~np.isnan(y)]
            if hit < N_SPLITS * 4 or len(np.unique(fin)) < min_distinct:
                continue
            name = f"ap{row['phenotype_id']}_{re.sub(r'[^A-Za-z0-9]+', '-', str(row['name']))[:28]}"
            out[name] = y
            seen_studies.add(study)
    if len(out) >= max_traits:
        return out

    for name, rel in PHENO.items():                   # committed-CSV fallback / top-up
        if len(out) >= max_traits:
            break
        fp = REPO / rel
        if not fp.exists():
            continue
        d = pd.read_csv(fp)
        y, hit = _vec(zip(d.accession_id, d.phenotype_value))
        fin = y[~np.isnan(y)]
        if hit >= N_SPLITS * 4 and len(np.unique(fin)) >= min_distinct:
            out[name] = y
    return out


def load_bloom(*, n_markers: int):
    """Bloom (bi-parental yeast cross) genotypes + traits, strided to a comparable marker budget.

    The committed `run_arm` loads bloom at stride 1 = 11,623 markers, while the positive arm runs on a
    ~1,500-marker subsample. That is both ~7h of extra compute AND an unfair asymmetry: a differential
    comparison between two arms should not also differ by an order of magnitude in marker count. The
    stride is computed to land near `n_markers`.
    """
    import importlib

    from scripts.f3_substrate_census import _unpack_loaded

    mod = importlib.import_module("scripts.yeast_bloom_gp_arm")
    probe, _ = _unpack_loaded(mod.load_genotype(str(mod._GENO), 1))
    stride = max(1, probe.shape[1] // max(1, n_markers))
    X, ids = _unpack_loaded(mod.load_genotype(str(mod._GENO), stride))
    ph = mod.load_phenotype(str(mod._PHENO))
    from scripts.f3_deconfound_calibration import _phenotype_traits
    traits = _phenotype_traits(ph, len(ids), 99)
    return X, ids, traits, {"stride": stride, "n_markers_used": int(X.shape[1]),
                            "n_markers_full": int(probe.shape[1])}


def write_crosscheck(out_path, ids, repo_pg, coarse_pg, rung_sweep,
                     repo_groups, kin_groups, kin_meta) -> float:
    """Pair-counting agreement between this repo's IBS grouping and the 1001G project's kinship.

    Adjusted Rand would pull in sklearn.metrics; pair counting needs no extra dep and is exactly
    interpretable. REFUSES on an empty common set rather than writing an all-zero artifact.
    """
    common = [i for i in ids if i in kin_groups and i in repo_groups]
    n = len(common)
    if n == 0:
        raise ValueError("no accession is in both groupings -- refusing to write a zeroed crosscheck")
    rg = np.array([repo_groups[i] for i in common])
    kg = np.array([kin_groups[i] for i in common])
    same_repo = rg[:, None] == rg[None, :]
    same_kin = kg[:, None] == kg[None, :]
    off = ~np.eye(n, dtype=bool)
    both = int((same_repo & same_kin & off).sum() // 2)
    only_repo = int((same_repo & ~same_kin & off).sum() // 2)
    only_kin = int((~same_repo & same_kin & off).sum() // 2)
    neither = int(((~same_repo) & (~same_kin) & off).sum() // 2)
    jacc = both / max(1, both + only_repo + only_kin)
    crosscheck = {
        "schema": "f3-ibs-crosscheck-v1", "generated": date.today().isoformat(),
        "n_common_accessions": n,
        "repo_ibs_grouping_finest": repo_pg.as_dict(),
        "repo_ibs_grouping_coarsest": coarse_pg.as_dict(),
        "admissible_rung_sweep": [{k: (round(v, 4) if isinstance(v, float) else v)
                                   for k, v in r.items() if k != "refusal"} for r in rung_sweep],
        "project_kinship_grouping": kin_meta | {"n_clusters": int(len(set(kg)))},
        "pair_agreement": {"same_in_both": both, "same_only_repo_ibs": only_repo,
                           "same_only_project_kinship": only_kin, "different_in_both": neither,
                           "jaccard_of_same_cluster_pairs": round(jacc, 4)},
        "honest_limits": [
            "The two groupings derive from DIFFERENT distances: this repo's IBS over a MAF-filtered "
            "marker SUBSAMPLE vs the 1001G project's own IBS kinship at MAC>=5 over all markers. "
            "They are closely related, not the same quantity.",
            "Pair-counting agreement is reported instead of adjusted Rand to avoid adding a sklearn "
            "metrics dependency; it is exactly interpretable but NOT chance-corrected.",
            "The calibration uses the REPO IBS grouping, so the verdict does not depend on the "
            "project's converted kinship.",
        ],
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(crosscheck, indent=2), encoding="utf-8")
    print(f"      repo IBS: {repo_pg.n_clusters} clusters p={repo_pg.percentile:g} "
          f"largest={repo_pg.largest_cluster_fraction:.3f} | kinship: {len(set(kg))} clusters "
          f"p={kin_meta['percentile']:g} | same-pair Jaccard={jacc:.4f}", flush=True)
    return jacc


PREAMBLE_DIR = Path("D:/dna_decode_cache/f3_calibration")


def cached_preamble(*, n_blocks: int, n_markers: int, maf_floor: float, max_largest: float):
    """(X, ids, meta, fine, coarse, sweep), reusing a disk cache when it matches exactly.

    The cache name carries every parameter that changes the matrix, and the stored ids are checked
    against the stored shape on load. A mismatch DISCARDS the cache -- it never adapts to it.
    """
    import pickle

    PREAMBLE_DIR.mkdir(parents=True, exist_ok=True)
    f = PREAMBLE_DIR / f"preamble_b{n_blocks}_m{n_markers}_maf{maf_floor}_ml{max_largest:.4f}.pkl"
    if f.exists():
        try:
            d = pickle.loads(f.read_bytes())
            X, ids = d["X"], d["ids"]
            if X.shape[0] == len(ids) and X.shape[1] == d["meta"]["n_markers_used"]:
                print(f"      [preamble cache hit] {f.name}", flush=True)
                return X, ids, d["meta"], d["fine"], d["coarse"], d["sweep"]
            print("      [preamble cache DISCARDED: shape/ids mismatch]", flush=True)
        except Exception as e:                        # noqa: BLE001
            print(f"      [preamble cache DISCARDED: {type(e).__name__}]", flush=True)

    X, ids, meta = load_1001g_genotypes(n_blocks=n_blocks, n_markers=n_markers, maf_floor=maf_floor)
    fine, coarse, sweep = admissible_rungs(X, ids, max_largest=max_largest)
    f.write_bytes(pickle.dumps({"X": X, "ids": ids, "meta": meta,
                                "fine": fine, "coarse": coarse, "sweep": sweep}))
    return X, ids, meta, fine, coarse, sweep


def admissible_rungs(X, ids, *, max_largest: float):
    """(finest, coarsest) PanelGroups over the rungs that clear BOTH guards at this fold share.

    `select_threshold` takes the SMALLEST passing rung on a named bias -- it under-merges, which biases
    AGAINST detecting a de-confounding response. The smoke showed that bias can be total (283 clusters
    over 1,135 individuals => grouped CV ~ random CV => a 0.0015 drop). Reporting BOTH endpoints of the
    admissible set brackets the response instead of resting on one end of it.
    """
    from dna_decode.eval.panel_structure import (IBS_PERCENTILE_LADDER, ibs_distance_matrix,
                                                 panel_groups, sweep_thresholds)

    dm = ibs_distance_matrix(X, ids)
    sweep = sweep_thresholds(dm, ids, IBS_PERCENTILE_LADDER, max_largest=max_largest)
    passing = [r["percentile"] for r in sweep if r["passes"]]
    if not passing:
        raise RuntimeError("no admissible rung at this fold share")
    fine = panel_groups(X, ids, percentile=min(passing), max_largest=max_largest)
    coarse = panel_groups(X, ids, percentile=max(passing), max_largest=max_largest)
    return fine, coarse, sweep


def kinship_grouping(ids, *, max_largest: float):
    """The 2026-10-09 grouping: the 1001G project's OWN kinship, converted to a distance."""
    import h5py

    from dna_decode.eval.panel_structure import select_threshold
    from dna_decode.eval.phylogeny import DistanceMatrix

    with h5py.File(KINSHIP_H5, "r") as f:
        k = f["kinship"][:]
        kacc = [s.decode() for s in f["accessions"][:]]
    d = (k.max() - k) / (k.max() - k.min())
    d = (d + d.T) / 2.0
    np.fill_diagonal(d, 0.0)
    dm = DistanceMatrix(strain_ids=kacc, matrix=d)
    pct, thr, _ = select_threshold(dm, kacc, max_largest=max_largest)
    from dna_decode.eval.clonality import greedy_representative_clusters_from_matrix
    assigned = greedy_representative_clusters_from_matrix(dm, thr)
    return {a: assigned[a] for a in kacc}, {"source": "1001G project kinship (IBS MAC>=5), converted",
                                            "percentile": pct, "threshold": thr}


def matched_rung(X, ids, *, target_clusters_per_individual: float, max_largest: float):
    """The admissible rung whose GRANULARITY most nearly matches a target, not its percentile.

    The two control panels have different structure, so the SAME percentile buys very different fold
    geometry: on 1001G p=10 gives 90 clusters over 1,135 accessions, while on the bloom cross p=20
    gives 22 over 1,008. Comparing those two as "the coarse end" would compare a 12.6-individual
    cluster against a 46-individual one, and the drop a de-confounder produces depends on exactly that.

    So the arms are matched on CLUSTERS PER INDIVIDUAL (equivalently mean cluster size), which is the
    quantity that drives the test. The rule is fixed before any drop is read, and the alternative
    (percentile matching) is rejected on the geometry above, not on the results it would give.
    """
    from dna_decode.eval.panel_structure import (IBS_PERCENTILE_LADDER, ibs_distance_matrix,
                                                 panel_groups, sweep_thresholds)

    dm = ibs_distance_matrix(X, ids)
    sweep = sweep_thresholds(dm, ids, IBS_PERCENTILE_LADDER, max_largest=max_largest)
    passing = [r for r in sweep if r["passes"]]
    if not passing:
        raise RuntimeError("no admissible rung at this fold share")
    n = len(ids)
    best = min(passing, key=lambda r: abs(r["n_clusters"] / n - target_clusters_per_individual))
    pg = panel_groups(X, ids, percentile=best["percentile"], max_largest=max_largest)
    return pg, {"target_clusters_per_individual": round(target_clusters_per_individual, 6),
                "chosen_percentile": best["percentile"], "n_clusters": pg.n_clusters,
                "clusters_per_individual": round(pg.n_clusters / n, 6),
                "mean_cluster_size": round(n / pg.n_clusters, 2),
                "candidates": [{"percentile": r["percentile"], "n_clusters": r["n_clusters"],
                                "clusters_per_individual": round(r["n_clusters"] / n, 6)}
                               for r in passing]}


def trait_rung_admissibility(groups_by_id, ids, traits, *, max_largest: float) -> dict:
    """Per-trait largest-cluster share on that trait's NON-MISSING subset, and admissibility.

    Grouped CV folds a trait's non-missing rows, so admissibility is a PER-TRAIT property. A rung can
    clear the one-fold share on the full panel and be inadmissible for a trait whose missingness
    concentrates one cluster -- which is exactly how a bracket endpoint got mis-specified here (the
    p=20 rung clears 0.1947 on all 1,135 accessions and reaches 0.2104 on FT10's 1,003).
    """
    g_all = np.array([groups_by_id[i] for i in ids])
    out = {}
    for nm, y in traits.items():
        m = ~np.isnan(y)
        _, cnt = np.unique(g_all[m], return_counts=True)
        share = float(cnt.max() / cnt.sum())
        out[nm] = {"n": int(m.sum()), "n_clusters": int(len(cnt)),
                   "largest_share": round(share, 4), "admissible": bool(share <= max_largest)}
    return out


def paired_rung_delta(fine_arm: dict, coarse_arm: dict) -> dict:
    """Fine-vs-coarse drop compared ONLY on traits scoreable at both rungs.

    The coarse rung refuses any trait whose non-missing subset pushes a cluster over the one-fold
    share. On this panel those are the highest-n traits, which are also the strongest-signal and
    smallest-drop ones -- so an unpaired median comparison would credit the coarse rung for ATTRITION.
    Reports the attrition explicitly rather than hiding it behind a restricted median.
    """
    if not fine_arm or not coarse_arm:
        return {"status": "an arm is absent this call"}
    f = {r["trait"]: r["drop"] for r in fine_arm["per_trait"] if "drop" in r}
    c = {r["trait"]: r["drop"] for r in coarse_arm["per_trait"] if "drop" in r}
    common = [t for t in f if t in c]
    only_f = sorted(set(f) - set(c))
    only_c = sorted(set(c) - set(f))
    errs = {r["trait"]: r["error"][:120] for r in coarse_arm["per_trait"] if "error" in r}
    if not common:
        return {"status": "no trait scoreable at both rungs", "scoreable_fine_only": only_f,
                "coarse_errors": errs}
    deltas = [c[t] - f[t] for t in common]
    return {
        "status": "ok", "n_common_traits": len(common), "common_traits": common,
        "median_drop_fine_on_common": round(float(np.median([f[t] for t in common])), 4),
        "median_drop_coarse_on_common": round(float(np.median([c[t] for t in common])), 4),
        "median_delta_coarse_minus_fine": round(float(np.median(deltas)), 4),
        "n_coarser_drops_more": int(sum(1 for d in deltas if d > 0)),
        "dropped_by_coarse_rung": only_f,
        "gained_by_coarse_rung": only_c,
        "coarse_rung_refusals": errs,
    }


CKPT = Path("D:/dna_decode_cache/f3_calibration/cells.jsonl")
SEP = bytes([31])  # hash field separator


def _cell_key(label, trait, X, groups, *, n_splits, n_perm, seed) -> str:
    """sha256 over EVERY input that changes the result.

    Omitting one would let a stale row pose as a fresh one. `groups` is included because the whole
    point of the run is to vary it (finest vs coarsest rung) while holding everything else fixed.
    """
    h = hashlib.sha256()
    for part in (label, trait, str(n_splits), str(n_perm), str(seed)):
        h.update(part.encode()); h.update(SEP)
    h.update(np.ascontiguousarray(X, dtype=np.float64).tobytes())
    h.update(SEP)
    h.update(np.asarray([str(g) for g in groups], dtype="U32").tobytes())
    return h.hexdigest()


RAND_CKPT = Path("D:/dna_decode_cache/f3_calibration/random_arm.jsonl")


def _random_key(trait, X, *, n_splits, n_perm, seed) -> str:
    """Key for the UNGROUPED arm, which by construction cannot depend on the grouping.

    Deliberately omits `groups` -- that is the point, not an oversight. `cv_ridge_gp(..., groups=None)`
    never sees a grouping, and the banked cells confirm it empirically: FT10 0.692, DTF1 0.7046,
    rSurvival 0.2412 and Raffinose 0.1817 are byte-identical at p=5 and p=10. Reusing it halves the
    cost of every coarse-rung cell. Everything else that changes the fit IS in the key.
    """
    h = hashlib.sha256()
    h.update(b"random-arm"); h.update(SEP)
    for part in (trait, str(n_splits), str(n_perm), str(seed)):
        h.update(part.encode()); h.update(SEP)
    h.update(np.ascontiguousarray(X, dtype=np.float64).tobytes())
    return h.hexdigest()


def _load_random() -> dict:
    if not RAND_CKPT.exists():
        return {}
    out = {}
    for line in RAND_CKPT.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        if "key" in d and "row" in d:
            out[d["key"]] = d["row"]
    return out


def _save_random(key: str, row: dict) -> None:
    RAND_CKPT.parent.mkdir(parents=True, exist_ok=True)
    with RAND_CKPT.open("a", encoding="utf-8") as fh:
        print(json.dumps({"key": key, "row": row}), file=fh)


def _load_cells() -> dict:
    if not CKPT.exists():
        return {}
    out = {}
    for line in CKPT.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError:                  # a torn final line -> drop it, never guess
            continue
        if "key" in d and "row" in d:
            out[d["key"]] = d["row"]
    return out


def _save_cell(key: str, row: dict) -> None:
    CKPT.parent.mkdir(parents=True, exist_ok=True)
    with CKPT.open("a", encoding="utf-8") as fh:
        print(json.dumps({"key": key, "row": row}), file=fh)


def run_1001g_arm(X, ids, traits, groups, *, n_splits=N_SPLITS, n_perm=N_PERM, seed=0,
                  label="1001g", resume: bool = True) -> dict:
    from dna_decode.eval.genomic_prediction import cv_ridge_gp

    from scripts.f3_deconfound_calibration import sign_test

    cells = _load_cells() if resume else {}
    rand = _load_random() if resume else {}
    rows = []
    for name, y in traits.items():
        key = _cell_key(label, name, X, groups, n_splits=n_splits, n_perm=n_perm, seed=seed)
        if key in cells:
            row = dict(cells[key]); row["from_checkpoint"] = True
            rows.append(row)
            print(f"        [cached] {row}", flush=True)
            continue
        rk = _random_key(name, X, n_splits=n_splits, n_perm=n_perm, seed=seed)
        try:
            cached_rand = rand.get(rk)
            if cached_rand is None:
                fit = cv_ridge_gp(X, y, trait=name, n_splits=n_splits, seed=seed, n_perm=n_perm)
                a = {"n": fit.n, "predictive_r": fit.predictive_r,
                     "beats_null": bool(fit.beats_null), "null_r_p95": fit.null_r_p95}
                _save_random(rk, a)
                rand[rk] = a
            else:
                a = cached_rand
            b = cv_ridge_gp(X, y, trait=name, n_splits=n_splits, seed=seed, n_perm=n_perm,
                            groups=groups)
        except Exception as e:                        # noqa: BLE001
            # PRINTED, never silent. An unprinted error row makes an arm that scored 6 of 8 traits
            # look like an arm that scored 8 -- the failure mode where a test that errors reads as a
            # finding. Error cells are deliberately NOT checkpointed: the cause may be fixable.
            msg = f"{type(e).__name__}: {e}"
            rows.append({"trait": name, "error": msg})
            print(f"        [ERROR] {name}: {msg[:200]}", flush=True)
            continue
        row = {
            "trait": name, "n": a["n"],
            "r_random": round(a["predictive_r"], 4), "r_grouped": round(b.predictive_r, 4),
            "drop": round(a["predictive_r"] - b.predictive_r, 4),
            "random_beats_own_null": bool(a["beats_null"]),
            "grouped_beats_own_null": bool(b.beats_null),
            "null_p95_random": round(a["null_r_p95"], 4),
            "null_p95_grouped": round(b.null_r_p95, 4),
            "n_groups": b.n_groups,
            "random_arm_reused": cached_rand is not None,
        }
        _save_cell(key, row)
        rows.append(row)
        print(f"        {row}", flush=True)
    scored = [r for r in rows if "drop" in r]
    drops = [r["drop"] for r in scored]
    return {
        "panel": label, "design": "natural_population",
        "n_traits": len(scored), "per_trait": rows,
        "median_drop": float(np.median(drops)) if drops else None,
        "n_dropped": sum(1 for d in drops if d > 0),
        "sign_test": sign_test(drops),
        "grouped_beats_own_null": (sum(1 for r in scored if r["grouped_beats_own_null"])
                                   > len(scored) / 2 if scored else False),
    }


def main(argv=None) -> int:
    from scripts.f3_deconfound_calibration import (MIN_TRAITS_FOR_SIGN_TEST, VERDICT_MEANING,
                                                   VERDICTS, compute_verdict, run_arm)
    from scripts.f3_substrate_census import PANELS
    from dna_decode.eval.panel_structure import panel_groups

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true",
                    help="fast correctness pass: few markers, few perms, 1 trait")
    ap.add_argument("--n-markers", type=int, default=N_MARKERS)
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--coarse-percentile", type=float, default=None,
                    help="IBS ladder rung for the COARSE end of the bracket. Default = the coarsest "
                         "rung the sweep passes, which is evaluated on the FULL panel and can be "
                         "inadmissible for a trait's non-missing subset; the run REFUSES such a rung. "
                         "See wiki/f3_rung_census_1001g_*.json for the per-trait census.")
    ap.add_argument("--only-arm", default=None,
                    choices=["pos-fine", "pos-coarse", "neg-fine", "neg-coarse"],
                    help="compute only this arm/rung this call (chunked driving); the final "
                         "artifact is written only when every arm has its cells")
    ap.add_argument("--only-trait", default=None,
                    help="substring; restrict to traits whose name contains it")
    ap.add_argument("--max-traits", type=int, default=8,
                    help="per control arm; must clear MIN_TRAITS_FOR_SIGN_TEST (frozen) for a verdict")
    ap.add_argument("--out", type=Path,
                    default=REPO / "wiki" / f"f3_deconfound_calibration_1001g_{date.today().isoformat()}.json")
    ap.add_argument("--crosscheck-out", type=Path,
                    default=REPO / "wiki" / f"f3_ibs_crosscheck_1001g_{date.today().isoformat()}.json")
    a = ap.parse_args(argv)

    n_markers = 600 if a.smoke else a.n_markers
    n_perm = 3 if a.smoke else a.n_perm
    n_blocks = 8 if a.smoke else N_BLOCKS
    max_traits = 1 if a.smoke else a.max_traits

    max_largest = 1.0 / N_SPLITS
    print(f"[1/5] loading 1001G genotypes (blocks={n_blocks}, target markers={n_markers}) ...",
          flush=True)
    X, ids, geno_meta, repo_pg, coarse_pg, rung_sweep = cached_preamble(
        n_blocks=n_blocks, n_markers=n_markers, maf_floor=MAF_FLOOR, max_largest=max_largest)
    print(f"      X={X.shape}  {geno_meta}", flush=True)

    print("[2/5] loading phenotypes ...", flush=True)
    traits = load_1001g_traits(ids, max_traits=max_traits)
    if a.only_trait:
        traits = {k: v for k, v in traits.items() if a.only_trait in k}
        print(f"      --only-trait {a.only_trait!r} -> {list(traits)}", flush=True)
    print(f"      traits={ {k: int(np.sum(~np.isnan(v))) for k, v in traits.items()} }", flush=True)
    if not traits:
        print("REFUSING: no trait aligned to the genotype accessions")
        return 2

    print("[3/5] grouping -- repo IBS over raw genotypes AND the project kinship ...", flush=True)
    if a.coarse_percentile is not None:
        from dna_decode.eval.panel_structure import panel_groups
        coarse_pg = panel_groups(X, ids, percentile=a.coarse_percentile, max_largest=max_largest)
        print(f"      coarse end OVERRIDDEN to p={a.coarse_percentile:g} "
              f"({coarse_pg.n_clusters} clusters)", flush=True)
    repo_groups = {i: g for i, g in zip(repo_pg.ids, repo_pg.groups)}
    coarse_groups = {i: g for i, g in zip(coarse_pg.ids, coarse_pg.groups)}

    # REFUSE a bracket endpoint that is inadmissible for any trait. The sweep's own `passes` is
    # computed on the FULL panel; grouped CV folds each trait's non-missing subset. Without this
    # check the coarse arm silently scores only the traits that happen to clear the share, and its
    # median is then an attrition artifact rather than a measurement.
    adm_fine = trait_rung_admissibility(repo_groups, ids, traits, max_largest=max_largest)
    adm_coarse = trait_rung_admissibility(coarse_groups, ids, traits, max_largest=max_largest)
    bad_coarse = [nm for nm, v in adm_coarse.items() if not v["admissible"]]
    bad_fine = [nm for nm, v in adm_fine.items() if not v["admissible"]]
    print(f"      per-trait admissibility: fine p={repo_pg.percentile:g} "
          f"{len(traits) - len(bad_fine)}/{len(traits)} | coarse p={coarse_pg.percentile:g} "
          f"{len(traits) - len(bad_coarse)}/{len(traits)}", flush=True)
    if bad_fine:
        print(f"REFUSING: fine rung p={repo_pg.percentile:g} is inadmissible for {bad_fine}")
        return 2
    if bad_coarse and a.coarse_percentile is None:
        print(f"REFUSING: the sweep's coarsest passing rung p={coarse_pg.percentile:g} is "
              f"inadmissible for {bad_coarse} (their non-missing subsets exceed the "
              f"{max_largest:.0%} one-fold share). Pass --coarse-percentile with a rung the "
              f"per-trait census admits for every trait.")
        return 2
    if bad_coarse:
        print(f"REFUSING: --coarse-percentile {a.coarse_percentile:g} is inadmissible for "
              f"{bad_coarse}")
        return 2
    if a.only_arm and a.crosscheck_out.exists():
        # The cross-check is a one-off diagnostic and its artifact is already written; recomputing an
        # 1135x1135 kinship IBS on every chunked call would spend the call's budget on it.
        print(f"      [crosscheck skipped -- {a.crosscheck_out.name} already written]", flush=True)
        kin_groups, kin_meta = {}, {"skipped_this_call": True}
    else:
        kin_groups, kin_meta = kinship_grouping(ids, max_largest=max_largest)

    # CROSS-CHECK: do the two groupings agree? Guarded AS A WHOLE -- with the kinship load skipped the
    # common-accession set is empty, and an unguarded run would overwrite the already-correct artifact
    # with an all-zero one. Skipping a diagnostic is fine; silently zeroing it is not.
    if kin_meta.get("skipped_this_call"):
        print(f"      repo IBS: {repo_pg.n_clusters} clusters p={repo_pg.percentile:g} "
              f"largest={repo_pg.largest_cluster_fraction:.3f} | kinship crosscheck not recomputed",
              flush=True)
    else:
        write_crosscheck(a.crosscheck_out, ids, repo_pg, coarse_pg, rung_sweep,
                         repo_groups, kin_groups, kin_meta)

    want = (lambda nm: a.only_arm is None or a.only_arm == nm)

    print(f"[4/5] POSITIVE control -- 1001G, {len(traits)} trait(s), n_perm={n_perm}", flush=True)
    pos = pos_coarse = None
    if want("pos-fine"):
        print(f"      FINEST admissible rung p={repo_pg.percentile:g} "
              f"({repo_pg.n_clusters} clusters) ...", flush=True)
        pos = run_1001g_arm(X, ids, traits, [repo_groups[i] for i in ids], n_perm=n_perm,
                            label="1001g_finest_rung")
    if want("pos-coarse"):
        print(f"      COARSEST admissible rung p={coarse_pg.percentile:g} "
              f"({coarse_pg.n_clusters} clusters) ...", flush=True)
        pos_coarse = run_1001g_arm(X, ids, traits, [coarse_groups[i] for i in ids], n_perm=n_perm,
                                   label="1001g_coarsest_rung")

    print("[5/5] NEGATIVE control -- bloom (bi-parental cross), SAME machinery + marker budget ...",
          flush=True)
    neg_spec = next(s for s in PANELS if s.design == "constructed_biparental")
    bX, bids, btraits, bmeta = load_bloom(n_markers=n_markers)
    btraits = dict(list(btraits.items())[:max_traits])
    print(f"      bloom X={bX.shape} {bmeta}  traits={list(btraits)}", flush=True)
    # GRANULARITY-MATCHED to the positive arm, not percentile-matched. See matched_rung() for why:
    # the same percentile buys different fold geometry on two differently-structured panels, and the
    # drop a de-confounder produces depends on cluster size rather than on the ladder rung's name.
    bfine, bfine_match = matched_rung(
        bX, bids, target_clusters_per_individual=repo_pg.n_clusters / len(ids),
        max_largest=max_largest)
    bcoarse, bcoarse_match = matched_rung(
        bX, bids, target_clusters_per_individual=coarse_pg.n_clusters / len(ids),
        max_largest=max_largest)
    bsweep = bfine_match["candidates"]
    bgroups = {i: g for i, g in zip(bfine.ids, bfine.groups)}
    print(f"      granularity-matched fine  p={bfine.percentile:g} "
          f"({bfine.n_clusters} clusters, mean size {bfine_match['mean_cluster_size']}) "
          f"vs positive mean size {len(ids) / repo_pg.n_clusters:.2f}", flush=True)
    print(f"      granularity-matched coarse p={bcoarse.percentile:g} "
          f"({bcoarse.n_clusters} clusters, mean size {bcoarse_match['mean_cluster_size']}) "
          f"vs positive mean size {len(ids) / coarse_pg.n_clusters:.2f}", flush=True)
    badm_f = trait_rung_admissibility({i: g for i, g in zip(bfine.ids, bfine.groups)},
                                      bids, btraits, max_largest=max_largest)
    badm_c = trait_rung_admissibility({i: g for i, g in zip(bcoarse.ids, bcoarse.groups)},
                                      bids, btraits, max_largest=max_largest)
    bad_b = [nm for nm, v in list(badm_f.items()) + list(badm_c.items()) if not v["admissible"]]
    if bad_b:
        print(f"REFUSING: a granularity-matched bloom rung is inadmissible for {sorted(set(bad_b))}")
        return 2
    neg = neg_coarse = None
    if want("neg-fine"):
        neg = run_1001g_arm(bX, bids, btraits, [bgroups[i] for i in bids], n_perm=n_perm,
                            label="bloom_finest_rung")
        neg["design"] = "constructed_biparental"
        neg["genotype_meta"] = bmeta
        neg["grouping"] = bfine.as_dict()
    # The bracket must be SYMMETRIC: comparing the positive control at its COARSEST rung against the
    # negative at its FINEST would be apples to oranges, and the coarsest-vs-coarsest pair is the only
    # fair reading of the coarse end.
    bcgroups = {i: g for i, g in zip(bcoarse.ids, bcoarse.groups)}
    print(f"      COARSEST admissible rung p={bcoarse.percentile:g} "
          f"({bcoarse.n_clusters} clusters) ...", flush=True)
    if want("neg-coarse"):
        neg_coarse = run_1001g_arm(bX, bids, btraits, [bcgroups[i] for i in bids], n_perm=n_perm,
                                   label="bloom_coarsest_rung")
        neg_coarse["design"] = "constructed_biparental"
        neg_coarse["grouping"] = bcoarse.as_dict()

    if a.only_arm:
        # A single-arm call BANKS cells and stops. Composing an artifact from one arm would publish a
        # verdict the run has no evidence for; the full no-selector run replays every cell from the
        # checkpoint in seconds and composes it then.
        done = [nm for nm, arm in (("pos-fine", pos), ("pos-coarse", pos_coarse),
                                   ("neg-fine", neg), ("neg-coarse", neg_coarse)) if arm]
        print(f"\n--only-arm {a.only_arm}: banked cells for {done}; NO artifact written "
              f"(run without --only-arm to compose it).")
        return 0

    print(f"      n_traits={neg['n_traits']} median_drop(finest)={neg['median_drop']} "
          f"median_drop(coarsest)={neg_coarse['median_drop']}", flush=True)

    verdict = compute_verdict(n_intact=2, has_negative=True, has_positive=True,
                              positive_fold_compatible=True, neg=neg, pos=pos)
    out = {
        "schema": "f3-deconfound-calibration-1001g-v1",
        "generated": date.today().isoformat(),
        "verdict": verdict, "verdict_meaning": VERDICT_MEANING[verdict],
        "all_possible_verdicts": list(VERDICTS),
        "smoke": bool(a.smoke),
        "genotype_meta": geno_meta,
        "n_perm": n_perm, "n_splits": N_SPLITS,
        "control_selection": {
            "rule": "BY CONSTRUCTION from each panel's `design` field; NOT by pc1",
            "negative_control": neg_spec.name, "positive_control": "1001g_arabidopsis",
            "positive_basis": "by sampling: 1,135 natural Arabidopsis accessions (1001 Genomes)",
        },
        "criterion": ("DIFFERENTIAL: the negative control keeps beating its own grouped null AND the "
                      "positive control's paired per-trait drop exceeds the negative control's. "
                      "Asserts no magnitude."),
        "negative_control_result": neg,
        "positive_control_result": pos,
        "positive_control_result_coarsest_rung": pos_coarse,
        "negative_control_result_coarsest_rung": neg_coarse,
        "grouping_bracket": {
            "why": ("select_threshold takes the SMALLEST passing rung on a named bias -- it "
                    "under-merges, biasing AGAINST detecting a de-confounding response. Running BOTH "
                    "endpoints of the admissible set separates 'the instrument is insensitive' from "
                    "'the grouping was too fine to de-confound anything', which the "
                    "INVALID_INSENSITIVE verdict text alone cannot."),
            "finest": {"percentile": repo_pg.percentile, "n_clusters": repo_pg.n_clusters,
                       "largest": round(repo_pg.largest_cluster_fraction, 4),
                       "median_drop": pos["median_drop"], "n_traits": pos["n_traits"]},
            "coarsest": {"percentile": coarse_pg.percentile, "n_clusters": coarse_pg.n_clusters,
                         "largest": round(coarse_pg.largest_cluster_fraction, 4),
                         "median_drop": pos_coarse["median_drop"],
                         "n_traits": pos_coarse["n_traits"]},
            "per_trait_admissibility": {
                "bar": max_largest, "fine": adm_fine, "coarse": adm_coarse,
                "coarse_percentile_overridden_to": a.coarse_percentile,
                "why_overridden": (
                    "The sweep's `passes` flag is computed on the FULL 1,135-accession panel, but "
                    "grouped CV folds each trait's NON-MISSING subset. The coarsest rung the sweep "
                    "passes (p=20, largest 0.1947 on the full panel) is inadmissible for 6 of 8 "
                    "traits once missingness is accounted for -- FT10 reaches 0.2104 on its 1,003 "
                    "rows. The coarsest rung admissible for EVERY trait is p=10 (90 clusters); see "
                    "wiki/f3_rung_census_1001g_*.json. The run now REFUSES an inadmissible endpoint "
                    "rather than scoring whichever traits happen to clear it."),
            },
            "paired_fine_vs_coarse": {
                "positive": paired_rung_delta(pos, pos_coarse),
                "negative": paired_rung_delta(neg, neg_coarse),
                "why": ("A rung's median over the traits IT could score is NOT comparable to the other "
                        "rung's median over a DIFFERENT trait set. The coarse rung refuses a trait "
                        "whose NON-MISSING subset pushes a cluster over the one-fold share, and the "
                        "traits it refuses are the HIGHEST-n ones -- which on this panel are also the "
                        "strongest-signal, smallest-drop ones. Dropping them inflates the coarse "
                        "median through ATTRITION alone, so the fine-vs-coarse reading is restricted "
                        "to the traits scoreable at BOTH rungs."),
            },
            "granularity_matching": {
                "rule": ("the negative arm's rungs are chosen to match the positive arm's CLUSTERS "
                         "PER INDIVIDUAL, not its ladder percentile; the rule was fixed before any "
                         "bloom drop was read"),
                "negative_fine": bfine_match, "negative_coarse": bcoarse_match,
                "positive_fine_clusters_per_individual": round(repo_pg.n_clusters / len(ids), 6),
                "positive_coarse_clusters_per_individual": round(coarse_pg.n_clusters / len(ids), 6),
                "why_not_percentile": ("on 1001G p=10 gives 90 clusters over 1,135 accessions while "
                                       "on the bloom cross p=20 gives 22 over 1,008 -- a 12.6-member "
                                       "cluster against a 46-member one. Percentile matching would "
                                       "compare different operations and call it the same end of the "
                                       "bracket."),
                "note_on_range_width": ("the negative control's admissible range is WIDE (p=1..20, "
                                        "260..22 clusters) and the positive's is NARROW (p=5..10, "
                                        "157..90): a bi-parental cross forms no dominant cluster, so "
                                        "few rungs are excluded. The WIDTH of the admissible range is "
                                        "itself a structure signal, observed here and not used as a "
                                        "test."),
            },
            "symmetric_comparison": {
                "finest_positive_minus_negative": (
                    None if pos["median_drop"] is None or neg["median_drop"] is None
                    else round(pos["median_drop"] - neg["median_drop"], 4)),
                "coarsest_positive_minus_negative": (
                    None if pos_coarse["median_drop"] is None or neg_coarse["median_drop"] is None
                    else round(pos_coarse["median_drop"] - neg_coarse["median_drop"], 4)),
                "note": ("both arms are compared AT THE SAME END of their own admissible ranges; "
                         "a positive-coarsest vs negative-finest comparison would be apples to "
                         "oranges and is deliberately not reported as the differential."),
            },
            "reading": ("if BOTH endpoints show a drop no larger than the negative control's, the "
                        "insensitivity is the INSTRUMENT's; if the coarsest shows a clear drop and the "
                        "finest does not, the fineness was the cause and select_threshold's bias is "
                        "the binding constraint."),
        },
        "crosscheck_artifact": a.crosscheck_out.name,
        "honest_limits": [
            "The two controls have DIFFERENT traits, so the between-control comparison is two-sample, "
            "not paired; only the within-control sign tests are paired.",
            "The positive control's genotypes are IMPUTED, which can flatten or sharpen the very "
            "relatedness structure a de-confounding calibration measures.",
            f"Markers are a MAF>={MAF_FLOOR} subsample read from chunk-aligned blocks, not an "
            "LD-pruned genome-wide set.",
            f"n_perm={n_perm} bounds the resolution of every beats_null decision.",
            "MIN_TRAITS_FOR_SIGN_TEST is asserted, not derived.",
            "The positive control's traits are ONE PER STUDY (163 usable AraPheno traits span only 8 "
            "studies, 121 of them one rainfall factorial), which is the best independence this "
            "substrate allows -- but FT10 and DTF1 are both flowering time, so they remain "
            "BIOLOGICALLY correlated across studies.",
            "Both arms run at a comparable marker budget by construction (bloom is strided to match), "
            "so the differential is not confounded by marker count.",
            "A PASS bounds the INSTRUMENT's differential response. It does not validate any "
            "downstream organism-level prediction.",
            "The VERDICT is computed from the FINEST admissible rung (what select_threshold would "
            "actually deploy). The coarsest-rung arm is a DIAGNOSTIC that brackets the admissible "
            "range; it does not feed the verdict, so the bar is not being chosen after the fact.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print()
    print(f"VERDICT: {verdict}")
    print(f"  {VERDICT_MEANING[verdict]}")
    print(f"  negative median_drop  finest={neg['median_drop']}  coarsest={neg_coarse['median_drop']}")
    print(f"  positive median_drop  finest(p{repo_pg.percentile:g})={pos['median_drop']}  "
          f"coarsest(p{coarse_pg.percentile:g})={pos_coarse['median_drop']}")
    print("  SYMMETRIC differential (pos-neg):  "
          f"finest={None if pos['median_drop'] is None or neg['median_drop'] is None else round(pos['median_drop']-neg['median_drop'],4)}  "
          f"coarsest={None if pos_coarse['median_drop'] is None or neg_coarse['median_drop'] is None else round(pos_coarse['median_drop']-neg_coarse['median_drop'],4)}")
    print(f"-> {a.out}")
    print(f"-> {a.crosscheck_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
