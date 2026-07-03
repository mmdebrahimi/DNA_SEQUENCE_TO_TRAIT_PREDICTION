"""Offline pin for the DGRP learned-path decoder: VCF/phenotype parsing + the de-confounding wiring."""
import gzip
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.dgrp_learned_decoder import (  # noqa: E402
    _gt_dosage,
    build_snp_matrix,
    load_phenotype,
    run,
)


def test_gt_dosage():
    assert _gt_dosage("0/0:x") == 0.0 and _gt_dosage("1/1") == 1.0
    assert _gt_dosage("0/1") == 0.5 and np.isnan(_gt_dosage("./."))


def test_load_phenotype(tmp_path):
    p = tmp_path / "ph.json"
    p.write_text(json.dumps({"original_data": {"DGRP_101": {"F": [40, 60], "M": [20, 40]},
                                                "line_2": {"F": [10], "M": [None]}}}), encoding="utf-8")
    ph = load_phenotype(p)
    assert ph["101"] == 40.0 and ph["2"] == 10.0


def _synth_vcf(path: Path, n_lines=60, n_snps=400, rng=None):
    rng = rng or np.random.default_rng(0)
    samples = [f"DGRP-{i:03d}" for i in range(n_lines)]
    hdr = ["##fileformat=VCFv4.2", "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                              "FILTER", "INFO", "FORMAT"] + samples)]
    rows = []
    for j in range(n_snps):
        af = rng.uniform(0.1, 0.9)
        gts = ["1/1" if rng.random() < af else "0/0" for _ in range(n_lines)]
        rows.append("\t".join(["2L", str(1000 + j), ".", "A", "T", ".", "PASS", ".", "GT"] + gts))
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write("\n".join(hdr + rows) + "\n")
    return samples


def test_build_matrix_and_run(tmp_path):
    rng = np.random.default_rng(1)
    vcf = tmp_path / "v.vcf.gz"
    samples = _synth_vcf(vcf, n_lines=60, n_snps=400, rng=rng)
    npz = tmp_path / "m.npz"
    info = build_snp_matrix(vcf, npz, keep_every=1, maf_min=0.05, callrate_min=0.5)
    assert info["n_lines"] == 60 and info["n_snps_kept"] > 100
    # pure-structure phenotype (random per line, no genotype signal) -> within-clade r2 ~ 0 (negative verdict)
    ph = {str(int(s.split("-")[1])): float(rng.normal()) for s in samples}
    pj = tmp_path / "ph.json"
    pj.write_text(json.dumps({"original_data": {f"DGRP_{k}": {"F": [v]} for k, v in ph.items()}}),
                  encoding="utf-8")
    res = run(npz, pj, k_clades=4)
    assert res["n_lines"] == 60 and "naive_cv_r2" in res
    assert res["within_clade_r2"] <= 0.3          # no real within-clade genotype signal
    assert res["verdict"] in ("LEARNED_PATH_NEGATIVE_UNDER_DECONFOUNDING", "WITHIN_CLADE_SIGNAL_PRESENT")
    assert "discrete_structure_degenerate" in res
