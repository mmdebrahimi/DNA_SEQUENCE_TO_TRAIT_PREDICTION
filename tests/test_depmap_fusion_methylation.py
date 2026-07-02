"""Offline end-to-end pin for the DepMap fusion + methylation modalities on synthetic CCLE-shaped files.

Builds minimal PRISM + fusion + RRBS-methylation fixtures with a PLANTED within-lineage fusion->sensitize and
methylation->sensitize signal, then runs the real loaders + run() — no network, no D: data.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.depmap_fusion_methylation import (  # noqa: E402
    fusion_presence,
    methylation_genes,
    run,
    _direction_ok,
)

LINEAGES = ["lung", "skin", "colon"]
N_PER = 16
CELLS = [(f"ACH-{i:06d}", f"CELL{i}_{LINEAGES[i % 3].upper()}", LINEAGES[i % 3])
         for i in range(3 * N_PER)]


def _write_fixture(d: Path, rng):
    d.mkdir(parents=True, exist_ok=True)
    dm = [c[0] for c in CELLS]; ccle = [c[1] for c in CELLS]; tis = [c[2] for c in CELLS]
    rows = [f"r_{i}" for i in range(len(CELLS))]
    # cells + treat
    pd.DataFrame({"row_name": rows, "depmap_id": dm, "ccle_name": ccle, "primary_tissue": tis}).to_csv(
        d / "prism_cells.csv", index=False)
    pd.DataFrame({"column_name": ["col_criz", "col_tmz"], "name": ["crizotinib", "temozolomide"]}).to_csv(
        d / "prism_treat.csv", index=False)
    # plant ALK fusion in ~4 cells per lineage; those get LOWER crizotinib LFC (sensitize)
    alk = set()
    for L in range(3):
        idx = [i for i in range(len(CELLS)) if i % 3 == L][:4]
        alk.update(idx)
    base = np.repeat(rng.normal(0, 2.0, 3), N_PER)                 # between-lineage structure
    criz = base + rng.normal(0, 0.2, len(CELLS))
    for i in alk:
        criz[i] -= 2.5                                             # fusion+ => sensitive
    # plant MGMT methylation continuous; within lineage high methyl -> low TMZ LFC (sensitize)
    methyl = np.clip(rng.uniform(0, 1, len(CELLS)), 0, 1)
    base2 = np.repeat(rng.normal(0, 2.0, 3), N_PER)
    tmz = base2 - 1.5 * methyl + rng.normal(0, 0.15, len(CELLS))
    pd.DataFrame({"col_criz": criz, "col_tmz": tmz}, index=rows).to_csv(d / "prism_lfc.csv")
    # mutations (empty-ish; a couple non-target so gene_matrix is non-trivial but ALK/MGMT absent)
    pd.DataFrame({"Hugo_Symbol": ["TP53", "KRAS"], "Tumor_Sample_Barcode": [dm[0], dm[1]],
                  "Variant_Classification": ["Missense", "Missense"],
                  "isDeleterious": ["True", "True"], "isCOSMIChotspot": ["False", "True"]}).to_csv(
        d / "ccle_mutations.csv", index=False)
    # fusions
    fu = pd.DataFrame({"Unnamed: 0": range(len(alk)), "DepMap_ID": [dm[i] for i in sorted(alk)],
                       "#FusionName": ["EML4--ALK"] * len(alk)})
    fu.to_csv(d / "ccle_fusions.csv", index=False)
    # RRBS methylation: TSS_id, gene, chr, fpos, tpos, strand, avg_coverage, then cell cols (ccle names)
    meta = ["TSS_id", "gene", "chr", "fpos", "tpos", "strand", "avg_coverage"]
    hdr = meta + ccle
    lines = ["\t".join(hdr)]
    mvals = ["MGMT_1", "MGMT", "chr10", "1", "2", "+", "30"] + [f"{v:.4f}" for v in methyl]
    other = ["FOO_1", "FOO", "chr1", "1", "2", "+", "30"] + ["0.5"] * len(CELLS)
    lines += ["\t".join(mvals), "\t".join(other)]
    (d / "ccle_rrbs_tss_1kb.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return alk


def test_direction_helper():
    assert _direction_ok(-0.5, "sensitize") and not _direction_ok(0.5, "sensitize")
    assert _direction_ok(0.5, "resist") and not _direction_ok(-0.5, "resist")


def test_fusion_presence_binary(tmp_path):
    rng = np.random.default_rng(0)
    alk = _write_fixture(tmp_path, rng)
    keep = {c[0] for c in CELLS}
    fus = fusion_presence(tmp_path, keep)
    assert fus["--ALK"].sum() == len(alk)                         # exactly the planted fusion+ cells
    assert set(fus.index[fus["--ALK"] == 1]) == {CELLS[i][0] for i in alk}


def test_methylation_bridge_and_average(tmp_path):
    rng = np.random.default_rng(1)
    _write_fixture(tmp_path, rng)
    keep = {c[0] for c in CELLS}
    ccle2dm = {c[1]: c[0] for c in CELLS}
    meth = methylation_genes(tmp_path, ccle2dm, keep)
    assert "MGMT" in meth.index
    # columns bridged to depmap_id, values in 0..1
    assert set(meth.columns) <= keep
    assert meth.loc["MGMT"].between(0, 1).all()


def test_run_recovers_both_planted_signals(tmp_path):
    rng = np.random.default_rng(2)
    _write_fixture(tmp_path, rng)
    res = run(tmp_path)
    assert res["n_cells"] == 3 * N_PER
    # fusion: crizotinib x ALK planted sensitize -> negative within-lineage t, feature_match True
    fc = [c for c in res["fusion_cases"] if c["drug"] == "crizotinib"]
    assert fc and fc[0]["fusion_within_lineage_t"] < -1.0 and fc[0]["feature_match"]
    # methylation: temozolomide x MGMT planted sensitize -> negative within-lineage rho, feature_match True
    mc = [c for c in res["methylation_cases"] if c["drug"] == "temozolomide"]
    assert mc and mc[0]["methylation_within_lineage_rho"] < -0.1 and mc[0]["feature_match"]
