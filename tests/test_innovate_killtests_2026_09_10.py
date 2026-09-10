"""Executable kill-tests for the 2026-09-10 /innovate sweep.

POLARITY IS INVERTED FROM NORMAL TESTS AND THAT IS DELIBERATE. Each test here is the DISPROOF of one
candidate claim: it PASSES (exit 0) exactly when disproving evidence IS found, which the falsification
engine reads as `killed`. A test that FAILS means the disproof was not found and the claim SURVIVED its
kill-test. Read every name as "disproof of ...".

They are pytest functions because the engine's gated runner allowlists `pytest` and fail-closes on
`python -c` / bare script invocation -- a `test-exit-0` kill-test that is not a pytest node comes back
`unfalsified` without ever executing, which looks like a considered verdict and is not one.

These are kept in-tree as permanent probes: each one encodes a question about the repo whose answer
could change, and a silent flip is exactly what nobody would otherwise notice.
"""
from __future__ import annotations

import glob
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def test_disproof_of_ar_bank_onesided_is_novel():
    """CLAIM: reading the single-class AR Bank cohorts one-sided is an unexploited move.
    DISPROOF: a committed memo already argues exactly that framing."""
    memo = ROOT / "wiki" / "ar_isolate_bank_external_validation_2026-07-18.md"
    assert memo.exists(), "memo absent -> no disproof"
    txt = memo.read_text(encoding="utf-8", errors="replace").lower()
    assert "one-sided" in txt and "one class is powered" in txt, (
        "the memo does not argue the one-sided framing -> claim survives")


def test_disproof_of_ar_bank_results_never_reach_the_trust_surface():
    """CLAIM: the AR Bank one-sided numbers exist but never reach the standing trust surface.
    DISPROOF: an AR Bank cell appears on the report card."""
    card = ROOT / "wiki" / "decoder_validation_report_card.json"
    assert card.exists(), "report card absent -> no disproof"
    assert "ar_bank" in card.read_text(encoding="utf-8", errors="replace").lower(), (
        "no ar_bank cell on the report card -> claim survives")


def test_disproof_of_doubt_layer_unreachable_for_bacterial_target_sites():
    """CLAIM: the L2 completeness screen is structurally unreachable for bacterial target-site cells,
    so a substitution at an uncatalogued codon gets a false clean bill.
    DISPROOF: a completeness signal is actually produced for a bacterial target-site drug."""
    from dna_decode.eval.doubt import target_site_doubt
    block = target_site_doubt("ciprofloxacin", {"gyrA": {"D86N"}})
    kinds = [s.kind for s in block.signals]
    assert any("complet" in k for k in kinds), (
        f"only {kinds} produced -> no completeness signal -> claim survives")


def test_disproof_of_source_concentration_never_applied_to_ar_bank():
    """CLAIM: the project's own source-concentration bar was applied to the NCBI-PD arm and never to
    the AR Bank arm, whose cohorts are single curated deposits.
    DISPROOF: some AR Bank artifact carries a source-diversity field."""
    paths = glob.glob(str(ROOT / "wiki" / "ar_bank_*validation*.json")) + \
        glob.glob(str(ROOT / "wiki" / "external_validation_ar_bank_*.json"))
    assert paths, "no AR Bank artifacts -> cannot disprove"
    fields = ("source_concentration", "bioprojects", "n_bioprojects", "largest_source_share")
    carriers = [p for p in paths
                if any(f in json.loads(Path(p).read_text(encoding="utf-8")) for f in fields)]
    assert carriers, f"0 of {len(paths)} AR Bank artifacts carry a source-diversity field -> claim survives"


def test_disproof_of_lineage_collapse_is_mechanism_dependent():
    """CLAIM (abduction): the real discriminator is causal-variant<->clade coupling, not population
    design -- so chromosomal-mechanism cipro should lose MORE sensitivity to lineage collapse than
    the horizontally-acquired drugs.
    DISPROOF: acquired-gene cells collapse at least as much as the chromosomal one."""
    d = json.loads((ROOT / "wiki" / "provdisjoint_lineage_metrics.json").read_text(encoding="utf-8"))
    cells = d["cells"]

    def drop(c):
        w = c.get("thresholds", {}).get("0.001", {}).get("cluster_weighted", {})
        if c.get("raw_sens") is None or w.get("sens") is None:
            return None
        return c["raw_sens"] - w["sens"]

    chrom = [v for c in cells if c.get("drug") == "ciprofloxacin" and (v := drop(c)) is not None]
    acq = [v for c in cells if c.get("drug") != "ciprofloxacin" and (v := drop(c)) is not None]
    assert chrom and acq, "insufficient cells to compare -> cannot disprove"
    assert sum(acq) / len(acq) >= sum(chrom) / len(chrom), (
        f"acquired mean drop {sum(acq)/len(acq):.3f} < chromosomal {sum(chrom)/len(chrom):.3f} "
        f"-> mechanism split holds -> claim survives")


def test_disproof_of_errors_are_isolate_clustered():
    """CLAIM: AR Bank decoder errors cluster on ISOLATES rather than being drug-independent, which
    would implicate assembly/sample quality rather than N independent catalog gaps.
    DISPROOF: no BioSample is mis-called under more than one cohort."""
    errs = defaultdict(set)
    for f in glob.glob(str(ROOT / "data" / "raw" / "ar_bank_*" / "predictions_strict.json")):
        cohort = Path(f).parent.name
        for r in json.loads(Path(f).read_text(encoding="utf-8")):
            if r.get("prediction") in ("R", "S") and r["prediction"] != r.get("label"):
                errs[r["biosample"]].add(cohort)
    assert errs, "no errors at all -> cannot disprove"
    assert not any(len(v) > 1 for v in errs.values()), (
        f"multi-cohort error isolates exist: "
        f"{ {k: sorted(v) for k, v in errs.items() if len(v) > 1} } -> claim survives")


def test_disproof_of_pooled_embedding_is_ancestry_dominated():
    """CLAIM (abduction): mean-pooling dilutes causal signal to ~1/N_genes so the pooled vector's
    leading direction is ancestry by construction -- locus resolution, not population design, is the
    necessary condition.
    DISPROOF: no single direction dominates the pooled representation."""
    h5py = __import__("h5py")
    import numpy as np
    cache = ROOT / "data" / "processed" / "mini_cipro_nt_cache.h5"
    assert cache.exists(), "cache absent -> cannot disprove"
    with h5py.File(cache, "r") as fh:
        grp = fh["strains"] if "strains" in fh else fh
        strains = list(grp)[:24]
        X = np.array([np.mean([grp[s][g][()] for g in grp[s]], axis=0) for s in strains])
    X = X - X.mean(0)
    var = np.linalg.svd(X, compute_uv=False) ** 2
    frac = float(var[0] / var.sum())
    assert frac < 0.5, f"PC1 explains {frac:.3f} of pooled variance -> ancestry-dominated -> claim survives"
