"""F1 falsifier: can the transfer benchmark return a POSITIVE and a NEGATIVE on known inputs?

A benchmark that cannot do both is not a benchmark. Three arms, and only two of them can earn a PASS:

  A1  scoring-core control      Bloom BYxRM yeast. POSITIVE expected. `is_transfer_evidence=False` --
                                a single-cross segregant panel has ONE group, so nothing
                                organism-level is held out. It proves the continuous scoring core
                                recovers real signal, NOT that transfer is measurable.
  B   synthetic negative        g groups with a per-group phenotype offset and ZERO within-group
                                signal, through folds -> leakage audit -> gauntlet. The ONLY arm that
                                exercises the live machinery on a known negative.
  C   verdict-replay regression The committed Arabidopsis per-seed scores through `verdict`.
                                `counts_toward_pass=False`: it checks the verdict FUNCTION, not the
                                machinery, and letting it satisfy the negative side would be a
                                loophole -- the benchmark could pass by replaying a historical
                                conclusion while the live path was only synthetic-tested.

THERE IS NO TRANSFER-POSITIVE ARM, and that is the honest state rather than an omission. Every
genotype-phenotype panel on disk is a single constructed population, so a held-out-organism fold has
nothing to hold out. `transfer_positive` is emitted as the first-class value
TRANSFER_POSITIVE_UNMEASURED on every run.

THE CROSS-PANEL CONTROL IS DELIBERATELY ABSENT. An earlier draft offered yeast/mouse/Arabidopsis-MAGIC
together, barred from PASS. It is not built: its result is knowable without running it (three
mutually unmappable identifier schemes, no orthology map anywhere here), so "~0" restates the input
rather than measuring anything, and a number with no information content that can be quoted out of
context is a liability. A test asserts the flag does not exist.

    uv run python scripts/transfer_benchmark_falsifier.py
    uv run python scripts/transfer_benchmark_falsifier.py --offline-only

Exit 0 = both directions discriminated · 1 = bar not met · 3 = refused / INDETERMINATE.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# ---------------------------------------------------------------------------------------------------
# THE BAR, FROZEN ABOVE ANY DATA LOADING. Pre-committing the verdict-conditional response before
# results land is the single largest discipline against motivated reasoning; a test asserts this
# constant appears before the first read of any substrate.
# ---------------------------------------------------------------------------------------------------
PREREGISTERED = {
    "pass_requires": ["A1_positive", "B_negative"],
    "A1_positive": "cv_ridge_gp predictive r beats its own label-permutation null p95",
    "B_negative": "gauntlet verdict is NOT BEATS_ALL_CONTROLS on a zero-within-group-signal input",
    "C_counts_toward_pass": False,
    "C_requirement": "every replayed seed must be NON-positive DESPITE the embedding winning global r2",
    "transfer_positive": "TRANSFER_POSITIVE_UNMEASURED unless a live held-out-organism substrate is supplied",
    "frozen_at": "plan v3, before any arm was run",
}

TRANSFER_POSITIVE_UNMEASURED = "TRANSFER_POSITIVE_UNMEASURED"
INDETERMINATE = "INDETERMINATE"
DISCRIMINATED = "DISCRIMINATED_BOTH_DIRECTIONS"
BAR_NOT_MET = "BAR_NOT_MET"

#: geno_v2.txt, NOT BYxRM_GenoData.txt. The latter is a TRUNCATED download -- exactly 1,048,576 bytes
#: (1 MiB), 511 of 11,623 marker rows, and its last line carries 626 fields instead of 1009. The plan
#: cited it; building A1 against it would have loaded a ragged fifth of the panel and produced a
#: number that looks like a result. `geno_direct.txt` in the same directory is a 199-byte HTTP 403
#: page. Same corruption class as the recorded GPS-assembly cohort.
BLOOM_GENO = Path("D:/dna_decode_cache/bloom/geno_v2.txt")
BLOOM_PHENO = Path("D:/dna_decode_cache/bloom/BYxRM_PhenoData.txt")
G2_FIXTURE = ROOT / "wiki" / "phase2_arabidopsis_result_2026-06-12.json"


class ArmRefused(RuntimeError):
    """An arm could not run. Never a silent skip -- it forces INDETERMINATE."""


def assert_genotype_file_intact(path: Path) -> dict:
    """Refuse a truncated genotype matrix rather than quietly loading a fifth of it.

    The check is RAGGEDNESS, not size: a truncated tab matrix has a final row with fewer fields than
    its header. A size threshold would be a guess; a ragged row is proof.
    """
    if not path.exists():
        raise ArmRefused(f"genotype file missing: {path}")
    with open(path, encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        n_cols, n_rows, last = len(header), 0, []
        for line in f:
            if line.strip():
                n_rows += 1
                last = line.rstrip("\n").split("\t")
    if len(last) != n_cols:
        raise ArmRefused(
            f"{path.name} is TRUNCATED: header has {n_cols} fields but the last row has "
            f"{len(last)}. Size is {path.stat().st_size} bytes"
            f"{' (exactly 1 MiB -- a classic cut download)' if path.stat().st_size == 1048576 else ''}. "
            f"Refusing rather than scoring a ragged matrix.")
    return {"path": str(path), "n_cols": n_cols, "n_rows": n_rows,
            "size_bytes": path.stat().st_size, "intact": True}


# --- Arm A1 ----------------------------------------------------------------------------------------

def arm_a1_scoring_core(trait: str = "Maltose", marker_stride: int = 16) -> dict:
    """Bloom BYxRM through cv_ridge_gp + its label-permutation null.

    NOT routed through the gauntlet, and the reason is structural rather than a shortcut: the gauntlet
    is group-based and a single-cross panel has ONE group, so `cv_r2(groups=...)` would skip it
    (train side empty) and return nan. The honest scoring-core control is K-fold with a permutation
    null, which is what `scripts/yeast_bloom_gp_arm.py` already does.
    """
    from dna_decode.eval.genomic_prediction import cv_ridge_gp
    from scripts.yeast_bloom_gp_arm import load_genotype, load_phenotype

    integrity = assert_genotype_file_intact(BLOOM_GENO)
    if not BLOOM_PHENO.exists():
        raise ArmRefused(f"phenotype file missing: {BLOOM_PHENO}")

    seg_g, _, G = load_genotype(str(BLOOM_GENO), marker_stride)
    seg_p, traits, P = load_phenotype(str(BLOOM_PHENO))
    if trait not in traits:
        raise ArmRefused(f"trait {trait!r} not in the phenotype file; available: {traits[:6]} ...")

    idx_p = {s: i for i, s in enumerate(seg_p)}
    shared = [s for s in seg_g if s in idx_p]
    if len(shared) < 100:
        raise ArmRefused(f"only {len(shared)} segregants align between genotype and phenotype")
    gi = {s: i for i, s in enumerate(seg_g)}
    X = G[[gi[s] for s in shared], :]
    y = P[[idx_p[s] for s in shared], traits.index(trait)]

    res = cv_ridge_gp(X, y, trait=trait)
    d = res.as_dict()
    return {"arm": "A1_scoring_core_control", "is_transfer_evidence": False,
            "positive": bool(res.beats_null), "trait": trait, "marker_stride": marker_stride,
            "n_aligned": len(shared), "genotype_integrity": integrity, "result": d,
            "note": "a single-cross panel holds out NO organism; this proves the scoring core, "
                    "not that transfer is measurable"}


# --- Arm B -----------------------------------------------------------------------------------------

def arm_b_synthetic_negative(n_groups: int = 4, per_group: int = 80, seed: int = 0) -> dict:
    """Pure population structure: per-group offsets, ZERO within-group signal.

    The fixture shape is `tests/test_deconfound_package.py:23` with the signal term removed -- the
    same construction the primitive's own test uses to prove de-confounding works.
    """
    import numpy as np

    from dna_decode.eval.transfer import held_out_group_folds, leakage_audit, require_clean_audit
    from dna_decode.eval.transfer_gauntlet import BEATS_ALL_CONTROLS, run_gauntlet, verdict

    rng = np.random.default_rng(seed)
    n = n_groups * per_group
    groups = np.repeat([f"g{i}" for i in range(n_groups)], per_group)
    ids = [f"syn{i}" for i in range(n)]
    X = rng.integers(0, 2, (n, 3)).astype(float)
    y = np.repeat(rng.normal(0, 5, n_groups), per_group) + rng.normal(0, 0.3, n)   # NO signal term

    folds = held_out_group_folds({ids[i]: groups[i] for i in range(n)})
    audits = []
    for f in folds:
        a = leakage_audit(f.train_ids, f.test_ids)
        audits.append(require_clean_audit(a))

    rep = run_gauntlet(X, y, groups, ids=ids, n_draws=100)
    rep.assert_all_controls_ran()
    rep.assert_same_evaluation_basis()
    v = verdict(rep.matrix, rep)
    return {"arm": "B_synthetic_structure_negative", "machinery_exercised": True,
            "negative": v != BEATS_ALL_CONTROLS, "verdict": v,
            "n_folds": len(folds), "n_leakage_audits_clean": len(audits),
            "report": rep.as_dict()}


# --- Arm C -----------------------------------------------------------------------------------------

def arm_c_verdict_replay() -> dict:
    """Replay the committed Arabidopsis scores through `verdict`. Cannot earn a PASS."""
    from dna_decode.eval.transfer_gauntlet import (BEATS_ALL_CONTROLS, CONTROL_NAMES,
                                                   GauntletReport, MetricMatrix, verdict)

    if not G2_FIXTURE.exists():
        raise ArmRefused(f"replay fixture missing: {G2_FIXTURE}")
    fx = json.loads(G2_FIXTURE.read_text(encoding="utf-8"))
    m = fx["metrics"]
    per_seed = []
    for i, seed in enumerate(fx["seeds"]):
        matrix = MetricMatrix(
            zero_shot_group_r2=m["global_r2"]["embedding"][i],
            kshot_pooled_r2=m["global_r2"]["embedding"][i],
            kshot_within_group_r2=m["within_group_r2"]["embedding"][i],
            kshot_group_offset_baseline=m["global_r2"]["structure_only"],
            n_groups_used_within=1,            # the run did produce a within-group number
        )
        rep = GauntletReport(controls={n: 0.0 for n in CONTROL_NAMES}
                             | {"structure_only_baseline": m["within_group_r2"]["structure_only"]})
        v = verdict(matrix, rep)
        per_seed.append({"seed": seed, "verdict": v, "non_positive": v != BEATS_ALL_CONTROLS,
                         "within_group_r2": matrix.kshot_within_group_r2,
                         "won_global_r2": matrix.zero_shot_group_r2 > m["global_r2"]["structure_only"]})
    return {"arm": "C_verdict_replay_regression", "mode": "verdict_replay",
            "machinery_exercised": False, "counts_toward_pass": False,
            "all_seeds_non_positive": all(s["non_positive"] for s in per_seed),
            "won_global_r2_on_every_seed": all(s["won_global_r2"] for s in per_seed),
            "per_seed": per_seed, "source": fx["source"]}


# --- driver ----------------------------------------------------------------------------------------

def run(offline_only: bool = False, trait: str = "Maltose", marker_stride: int = 16) -> dict:
    arms, refusals = {}, []

    if offline_only:
        refusals.append("A1_scoring_core_control: skipped (--offline-only)")
    else:
        try:
            arms["A1"] = arm_a1_scoring_core(trait, marker_stride)
        except ArmRefused as exc:
            refusals.append(f"A1_scoring_core_control: {exc}")
    for key, fn in (("B", arm_b_synthetic_negative), ("C", arm_c_verdict_replay)):
        try:
            arms[key] = fn()
        except ArmRefused as exc:
            refusals.append(f"{key}: {exc}")

    # PASS depends on A1 and B ONLY. Arm C must not be load-bearing for a pass -- otherwise the
    # benchmark could be blocked (or carried) by a regression fixture rather than by the machinery.
    #
    # The asymmetry is deliberate and is NOT the same thing twice:
    #   C could not RUN (fixture missing)  -> recorded in `refusals`, does NOT block a pass
    #   C RAN and came out positive        -> a verdict-function regression, forces BAR_NOT_MET
    # "cannot run" and "ran and broke" are different facts and must not share an outcome.
    a1_ok = bool(arms.get("A1", {}).get("positive"))
    b_ok = bool(arms.get("B", {}).get("negative"))
    c_regressed = "C" in arms and not arms["C"].get("all_seeds_non_positive")
    blocking = [r for r in refusals if not r.startswith("C:")]

    if blocking or "A1" not in arms or "B" not in arms:
        status, exit_code = INDETERMINATE, 3
    elif a1_ok and b_ok and not c_regressed:
        status, exit_code = DISCRIMINATED, 0
    else:
        status, exit_code = BAR_NOT_MET, 1

    return {
        "schema": "transfer-benchmark-falsifier-v1",
        "generated": date.today().isoformat(),
        "preregistered": PREREGISTERED,
        "status": status,
        "exit_code": exit_code,
        "transfer_positive": TRANSFER_POSITIVE_UNMEASURED,
        "transfer_positive_reason": ("every genotype-phenotype panel on disk is a single constructed "
                                     "population, so a held-out-organism fold has nothing to hold "
                                     "out. This is a measured absence, not an omission."),
        "arms": arms,
        "refusals": refusals,
    }


def render_md(res: dict) -> str:
    L = [f"# Transfer-benchmark falsifier — {res['generated']}", "",
         f"**Status:** `{res['status']}` (exit {res['exit_code']})", "",
         f"**transfer_positive:** `{res['transfer_positive']}`", "",
         res["transfer_positive_reason"], "", "## Arms", ""]
    for k, a in res["arms"].items():
        L.append(f"- **{k} — {a['arm']}**")
        for field in ("positive", "negative", "verdict", "is_transfer_evidence",
                      "machinery_exercised", "counts_toward_pass", "all_seeds_non_positive",
                      "won_global_r2_on_every_seed"):
            if field in a:
                L.append(f"  - `{field}`: {a[field]}")
    if res["refusals"]:
        L += ["", "## Refusals", ""] + [f"- {r}" for r in res["refusals"]]
    L += ["", "## Pre-registered bar", "", "```json",
          json.dumps(res["preregistered"], indent=2), "```", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="F1 transfer-benchmark falsifier.")
    ap.add_argument("--offline-only", action="store_true",
                    help="skip Arm A1 (needs D:); forces INDETERMINATE")
    ap.add_argument("--trait", default="Maltose")
    ap.add_argument("--marker-stride", type=int, default=16)
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args(argv)

    res = run(args.offline_only, args.trait, args.marker_stride)
    out = Path(args.out_dir) if args.out_dir else (ROOT / "wiki")
    stem = out / f"transfer_benchmark_falsifier_{res['generated']}"
    out.mkdir(parents=True, exist_ok=True)
    stem.with_suffix(".json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    stem.with_suffix(".md").write_text(render_md(res), encoding="utf-8")

    print(f"status: {res['status']}  (exit {res['exit_code']})")
    print(f"transfer_positive: {res['transfer_positive']}")
    for k, a in res["arms"].items():
        flag = a.get("positive", a.get("negative", a.get("all_seeds_non_positive")))
        print(f"  {k:3s} {a['arm']:34s} -> {flag}  {a.get('verdict', '')}")
    for r in res["refusals"]:
        print(f"  REFUSED {r}")
    print(f"wrote {stem.with_suffix('.json').name} + .md")
    return res["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())
