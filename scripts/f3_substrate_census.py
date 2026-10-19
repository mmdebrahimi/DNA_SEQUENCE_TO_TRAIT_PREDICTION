"""F3 Step 1 (W0 gate): substrate INTEGRITY + structure census over the cached panels.

Read-only. Exit 0 ALWAYS -- this is a census, not a gate on the suite. The consumer
(`f3_deconfound_calibration.py`) is what refuses on a census with no intact panel.

WHY INTEGRITY AND NOT PRESENCE. `file exists` is not `file intact`, measured in this repo:
23 cached GPS assemblies passed an existence check while 11 were 199-byte HTTP 403 pages and
10 were TRUNCATED gzip, all under genome filenames. Fed to a tool they give a generic failure
indistinguishable from a real negative result. A HEAD-ONLY read is the WRONG guard and was
written first that time: a truncated gzip decompresses its opening blocks perfectly and raises
only at the END, so reading 2048 bytes passes exactly the files that later fail.

Paths come from the OWNING ARM MODULE (`_GENO` / `_PHENO`), and the genotype is parsed by that
arm's OWN `load_genotype`, so the matrix this census describes is the matrix the arm would use.
`tests/test_f3_substrate_census.py` pins the dgrp paths against its module source, since that
one exposes them only as argparse defaults.

Run: uv run python scripts/f3_substrate_census.py
"""
from __future__ import annotations

import argparse
import gzip
import importlib
import json
import sys
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

HEAD_ONLY_BYTES = 2048          # the WRONG guard, kept to pin that it is wrong


@dataclass
class PanelSpec:
    name: str
    arm_module: str
    kind: str                                  # "arm_loader" | "npz"
    design: str                                # "constructed_biparental" | "constructed_multiparent" | "natural_population"
    design_basis: str                          # how `design` is known -- measured or by construction
    npz_geno: str | None = None
    npz_pheno: str | None = None
    notes: list[str] = field(default_factory=list)


# `design` is NOT a measurement: it is a property of how each panel was built, which is what
# makes the F3 negative control's claim ("ancestry randomized BY CONSTRUCTION") stronger than
# the positive control's ("structure MEASURED"). Recorded here so the asymmetry is explicit.
PANELS: list[PanelSpec] = [
    PanelSpec("bloom", "scripts.yeast_bloom_gp_arm", "arm_loader",
              design="constructed_biparental",
              design_basis="by construction: a single BY x RM cross, recombination randomizes "
                           "ancestry; this is the F3 NEGATIVE control",
              notes=["markers x individuals in the file; the arm's loader transposes"]),
    PanelSpec("arabmagic", "scripts.arabmagic_gp_arm", "arm_loader",
              design="constructed_multiparent",
              design_basis="by construction: a MAGIC multi-parent intercross. Residual FAMILY "
                           "structure, NOT natural population structure -- any positive-control "
                           "role must be scoped accordingly"),
    PanelSpec("bxd", "scripts.bxd_gp_arm", "arm_loader",
              design="constructed_multiparent",
              design_basis="by construction: BXD recombinant inbred lines from two founders, "
                           "with known family structure"),
    PanelSpec("dgrp", "scripts.dgrp_learned_decoder", "npz",
              design="natural_population",
              design_basis="by sampling: inbred lines drawn from ONE natural Drosophila "
                           "population (Raleigh), with documented cryptic relatedness",
              npz_geno="D:/dna_decode_cache/dgrp/dgrp_snp_matrix.npz",
              npz_pheno="D:/dna_decode_cache/dgrp/pheno_2798.json",
              notes=["paths are argparse defaults in the arm, not module constants; pinned by test"]),
]


# ----------------------------------------------------------------- raw integrity (pure-ish)

def verify_gzip_fully(path: Path) -> tuple[bool, str | None, int]:
    """Decompress the WHOLE member. Returns (ok, reason, n_bytes).

    A truncated gzip raises only at the end, so this must read to EOF -- see the module
    docstring. `HEAD_ONLY_BYTES` exists solely so a test can pin that the cheap check passes
    files this one rejects.
    """
    try:
        n = 0
        with gzip.open(path, "rb") as fh:
            while True:
                chunk = fh.read(1 << 20)
                if not chunk:
                    break
                n += len(chunk)
        return True, None, n
    except Exception as e:                            # noqa: BLE001 - any decompress failure
        return False, f"{type(e).__name__}: {e}", 0


def head_only_gzip_ok(path: Path, n_bytes: int = HEAD_ONLY_BYTES) -> bool:
    """The WRONG guard. Present to be refuted by test, never called by the census."""
    try:
        with gzip.open(path, "rb") as fh:
            fh.read(n_bytes)
        return True
    except Exception:                                 # noqa: BLE001
        return False


def raw_text_integrity(path: Path, *, comment: str = "#") -> dict:
    """Parse a delimited text file IN FULL and describe its shape.

    Records the field-width HISTOGRAM rather than just "ragged or not", so a single stray row
    is distinguishable from a file with two different shapes.
    """
    out: dict = {"path": str(path), "exists": path.exists()}
    if not path.exists():
        out["verdict"] = "unusable:absent"
        return out
    out["size_bytes"] = path.stat().st_size
    if out["size_bytes"] == 0:
        out["verdict"] = "unusable:empty_file"
        return out

    if path.suffix == ".gz":
        ok, reason, n = verify_gzip_fully(path)
        out["gzip_fully_decompressed"] = ok
        out["decompressed_bytes"] = n
        if not ok:
            out["verdict"] = f"unusable:truncated_or_corrupt_gzip ({reason})"
            return out
        text = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    else:
        text = path.read_text(encoding="utf-8", errors="replace")

    lines = text.splitlines()
    out["n_lines_total"] = len(lines)
    data = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith(comment)]
    out["n_comment_lines"] = sum(1 for ln in lines if ln.lstrip().startswith(comment))
    out["n_data_lines"] = len(data)
    if not data:
        out["verdict"] = "unusable:no_data_rows"
        return out

    delim = "\t" if data[0].count("\t") >= data[0].count(",") else ","
    out["delimiter"] = "tab" if delim == "\t" else "comma"
    widths: dict[int, int] = {}
    for ln in data:
        w = len(ln.split(delim))
        widths[w] = widths.get(w, 0) + 1
    out["field_width_histogram"] = {str(k): v for k, v in sorted(widths.items())}
    out["field_widths_agree"] = len(widths) == 1
    if len(widths) > 1:
        out["verdict"] = f"unusable:ragged_rows (widths {sorted(widths)})"
        return out

    alphabet = sorted({tok for ln in data[1:200] for tok in ln.split(delim)[1:]})
    out["value_alphabet_sample"] = alphabet[:12]
    out["n_distinct_values_sampled"] = len(alphabet)
    out["verdict"] = "intact"
    return out


def matrix_integrity(X: np.ndarray) -> dict:
    """Shape + missingness of the loaded matrix, and the two refusal conditions a census can
    see without clustering."""
    out = {
        "n_individuals": int(X.shape[0]),
        "n_markers": int(X.shape[1]),
        "missing_fraction": float(np.mean(~np.isfinite(X))),
    }
    col_all_missing = int(np.sum(np.all(~np.isfinite(X), axis=0)))
    out["n_markers_entirely_missing"] = col_all_missing
    row_all_missing = int(np.sum(np.all(~np.isfinite(X), axis=1)))
    out["n_individuals_entirely_missing"] = row_all_missing
    if out["n_individuals"] < 2:
        out["verdict"] = "unusable:fewer_than_2_individuals"
    elif out["n_markers"] < 1:
        out["verdict"] = "unusable:no_markers"
    elif row_all_missing == out["n_individuals"]:
        out["verdict"] = "unusable:every_individual_entirely_missing"
    else:
        out["verdict"] = "intact"
    return out


# ----------------------------------------------------------------- per-panel census

def census_panel(spec: PanelSpec) -> dict:
    """Never raises: a census must describe every panel, including the broken ones."""
    from dna_decode.eval.panel_structure import structure_census

    row: dict = {"panel": spec.name, "arm_module": spec.arm_module, "kind": spec.kind,
                 "design": spec.design, "design_basis": spec.design_basis,
                 "notes": list(spec.notes)}
    try:
        mod = importlib.import_module(spec.arm_module)
    except Exception as e:                            # noqa: BLE001
        row["integrity"] = f"unusable:arm_import_failed ({type(e).__name__}: {e})"
        return row

    if spec.kind == "arm_loader":
        geno_path, pheno_path = Path(mod._GENO), Path(mod._PHENO)
        row["path_provenance"] = "module constant _GENO/_PHENO"
    else:
        geno_path, pheno_path = Path(spec.npz_geno or ""), Path(spec.npz_pheno or "")
        row["path_provenance"] = "argparse default in the arm (pinned by test)"
    row["genotype_path"] = str(geno_path)
    row["phenotype_path"] = str(pheno_path)

    # --- raw pass
    if spec.kind == "npz":
        row["raw"] = {"path": str(geno_path), "exists": geno_path.exists(),
                      "verdict": "intact" if geno_path.exists() else "unusable:absent"}
        if geno_path.exists():
            row["raw"]["size_bytes"] = geno_path.stat().st_size
    else:
        row["raw"] = raw_text_integrity(geno_path)
    row["phenotype_present"] = pheno_path.exists()
    if row["raw"]["verdict"] != "intact":
        row["integrity"] = row["raw"]["verdict"]
        return row

    # --- load through the arm's OWN loader, so this describes the matrix the arm would use
    try:
        if spec.kind == "npz":
            d = np.load(geno_path, allow_pickle=True)
            X = np.asarray(d["X"], dtype=float)
            ids = [str(s) for s in d["samples"]]
            row["loader"] = "np.load(npz)['X']"
        else:
            loaded = mod.load_genotype(str(geno_path), 1)
            X, ids = _unpack_loaded(loaded)
            row["loader"] = f"{spec.arm_module}.load_genotype"
    except Exception as e:                            # noqa: BLE001
        row["integrity"] = f"unusable:load_failed ({type(e).__name__}: {e})"
        return row

    row["matrix"] = matrix_integrity(X)
    if row["matrix"]["verdict"] != "intact":
        row["integrity"] = row["matrix"]["verdict"]
        return row
    row["integrity"] = "intact"

    # --- structure census (never stops on a refusal)
    row["structure"] = structure_census(X, ids)
    try:
        row["structure"]["pc1_null_comparison"] = pc1_excess(X)
    except Exception as e:                            # noqa: BLE001
        row["structure"]["pc1_null_comparison"] = {"error": f"{type(e).__name__}: {e}"}
    # A grouping can clear the panel-level 60% bound and still be unusable for k-fold CV.
    # Record fold-compatibility at k=5 so the consumer does not discover it as a raise.
    row["fold_compatible_k5"] = _fold_compatible(X, ids, k=5)
    return row


def _fold_compatible(X, ids, k: int) -> dict:
    """Is there a rung that is BOTH fold-compatible (largest <= 1/k) and not mostly singletons?"""
    from dna_decode.eval.panel_structure import PanelStructureError, panel_groups

    try:
        pg = panel_groups(X, ids, max_largest=1.0 / k)
        return {"k": k, "compatible": True, "percentile": pg.percentile,
                "n_clusters": pg.n_clusters,
                "largest_cluster_fraction": round(pg.largest_cluster_fraction, 4),
                "singleton_fraction": round(pg.singleton_fraction, 4)}
    except PanelStructureError as e:
        return {"k": k, "compatible": False, "reason": str(e)}


def _unpack_loaded(loaded):
    """The arms' `load_genotype` return shapes differ; accept the ones that occur."""
    if isinstance(loaded, tuple):
        arr = None
        ids = None
        for part in loaded:
            a = np.asarray(part, dtype=object)
            if a.ndim == 2 and arr is None:
                arr = np.asarray(part, dtype=float)
            elif a.ndim == 1 and ids is None and arr is not None and a.shape[0] == arr.shape[0]:
                ids = [str(x) for x in part]
        if arr is None:
            raise ValueError(f"no 2-D matrix in loader return of {len(loaded)} parts")
        if ids is None:
            ids = [f"i{j}" for j in range(arr.shape[0])]
        return arr, ids
    arr = np.asarray(loaded, dtype=float)
    return arr, [f"i{j}" for j in range(arr.shape[0])]


def pc1_excess(X: np.ndarray, *, n_markers: int = 1200, n_perm: int = 3,
               seed: int = 0) -> dict:
    """pc1 against the panel's OWN permutation null, at a MATCHED marker count.

    Why this exists: raw pc1 is NOT comparable across panels, because its no-structure
    expectation scales as ~1/n. Within-marker shuffling destroys inter-individual correlation
    while preserving allele frequencies, so `excess = observed / null` is panel-native.

    It does NOT rescue pc1 as a SELECTION statistic -- see `PC1_SELECTION_VERDICT`. pc1
    measures TOTAL genotypic correlation (linkage disequilibrium PLUS population structure),
    and a bi-parental cross maximises LD by construction while having zero structure.
    """
    from dna_decode.eval.panel_structure import pc1_variance_fraction

    rng = np.random.default_rng(seed)
    cols = rng.choice(X.shape[1], size=min(n_markers, X.shape[1]), replace=False)
    Xs = X[:, cols]
    obs = pc1_variance_fraction(Xs)
    nulls = []
    for _ in range(n_perm):
        P = Xs.copy()
        for j in range(P.shape[1]):
            P[:, j] = P[rng.permutation(P.shape[0]), j]
        nulls.append(pc1_variance_fraction(P))
    nm = float(np.mean(nulls))
    return {"pc1_observed": obs, "pc1_null_mean": nm,
            "pc1_excess": (obs / nm if nm > 0 else None),
            "n_markers_used": int(len(cols)), "n_perm": n_perm}


PC1_SELECTION_VERDICT = (
    "PC1_UNFIT_FOR_SELECTING_A_STRUCTURED_PANEL. The plan pre-registered 'the intact panel "
    "with the highest measured pc1_variance_fraction' as the positive control. Measured, that "
    "rule is self-contradictory on this panel set: by raw pc1 the ordering is "
    "bxd > bloom > dgrp > arabmagic, and by excess-over-its-own-null it is "
    "bloom > bxd > arabmagic > dgrp -- so under BOTH normalisations the one NATURAL population "
    "(dgrp) ranks last or third, and under excess the top panel is BLOOM, the bi-parental cross "
    "that must be the NEGATIVE control. Mechanism, which is structural rather than statistical: "
    "pc1 measures total genotypic correlation = LD + structure, and a bi-parental cross "
    "maximises LD by construction while carrying no population structure. The control pair is "
    "therefore selected BY CONSTRUCTION from each panel's `design` field (declared before any "
    "number was read), with cluster SHAPE as independent corroboration -- not by pc1."
)


def build_census() -> dict:
    rows = [census_panel(s) for s in PANELS]
    intact = [r for r in rows if r.get("integrity") == "intact"]
    return {
        "schema": "f3-substrate-census-v1",
        "generated": date.today().isoformat(),
        "read_only": True,
        "n_panels": len(rows),
        "n_intact": len(intact),
        "intact_panels": [r["panel"] for r in intact],
        "panels": rows,
        "pc1_selection_verdict": PC1_SELECTION_VERDICT,
        "fold_compatible_panels": [r["panel"] for r in rows
                                   if (r.get("fold_compatible_k5") or {}).get("compatible")],
        "honest_limits": [
            "A census bounds INTEGRITY and STRUCTURE only. It says nothing about whether a "
            "panel's phenotypes are usable, nor whether its grouping is the right one.",
            "`design` is recorded by CONSTRUCTION (how the panel was built), not measured. "
            "pc1_variance_fraction is the measured statistic.",
            "A clean cluster shape is NOT evidence of structure: a percentile threshold "
            "manufactures well-shaped clusters on an unstructured panel (measured in "
            "eval/panel_structure.py). Read pc1_variance_fraction for existence.",
            "pc1 cannot SELECT a structured panel -- see pc1_selection_verdict. It conflates "
            "linkage disequilibrium with population structure, and a bi-parental cross "
            "maximises the former while having none of the latter.",
            "pc1_excess uses 3 permutations at 1200 subsampled markers: enough to establish "
            "an ordering, not a precise effect size.",
        ],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path,
                    default=REPO / "wiki" / f"f3_substrate_census_{date.today().isoformat()}.json")
    ap.add_argument("--self-check", action="store_true",
                    help="run the offline fixtures only; touch no cached panel")
    a = ap.parse_args(argv)

    if a.self_check:
        return _self_check()

    cen = build_census()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(cen, indent=2), encoding="utf-8")

    print(f"F3 substrate census -> {a.out}")
    print(f"  {cen['n_intact']} of {cen['n_panels']} panels intact")
    for r in cen["panels"]:
        st = r.get("structure") or {}
        pc1 = st.get("pc1_variance_fraction")
        g = st.get("groups") or {}
        line = f"  {r['panel']:<10} {r.get('integrity','?'):<40} design={r['design']}"
        if pc1 is not None:
            line += f" pc1={pc1:.4f}"
        if g:
            line += (f" clusters={g['n_clusters']} largest={g['largest_cluster_fraction']:.2f}"
                     f" p={g['percentile']:g}")
        elif st.get("groups_refused"):
            line += "  groups=REFUSED"
        print(line)
    return 0                                           # a census never gates


def _self_check() -> int:
    import tempfile
    ok = True
    with tempfile.TemporaryDirectory() as td:
        t = Path(td)
        # a truncated gzip: valid header, body cut
        good = t / "good.gz"
        with gzip.open(good, "wb") as fh:
            fh.write(b"x" * (1 << 20))
        trunc = t / "trunc.gz"
        trunc.write_bytes(good.read_bytes()[: good.stat().st_size // 2])
        full_ok, _, _ = verify_gzip_fully(trunc)
        head_ok = head_only_gzip_ok(trunc)
        print(f"  truncated gzip: full-decompress ok={full_ok} (want False), "
              f"head-only ok={head_ok} (want True -- the wrong guard)")
        ok &= (full_ok is False)
        ragged = t / "ragged.csv"
        ragged.write_text("id,a,b\nx,1,2\ny,1\n")
        v = raw_text_integrity(ragged)["verdict"]
        print(f"  ragged rows -> {v}")
        ok &= v.startswith("unusable:ragged_rows")
    print("SELF-CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
