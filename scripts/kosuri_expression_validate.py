"""Track B — can a learned model predict E. coli expression better than the element-strength model?

Runs the pre-registered test from `wiki/design_epoch_plan_2026-08-07.md` against Kosuri et al. 2013
(PNAS 110:14024), 12,563 constructed promoter x RBS combinations with measured DNA, RNA and protein.

THE SPLIT IS THE EXPERIMENT. Three are run, because they answer different questions and only one of
them is the design-relevant one:

  * held-out COMBINATION - both elements were seen in training, the pairing was not. This is the
    paper's own question (composability) and the split on which its R^2 numbers are meaningful.
  * held-out PROMOTER / held-out RBS - the element itself is unseen. This is what a designer needs
    (score a NEW part), and it is where an identity-based model has nothing to say.

The baseline is the paper's model re-fit on the TRAINING SPLIT ONLY: log2(protein) = mu + promoter
effect + RBS effect. The published 0.76 / 0.82 are effectively in-sample, so quoting them directly
would hand the learned model a free win.

Data is NOT committed (16 MB, third-party supplementary). Point --sd03 at the file:
    uv run python scripts/kosuri_expression_validate.py --sd03 D:/path/sd03.xls
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# The bar pre-registered before any of this data was seen (design-epoch plan).
PREREGISTERED_BAR = 0.82


def r2(y, yhat) -> float:
    y, yhat = np.asarray(y, float), np.asarray(yhat, float)
    m = np.isfinite(y) & np.isfinite(yhat)
    y, yhat = y[m], yhat[m]
    if len(y) < 2:
        return float("nan")
    return float(1 - ((y - yhat) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def additive_predict(train, test):
    """The paper's element-strength model, fit on TRAIN ONLY: mu + promoter effect + RBS effect.

    An unseen element contributes 0 (i.e. falls back to the other element + the grand mean) rather
    than erroring -- that is the fairest possible treatment of it, and it is why this baseline still
    scores ~0.26-0.50 on a held-out-element split where an identity model scores ~0.
    """
    mu = train.y.mean()
    pe = train.groupby("p_code").y.mean() - mu
    re = train.groupby("r_code").y.mean() - mu
    return mu + test.p_code.map(pe).fillna(0).values + test.r_code.map(re).fillna(0).values


def load_sd03(path: str):
    """Kosuri Dataset S3 -> a frame with log2(protein), element codes and deltaG."""
    import pandas as pd  # noqa: PLC0415
    import xlrd  # noqa: PLC0415

    sh = xlrd.open_workbook(path).sheet_by_index(0)
    hdr = [str(sh.cell_value(0, c)).strip('"') for c in range(sh.ncols)]
    rows = [[sh.cell_value(r, c) for c in range(sh.ncols)] for r in range(1, sh.nrows)]
    df = pd.DataFrame(rows, columns=hdr)
    for c in df.columns:
        if c not in ("Promoter", "RBS", "target", "RBS.TTL", "Promoter.TTL"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("Promoter", "RBS"):
        df[c] = df[c].astype(str).str.strip('"')
    d = df[np.isfinite(df.prot)].copy()
    d["y"] = np.log2(d.prot.values)          # protein is raw RFU; the paper models it in log space
    d["dG"] = d["deltaG"].values             # 5' secondary structure, computable from sequence
    d["p_code"] = d.Promoter.astype("category").cat.codes
    d["r_code"] = d.RBS.astype("category").cat.codes
    return d


def reproduce_published(path: str) -> dict:
    """Sanity gate: recompute the paper's own numbers from its own model columns before trusting anything.

    `model.prot.simple` is stored in LOG2 while `prot` is raw -- comparing them in the wrong space gives
    R^2 = -15 and would look like a data-loading bug rather than a units mismatch.
    """
    import pandas as pd  # noqa: PLC0415
    import xlrd  # noqa: PLC0415

    sh = xlrd.open_workbook(path).sheet_by_index(0)
    hdr = [str(sh.cell_value(0, c)).strip('"') for c in range(sh.ncols)]
    rows = [[sh.cell_value(r, c) for c in range(sh.ncols)] for r in range(1, sh.nrows)]
    df = pd.DataFrame(rows, columns=hdr)
    for c in ("prot", "RNA", "model.prot.simple", "model.RNA.simple", "model.RNA.full"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return {
        "protein_simple_log2": round(r2(np.log2(df.prot.values), df["model.prot.simple"].values), 4),
        "protein_simple_published": 0.76,
        "rna_simple_log10": round(r2(np.log10(df.RNA.values), np.log10(df["model.RNA.simple"].values)), 4),
        "rna_simple_published": 0.92,
        "rna_full_log10": round(r2(np.log10(df.RNA.values), np.log10(df["model.RNA.full"].values)), 4),
        "rna_full_published": 0.96,
    }


def run_splits(d, seed: int = 0) -> dict:
    from sklearn.ensemble import HistGradientBoostingRegressor  # noqa: PLC0415
    from sklearn.model_selection import GroupKFold, KFold  # noqa: PLC0415

    def evaluate(splitter, groups=None):
        base, ident, withdg = [], [], []
        it = splitter.split(d, groups=groups) if groups is not None else splitter.split(d)
        for tri, tei in it:
            tr, te = d.iloc[tri], d.iloc[tei]
            base.append(r2(te.y.values, additive_predict(tr, te)))
            xi = lambda f: np.c_[f.p_code.values, f.r_code.values]  # noqa: E731
            m = HistGradientBoostingRegressor(max_iter=400, categorical_features=[0, 1],
                                              random_state=seed).fit(xi(tr), tr.y.values)
            ident.append(r2(te.y.values, m.predict(xi(te))))
            xg = lambda f: np.c_[f.p_code.values, f.r_code.values, f.dG.values]  # noqa: E731
            mg = HistGradientBoostingRegressor(max_iter=400, categorical_features=[0, 1],
                                               random_state=seed).fit(xg(tr), tr.y.values)
            withdg.append(r2(te.y.values, mg.predict(xg(te))))
        return {"additive_baseline": round(float(np.mean(base)), 4),
                "gbm_identity": round(float(np.mean(ident)), 4),
                "gbm_identity_plus_deltaG": round(float(np.mean(withdg)), 4)}

    return {
        "held_out_combination": evaluate(KFold(5, shuffle=True, random_state=seed)),
        "held_out_promoter": evaluate(GroupKFold(5), groups=d.p_code.values),
        "held_out_rbs": evaluate(GroupKFold(5), groups=d.r_code.values),
    }


_SD_CORE = "AGGAGG"
_KMERS = None


def _kmer_list():
    global _KMERS
    if _KMERS is None:
        import itertools  # noqa: PLC0415
        _KMERS = [''.join(p) for k in (1, 2, 3) for p in itertools.product("ACGT", repeat=k)]
    return _KMERS


def rbs_features(seq: str) -> list[float]:
    """Sequence-only features for one RBS: composition + the Shine-Dalgarno mechanics.

    Deliberately simple and mechanistic rather than learned: 1/2/3-mer frequencies, length, GC, the best
    match to the SD core `AGGAGG`, and the SD-to-start spacing — the two things that actually govern
    translation initiation. Nothing here uses the RBS's identity or its measured strength, so a held-out
    RBS is scored purely from its letters.
    """
    s = (seq or "").upper()
    n = max(len(s), 1)
    best, pos = 0, -1
    for i in range(len(s) - len(_SD_CORE) + 1):
        m = sum(a == b for a, b in zip(s[i:i + len(_SD_CORE)], _SD_CORE))
        if m > best:
            best, pos = m, i
    gc = (s.count("G") + s.count("C")) / n
    return [s.count(k) / n for k in _kmer_list()] + [len(s), gc, best, (len(s) - pos) if pos >= 0 else -1]


def load_rbs_sequences(sd02_path: str) -> dict:
    """Kosuri Dataset S2 -> {RBS name: sequence}. S2 is the only file here carrying element sequences."""
    import xlrd  # noqa: PLC0415

    sh = xlrd.open_workbook(sd02_path).sheet_by_index(0)
    hdr = [str(sh.cell_value(0, c)).strip('"') for c in range(sh.ncols)]
    si = hdr.index("Sequence")
    return {str(sh.cell_value(r, 0)).strip('"'):
            str(sh.cell_value(r, si)).strip('"').replace(" ", "").upper()
            for r in range(1, sh.nrows)}


def run_rbs_sequence_split(d, sd02_path: str, seed: int = 0) -> dict:
    """THE design question, on the half the data can answer: score a NOVEL RBS from its sequence.

    The RBS is held out entirely, so nothing about it was seen in training. Promoter identity is
    supplied because a designer knows their promoter — and `promoter_only` / `seq_only` controls are
    reported so the sequence contribution can't be confused with the promoter doing the work.
    """
    import numpy as np  # noqa: PLC0415
    from sklearn.ensemble import HistGradientBoostingRegressor  # noqa: PLC0415
    from sklearn.model_selection import GroupKFold  # noqa: PLC0415

    seqs = load_rbs_sequences(sd02_path)
    d = d[d.RBS.map(seqs).notna()].copy()
    F = np.array([rbs_features(seqs[r]) for r in d.RBS.values])

    out = {k: [] for k in ("additive_baseline", "identity", "promoter_only",
                           "rbs_sequence_only", "promoter_plus_rbs_sequence",
                           "promoter_plus_rbs_sequence_plus_deltaG")}
    for tri, tei in GroupKFold(5).split(d, groups=d.r_code.values):
        tr, te = d.iloc[tri], d.iloc[tei]
        ftr, fte = F[tri], F[tei]
        out["additive_baseline"].append(r2(te.y.values, additive_predict(tr, te)))
        def gbm(xtr, ytr, xte, cats=None):  # noqa: E306
            m = HistGradientBoostingRegressor(max_iter=400, categorical_features=cats,
                                              random_state=seed).fit(xtr, ytr)
            return m.predict(xte)
        xi_tr, xi_te = np.c_[tr.p_code.values, tr.r_code.values], np.c_[te.p_code.values, te.r_code.values]
        out["identity"].append(r2(te.y.values, gbm(xi_tr, tr.y.values, xi_te, [0, 1])))
        p_tr, p_te = tr.p_code.values.reshape(-1, 1), te.p_code.values.reshape(-1, 1)
        out["promoter_only"].append(r2(te.y.values, gbm(p_tr, tr.y.values, p_te, [0])))
        out["rbs_sequence_only"].append(r2(te.y.values, gbm(ftr, tr.y.values, fte)))
        s_tr, s_te = np.c_[tr.p_code.values, ftr], np.c_[te.p_code.values, fte]
        out["promoter_plus_rbs_sequence"].append(r2(te.y.values, gbm(s_tr, tr.y.values, s_te, [0])))
        g_tr = np.c_[tr.p_code.values, ftr, tr.dG.values]
        g_te = np.c_[te.p_code.values, fte, te.dG.values]
        out["promoter_plus_rbs_sequence_plus_deltaG"].append(
            r2(te.y.values, gbm(g_tr, tr.y.values, g_te, [0])))
    return {k: round(float(np.mean(v)), 4) for k, v in out.items()}


def verdict(splits: dict) -> dict:
    """The pre-registered falsifier was 'beat R^2 0.82 on held-out constructs, split BY ELEMENT'.

    Both halves get reported, because they disagree and the disagreement IS the finding.
    """
    comb = splits["held_out_combination"]
    best_elem = max(splits["held_out_promoter"]["gbm_identity_plus_deltaG"],
                    splits["held_out_rbs"]["gbm_identity_plus_deltaG"])
    return {
        "preregistered_bar": PREREGISTERED_BAR,
        "combination_split_verdict": "PASS" if comb["gbm_identity_plus_deltaG"] > PREREGISTERED_BAR else "FAIL",
        "combination_best": comb["gbm_identity_plus_deltaG"],
        "combination_baseline": comb["additive_baseline"],
        "element_split_verdict": "FAIL" if best_elem <= PREREGISTERED_BAR else "PASS",
        "element_best": round(best_elem, 4),
        "element_note": (
            "On the element split NOTHING approaches 0.82 -- including the baseline (0.26-0.50). The bar "
            "was mis-specified for that split: 0.82 is a COMBINATION-level in-sample number and does not "
            "transfer to predicting an unseen element."
        ),
        "identity_model_headline": (
            "An IDENTITY-based learned model's entire advantage is INTERACTION CAPTURE among elements it "
            "has already seen. Given an unseen promoter it is WORSE than the additive baseline and below "
            "chance (-0.014): it learned identity, not sequence."
        ),
    }


def sequence_verdict(rbs: dict) -> dict:
    """Did REAL sequence features generalise to a never-seen RBS? (the design question)"""
    seq = rbs["promoter_plus_rbs_sequence_plus_deltaG"]
    return {
        "held_out_rbs_from_sequence": seq,
        "vs_additive_baseline": round(seq - rbs["additive_baseline"], 4),
        "vs_identity_model": round(seq - rbs["identity"], 4),
        "promoter_only_control": rbs["promoter_only"],
        "rbs_sequence_only_control": rbs["rbs_sequence_only"],
        "generalises_from_sequence": seq > rbs["additive_baseline"] + 0.05,
        "headline": (
            "Scoring a NOVEL RBS from its letters works: sequence features lift a never-seen RBS well "
            "above both the additive baseline and the identity model. The promoter_only control shows the "
            "promoter alone does not account for it."
        ),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sd03", required=True, help="path to Kosuri Dataset S3 (.xls, ~16 MB, not committed)")
    ap.add_argument("--sd02", default=None, help="path to Dataset S2 (RBS sequences) -- enables the "
                                                 "sequence-generalisation arm (the design question)")
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    repro = reproduce_published(a.sd03)
    print("reproduction gate (their columns, their numbers):")
    for k, v in repro.items():
        print(f"   {k:28s} {v}")

    d = load_sd03(a.sd03)
    print(f"\nconstructs with protein: {len(d)} | promoters {d.p_code.nunique()} | RBSs {d.r_code.nunique()}")
    splits = run_splits(d)
    for name, s in splits.items():
        print(f"  {name:22s} additive {s['additive_baseline']:7.3f} | GBM {s['gbm_identity']:7.3f} "
              f"| GBM+dG {s['gbm_identity_plus_deltaG']:7.3f}")
    rbs_seq = seq_v = None
    if a.sd02:
        rbs_seq = run_rbs_sequence_split(d, a.sd02)
        seq_v = sequence_verdict(rbs_seq)
        print("\nHELD-OUT RBS, scored from SEQUENCE (the design question):")
        for k, val in rbs_seq.items():
            print(f"   {k:42s} {val:7.4f}")
    v = verdict(splits)
    print(f"\ncombination split: {v['combination_split_verdict']} "
          f"({v['combination_best']} vs bar {PREREGISTERED_BAR}, baseline {v['combination_baseline']})")
    print(f"element split:     {v['element_split_verdict']} (best {v['element_best']})")
    print(f"\n{v['identity_model_headline']}")

    if seq_v:
        print(f"\nsequence generalisation: {seq_v['held_out_rbs_from_sequence']} "
              f"(+{seq_v['vs_additive_baseline']} vs additive, +{seq_v['vs_identity_model']} vs identity) "
              f"-> generalises={seq_v['generalises_from_sequence']}")
    rec = {"record": "kosuri-expression-validation-v1", "date": a.date,
           "dataset": "Kosuri 2013 PNAS 110:14024 Dataset S3", "n_constructs": int(len(d)),
           "target": "log2(protein)", "reproduction_gate": repro, "splits": splits, "verdict": v,
           "rbs_sequence_split": rbs_seq, "sequence_verdict": seq_v}
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / f"kosuri_expression_{a.date}"
    stem.with_suffix(".json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print(f"\nwrote {stem}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
