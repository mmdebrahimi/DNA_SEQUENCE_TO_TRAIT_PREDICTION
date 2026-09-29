"""Does a LEARNED genotype alphabet beat the zero-shot incumbents on held-out POSITIONS?

THE QUESTION. The user's proposal is a genotype alphabet for a Genome LLM. This repo's only working
learned regime (constructed x molecular) is served by ZERO-SHOT protein-LM scoring -- nothing in it is
trained on measured phenotype. So a supervised genotype-token model would be the FIRST supervised entrant
there, and PEAR is the cheapest decisive substrate: 2,114 constructed CTX-M-14 variants with measured
cefotaxime fitness, an independent lab, and it clears all ten rejection gates.

THE SPLIT IS BY PROTEIN POSITION, AND THAT IS THE WHOLE EXPERIMENT. Each position carries ~19 variants, so
a random VARIANT split leaks the position mean and any model would look good by memorising "position 104 is
tolerant". Holding out whole positions forces a model to generalise across the alphabet, which is exactly
the claim a genotype alphabet makes.

WHAT COUNTS AS THE ALPHABET. A position-free map from (wt_aa -> mut_aa) to effect. Note what that means:
**BLOSUM62 already IS a learned genotype alphabet** -- a 20x20 substitution matrix learned from alignments
in 1992. So "learn an alphabet" is not a new idea here; the testable question is whether learning it from
THIS protein's own measured DMS beats the generic 1992 one, and whether either beats a PLM.

PRE-REGISTERED BARS (frozen before any number was seen; see also PREREGISTERED below):
  B1  the learned alphabet must beat BLOSUM62 on held-out positions.
      Rationale: if learning from this protein's own DMS cannot beat a 34-year-old generic matrix, a
      learned alphabet adds nothing and the GLM premise loses its cheapest support.
  B2  to "earn the next step" for a GLM, a supervised model must beat ESM2 zero-shot on the SAME held-out
      positions. Beating BLOSUM but losing to ESM2 is a real and informative outcome: it would say a
      POSITION-FREE alphabet is insufficient and the GLM needs position/structure context.
  B3  the decisive variant is E (supervised head ON TOP of the ESM2 score). E > B is the cornerstone
      question -- does supervision ADD to a pretrained representation, or merely re-derive it?

HONESTY RAILS BUILT IN:
  * Never compared to the published POOLED 0.352 -- a pooled number and a held-out-position number are
    different quantities. Every number here is recomputed on the identical held-out set.
  * `n_pair_unseen` is reported: a (wt,mut) pair absent from train must fall back, and if most fall back
    then "the learned alphabet" is really just its fallback. A high fallback rate invalidates B1.
  * Mid-ranks via scipy's default; BLOSUM62 ties heavily and tie-breaking by sort order has already
    shifted a median in this repo.
  * Sign convention asserted, not assumed: silent > missense > nonsense on measured fitness, checked
    before any correlation is read. A frame/coordinate error fails loudly here.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

REF = Path("D:/dna_decode_cache/pear/CTXM-14/Genotype_barcode_calling/Ref_CTXM.fasta")
TABLE = Path("D:/dna_decode_cache/pear/extracted/Figure3.A__data.tsv")
ESM_CACHE = Path("D:/dna_decode_cache/pear/ctxm14_esm2_650M_marginals.json")
OUT = ROOT / "wiki" / "pear_genotype_alphabet_2026-09-29.json"

PREREGISTERED = {
    "frozen_before_any_number_was_seen": True,
    "B1_learned_alphabet_beats_blosum62": "learned alphabet Spearman > BLOSUM62 Spearman on held-out positions",
    "B2_supervised_beats_esm2": "best supervised Spearman > ESM2 Spearman on the SAME held-out positions",
    "B3_decisive": "variant E (supervised head on top of ESM2) > variant B (ESM2 alone)",
    "invalidator": "if n_pair_unseen / n_test > 0.50 the learned alphabet is mostly its fallback and B1 is void",
    "not_compared_to": "the published POOLED 0.352 -- different quantity from a held-out-position number",
}

CODON = {  # standard table
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
AAS = "ACDEFGHIKLMNPQRSTVWY"


def load_ref() -> str:
    return "".join(l.strip() for l in REF.read_text().splitlines() if not l.startswith(">"))


def parse_variants(cds: str) -> list[dict]:
    """Nucleotide `C648T` -> protein consequence. REF base must match the CDS or we refuse."""
    rows, mismatched = [], 0
    for line in TABLE.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        gt, ctx = parts[0].strip(), parts[1].strip()
        if not gt or gt[0] not in "ACGT" or not gt[-1] in "ACGT":
            continue
        try:
            nt_pos, fit = int(gt[1:-1]), float(ctx)
        except ValueError:
            continue
        wt_nt, mut_nt = gt[0], gt[-1]
        if not (1 <= nt_pos <= len(cds)) or cds[nt_pos - 1] != wt_nt:
            mismatched += 1
            continue
        cpos = (nt_pos - 1) // 3                      # 0-based codon index
        off = (nt_pos - 1) % 3
        wt_cod = cds[cpos * 3: cpos * 3 + 3]
        mut_cod = wt_cod[:off] + mut_nt + wt_cod[off + 1:]
        wt_aa, mut_aa = CODON.get(wt_cod), CODON.get(mut_cod)
        if wt_aa is None or mut_aa is None:
            continue
        kind = "silent" if wt_aa == mut_aa else ("nonsense" if mut_aa == "*" else "missense")
        rows.append({"gt": gt, "pos": cpos + 1, "wt": wt_aa, "mut": mut_aa,
                     "kind": kind, "fitness": fit})
    if mismatched:
        raise SystemExit(f"REFUSING: {mismatched} variants whose REF base does not match the CDS -- "
                         "a coordinate/frame error, not a data quirk.")
    return rows


def assert_sign_convention(rows: list[dict]) -> dict:
    """Silent > missense > nonsense on MEASURED fitness. No frame-shifted mapping gives this by accident."""
    med = {k: float(np.median([r["fitness"] for r in rows if r["kind"] == k]))
           for k in ("silent", "missense", "nonsense")}
    ok = med["silent"] > med["missense"] > med["nonsense"]
    if not ok:
        raise SystemExit(f"REFUSING: consequence ordering is wrong ({med}) -- the protein mapping is broken.")
    return med


def esm2_marginals(protein: str) -> dict:
    """Masked-marginal log-prob per (position, aa). Cached -- ~291 forward passes."""
    if ESM_CACHE.exists():
        return json.loads(ESM_CACHE.read_text(encoding="utf-8"))
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    name = "facebook/esm2_t33_650M_UR50D"
    tok = AutoTokenizer.from_pretrained(name)
    model = AutoModelForMaskedLM.from_pretrained(name)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(dev).eval()
    enc = tok(protein, return_tensors="pt")
    ids = enc["input_ids"].to(dev)
    out: dict[str, dict] = {}
    with torch.no_grad():
        for i in range(len(protein)):
            masked = ids.clone()
            masked[0, i + 1] = tok.mask_token_id          # +1 for the CLS token
            logits = model(input_ids=masked).logits[0, i + 1]
            lp = torch.log_softmax(logits.float(), dim=-1)
            out[str(i + 1)] = {aa: float(lp[tok.convert_tokens_to_ids(aa)]) for aa in AAS}
            if (i + 1) % 50 == 0:
                print(f"    esm2 marginals {i+1}/{len(protein)}", flush=True)
    ESM_CACHE.parent.mkdir(parents=True, exist_ok=True)
    ESM_CACHE.write_text(json.dumps(out), encoding="utf-8")
    return out


def blosum62_fn():
    """The SHIPPED scorer, not a hand-typed matrix -- so this test uses the same BLOSUM62 the CLI ships."""
    from dna_decode.forward.variant_effect import blosum62_score
    # non-vacuity: a conservative substitution must score above a radical one, else the accessor is wrong
    if not blosum62_score("I", "V") > blosum62_score("W", "P"):
        raise SystemExit("REFUSING: blosum62_score does not order I>V above W>P -- wrong accessor.")
    return blosum62_score


def ridge(X: np.ndarray, y: np.ndarray, lam: float = 1.0) -> np.ndarray:
    X = np.column_stack([X, np.ones(len(X))])
    A = X.T @ X + lam * np.eye(X.shape[1])
    A[-1, -1] -= lam                                    # do not penalise the intercept
    return np.linalg.solve(A, X.T @ y)


def predict(w: np.ndarray, X: np.ndarray) -> np.ndarray:
    return np.column_stack([X, np.ones(len(X))]) @ w


def onehot(rows: list[dict], esm: dict | None, with_esm: bool) -> np.ndarray:
    cols = []
    for r in rows:
        wt = [1.0 if a == r["wt"] else 0.0 for a in AAS]
        mu = [1.0 if a == r["mut"] else 0.0 for a in AAS]
        feats = wt + mu
        if with_esm:
            feats = feats + [esm[str(r["pos"])][r["mut"]] - esm[str(r["pos"])][r["wt"]]]
        cols.append(feats)
    return np.asarray(cols, dtype=float)


def spearman(a, b) -> float:
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--holdout-stride", type=int, default=3,
                    help="every Nth unique POSITION (sorted) is held out; deterministic, no seed")
    args = ap.parse_args(argv)

    cds = load_ref()
    rows = parse_variants(cds)
    medians = assert_sign_convention(rows)
    mis = [r for r in rows if r["kind"] == "missense"]
    print(f"variants {len(rows)} | missense {len(mis)} | consequence medians {medians}")

    protein = "".join(CODON[cds[i:i + 3]] for i in range(0, len(cds) - len(cds) % 3, 3)).rstrip("*")
    print(f"protein length {len(protein)}")
    print("computing/loading ESM2-650M masked marginals ...", flush=True)
    esm = esm2_marginals(protein)
    bl = blosum62_fn()

    # ---- SPLIT BY POSITION (the experiment) ----
    positions = sorted({r["pos"] for r in mis})
    test_pos = set(positions[::args.holdout_stride])
    tr = [r for r in mis if r["pos"] not in test_pos]
    te = [r for r in mis if r["pos"] in test_pos]
    assert not ({r["pos"] for r in tr} & {r["pos"] for r in te}), "position leak"
    print(f"positions {len(positions)} -> train {len(positions)-len(test_pos)} / test {len(test_pos)}")
    print(f"variants  train {len(tr)} / test {len(te)}")

    y_te = np.array([r["fitness"] for r in te])
    y_tr = np.array([r["fitness"] for r in tr])

    # ---- A: BLOSUM62 (the 1992 generic alphabet) ----
    sA = np.array([float(bl(r["wt"], r["mut"])) for r in te])
    # ---- B: ESM2 zero-shot ----
    sB = np.array([esm[str(r["pos"])][r["mut"]] - esm[str(r["pos"])][r["wt"]] for r in te])
    # ---- C: learned alphabet, pair means from TRAIN only ----
    pair, permut = {}, {}
    for r in tr:
        pair.setdefault((r["wt"], r["mut"]), []).append(r["fitness"])
        permut.setdefault(r["mut"], []).append(r["fitness"])
    pair_mean = {k: float(np.mean(v)) for k, v in pair.items()}
    mut_mean = {k: float(np.mean(v)) for k, v in permut.items()}
    gmean = float(np.mean(y_tr))
    unseen = 0
    sC = []
    for r in te:
        k = (r["wt"], r["mut"])
        if k in pair_mean:
            sC.append(pair_mean[k])
        else:
            unseen += 1
            sC.append(mut_mean.get(r["mut"], gmean))
    sC = np.array(sC)
    # ---- D: ridge on one-hot(wt)+one-hot(mut) ----
    wD = ridge(onehot(tr, esm, False), y_tr)
    sD = predict(wD, onehot(te, esm, False))
    # ---- E: ridge on one-hot + the ESM2 score (the cornerstone variant) ----
    wE = ridge(onehot(tr, esm, True), y_tr)
    sE = predict(wE, onehot(te, esm, True))

    res = {
        "A_blosum62_zeroshot": spearman(sA, y_te),
        "B_esm2_zeroshot": spearman(sB, y_te),
        "C_learned_alphabet_pairmeans": spearman(sC, y_te),
        "D_ridge_onehot_wt_mut": spearman(sD, y_te),
        "E_ridge_onehot_PLUS_esm2": spearman(sE, y_te),
    }
    n_unseen_frac = unseen / len(te)
    verdicts = {
        "B1_learned_alphabet_beats_blosum62": res["C_learned_alphabet_pairmeans"] > res["A_blosum62_zeroshot"]
        or res["D_ridge_onehot_wt_mut"] > res["A_blosum62_zeroshot"],
        "B2_supervised_beats_esm2": max(res["C_learned_alphabet_pairmeans"], res["D_ridge_onehot_wt_mut"],
                                        res["E_ridge_onehot_PLUS_esm2"]) > res["B_esm2_zeroshot"],
        "B3_esm2_plus_head_beats_esm2_alone": res["E_ridge_onehot_PLUS_esm2"] > res["B_esm2_zeroshot"],
        "B1_VOID_fallback_dominates": n_unseen_frac > 0.50,
    }

    print("\n=== held-out-POSITION Spearman (higher = better; NOT comparable to the pooled 0.352) ===")
    for k, v in sorted(res.items(), key=lambda kv: -kv[1]):
        print(f"  {k:34s} {v:+.4f}")
    print(f"\n  n_pair_unseen {unseen}/{len(te)} = {n_unseen_frac:.3f}")
    print("\n=== PRE-REGISTERED VERDICTS ===")
    for k, v in verdicts.items():
        print(f"  {k:38s} {v}")

    OUT.write_text(json.dumps({
        "schema": "pear-genotype-alphabet-v1",
        "analysis_date": "2026-09-29",
        "substrate": "PEAR CTX-M-14 cefotaxime, constructed variation, independent lab",
        "split": f"by PROTEIN POSITION, every {args.holdout_stride}rd sorted unique position held out",
        "n_missense": len(mis), "n_train": len(tr), "n_test": len(te),
        "n_positions_total": len(positions), "n_positions_test": len(test_pos),
        "consequence_medians": medians,
        "spearman_heldout": res, "n_pair_unseen": unseen, "pair_unseen_frac": n_unseen_frac,
        "prereg": PREREGISTERED, "verdicts": verdicts,
        "honest_limits": [
            "ONE protein, ONE assay, ONE lab -- a position-held-out result on CTX-M-14, not a general law.",
            "The alphabet variants are POSITION-FREE by construction; that is the hypothesis under test, "
            "not a limitation of the implementation.",
            "ESM2 is zero-shot here and sees no measured fitness; the supervised variants see 2/3 of "
            "positions. That asymmetry FAVOURS the supervised side, so a supervised loss is decisive "
            "while a supervised win is only suggestive.",
            "Held out by position, not by structural neighbourhood -- adjacent positions can share local "
            "context, so B2/B3 may be optimistic for the supervised side.",
        ],
    }, indent=1), encoding="utf-8")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
