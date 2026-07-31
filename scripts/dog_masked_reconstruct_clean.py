"""F1' CLEAN re-run: NT-v2-100M masked reconstruction on canFam4, leakage-free comparison.

Fixes the two biases the 2026-07-31 adversarial review found in the deprecated F1 smoke:
  (1) Markov no longer fit-and-scored on the same slice -- here it is DISJOINT-fit (train region
      does not overlap the eval window) AND separately LEAVE-ONE-OUT scored; and
  (2) NT is scored by its per-base MARGINAL distribution (not the single argmax 6-mer), with per-base
      NLL as the PRIMARY endpoint.
The Markov baseline is given its BEST shot: a k-sweep picks the k with the LOWEST Markov NLL (the
hardest baseline for NT to beat) as the headline. Real NT weights; runs in the isolated
transformers==4.30.2 env. Writes wiki/dog_masked_reconstruct_clean_<date>.{json,md}.

    HF_HOME=D:/hf_cache uv run --isolated --with transformers==4.30.2 ... \
        python scripts/dog_masked_reconstruct_clean.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.interp.masked_reconstruct import score_from_predictions  # noqa: E402
from dna_decode.models.foundation import model_factory  # noqa: E402

_DEFAULT_SEQ = Path("D:/dna_decode_cache/dog_ref/NC_051804.1_20000000_20006000.txt")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Clean (leakage-free) NT masked-reconstruction on canFam4.")
    ap.add_argument("--seq-file", default=str(_DEFAULT_SEQ))
    ap.add_argument("--eval-bp", type=int, default=600, help="eval window (<= max_context)")
    ap.add_argument("--train-start", type=int, default=1200, help="disjoint Markov-train start (bp)")
    ap.add_argument("--n-tokens", type=int, default=40)
    ap.add_argument("--k-sweep", default="1,2,3,4,5,6,7,8")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--date", default="2026-07-31")
    args = ap.parse_args(argv)

    full = Path(args.seq_file).read_text().strip().upper()
    eval_window = full[: args.eval_bp]
    disjoint_train = full[args.train_start:]        # strictly after the eval window -> no overlap
    assert args.train_start >= args.eval_bp, "train region must start after the eval window"
    ks = [int(x) for x in args.k_sweep.split(",")]

    print(f"loading NT-v2-100M ({args.device}) ...", flush=True)
    model = model_factory("nucleotide_transformer", device=args.device)
    positions = list(range(args.n_tokens))
    print(f"eval {len(eval_window)} bp / {args.n_tokens} tokens; disjoint train {len(disjoint_train)} bp "
          f"(from bp {args.train_start}); k-sweep {ks}", flush=True)

    # Compute the expensive NT forward ONCE (strict: fail loudly if any position is dropped).
    preds = model.masked_token_predictions(eval_window, positions=positions, strict=True)
    print(f"NT preds computed for {len(preds)} tokens; sweeping Markov ...", flush=True)

    # Disjoint-fit Markov, k-swept; headline = the k giving the LOWEST Markov NLL (hardest baseline).
    disjoint_by_k = {
        k: score_from_predictions(preds, eval_window, model_name=model.name, markov_k=k,
                                  markov_fit_sequences=[disjoint_train], region_label="canfam4_disjoint")
        for k in ks
    }
    best_k = min(ks, key=lambda k: disjoint_by_k[k].markov_per_base_nll)
    headline = disjoint_by_k[best_k]

    # Secondary: leave-one-out Markov on the eval window itself (the other leakage control).
    loo = score_from_predictions(preds, eval_window, model_name=model.name, markov_k=best_k,
                                 region_label="canfam4_loo")

    out = {
        "headline": headline.as_dict(),
        "headline_markov_k": best_k,
        "markov_mode": "disjoint-fit (k chosen = lowest Markov NLL = hardest baseline)",
        "disjoint_nll_delta_by_k": {k: round(disjoint_by_k[k].nll_delta, 4) for k in ks},
        "disjoint_markov_nll_by_k": {k: round(disjoint_by_k[k].markov_per_base_nll, 4) for k in ks},
        "leave_one_out_secondary": loo.as_dict(),
        "nt_per_base_nll": round(headline.nt_per_base_nll, 4),
        "substrate": args.seq_file,
        "supersedes": "wiki/dog_masked_reconstruct_smoke_2026-07-31 (deprecated: leaky Markov + argmax NT)",
        "real_surface": True,
    }
    Path(f"wiki/dog_masked_reconstruct_clean_{args.date}.json").write_text(json.dumps(out, indent=2))

    verdict = ("NT BEATS baseline" if headline.nll_delta > 0 else "NT LOSES to baseline")
    print("\n=== NT-v2-100M CLEAN masked reconstruction on canFam4 (F1') ===")
    print(f"  PRIMARY per-base NLL delta (Markov - NT), disjoint, best-k={best_k}: "
          f"{headline.nll_delta:+.4f}  -> {verdict}")
    print(f"  NT per-base NLL        : {headline.nt_per_base_nll:.4f}")
    print(f"  Markov-{best_k} per-base NLL  : {headline.markov_per_base_nll:.4f} (disjoint-fit, hardest k)")
    print(f"  LOO-Markov NLL delta   : {loo.nll_delta:+.4f} (secondary leakage control)")
    print(f"  per-base marginal acc  : NT {headline.nt_per_base_marginal_accuracy:.4f} vs "
          f"Markov {headline.markov_per_base_accuracy:.4f} (delta {headline.accuracy_delta:+.4f})")
    print(f"  disjoint NLL delta by k: {out['disjoint_nll_delta_by_k']}")
    print(f"  n_bases={headline.n_bases_scored}  -> wiki/dog_masked_reconstruct_clean_{args.date}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
