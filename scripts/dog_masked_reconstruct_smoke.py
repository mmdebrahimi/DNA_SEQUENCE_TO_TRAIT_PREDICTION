"""F1 real-surface smoke: NT-v2-100M masked reconstruction on a real canFam4 slice.

Loads REAL Nucleotide-Transformer-v2-100M weights (cached at D:/hf_cache) and reconstructs masked
6-mer tokens of a real dog (canFam4 chr1) sequence, scoring the DELTA vs an order-5 Markov baseline
over the identical masked base set. This is the R3 real-surface integration test for the F1 engine:
"ran once on the real boundary + produced the real artifact". CPU by default (avoids the GTX 860M
display-GPU TDR); small token budget so it finishes in ~1-2 min.

    HF_HOME=D:/hf_cache uv run python scripts/dog_masked_reconstruct_smoke.py

Writes wiki/dog_masked_reconstruct_smoke_<date>.json. The headline is the delta, NOT raw accuracy.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.interp.masked_reconstruct import score_reconstruction  # noqa: E402
from dna_decode.models.foundation import model_factory  # noqa: E402

_DEFAULT_SEQ = Path("D:/dna_decode_cache/dog_ref/NC_051804.1_20000000_20006000.txt")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="NT masked-reconstruction smoke on canFam4.")
    ap.add_argument("--seq-file", default=str(_DEFAULT_SEQ))
    ap.add_argument("--window-bp", type=int, default=600, help="window fed to NT (<= max_context)")
    ap.add_argument("--n-tokens", type=int, default=40, help="how many 6-mer tokens to mask")
    ap.add_argument("--markov-k", type=int, default=5)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--date", default="2026-07-31")
    args = ap.parse_args(argv)

    full = Path(args.seq_file).read_text().strip().upper()
    window = full[: args.window_bp]
    if len(window) < args.window_bp:
        print(f"warn: sequence shorter than window ({len(window)} bp)", file=sys.stderr)

    print(f"loading NT-v2-100M ({args.device}) ...", flush=True)
    model = model_factory("nucleotide_transformer", device=args.device)

    positions = list(range(args.n_tokens))
    print(f"masking {args.n_tokens} tokens of a {len(window)} bp window; Markov fit on {len(full)} bp",
          flush=True)
    score = score_reconstruction(
        model, window, markov_k=args.markov_k, positions=positions,
        region_label="canfam4_chr1_NC_051804.1_smoke", markov_fit_sequences=[full],
    )
    d = score.as_dict()
    d["null_uniform_base_accuracy"] = 0.25
    d["substrate"] = args.seq_file
    d["real_surface"] = True

    out = Path(f"wiki/dog_masked_reconstruct_smoke_{args.date}.json")
    out.write_text(json.dumps(d, indent=2))

    print("\n=== NT-v2-100M masked reconstruction on canFam4 (smoke) ===")
    print(f"  HEADLINE delta (NT - Markov-{args.markov_k}): {d['headline_delta_lm_minus_markov']:+.4f}")
    print(f"  NT   base accuracy : {d['lm_base_accuracy']:.4f}   (token {d['lm_token_accuracy']:.4f}, "
          f"mean true-prob {d['lm_mean_true_prob']:.4f})")
    print(f"  Markov-{args.markov_k} base acc: {d['markov_base_accuracy']:.4f}")
    print(f"  null (uniform)     : 0.2500")
    print(f"  n_bases_scored={d['n_bases_scored']}  ->  {out}")
    print("  [raw accuracy is NOT the claim; the delta vs the cheap baseline is]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
