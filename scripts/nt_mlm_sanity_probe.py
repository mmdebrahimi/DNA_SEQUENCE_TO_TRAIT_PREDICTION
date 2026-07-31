"""Degeneracy guard for the NT masked-LM path (trust_remote_code version-drift check).

A trust_remote_code model can LOAD cleanly yet return near-random logits when the host transformers
differs from the one the vendored modeling code targets. This probe runs controls a WORKING DNA
masked-LM must pass, so a degenerate load fails loudly BEFORE any reconstruction number is trusted:

  (A) poly-A: mask a 6-mer inside "AAAA..." -> the true token "AAAAAA" must dominate (high prob).
  (B) poly-AT: mask inside "ATAT..." -> a periodic token must be strongly favored.
  (C) real canFam4 window token -> report true-token prob + entropy vs uniform.

VERDICT: DEGENERATE if the poly-A control's true-token prob < 0.5 (a working NT nails it near 1.0),
else OK. Runs in whatever env provides a compatible transformers (isolated NT-era env locally, or
Kaggle). Also logs library versions (the drift-diagnostic the lesson mandates).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.models.foundation import model_factory  # noqa: E402


def _probe_one(model, seq, token_index):
    preds = model.masked_token_predictions(seq, positions=[token_index])
    p = preds[0]
    return p


def main(argv=None) -> int:
    import transformers
    try:
        import huggingface_hub
        hub_v = huggingface_hub.__version__
    except Exception:  # noqa: BLE001
        hub_v = "?"
    print(f"transformers {transformers.__version__}  huggingface_hub {hub_v}", flush=True)

    print("loading NT-v2-100M (cpu) ...", flush=True)
    model = model_factory("nucleotide_transformer", device="cpu")
    vocab = None

    controls = [
        ("poly-A", "A" * 120, 5),
        ("poly-AT", "AT" * 60, 5),
    ]
    results = []
    for label, seq, ti in controls:
        p = _probe_one(model, seq, ti)
        results.append((label, p))
        print(f"  {label:8s} true={p.true_kmer} pred={p.pred_kmer} "
              f"true_prob={p.true_prob:.4f} pred_prob={p.pred_prob:.4f} "
              f"{'HIT' if p.true_kmer == p.pred_kmer else 'miss'}", flush=True)

    # real window (context) sanity
    seqfile = Path("D:/dna_decode_cache/dog_ref/NC_051804.1_20000000_20006000.txt")
    if seqfile.exists():
        real = seqfile.read_text().strip().upper()[:300]
        rp = _probe_one(model, real, 10)
        print(f"  real     true={rp.true_kmer} pred={rp.pred_kmer} "
              f"true_prob={rp.true_prob:.4f} pred_prob={rp.pred_prob:.4f}", flush=True)

    polya_true_prob = results[0][1].true_prob
    verdict = "OK" if polya_true_prob >= 0.5 else "DEGENERATE"
    print(f"\nVERDICT: {verdict}  (poly-A true-token prob = {polya_true_prob:.4f}; "
          f"a working NT nails poly-A near 1.0, degenerate ~ uniform {1/4104:.5f})")
    return 0 if verdict == "OK" else 3


if __name__ == "__main__":
    raise SystemExit(main())
