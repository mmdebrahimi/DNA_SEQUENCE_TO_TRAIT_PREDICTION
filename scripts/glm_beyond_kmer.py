"""Did the generator learn anything BEYOND 4-mer statistics? Diffusion vs an order-3 Markov chain.

THE WHOLE POINT. An order-3 Markov chain's sufficient statistic IS the 4-mer conditional
distribution. So the delta between the diffusion model's per-position cross-entropy and that chain's,
on IDENTICAL positions of identical held-out sequence, measures exactly one thing: whether the model
captured structure a 4-mer view cannot represent. That is the question the committed falsifier cannot
answer, because its k=4 features are saturated by the same null by construction (measured k-mer range
across Markov orders 1-5 is 0.3084 against positional 0.0344).

It needs no discriminator, so no power gate, no feature-mode choice, and none of `both` mode's 73%
dead weight. It is PAIRED per window, so it has power at small n. And it scores the MODEL rather than
a sample of 120 sequences, so sampling noise does not bound it.

TWO ARMS, BECAUSE THE OBVIOUS DESIGN IS INFORMATION-ASYMMETRIC. A masked-diffusion model is
BIDIRECTIONAL; a Markov chain sees `order` bases of LEFT context. Masking a contiguous suffix and
teacher-forcing the chain hands the chain true bases the diffusion model cannot see -- the mirror of
the asymmetry this project just caught in its own committed path, where a prompted generator was
compared against an unprompted null.

* ARM B (`left_only`, PRIMARY) -- mask position i AND everything right of it. Both models then
  predict position i from the same unmasked left prefix. A win here means the model uses longer LEFT
  context than three bases. Strict, and the information budget is matched by construction.
* ARM A (`bidirectional`) -- mask position i alone. The model sees both flanks, the chain sees three
  left bases. The asymmetry is REAL here and is not a confound: "uses more context than a 3-mer" IS
  the hypothesis. Reported second because it is the more favourable framing.

WHY THE COMPARISON IS INTERPRETABLE AT ALL. `wiki/glm_copying_penalty_2026-11-22.json` cleared this
checkpoint of copying (0 of 120 generated sequences share even 30 bp with the corpus), so a win here
cannot be recall. And the eval set is ORTHOLOGY-FREE -- MG1655 windows whose gene has no counterpart
in the draft assembly that is actually in the corpus -- so a win cannot be a memorised orthologue
either. Both of those had to land before this number meant anything.

HONEST LIMITS, stated here rather than discovered later.
* The eval set is accessory-gene biased (prophage / plasmid / HGT-derived), because strain-specific
  genes are. So a POSITIVE is conservative (the bias works against a model trained mostly on core
  promoters) while a NEGATIVE is confounded with out-of-distribution failure. Those are not
  symmetric and the verdict says so.
* Infilling and continuation skill are NOT generation quality. A model can score well here and
  sample badly. This is necessary and diagnostic, not sufficient.
* The chain is fit on the same corpus the model trained on, which is the fair pairing, and the eval
  windows are absent from it, so no leave-one-out correction applies.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dna_decode.glm.baselines import MarkovGenerator  # noqa: E402
from dna_decode.glm.corpus import build_natural_set  # noqa: E402
from dna_decode.glm.diffusion import MASK_ID, DiffusionConfig, build_model, encode  # noqa: E402
from scripts.glm_orthology_free_baseline import id_to_symbol, symbols  # noqa: E402

MG1655 = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
DRAFT_GFF = ROOT / "data" / "refseq_cache" / "GCA_059754755.1" / "annotations.gff3"
CORPUS = ROOT / "wiki" / "glm_training_corpus_2026-11-21.txt"
CKPT = ROOT / "data" / "glm_diffusion_ckpt" / "latest.pt"
OUT = ROOT / "wiki" / "glm_beyond_kmer_2026-11-22.json"

BASES = "ACGT"
ORDER = 3  # the committed markov_order; the whole comparison is defined against it


def markov_nll(chain: MarkovGenerator, seq: str, pos: int) -> float:
    """-log P(seq[pos] | seq[pos-ORDER:pos]) in NATS, using the chain's own fitted counts.

    Reads `_counts` directly rather than sampling from it: `_next` draws a base, which would measure
    sampling noise instead of the chain's assigned probability. Add-one smoothing is already in the
    table, so an unseen context backs off to uniform rather than to zero.
    """
    ctx = seq[pos - ORDER:pos]
    w = chain._counts.get(ctx) or [1, 1, 1, 1]
    tot = float(sum(w))
    p = w[BASES.index(seq[pos])] / tot
    return -math.log(p)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-windows", type=int, default=120)
    ap.add_argument("--positions", type=int, default=30, help="target positions sampled per window")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()

    import torch

    for p in (MG1655 / "annotations.gff3", DRAFT_GFF, CORPUS, CKPT):
        if not p.exists():
            print(f"REFUSING: missing {p}", file=sys.stderr)
            return 2

    # ---- the model ----
    blob = torch.load(CKPT, map_location="cpu", weights_only=False)
    cfg_d = dict(blob["config"])
    step = blob.get("step")
    cfg = DiffusionConfig(**{k: v for k, v in cfg_d.items()
                            if k in DiffusionConfig.__dataclass_fields__})
    model = build_model(cfg)
    model.load_state_dict(blob["model"])
    model.eval()
    print(f"checkpoint step {step}, length {cfg.length}, d_model {cfg.d_model}")

    # ---- the chain, fit on exactly what the model trained on ----
    corpus = [ln.strip().upper() for ln in CORPUS.open(encoding="utf-8")
              if ln.strip() and not ln.startswith("#")]
    cset = set(corpus)
    chain = MarkovGenerator(order=ORDER, seed=a.seed).fit(corpus)
    print(f"order-{ORDER} chain (POOLED, same data as the model) fit on {len(corpus)} windows")

    # THE STRONG NULL, and it is the one the hypothesis actually needs. Measured 2026-11-22: the
    # POOLED chain barely beats uniform (in-sample 1.3522 against ln(4)=1.3863, held-out MG1655
    # 1.3786) because the corpus spans 18 lineages at GC 0.24-0.72 and pooling flattens the
    # conditionals. "Beating an order-3 chain means beating 4-mers" only holds if the chain models
    # the distribution being SCORED -- otherwise the win is distribution shift, not order. So a
    # second chain is fit on exactly the E. COLI sequence the model saw.
    draft_wins, _ = build_natural_set(ROOT / "data" / "refseq_cache" / "GCA_059754755.1",
                                      length=cfg.length, limit=None, seed=0, kind="promoter")
    ecoli_fit = [w.seq.upper() for w in draft_wins if w.seq.upper() in cset]
    if len(ecoli_fit) < 500:
        print(f"REFUSING: only {len(ecoli_fit)} draft E. coli windows found in the corpus; the "
              "strong null would be fit on too little data to be the stronger of the two.",
              file=sys.stderr)
        return 5
    chain_ec = MarkovGenerator(order=ORDER, seed=a.seed).fit(ecoli_fit)
    print(f"order-{ORDER} chain (E. COLI-fit, the STRONG 4-mer null) fit on {len(ecoli_fit)} windows")

    # ---- the orthology-free held-out eval set ----
    only = symbols(MG1655 / "annotations.gff3") - symbols(DRAFT_GFF)
    wins, _ = build_natural_set(MG1655, length=cfg.length, limit=None, seed=0, kind="promoter")
    m = id_to_symbol(MG1655 / "annotations.gff3")
    held = [w for w in wins if w.seq.upper() not in cset]
    of = [w.seq.upper() for w in held if (m.get(w.gene_id) or "") in only]
    if len(of) < a.n_windows:
        print(f"REFUSING: {len(of)} orthology-free windows < requested {a.n_windows}", file=sys.stderr)
        return 3
    eval_seqs = of[:a.n_windows]
    print(f"eval: {len(eval_seqs)} orthology-free held-out MG1655 windows "
          f"(of {len(held)} held-out, {len(only)} MG1655-only symbols)")

    rng = random.Random(a.seed)
    per_window: list[dict] = []
    t0 = time.time()

    for wi, seq in enumerate(eval_seqs):
        # Positions need a full order-3 left context, so they start at ORDER. The SAME positions
        # serve every arm and both nulls -- that identity is the whole basis of the pairing.
        pos = sorted(rng.sample(range(ORDER, cfg.length), a.positions))
        ids = encode([seq], cfg.length)[0]

        row = {
            "markov_pooled": statistics.fmean(markov_nll(chain, seq, p) for p in pos),
            "markov_ecoli": statistics.fmean(markov_nll(chain_ec, seq, p) for p in pos),
            "uniform": math.log(4),
        }
        for arm in ("left_only", "bidirectional"):
            batch = ids.repeat(len(pos), 1).clone()
            for j, p in enumerate(pos):
                if arm == "left_only":
                    batch[j, p:] = MASK_ID      # nothing to the right is visible
                else:
                    batch[j, p] = MASK_ID       # both flanks visible
            with torch.no_grad():
                logits = model(batch)
            lp = torch.log_softmax(logits.float(), dim=-1)
            row[arm] = statistics.fmean(-float(lp[j, p, int(ids[p])]) for j, p in enumerate(pos))
        per_window.append(row)

        if (wi + 1) % 20 == 0:
            print(f"  {wi + 1}/{len(eval_seqs)} ({time.time() - t0:.0f}s)", flush=True)

    # ---- paired statistics ----
    def summarise(arm: str, null: str):
        pairs = [(w[arm], w[null]) for w in per_window]
        d = [mk - ce for ce, mk in pairs]          # positive = diffusion BETTER (lower CE)
        boot = []
        r = random.Random(a.seed + 1)
        for _ in range(a.boot):
            s = [d[r.randrange(len(d))] for _ in range(len(d))]
            boot.append(statistics.fmean(s))
        boot.sort()
        lo, hi = boot[int(0.025 * len(boot))], boot[int(0.975 * len(boot))]
        return {
            "arm": arm, "null": null,
            "n_windows": len(pairs),
            "diffusion_ce_nats": round(statistics.fmean(c for c, _ in pairs), 4),
            "null_ce_nats": round(statistics.fmean(mk for _, mk in pairs), 4),
            "paired_delta_nats": round(statistics.fmean(d), 4),
            "ci95": [round(lo, 4), round(hi, 4)],
            "resolved": bool(lo > 0 or hi < 0),
            "n_windows_diffusion_better": sum(1 for x in d if x > 0),
        }

    CELLS = [(arm, null) for arm in ("left_only", "bidirectional")
             for null in ("markov_ecoli", "markov_pooled", "uniform")]

    res = {
        "record": "glm-beyond-kmer-v1",
        "date": "2026-11-22",
        "question": "does the diffusion model capture structure an order-3 Markov chain cannot?",
        "why_order3": "an order-3 chain's sufficient statistic IS the 4-mer conditional distribution",
        "checkpoint_step": step,
        "eval": {"set": "orthology-free held-out MG1655 promoter windows",
                 "n_windows": len(eval_seqs), "positions_per_window": a.positions,
                 "strain_note": "corpus E. coli is GCA_059754755.1 (draft); eval is "
                                "GCF_000005845.2 (MG1655), and eval genes have NO draft counterpart",
                 "copying_cleared_by": "wiki/glm_copying_penalty_2026-11-22.json"},
        "cells": {f"{arm}|{null}": summarise(arm, null) for arm, null in CELLS},
        "primary_cell": "left_only|markov_ecoli",
        "why_primary": ("left_only matches the information budget (both see only the unmasked left "
                        "prefix) and markov_ecoli is the STRONG 4-mer null, fit on the E. coli the "
                        "model actually trained on. The pooled chain is reported because it is the "
                        "same-data pairing, but it barely beats uniform and is the WEAK null."),
        "seconds": round(time.time() - t0, 1),
    }

    pri = res["cells"]["left_only|markov_ecoli"]

    # UNIFORM SANITY GATE, evaluated BEFORE any delta is interpreted. Measured 2026-11-22: on the
    # orthology-free (accessory-gene) eval set the model AND both chains all score WORSE than
    # ln(4) -- a model confident about the wrong base loses to one that commits to nothing. A
    # positive delta between two sub-chance models is not evidence of learning, so reporting it
    # would be the same error class as a broken detector reading as a clean result.
    uni_cell = res["cells"]["left_only|uniform"]
    UNI = math.log(4)
    null_beats_uniform = pri["null_ce_nats"] < UNI
    model_beats_uniform = uni_cell["paired_delta_nats"] > 0
    if not (model_beats_uniform and null_beats_uniform):
        res["verdict"] = "EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM"
        res["verdict_detail"] = (
            f"On this eval set the model scores {uni_cell['diffusion_ce_nats']} nats and the strong "
            f"null {pri['null_ce_nats']}, both against uniform {UNI:.4f} -- so "
            f"{'neither' if not model_beats_uniform and not null_beats_uniform else 'at least one'} "
            "beats a random guess. The model-vs-null delta is therefore a comparison between "
            "sub-chance predictors and carries NO information about structure beyond 4-mers. "
            "REFUSING a verdict rather than reporting the delta. This is the accessory-gene "
            "out-of-distribution failure the staged plan predicted: filtering orthology by gene "
            "necessarily selects strain-specific (prophage / plasmid / HGT) loci, so no eval-side "
            "shortcut can be both orthology-free and representative. The representative "
            "re-pretrain -- excluding draft E. coli so that ALL MG1655 windows are orthology-free "
            "without any gene filtering -- is the only route to an informative measurement.")
    elif pri["resolved"] and pri["paired_delta_nats"] > 0:
        res["verdict"] = "BEYOND_KMER_CONSERVATIVE_POSITIVE"
        res["verdict_detail"] = (
            f"On matched left-only context the model beats the order-3 chain by "
            f"{pri['paired_delta_nats']} nats, CI {pri['ci95']}, on {pri['n_windows_diffusion_better']} "
            f"of {pri['n_windows']} windows. It captured structure a 4-mer view cannot represent. "
            "CONSERVATIVE: the eval set is accessory-gene biased, which works AGAINST a model "
            "trained mostly on core promoters, and copying plus orthology are both excluded.")
    elif pri["resolved"]:
        res["verdict"] = "WORSE_THAN_ORDER3_ON_THIS_SET"
        res["verdict_detail"] = (
            f"The chain beats the model by {-pri['paired_delta_nats']} nats, CI {pri['ci95']}. Do NOT "
            "close on this alone: the eval set is accessory-gene biased, so this is confounded with "
            "out-of-distribution failure rather than absence of signal beyond 4-mers. The "
            "representative re-pretrain is what settles it.")
    else:
        res["verdict"] = "UNRESOLVED"
        res["verdict_detail"] = (
            f"Paired delta {pri['paired_delta_nats']} nats with CI {pri['ci95']} spanning zero. "
            "Underpowered or genuinely equivalent; neither reading is supported over the other.")

    print("\n" + "=" * 86)
    print(f"{'arm vs null':<34}{'model':>8}{'null':>8}{'delta':>9}{'CI95':>20}{'better':>8}")
    print("-" * 86)
    for arm, null in CELLS:
        s = res["cells"][f"{arm}|{null}"]
        tag = " *" if f"{arm}|{null}" == res["primary_cell"] else "  "
        print(f"{arm}|{null}{tag:<4}"[:34].ljust(34)
              + f"{s['diffusion_ce_nats']:>8}{s['null_ce_nats']:>8}"
              + f"{s['paired_delta_nats']:>9}"
              + f"  [{s['ci95'][0]:>7}, {s['ci95'][1]:>7}]"
              + f"{s['n_windows_diffusion_better']:>5}/{s['n_windows']}")
    print("-" * 86)
    print("* = primary. positive delta = model BETTER (lower cross-entropy).")
    print("=" * 86)
    print(f"VERDICT: {res['verdict']}\n  {res['verdict_detail']}")
    OUT.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
