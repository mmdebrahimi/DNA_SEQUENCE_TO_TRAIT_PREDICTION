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
from dna_decode.glm.eval_admissibility import (  # noqa: E402
    REFUSED_BOTH_BELOW_UNIFORM,
    admissibility,
)
from scripts.glm_orthology_free_baseline import id_to_symbol, symbols  # noqa: E402

MG1655 = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
DRAFT_DIR = ROOT / "data" / "refseq_cache" / "GCA_059754755.1"
DRAFT_GFF = DRAFT_DIR / "annotations.gff3"
CORPUS = ROOT / "wiki" / "glm_training_corpus_2026-11-21.txt"
CKPT = ROOT / "data" / "glm_diffusion_ckpt" / "latest.pt"
OUT = ROOT / "wiki" / "glm_beyond_kmer_2026-11-22.json"

BASES = "ACGT"
ORDER = 3  # the committed markov_order; the whole comparison is defined against it

# The geometry that DEFINES the designed withhold (set A), from the `geometry` block of
# wiki/glm_training_corpus_manifest_2026-11-21.json. Only read on a `--eval-set` other than the
# default, because recomputing set A costs ten natural-set builds.
SCORING_LIMIT = 1200
N_SEEDS = 10

#: The selectable eval sets. `orthology_free` is the committed path and is deliberately FIRST so it
#: stays the argparse default, leaving the original invocation byte-unchanged.
EVAL_SETS = ("orthology_free", "sequence_withheld", "ortholog_bearing", "withheld_orthology_free")

SET_LABELS = {
    "orthology_free": "orthology-free held-out MG1655 promoter windows",
    "sequence_withheld": "sequence-withheld MG1655 promoter windows (the designed withhold)",
    "ortholog_bearing": "withheld MG1655 promoter windows WHOSE GENE HAS a draft counterpart",
    "withheld_orthology_free": "withheld MG1655 promoter windows whose gene has NO draft counterpart",
}

# The committed artifact hardcoded the orthology-free note, which is FALSE for two of these sets --
# `ortholog_bearing` is its exact opposite. A note that does not describe the set it ships with is
# worse than no note, so each set carries its own.
SET_STRAIN_NOTES = {
    "orthology_free": ("corpus E. coli is GCA_059754755.1 (draft); eval is GCF_000005845.2 "
                       "(MG1655), and eval genes have NO draft counterpart"),
    "sequence_withheld": ("corpus E. coli is GCA_059754755.1 (draft); eval is GCF_000005845.2 "
                          "(MG1655). Withheld BY SEQUENCE, so a 95%-identical draft ortholog may "
                          "remain in training -- this set is NOT orthology-free and MIXES windows "
                          "that have a draft counterpart with windows that do not"),
    "ortholog_bearing": ("corpus E. coli is GCA_059754755.1 (draft); eval is GCF_000005845.2 "
                         "(MG1655). Every window here is withheld by sequence BUT its gene DOES "
                         "have a draft counterpart, so a near-ortholog is in training by "
                         "construction. This is the EXPOSED arm of the leak contrast"),
    "withheld_orthology_free": ("corpus E. coli is GCA_059754755.1 (draft); eval is "
                                "GCF_000005845.2 (MG1655). Withheld by sequence AND the gene has no "
                                "draft counterpart. The CLEAN arm of the leak contrast"),
}


def load_corpus_windows(path: Path = CORPUS) -> list[str]:
    """The training corpus windows, uppercased -- the ONE way every consumer reads this file.

    Extracted from `main` so the set-inventory probe reads the corpus the scorer scores against
    rather than a second description of it. A probe that strips differently sees a different
    `held` set, and `held` is the basis of every candidate set below.
    """
    return [ln.strip().upper() for ln in path.open(encoding="utf-8")
            if ln.strip() and not ln.startswith("#")]


def ecoli_fit_windows(corpus_set: set[str], *, length: int) -> list[str]:
    """Draft E. coli promoter windows that are IN the corpus: the STRONG null's fitting set.

    NOTE the asymmetry this creates, because it is load-bearing for reading any two-set comparison:
    these windows come from the DRAFT assembly and are IN the corpus, while every candidate eval set
    below is MG1655 windows ABSENT from it. So the strong null is never fit on the windows it scores
    -- but it HAS seen the draft near-orthologs of any eval window that has one, which means a set
    retaining near-orthologs strengthens the null as well as the model.
    """
    draft_wins, _ = build_natural_set(DRAFT_DIR, length=length, limit=None, seed=0, kind="promoter")
    return [w.seq.upper() for w in draft_wins if w.seq.upper() in corpus_set]


def mg1655_candidate_sets(corpus_set: set[str], *, length: int) -> dict:
    """Set B (`held`) and set C (`of`) in ONE computation, because that is how they are defined.

    `held` = MG1655 promoter windows whose sequence is absent from the corpus TEXT. `of` = those
    whose gene symbol has no counterpart in the draft assembly. Returned together with the symbol
    inputs so a caller can report the filter's cost without recomputing anything.
    """
    only = symbols(MG1655 / "annotations.gff3") - symbols(DRAFT_GFF)
    wins, _ = build_natural_set(MG1655, length=length, limit=None, seed=0, kind="promoter")
    m = id_to_symbol(MG1655 / "annotations.gff3")
    held = [w for w in wins if w.seq.upper() not in corpus_set]
    of = [w.seq.upper() for w in held if (m.get(w.gene_id) or "") in only]
    return {"only": only, "all_windows": wins, "held": held, "orthology_free": of}


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


class EvalSetUnavailable(RuntimeError):
    """The requested set holds fewer windows than were asked for."""


def select_eval_windows(eval_set: str, *, of_seqs, withheld_seqs, n_windows: int,
                        selection_seed: int) -> tuple[list[str], dict]:
    """Pick the windows to score, and say exactly how they were picked.

    `orthology_free` keeps `of_seqs[:n]` -- MG1655 coordinate order, the committed slice -- because
    the reproduction pin depends on that exact selection and must not move.

    Every OTHER set is a SEEDED RANDOM SAMPLE of its sorted distinct sequences. Taking the first n
    of a sorted order would be a biased sample of DNA (a lexicographic prefix over-selects A-rich
    windows), and the sets differ in size, so an unseeded or order-dependent pick would make a
    between-set difference partly a sampling difference. Deduplicated first: the orthology-free LIST
    carries 14 repeated sequences and scoring one twice would inflate n and break the paired
    bootstrap's independence.
    """
    if eval_set not in EVAL_SETS:
        raise EvalSetUnavailable(f"unknown eval set {eval_set!r}; expected one of {EVAL_SETS}")
    # A non-positive request must REFUSE, not slice. `n_windows=0` returned an empty eval set and
    # `-1` returned `pool[:-1]` -- every window but the last -- both silently, and a downstream
    # artifact would then report whatever it happened to score as the run.
    if n_windows <= 0:
        raise EvalSetUnavailable(f"n_windows must be positive, got {n_windows}")

    if eval_set == "orthology_free":
        pool: list[str] = list(of_seqs)
        how = "the first n in MG1655 coordinate order (the committed slice, preserved)"
        sampled = False
    else:
        a = {s.upper() for s in withheld_seqs}
        c = {s.upper() for s in of_seqs}
        chosen = {"sequence_withheld": a,
                  "ortholog_bearing": a - c,
                  "withheld_orthology_free": a & c}[eval_set]
        pool = sorted(chosen)
        how = f"a seeded random sample (seed {selection_seed}) of the sorted distinct sequences"
        sampled = True

    if len(pool) < n_windows:
        raise EvalSetUnavailable(
            f"{eval_set}: {len(pool)} windows available < requested {n_windows}")

    seqs = (random.Random(selection_seed).sample(pool, n_windows) if sampled
            else pool[:n_windows])
    return seqs, {
        "set": eval_set,
        "n_available": len(pool),
        "n_scored": len(seqs),
        "selection": how,
        "selection_seed": selection_seed if sampled else None,
        "deduplicated": sampled,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-windows", type=int, default=120)
    ap.add_argument("--positions", type=int, default=30, help="target positions sampled per window")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--eval-set", choices=EVAL_SETS, default="orthology_free",
                    help="which candidate set to score; the default is the committed path")
    ap.add_argument("--selection-seed", type=int, default=0,
                    help="seeds the window sample for every set except the committed default")
    ap.add_argument("--out", type=Path, default=OUT,
                    help="artifact path; the default OVERWRITES the committed 2026-11-22 artifact")
    a = ap.parse_args()
    RUN_DATE = time.strftime("%Y-%m-%d")

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
    corpus = load_corpus_windows()
    cset = set(corpus)
    chain = MarkovGenerator(order=ORDER, seed=a.seed).fit(corpus)
    print(f"order-{ORDER} chain (POOLED, same data as the model) fit on {len(corpus)} windows")

    # THE STRONG NULL, and it is the one the hypothesis actually needs. Measured 2026-11-22: the
    # POOLED chain barely beats uniform (in-sample 1.3522 against ln(4)=1.3863, held-out MG1655
    # 1.3786) because the corpus spans 18 lineages at GC 0.24-0.72 and pooling flattens the
    # conditionals. "Beating an order-3 chain means beating 4-mers" only holds if the chain models
    # the distribution being SCORED -- otherwise the win is distribution shift, not order. So a
    # second chain is fit on exactly the E. COLI sequence the model saw.
    ecoli_fit = ecoli_fit_windows(cset, length=cfg.length)
    if len(ecoli_fit) < 500:
        print(f"REFUSING: only {len(ecoli_fit)} draft E. coli windows found in the corpus; the "
              "strong null would be fit on too little data to be the stronger of the two.",
              file=sys.stderr)
        return 5
    chain_ec = MarkovGenerator(order=ORDER, seed=a.seed).fit(ecoli_fit)
    print(f"order-{ORDER} chain (E. COLI-fit, the STRONG 4-mer null) fit on {len(ecoli_fit)} windows")

    # ---- the eval set ----
    _sets = mg1655_candidate_sets(cset, length=cfg.length)
    only, held, of = _sets["only"], _sets["held"], _sets["orthology_free"]

    # Set A costs ten natural-set builds, so it is only computed when a set that needs it is asked
    # for. The default path never touches it and so is unchanged in cost as well as in behaviour.
    withheld: set[str] = set()
    if a.eval_set != "orthology_free":
        from dna_decode.glm.training_corpus import falsifier_eval_windows
        withheld = falsifier_eval_windows(MG1655, length=cfg.length, limit=SCORING_LIMIT,
                                          seeds=range(N_SEEDS), kind="promoter")
        print(f"designed withhold (set A): {len(withheld)} windows")

    try:
        eval_seqs, eval_prov = select_eval_windows(
            a.eval_set, of_seqs=of, withheld_seqs=withheld,
            n_windows=a.n_windows, selection_seed=a.selection_seed)
    except EvalSetUnavailable as exc:
        print(f"REFUSING: {exc}", file=sys.stderr)
        return 3
    print(f"eval: {len(eval_seqs)} windows from {a.eval_set} "
          f"({eval_prov['n_available']} available; {eval_prov['selection']})")
    print(f"      context: {len(held)} held-out, {len(only)} MG1655-only symbols")

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
        "date": RUN_DATE,
        "question": "does the diffusion model capture structure an order-3 Markov chain cannot?",
        "why_order3": "an order-3 chain's sufficient statistic IS the 4-mer conditional distribution",
        "checkpoint_step": step,
        "eval": {"set": SET_LABELS[a.eval_set],
                 "set_key": a.eval_set,
                 "n_windows": len(eval_seqs), "positions_per_window": a.positions,
                 "n_predictions": len(eval_seqs) * a.positions,
                 "strain_note": SET_STRAIN_NOTES[a.eval_set],
                 "provenance": eval_prov,
                 "copying_cleared_by": "wiki/glm_copying_penalty_2026-11-22.json"},
        # The committed artifact recorded NO seed and NO boot, which is why its confidence intervals
        # cannot be reproduced from it while its point estimates can. Recorded here so this one can.
        "run_args": {"n_windows": a.n_windows, "positions": a.positions, "seed": a.seed,
                     "boot": a.boot, "eval_set": a.eval_set,
                     "selection_seed": a.selection_seed},
        "cells": {f"{arm}|{null}": summarise(arm, null) for arm, null in CELLS},
        # PER-WINDOW rows, persisted so a BETWEEN-SET contrast is computable without re-running.
        # The committed artifact kept only cell summaries, which makes a two-sample bootstrap on a
        # difference of paired deltas impossible after the fact -- and two cell CIs overlapping is
        # not that test: two intervals can both span zero while their difference does not.
        "per_window": [{k: round(v, 6) for k, v in row.items()} for row in per_window],
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
    # The gate now has exactly ONE definition, in `dna_decode.glm.eval_admissibility`. The model
    # side is passed as the uniform cell's PAIRED DELTA rather than its cross-entropy, because that
    # is the value this block originally compared and it is rounded to 4dp BEFORE comparison --
    # passing the delta makes the extraction bit-identical rather than merely equivalent.
    uni_cell = res["cells"]["left_only|uniform"]
    UNI = math.log(4)
    gate = admissibility(null_ce_nats=pri["null_ce_nats"],
                         model_ce_nats=uni_cell["diffusion_ce_nats"],
                         model_uniform_delta_nats=uni_cell["paired_delta_nats"])
    res["admissibility"] = gate.to_dict()
    if not gate.admissible:
        res["verdict"] = "EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM"
        res["verdict_detail"] = (
            f"On this eval set the model scores {uni_cell['diffusion_ce_nats']} nats and the strong "
            f"null {pri['null_ce_nats']}, both against uniform {UNI:.4f} -- so "
            f"{'neither' if gate.state == REFUSED_BOTH_BELOW_UNIFORM else 'at least one'} "
            "beats a random guess. The model-vs-null delta is therefore a comparison between "
            "sub-chance predictors and carries NO information about structure beyond 4-mers. "
            "REFUSING a verdict rather than reporting the delta. This is the accessory-gene "
            "out-of-distribution failure the staged plan predicted: filtering orthology by gene "
            "necessarily selects strain-specific (prophage / plasmid / HGT) loci, so no eval-side "
            "shortcut can be both orthology-free and representative. The representative "
            "re-pretrain -- excluding draft E. coli so that ALL MG1655 windows are orthology-free "
            "without any gene filtering -- is the only route to an informative measurement.")
        # The paragraph above is preserved VERBATIM because the committed artifact carries it and
        # the reproduction pin depends on it. But it attributes the refusal to the EVAL SET, and
        # that attribution is only available when the NULL also failed. When the null clears
        # uniform and only the model does not, the set is demonstrably predictable and the
        # accessory-gene story does not explain the result -- so say so beside the text rather than
        # editing it. Measured 2026-11-24: the strong null reads 1.3571 on all 648 windows of this
        # same set, so this branch is reachable and is not hypothetical.
        if gate.state != REFUSED_BOTH_BELOW_UNIFORM:
            res["verdict_detail_scope_note"] = (
                f"The paragraph above is the committed 2026-11-22 wording, kept verbatim. It does "
                f"NOT account for THIS run: the gate state is {gate.state}, i.e. the strong null "
                f"CLEARED uniform on this set and only the model did not. That is a MODEL negative "
                f"on a set shown to be predictable, not an uninformative eval set, so the "
                f"accessory-gene out-of-distribution explanation and the claim that a re-pretrain "
                f"is 'the only route' do not follow from it. See "
                f"wiki/glm_eval_set_inventory_2026-11-23.json.")
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
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\n-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
