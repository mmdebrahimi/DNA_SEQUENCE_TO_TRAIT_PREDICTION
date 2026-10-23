"""The LoRA fine-tune's train corpus, PROVEN disjoint from the falsifier's eval half.

    uv run python scripts/glm_finetune_window_split.py            # emit the manifest
    uv run python scripts/glm_finetune_window_split.py --check     # verify only, write nothing

Family G-D's refined goal requires fine-tuning to bring distinguishability below the 3-mer Markov
null's 0.6679 **"on the frozen natural comparison set, with the corpus fingerprint stamped so the
comparison is auditable."** Both halves of that sentence are load-bearing, and they are what this
script exists to make true:

  * **FROZEN comparison set.** The falsifier scores every candidate on its EVAL half and states the
    contract in its own output: *"no generator is scored against sequences it was trained on."* A
    LoRA run that touched even one eval window would make the post-fine-tune number meaningless --
    and it would look like a WIN, because a generator reproduces what it memorised. So the train set
    is derived by EXCLUSION from eval, and the disjointness is asserted on content, not assumed from
    index arithmetic.
  * **Stamped fingerprint.** The corpus is rebuilt here and its fingerprint compared against the one
    recorded in the baseline artifact. A silent corpus drift (a different GFF3, a different window
    rule) would move the eval set underneath the comparison, so a mismatch REFUSES rather than
    warning.

**This script does NOT build its own window loader, deliberately.** `glm_generator_falsifier`
already owns this resource, and a second reader of the same files is the recorded way a solved trap
gets re-derived -- the guard lives inside the function that was not called. So `build_natural` is
IMPORTED and its split is REPRODUCED, which also means any future change to the window rule moves
both sides together instead of silently separating them.

**The finding that makes this worth running.** The falsifier takes `--limit-natural 1200` of the
corpus's 4,340 upstream windows -- 600 fit + 600 eval -- and the remaining ~3,140 are UNUSED by it.
Those are disjoint from eval by construction, so the fine-tune corpus is not limited to the 600-window
fit half: it is `fit + unused`, roughly 6x more data, at zero leakage risk. That is a property of the
existing split rather than a new design, which is exactly why it was worth measuring before writing a
training script against the smaller assumption.

Stdlib + the project's own corpus loader only. No torch, no GPU, no network, and nothing on an
external volume: the MG1655 reference is local at `data/cache/refseq/GCF_000005845.2/`.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
BASELINE = ROOT / "wiki" / "glm_generator_falsifier_REAL_2026-10-07.json"
OUT = ROOT / "wiki" / "glm_finetune_window_split.json"

#: The falsifier's live geometry, as the LAUNCHER invokes it -- not its argparse defaults.
#: `run_glm_anchor_detached.bat` passes `--limit-natural 1200`; the parser's own default is 2000.
#: Reproducing the defaults instead of the invocation would build a DIFFERENT split and the
#: disjointness proof would be about the wrong set.
LENGTH_BP = 150
LIMIT_NATURAL = 1200
SEED = 7
KIND = "promoter"


def _sha(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii", "replace")).hexdigest()[:16]


def _as_dict(stats) -> dict:
    """`build_natural` returns a `CorpusStats` DATACLASS, not a dict.

    Normalised here rather than reached through a getattr chain at each use, so the comparison
    against the recorded JSON fingerprint is dict-to-dict and a renamed field surfaces as a missing
    key instead of silently reading `None` and comparing equal to nothing.
    """
    import dataclasses

    if dataclasses.is_dataclass(stats):
        return dataclasses.asdict(stats)
    if isinstance(stats, dict):
        return stats
    return dict(vars(stats))


def load_windows() -> tuple[list, dict]:
    """Every upstream window, via the module that OWNS this resource."""
    sys.path.insert(0, str(SCRIPTS))
    import glm_generator_falsifier as gf  # lazy third-party imports; safe without torch

    # build_natural shuffles by seed then truncates to `limit`. Pass limit=None to get the FULL
    # corpus in the falsifier's own post-shuffle order, so a prefix of it IS its split.
    wins, stats = gf.build_natural(LENGTH_BP, None, SEED, kind=KIND)
    return wins, _as_dict(stats)


def _seq(w) -> str:
    """A window may be a bare string or a record; take the sequence either way."""
    if isinstance(w, str):
        return w
    for attr in ("seq", "sequence"):
        v = getattr(w, attr, None)
        if isinstance(v, str):
            return v
    if isinstance(w, dict):
        for k in ("seq", "sequence"):
            if isinstance(w.get(k), str):
                return w[k]
    raise TypeError(f"cannot read a sequence from {type(w).__name__}")


def split(wins: list) -> dict:
    """Reproduce the falsifier's split, then derive the train set BY EXCLUSION from eval."""
    used = wins[:LIMIT_NATURAL]
    half = len(used) // 2
    fit, evl = used[:half], used[half:]
    unused = wins[LIMIT_NATURAL:]

    eval_seqs = {_seq(w) for w in evl}
    train = [w for w in (fit + unused) if _seq(w) not in eval_seqs]

    return {
        "fit": fit, "eval": evl, "unused": unused, "train": train,
        "eval_seqs": eval_seqs,
    }


def verify(s: dict, stats: dict, baseline: dict) -> tuple[list[str], list[str]]:
    """`(failures, notes)`. A failure REFUSES; the comparison would not be auditable."""
    fail: list[str] = []
    notes: list[str] = []

    # --- the disjointness proof, on CONTENT and not on index arithmetic -------------------------
    train_seqs = {_seq(w) for w in s["train"]}
    overlap = train_seqs & s["eval_seqs"]
    if overlap:
        fail.append(
            f"LEAKAGE: {len(overlap)} window(s) appear in BOTH the fine-tune train set and the "
            f"falsifier's eval half. A LoRA run on these would look like a win because a generator "
            f"reproduces what it memorised."
        )
    else:
        notes.append(f"disjoint: 0 of {len(train_seqs)} train windows appear in the {len(s['eval_seqs'])} eval windows")

    # --- the fingerprint, against the recorded baseline ------------------------------------------
    rec = (baseline.get("natural_sets") or {}).get(KIND) or {}
    for k in ("n_cds_total", "n_windows", "n_overlapping_cds"):
        if k in rec and stats.get(k) != rec[k]:
            fail.append(f"corpus drift: {k} is {stats.get(k)}, baseline recorded {rec[k]}")
    rs, bs = (stats.get("strand_counts") or {}), (rec.get("strand_counts") or {})
    if bs and rs != bs:
        fail.append(f"corpus drift: strand_counts {rs} != baseline {bs}")
    if not fail:
        notes.append(f"fingerprint matches the baseline: n_windows={stats.get('n_windows')}, "
                     f"strand={stats.get('strand_counts')}, overlapping={stats.get('n_overlapping_cds')}")

    # --- the split the falsifier itself recorded ---------------------------------------------------
    fes = baseline.get("fit_eval_split") or {}
    if fes.get("n_fit") is not None and len(s["fit"]) != fes["n_fit"]:
        fail.append(f"fit half is {len(s['fit'])}, baseline recorded n_fit={fes['n_fit']}")
    if fes.get("n_eval") is not None and len(s["eval"]) != fes["n_eval"]:
        fail.append(f"eval half is {len(s['eval'])}, baseline recorded n_eval={fes['n_eval']}")

    if not s["unused"]:
        notes.append("no unused tail -- the train set is the fit half alone")
    else:
        notes.append(f"train set is fit({len(s['fit'])}) + unused tail({len(s['unused'])}) = "
                     f"{len(s['train'])}, i.e. {len(s['train']) / max(1, len(s['fit'])):.1f}x the fit "
                     f"half, all disjoint from eval")
    return fail, notes


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="verify only; write no manifest")
    a = ap.parse_args(argv)

    if not BASELINE.exists():
        print(f"REFUSE: baseline {BASELINE.name} absent; the fingerprint cannot be checked")
        return 3
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))

    wins, stats = load_windows()
    s = split(wins)
    fail, notes = verify(s, stats, baseline)

    print(f"corpus: {len(wins)} upstream windows of {LENGTH_BP} bp (seed {SEED}, kind {KIND})")
    print(f"falsifier split: fit {len(s['fit'])} / eval {len(s['eval'])} "
          f"(from --limit-natural {LIMIT_NATURAL})")
    print(f"fine-tune train: {len(s['train'])}  |  HELD OUT (never trained on): {len(s['eval'])}")
    for n in notes:
        print(f"  note   : {n}")
    for f in fail:
        print(f"  REFUSE : {f}")
    if fail:
        return 2

    if a.check:
        print("--check: verified, nothing written")
        return 0

    OUT.write_text(json.dumps({
        "schema": "glm-finetune-window-split-v1",
        "date": str(date.today()),
        "purpose": ("the LoRA fine-tune's train corpus for family G-D, derived BY EXCLUSION from the "
                    "falsifier's eval half so the frozen comparison set stays frozen"),
        "geometry": {"length_bp": LENGTH_BP, "limit_natural": LIMIT_NATURAL, "seed": SEED,
                     "kind": KIND,
                     "note": "limit_natural is the LAUNCHER's 1200, not the parser's 2000 default"},
        "owner_module": "scripts/glm_generator_falsifier.py::build_natural (imported, not reimplemented)",
        "corpus_fingerprint": stats,
        "fingerprint_matches_baseline": True,
        "baseline": BASELINE.name,
        "counts": {"corpus_total": len(wins), "falsifier_fit": len(s["fit"]),
                   "falsifier_eval": len(s["eval"]), "unused_tail": len(s["unused"]),
                   "finetune_train": len(s["train"])},
        "disjointness": {"method": "set intersection on window SEQUENCE content, not on indices",
                         "train_vs_eval_overlap": 0,
                         "why": ("a LoRA run that touched an eval window would look like a WIN, "
                                 "because a generator reproduces what it memorised")},
        "held_out_eval_sha16": sorted(_sha(_seq(w)) for w in s["eval"]),
        "finetune_train_sha16": sorted(_sha(_seq(w)) for w in s["train"]),
        "honest_limits": [
            "This is a SPLIT, not a fine-tune. No model was loaded and no training was run.",
            "Disjointness is proven for the CONTENT of these windows. Two distinct genomic loci "
            "with identical 150 bp sequence would collapse to one key -- that is the conservative "
            "direction (it can only move a window OUT of train), but it means the train count is a "
            "lower bound rather than an exact partition.",
            "The bar it serves (below 0.6679, mode `both`) is the falsifier's, unchanged.",
        ],
    }, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
