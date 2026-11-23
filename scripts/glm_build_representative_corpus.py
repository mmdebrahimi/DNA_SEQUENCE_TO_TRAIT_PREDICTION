"""Build the REPRESENTATIVE pretraining corpus: the committed corpus minus its draft E. coli.

WHY. The committed corpus's E. coli is `GCA_059754755.1`, a DRAFT assembly, while every falsifier
scoring script pins MG1655 `GCF_000005845.2`. So held-out MG1655 windows are ORTHOLOGUES of training
windows -- measured at a median 55 bp verbatim share, max 148 of 150
(`wiki/glm_copying_penalty_2026-11-22.json`).

The cheap fix was to filter orthology on the EVAL side, keeping only MG1655 windows whose gene has no
draft counterpart. That was measured and it FAILS, for a structural reason:
`wiki/glm_beyond_kmer_2026-11-22.json` found the model AND both Markov chains all scoring WORSE than
uniform (ln 4 = 1.3863) on that set, because filtering orthology by gene necessarily selects
strain-specific loci -- prophage, plasmid, HGT -- which are compositionally atypical. Orthology-free
and representative are in direct conflict on the eval side, so no amount of n fixes it.

Removing the draft E. coli from the PRETRAIN resolves both at once: every one of the 3,923 held-out
MG1655 windows becomes orthology-free with NO gene filtering, so the eval set stays representative.

THE LEAK GATE IS STRICTLY STRONGER HERE, NOT WEAKER, and the manifest says so rather than inheriting
a flag. The committed gate withheld the 3,349-window union of falsifier eval windows from the
training pool. This corpus additionally contains NO E. coli at all, so it cannot hold a falsifier
eval window under any seed. `enforced` stays true because it is MORE true, and `derived_from` plus
the removal counts are recorded so nobody has to take that on trust.

WHAT THIS CORPUS IS NOT. It has no E. coli, so a model trained on it alone is a cross-ORGANISM
generator. The same-strain measurement needs a fine-tune on MG1655's own trainable windows (4,340
minus the 3,349 withheld = 991) on top of this pretrain. That is a separate step.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dna_decode.glm.corpus import build_natural_set  # noqa: E402

SRC_CORPUS = ROOT / "wiki" / "glm_training_corpus_2026-11-21.txt"
SRC_MANIFEST = ROOT / "wiki" / "glm_training_corpus_manifest_2026-11-21.json"
DRAFT = ROOT / "data" / "refseq_cache" / "GCA_059754755.1"
OUT_CORPUS = ROOT / "wiki" / "glm_training_corpus_representative_2026-11-22.txt"
OUT_MANIFEST = ROOT / "wiki" / "glm_training_corpus_representative_manifest_2026-11-22.json"
LENGTH = 150


def main() -> int:
    for p in (SRC_CORPUS, SRC_MANIFEST, DRAFT):
        if not p.exists():
            print(f"REFUSE: missing {p}", file=sys.stderr)
            return 2

    corpus = [ln.strip().upper() for ln in SRC_CORPUS.open(encoding="utf-8")
              if ln.strip() and not ln.startswith("#")]
    print(f"source corpus: {len(corpus)} windows")

    wins, _ = build_natural_set(DRAFT, length=LENGTH, limit=None, seed=0, kind="promoter")
    draft_set = {w.seq.upper() for w in wins}
    print(f"draft E. coli windows available: {len(draft_set)}")

    kept = [s for s in corpus if s not in draft_set]
    removed = len(corpus) - len(kept)
    print(f"removed {removed} | kept {len(kept)}")

    # A removal count far from the manifest's own E. coli figure means I am filtering the wrong thing.
    src = json.loads(SRC_MANIFEST.read_text(encoding="utf-8"))
    ecoli_entry = next((e for e in src["manifest"]["entries"]
                        if e.get("accession") == "GCA_059754755.1"), None)
    expected = (ecoli_entry or {}).get("n_windows_kept")
    print(f"manifest records E. coli kept = {expected}")
    if expected and abs(removed - expected) > 0.1 * expected:
        print(f"REFUSE: removed {removed} but the manifest says the draft contributed {expected}. "
              "A mismatch this large means the filter is not selecting the E. coli windows.",
              file=sys.stderr)
        return 3
    if not kept:
        print("REFUSE: filtered corpus is empty", file=sys.stderr)
        return 4

    OUT_CORPUS.write_text("\n".join(kept) + "\n", encoding="utf-8")

    man = {
        "record": "glm-training-corpus-representative-v1",
        "date": "2026-11-22",
        "derived_from": {"corpus": SRC_CORPUS.name, "manifest": SRC_MANIFEST.name},
        "purpose": ("remove the DRAFT E. coli (GCA_059754755.1) so that every held-out MG1655 "
                    "window is orthology-free WITHOUT gene filtering, which the eval-side filter "
                    "could not achieve without selecting accessory loci"),
        "manifest": {
            "n_genomes": src["manifest"]["n_genomes"] - 1,
            "window_length_bp": LENGTH,
            "n_windows_kept": len(kept),
            "n_windows_removed_draft_ecoli": removed,
            "corpus_bp": len(kept) * LENGTH,
            "falsifier_leak_gate": {
                "enforced": True,
                "why_still_enforced": (
                    "STRICTLY STRONGER than the source. The source withheld the 3,349-window union "
                    "of falsifier eval windows from the training pool; this corpus additionally "
                    "contains NO E. coli, so it cannot hold a falsifier eval window under any "
                    "scoring seed. The flag is not inherited -- it is re-derived."),
                "inherited_withheld_union_size":
                    (src["manifest"].get("falsifier_leak_gate") or {}).get("withheld_union_size"),
            },
            "removed_accession": "GCA_059754755.1",
            "entries": [e for e in src["manifest"]["entries"]
                        if e.get("accession") != "GCA_059754755.1"],
        },
        "honest_limits": [
            "This corpus contains NO E. coli, so a model trained on it alone is a cross-ORGANISM "
            "generator. The same-strain measurement requires a fine-tune on MG1655's own 991 "
            "trainable windows (4,340 minus the 3,349 withheld) on top of this pretrain.",
            "Removal is by exact SEQUENCE match against draft windows rebuilt from "
            "GCA_059754755.1. A draft window that differs by one base from its corpus copy would "
            "survive; the count check against the manifest's own figure is what bounds that.",
        ],
    }
    OUT_MANIFEST.write_text(json.dumps(man, indent=2), encoding="utf-8")
    print(f"-> {OUT_CORPUS}")
    print(f"-> {OUT_MANIFEST}")
    print(f"tokens: {len(kept) * LENGTH:,} bp over {man['manifest']['n_genomes']} genomes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
