"""Extract (claim, cited artifact) pairs from CLAUDE.md for the semantic staleness audit.

WHY CLAUDE.md FIRST. It is loaded into EVERY session, so one stale line there misdirects every future run
-- which is not hypothetical: five separate stale claims in it have been caught and pinned as regression
tests, and each was repeated to readers until someone happened to notice. The 542-memo `wiki/` corpus is
larger but each memo is read rarely and is DATED (a memo records a moment, so "stale" barely applies).
CLAUDE.md is the live surface. Auditing it first is the highest value per token.

WHAT A "CLAIM" IS HERE. A bullet in CLAUDE.md that cites at least one repo path. The bullet is the claim;
the path is the artifact it is about. A bullet citing several paths yields several pairs -- each is judged
separately, because a bullet can be accurate about one artifact and stale about another (the ProSST case
was exactly that shape).

WHAT THIS DELIBERATELY DOES NOT DO. It does not judge anything. Extraction is mechanical and PURE so the
model's input is reproducible and the pairing logic is testable without a GPU. The judging happens in
`scripts/kaggle_staleness_auditor.py`; the adjudication happens by hand afterwards.

Run: uv run python scripts/staleness_corpus.py [--json] [--limit N]
"""
from __future__ import annotations

import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CLAUDE_MD = ROOT / "CLAUDE.md"

_PATH = re.compile(r"`([A-Za-z0-9_./-]+/[A-Za-z0-9_./-]+\.(?:py|md|json|yaml|yml|tsv|fna|csv|bed|lock|toml))`")

# A bullet only makes a checkable STATUS claim if it says something about state. A bullet that merely
# points at a file ("see X for the schema") has nothing to be stale about, and auditing it would spend
# GPU to learn nothing. Keeping the filter EXPLICIT and narrow beats a vague relevance heuristic.
_STATUS = re.compile(
    r"\b(?:deferred|pending|blocked|not yet|unrun|never run|still to|TODO|"
    r"shipped|landed|closed|complete|done|ran|runs|validated|measured|"
    r"infeasible|unavailable|exhausted|no free|absent|missing|unreachable)\b", re.I)

ARTIFACT_HEAD_CHARS = 6000
MAX_CLAIM_CHARS = 1200          # a whole CLAUDE.md bullet can be enormous; the head carries the claim


@dataclass(frozen=True)
class Pair:
    pair_id: str
    claim: str
    artifact: str
    bullet_index: int


def split_bullets(text: str) -> list[str]:
    """CLAUDE.md's claims live in top-level `- ` bullets. PURE.

    Splits on a newline followed by "- " at column 0, which is how this file's bullets are written.
    Sub-bullets stay attached to their parent, deliberately: a nested clause is usually a qualification
    OF the parent claim, and separating them would strip the qualifier that makes the claim true.
    """
    parts = re.split(r"\n(?=- )", text)
    return [p.strip() for p in parts if p.strip().startswith("- ")]


# How close a status word must be to a path for the claim to be ABOUT that path. Measured, not guessed:
# CLAUDE.md bullets are long and cite many files, so pairing a bullet with every path it mentions produced
# obvious mismatches (a claim about `viz/browser.py` paired with `config/datasources.yaml` because both
# appeared in one bullet). A local window keeps the pair honest; the price is recall on claims whose
# status word sits far from the path, which is the safe direction -- a missed pair costs nothing, a
# mismatched pair sends the model to judge a claim against the wrong artifact.
_LOCAL_WINDOW = 400


def _local_claim(bullet: str, match: re.Match) -> str | None:
    """The sentence-ish region around a cited path, IF it carries a status word. PURE.

    Returns None when nothing near the path makes a status claim -- that citation is a pointer
    ("see X for the schema"), not an assertion about X, and has nothing to be stale about.
    """
    lo = max(0, match.start() - _LOCAL_WINDOW)
    hi = min(len(bullet), match.end() + _LOCAL_WINDOW)
    region = bullet[lo:hi]
    return region if _STATUS.search(region) else None


def extract_pairs(text: str, *, require_status: bool = True) -> list[Pair]:
    """(claim, artifact) pairs. PURE -- no I/O, so the pairing rule is testable without a GPU.

    The claim shipped is the LOCAL region around the citation, not the whole bullet: a CLAUDE.md bullet
    can cite a dozen files and only assert something about one of them.
    """
    out = []
    for bi, bullet in enumerate(split_bullets(text)):
        if require_status and not _STATUS.search(bullet):
            continue
        seen = set()
        for m in _PATH.finditer(bullet):
            path = m.group(1)
            if path in seen:
                continue
            claim = _local_claim(bullet, m) if require_status else bullet[:MAX_CLAIM_CHARS]
            if claim is None:
                continue
            seen.add(path)
            out.append(Pair(pair_id=f"b{bi:03d}_{len(seen)}_{Path(path).stem[:28]}",
                            claim=claim[:MAX_CLAIM_CHARS], artifact=path, bullet_index=bi))
    return out


def resolvable(pairs: list[Pair]) -> list[Pair]:
    """Drop pairs whose artifact does not exist -- `test_claude_md_citations` already guards that class,
    and an absent artifact gives the model nothing to judge against."""
    return [p for p in pairs if (ROOT / p.artifact).exists()]


def build_corpus(limit: int | None = None) -> list[dict]:
    """The shippable payload: pair + artifact head + structural facts."""
    from kaggle_staleness_auditor import structural_facts
    pairs = resolvable(extract_pairs(CLAUDE_MD.read_text(encoding="utf-8", errors="replace")))
    if limit:
        pairs = pairs[:limit]
    out = []
    for p in pairs:
        text = (ROOT / p.artifact).read_text(encoding="utf-8", errors="replace")
        # A DATA file is summarised by SHAPE, never by a raw head. The claim about a .tsv is about its
        # row count and columns; 6000 characters of accession IDs answer nothing and tokenize so badly
        # they crashed two full-corpus runs at the same item. See tabular_digest for the measurement.
        from kaggle_staleness_auditor import TABULAR_SUFFIXES, tabular_digest
        if p.artifact.endswith(TABULAR_SUFFIXES):
            body = tabular_digest(text[:ARTIFACT_HEAD_CHARS])
        else:
            body = text[:ARTIFACT_HEAD_CHARS]
        out.append({"item_id": p.pair_id, "claim": p.claim, "artifact": p.artifact,
                    "facts": structural_facts(p.artifact, text),
                    "artifact_text": body})
    return out


def main() -> int:
    sys.path.insert(0, str(ROOT / "scripts"))
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])
    corpus = build_corpus(limit)

    if "--json" in sys.argv:
        print(json.dumps(corpus, indent=2))
        return 0

    text = CLAUDE_MD.read_text(encoding="utf-8", errors="replace")
    all_bullets = split_bullets(text)
    with_status = [b for b in all_bullets if _STATUS.search(b)]
    print(f"CLAUDE.md: {len(all_bullets)} bullets, {len(with_status)} make a status claim")
    print(f"-> {len(corpus)} (claim, artifact) pairs with a resolvable artifact\n")
    by_art: dict[str, int] = {}
    for c in corpus:
        by_art[c["artifact"]] = by_art.get(c["artifact"], 0) + 1
    print(f"distinct artifacts: {len(by_art)}")
    for a, n in sorted(by_art.items(), key=lambda kv: -kv[1])[:5]:
        print(f"  {n}x {a}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
