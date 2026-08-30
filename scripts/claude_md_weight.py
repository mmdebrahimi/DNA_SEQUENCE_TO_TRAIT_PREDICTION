"""What does CLAUDE.md cost every session, and which bullets could be a pointer instead of prose?

CLAUDE.md is the only project text auto-loaded into EVERY session, so its whole word count is a fixed
per-session tax on context -- and its structure decides whether a session gets the project right. That is
not theoretical here: a stale opening paragraph caused a ~4x scope understatement that survived a full day
(see `## READ THIS FIRST`).

THE DISCIPLINE THIS MEASURES, AND WHY IT IS NOT "SHORTER IS BETTER". The body's value is hard-won gotchas
that nothing else records; deleting those to save tokens would be strictly destructive. The question is
narrower and answerable: **does this bullet's detail exist anywhere else?** A bullet citing a `wiki/` memo
that really is on disk can hold the rule and the headline and point at the memo for the derivation. A
bullet citing nothing IS the only store, and must stay whole regardless of length.

So this reports two things per bullet -- its size, and whether it has an external store -- and REFUSES to
recommend compressing anything it cannot prove is stored elsewhere. It never edits; it is a report.

Run: uv run python scripts/claude_md_weight.py [--top 15]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CLAUDE_MD = ROOT / "CLAUDE.md"

# A bullet this long is worth looking at; below it the tax is not worth the risk of touching it.
LONG_BULLET_WORDS = 250
CHARS_PER_TOKEN = 4  # rough, and stated as rough -- the ratio is for orientation, not billing

CITATION = re.compile(r"`?((?:wiki|plans|features|executed_plans)/[A-Za-z0-9_./*{}<>\-,]+)`?")


def bullets(text: str) -> list[str]:
    return [b for b in re.split(r"\n(?=- \*\*)", text) if b.startswith("- ")]


def expand_braces(path: str) -> list[str]:
    """Expand EVERY brace group, not just the first.

    The first version handled one, so a real citation like
    `wiki/forward_inverse_roundtrip_2026-07-1{6,7}.{md,json}` came out as a literal path containing
    `.{md,json}` and was reported BROKEN. Three of the six "broken citations" in the first run were this
    bug -- a checker whose false positives look exactly like the defect it hunts is worse than none.
    """
    if "{" not in path:
        return [path]
    pre, rest = path.split("{", 1)
    if "}" not in rest:
        return [path]
    opts, post = rest.split("}", 1)
    return [p for o in opts.split(",") for p in expand_braces(f"{pre}{o}{post}")]


def cited_groups(bullet: str) -> list[tuple[str, list[str]]]:
    """One entry per CITATION AS WRITTEN, with its expansions.

    Grouping is the point. `wiki/forward_inverse_{a,b}_2026-07-1{6,7}.{md,json}` means "these files, each
    on one of those dates" -- it is satisfied when the referenced work is findable, not when all twelve
    expansions exist. Flattening first made a satisfied citation report six broken siblings, which is the
    same false-positive class the fix above removed.
    """
    return [(raw.rstrip(".,;)"), expand_braces(raw.rstrip(".,;)")))
            for raw in CITATION.findall(bullet)]


def resolves(path: str) -> bool:
    """Does this citation point at something real?

    A `<placeholder>` (e.g. `wiki/genome_map_spike_verdict_<date>.json`) is a TEMPLATE, not a path -- the
    honest check is whether any file matches its literal prefix. Treating it as a literal filename was the
    other half of the first run's false positives.
    """
    p = ROOT / path
    if p.exists():
        return True
    if "<" in path:                       # template citation -> satisfied by any file sharing its prefix
        prefix = Path(path).name.split("<", 1)[0]
        parent = (ROOT / path).parent
        return bool(prefix) and parent.exists() and any(parent.glob(prefix + "*"))
    if "*" in path:                       # a glob citation counts if anything matches
        parent = (ROOT / path).parent
        return parent.exists() and any(parent.glob(Path(path).name))
    return False


def analyse(bullet: str) -> dict:
    groups = cited_groups(bullet)
    stored = [raw for raw, exp in groups if any(resolves(e) for e in exp)]
    missing = [raw for raw, exp in groups if not any(resolves(e) for e in exp)]
    cites = [raw for raw, _ in groups]
    head = bullet.split("**")[1][:70] if "**" in bullet else bullet[2:72]
    return {"words": len(bullet.split()), "head": head.strip(),
            "n_cited": len(cites), "n_stored": len(stored), "missing": missing,
            # Only a bullet with a REAL external store is a compression candidate. No store -> keep whole.
            "compressible": bool(stored) and len(bullet.split()) >= LONG_BULLET_WORDS}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    text = CLAUDE_MD.read_text(encoding="utf-8")
    bl = bullets(text)
    rows = sorted((analyse(b) for b in bl), key=lambda r: -r["words"])

    total_words = len(text.split())
    approx_tokens = len(text) // CHARS_PER_TOKEN
    print(f"CLAUDE.md: {total_words:,} words across {len(bl)} bullets "
          f"(~{approx_tokens:,} tokens, loaded EVERY session)")
    print(f"\n  {'words':>6}  {'cited':>5} {'stored':>6}  cand  bullet")
    for r in rows[:args.top]:
        print(f"  {r['words']:>6}  {r['n_cited']:>5} {r['n_stored']:>6}  "
              f"{'YES ' if r['compressible'] else '  - '}  {r['head']}")

    cand = [r for r in rows if r["compressible"]]
    keep = [r for r in rows if r["words"] >= LONG_BULLET_WORDS and not r["compressible"]]
    print(f"\n  compression CANDIDATES (long AND provably stored elsewhere): {len(cand)}"
          f"  = {sum(r['words'] for r in cand):,} words")
    print(f"  long but NO external store -- MUST STAY WHOLE:               {len(keep)}"
          f"  = {sum(r['words'] for r in keep):,} words")

    broken = [(r["head"], m) for r in rows for m in r["missing"]]
    if broken:
        print(f"\n  BROKEN CITATIONS ({len(broken)}) -- a bullet promising a memo that is not on disk is"
              f" worse than no pointer:")
        for head, m in broken[:12]:
            print(f"    {m}   <- {head[:52]}")
    else:
        print("\n  every cited wiki/plans path resolves.")
    print("\n  This is a REPORT. It never edits, and it refuses to call a bullet compressible")
    print("  unless a cited file really exists -- an uncited bullet is the only copy of what it knows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
