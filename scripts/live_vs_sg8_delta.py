"""READ-ONLY diagnostic: every path where the live tree and the `sg8` clone differ.

WHY
---
`C:/Users/Laptop Portable/sg8` is a real clone of the project's GitHub remote; the
live working directory here is NOT a git repo at all, so there is no `git status` that
can show the drift. Before anything is staged or committed, the reconcile has to be
done against a LIST rather than from memory -- so this enumerates the delta and writes
a manifest.

It is strictly read-only with respect to BOTH trees. The only thing it writes is its
own artifact.

The classifier core (`classify`) is a pure function over two path->fingerprint maps,
so it is unit-testable on synthetic trees with no reference to either real location.

TWO DELTAS, AND THE NORMALISED ONE IS THE HEADLINE
--------------------------------------------------
Measured on the first real run: a raw byte-hash comparison reported **1,934**
content-differs, and a 60-file sample showed **59 of 60 were LINE-ENDING-ONLY**. The
cause is structural, not incidental -- `.gitattributes` forces LF for the
integrity-hashed surfaces so a pinned sha256 is stable across platforms, so the `sg8`
checkout holds LF while the live working copy holds CRLF. Reconciling against a
1,934-line manifest that is ~98% noise defeats the point of writing a manifest, so
this reports BOTH and leads with the EOL-normalised delta.

Normalisation is applied ONLY to files that look like text (no NUL byte in the first
8 KiB). Collapsing CRLF inside a genuine binary could manufacture false identity,
which is the opposite failure and worse.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

LIVE_ROOT = Path(r"C:/Users/Laptop Portable/pythonprojects/dna_decode")
SG8_ROOT = Path(r"C:/Users/Laptop Portable/sg8")
OUT_PATH = Path("wiki/live_vs_sg8_delta_2026-10-19.json")
SCHEMA = "live-vs-sg8-delta-v1"

# Skipped wholesale. `data` matters most, but NOT for the reason this comment first
# gave. It used to say "a directory junction to the offline D: drive" -- inherited from
# CLAUDE.md and FALSE as of 2026-10-19: `dir /AL` finds no reparse point anywhere in the
# project, so `data/` is a REAL directory on C: holding only a sparse remnant (52 of 59
# files tracked under it in origin/main are absent locally). Skipping it is still right
# -- a partially-populated tree would report absent-on-C: files as live-vs-sg8
# differences, which is a drive-state artifact, not drift. The caches are regenerable
# and would swamp the signal.
SKIP_DIRS = frozenset({
    ".git", ".venv", "venv", "__pycache__", ".pytest_cache", ".ruff_cache",
    ".mypy_cache", "node_modules", "data", ".idea", ".vscode", "htmlcov",
    ".repo-index", "egg-info",
})

LIVE_ONLY = "live-only"
SG8_ONLY = "sg8-only"
CONTENT_DIFFERS = "content-differs"
IDENTICAL = "identical"


@dataclass
class Delta:
    live_only: list[str] = field(default_factory=list)
    sg8_only: list[str] = field(default_factory=list)
    content_differs: list[str] = field(default_factory=list)
    identical: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            LIVE_ONLY: sorted(self.live_only),
            SG8_ONLY: sorted(self.sg8_only),
            CONTENT_DIFFERS: sorted(self.content_differs),
            IDENTICAL: sorted(self.identical),
        }

    def counts(self) -> dict:
        return {
            LIVE_ONLY: len(self.live_only),
            SG8_ONLY: len(self.sg8_only),
            CONTENT_DIFFERS: len(self.content_differs),
            IDENTICAL: len(self.identical),
        }


def classify(live_index: dict[str, str], sg8_index: dict[str, str]) -> Delta:
    """Pure classifier over two `{relpath: fingerprint}` maps."""
    delta = Delta()
    for path in set(live_index) | set(sg8_index):
        in_live, in_sg8 = path in live_index, path in sg8_index
        if in_live and not in_sg8:
            delta.live_only.append(path)
        elif in_sg8 and not in_live:
            delta.sg8_only.append(path)
        elif live_index[path] == sg8_index[path]:
            delta.identical.append(path)
        else:
            delta.content_differs.append(path)
    return delta


def looks_like_text(blob: bytes) -> bool:
    """No NUL byte in the first 8 KiB. Deliberately conservative."""
    return b"\x00" not in blob[:8192]


def _fingerprints(path: Path) -> tuple[str, str]:
    """`(raw_sha256, eol_normalised_sha256)`.

    For a binary the two are identical -- normalisation is NOT applied, because
    collapsing CRLF inside real binary content could make two different files hash
    the same, which is a far worse error than over-reporting a text diff.
    """
    blob = path.read_bytes()
    raw = hashlib.sha256(blob).hexdigest()
    if not looks_like_text(blob):
        return raw, raw
    norm = hashlib.sha256(blob.replace(b"\r\n", b"\n")).hexdigest()
    return raw, norm


def index_tree(
    root: Path, *, skip_dirs: frozenset[str] = SKIP_DIRS
) -> tuple[dict[str, str], dict[str, str]]:
    """`({relpath: raw_sha}, {relpath: normalised_sha})` for every readable file."""
    root = Path(root)
    raw_index: dict[str, str] = {}
    norm_index: dict[str, str] = {}
    if not root.is_dir():
        return raw_index, norm_index
    stack = [root]
    while stack:
        current = stack.pop()
        try:
            entries = list(current.iterdir())
        except OSError:
            continue
        for entry in entries:
            try:
                if entry.is_dir():
                    if entry.name in skip_dirs or entry.name.endswith(".egg-info"):
                        continue
                    stack.append(entry)
                elif entry.is_file():
                    rel = entry.relative_to(root).as_posix()
                    raw_index[rel], norm_index[rel] = _fingerprints(entry)
            except OSError:
                continue
    return raw_index, norm_index


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--live", default=str(LIVE_ROOT))
    ap.add_argument("--sg8", default=str(SG8_ROOT))
    ap.add_argument("--out", default=str(OUT_PATH))
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    live_root, sg8_root = Path(args.live), Path(args.sg8)
    if not sg8_root.is_dir():
        print(f"SG8 ABSENT: {sg8_root} is not a directory -- nothing to compare.")
        return 2

    live_raw, live_norm = index_tree(live_root)
    sg8_raw, sg8_norm = index_tree(sg8_root)
    delta = classify(live_norm, sg8_norm)        # HEADLINE: EOL-normalised
    raw_delta = classify(live_raw, sg8_raw)      # secondary: raw bytes

    eol_only = sorted(set(raw_delta.content_differs) - set(delta.content_differs))

    out = live_root / args.out if not Path(args.out).is_absolute() else Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "_schema": SCHEMA,
        "generated": "2026-10-19",
        "live_root": str(live_root).replace("\\", "/"),
        "sg8_root": str(sg8_root).replace("\\", "/"),
        "skipped_dirs": sorted(SKIP_DIRS),
        "note": (
            "READ-ONLY. `data/` is skipped because it holds only a sparse remnant: 52 of "
            "the 59 files tracked under it in origin/main are absent locally. Walking it "
            "would report a drive-state artifact as live-vs-sg8 drift. (An earlier "
            "version of this note said `data/` is a directory junction to D:. That was "
            "inherited from CLAUDE.md and is FALSE -- there is no reparse point anywhere "
            "in the project as of 2026-10-19.)"
        ),
        "which_delta_to_reconcile_against": (
            "`delta` (EOL-normalised). `raw_delta` is kept only to show how much of the "
            "byte-level difference is line endings: .gitattributes forces LF in the "
            "repo blob so the integrity-hashed surfaces keep a stable sha256 across "
            "platforms, while the live working copy is CRLF. Reconciling against "
            "`raw_delta` would mean reading a list that is overwhelmingly noise."
        ),
        "counts": delta.counts(),
        "raw_counts": raw_delta.counts(),
        "n_eol_only_differences": len(eol_only),
        "n_live_files": len(live_raw),
        "n_sg8_files": len(sg8_raw),
        "delta": delta.as_dict(),
        "raw_delta": raw_delta.as_dict(),
        "eol_only_differences": eol_only,
    }
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    if not args.quiet:
        counts, raw_counts = delta.counts(), raw_delta.counts()
        print(f"WROTE {out}")
        print(f"  live files      : {len(live_raw)}")
        print(f"  sg8 files       : {len(sg8_raw)}")
        print("  -- EOL-NORMALISED (reconcile against this) --")
        for key in (LIVE_ONLY, SG8_ONLY, CONTENT_DIFFERS, IDENTICAL):
            print(f"  {key:16s}: {counts[key]}")
        print(f"  -- raw bytes: content-differs {raw_counts[CONTENT_DIFFERS]} "
              f"({len(eol_only)} of them line-endings only) --")
        for key in (LIVE_ONLY, SG8_ONLY, CONTENT_DIFFERS):
            sample = sorted(getattr(delta, key.replace("-", "_")))[:10]
            if sample:
                print(f"  sample {key}:")
                for s in sample:
                    print(f"      {s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
