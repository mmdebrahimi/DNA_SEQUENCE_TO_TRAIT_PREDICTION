"""What DATA this project has, derived live -- and, for each dataset, what already consumed it.

WHY THIS EXISTS. `scripts/project_status.py` fixed the same failure for project SCOPE: a written figure
goes stale silently, so every number there is computed at run time. This is its data-side sibling, and it
exists because the rediscovery failure kept happening on the DATA surface instead:

  * 2026-09-09: a session called generalizing a method "offline and cheap", not knowing that the two
    per-isolate checkpoints it needed were already sitting on the D: cache.
  * The same session called a whole track "no agent lever" when the script that advances it had existed
    for months.
  * The same session proposed a measurement that had been completed two days earlier and was sitting in
    a committed artifact.

(Deliberately described rather than named: naming the datasets here would make this module a reference
to them, which is defect 3 below -- the bug this file spent four rounds learning.)

Three rediscoveries in one conversation, on a surface of ~83 cohort dirs + ~42 processed artifacts + ~84
off-repo caches + ~1,100 wiki artifacts. That is far past what any context window holds, so the fix cannot
be "remember better" or "write a good summary" -- a summary is a written number.

THE LOAD-BEARING FIELDS ARE DERIVED, NOT WRITTEN. For every dataset:

  consumed_by  -- which scripts/modules reference it   } computed by a reverse token index over the real
  reported_in  -- which wiki artifacts report on it    } tree, so they cannot drift from what ships

Those two answer the question that was actually missed each time: "has anything already used this?"
A hand-written map would have answered it wrongly within a week.

WHAT IS HAND-CURATED, AND WHY THAT PART IS SMALL. Semantics -- what a dataset IS, how it was used, what
it is blocked on, what it could still support -- is genuinely not derivable. It lives in
`dna_decode/data/inventory_notes.py`, and it is guarded in BOTH directions: a note naming a path that no
longer exists FAILS `--self-check` (the staleness guard), and a substantial dataset with no note is
reported as UNDOCUMENTED rather than silently omitted (the coverage guard). The map reports its own gaps.

READ-ONLY. No network, no Docker, no model load. Exit 0 always except `--self-check`, which is a gate.

Run: uv run python scripts/data_inventory.py            (writes wiki/data_inventory.{md,json})
     uv run python scripts/data_inventory.py --self-check
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from dataclasses import dataclass, field
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Where datasets live. `kind` is the shape of the thing, not a judgement about its value.
SCAN_ROOTS: list[tuple[str, Path, str]] = [
    ("raw", ROOT / "data" / "raw", "cohort"),
    ("processed", ROOT / "data" / "processed", "processed"),
    ("ref", ROOT / "data", "reference"),
    ("cache", Path("D:/dna_decode_cache"), "external_cache"),
]

# Where to look for USE of a dataset. Code first, then the artifacts that report results.
CODE_ROOTS = [ROOT / "scripts", ROOT / "dna_decode", ROOT / "tests"]
ARTIFACT_ROOTS = [ROOT / "wiki", ROOT / "project_state", ROOT / "plans"]

CODE_SUFFIXES = {".py", ".sh", ".bat", ".R"}
ARTIFACT_SUFFIXES = {".md", ".json"}

# THE INVENTORY'S OWN FILES MUST NOT COUNT AS USE. The invariant is: DESCRIBING a dataset must never
# make it look CONSUMED.
#
# It bit FOUR times, each worse than the last, and every one was caught by RE-RUNNING rather than by
# reading the code -- each bad run printed a clean success:
#   1. wiki/data_inventory.{md,json} land in an ARTIFACT_ROOT and name every dataset by construction,
#      so run #2 showed every dataset `reported_in >= 1` -- orphans collapsed 59 -> 8, destroying 86%
#      of the signal while still printing a successful run.
#   2. dna_decode/data/inventory_notes.py lives in a CODE_ROOT, so writing a note about a forgotten
#      dataset marked it `consumed_by` the notes file -- i.e. documenting a dataset as UNUSED was
#      itself enough to report it as used. That is the worse failure: it would tell the next reader
#      the forgotten asset is in active use.
#   3. THIS FILE, because the comment above once NAMED an orphan while explaining the bug, and that
#      mention alone un-orphaned it. Explanatory prose about a thing is lexically indistinguishable
#      from a reference to it -- the same defect class that made a prose-threshold guard read its own
#      explanatory comment as a live claim.
#   4. tests/test_data_inventory.py -- the test asserting that (3) does not recur has to NAME the
#      datasets it checks for, so the guard against naming datasets named two of them.
#
# Four instances killed the hand-listed set: a literal allow/deny list beside the data that defines it
# always under-covers, and each new inventory-system file silently re-opened the hole. The exclusion is
# therefore DERIVED from a filename pattern, so a file added later is covered without being remembered.
# No dataset key is written literally anywhere in this module.
SELF_REFERENCE_PATTERNS = ("data_inventory", "inventory_notes")


def is_self_reference(rel_path: str) -> bool:
    """True for any file belonging to the inventory system itself.

    Pattern-matched rather than enumerated: this bug recurred four times, each time because a NEW
    inventory-owned file was not in the list. Describing a dataset must never make it look consumed.
    """
    return any(pat in rel_path for pat in SELF_REFERENCE_PATTERNS)

# A file bigger than this is a data blob, not a reference to one; reading it would be slow and would
# also produce spurious "mentions" from raw content.
MAX_TEXT_BYTES = 2_000_000

# Directories under data/ that are NOT datasets (they are scanned as their own roots, or are scratch).
REF_SKIP = {"raw", "processed", "cache", "__pycache__"}

_TOKEN = re.compile(r"[A-Za-z0-9_.\-]{4,}")


@dataclass
class Dataset:
    key: str
    root: str
    kind: str
    path: str
    exists: bool
    n_entries: int = 0
    size_bytes: int = 0
    size_truncated: bool = False
    consumed_by: list[str] = field(default_factory=list)
    reported_in: list[str] = field(default_factory=list)
    note: dict | None = None

    @property
    def documented(self) -> bool:
        return self.note is not None

    @property
    def orphan(self) -> bool:
        """Nothing in the tree references it. Either dead weight, or a forgotten asset -- the exact
        shape of the thing that gets rediscovered."""
        return not self.consumed_by and not self.reported_in


def _dir_stats(p: Path, file_cap: int = 4000) -> tuple[int, int, bool]:
    """(n_entries, size_bytes, truncated). Capped: a genome cache can hold 100k files and this is an
    orientation tool, not a disk auditor."""
    if p.is_file():
        return 1, p.stat().st_size, False
    n_entries = 0
    total = 0
    seen = 0
    truncated = False
    try:
        n_entries = sum(1 for _ in os.scandir(p))
    except OSError:
        return 0, 0, False
    for dirpath, _dirnames, filenames in os.walk(p):
        for fn in filenames:
            seen += 1
            if seen > file_cap:
                truncated = True
                break
            try:
                total += os.path.getsize(os.path.join(dirpath, fn))
            except OSError:
                pass
        if truncated:
            break
    return n_entries, total, truncated


def scan_datasets(scan_roots=None) -> list[Dataset]:
    out: list[Dataset] = []
    for root_name, root_path, kind in (scan_roots or SCAN_ROOTS):
        if not root_path.exists():
            # A missing root is INFORMATION (D: unplugged), not an error -- record it and move on.
            out.append(Dataset(key=f"<{root_name}-root-missing>", root=root_name, kind=kind,
                               path=str(root_path), exists=False))
            continue
        for entry in sorted(root_path.iterdir()):
            if entry.name.startswith(".") or entry.name == "__pycache__":
                continue
            if root_name == "ref":
                # data/ holds the other scan roots plus loose scratch files; only directories that are
                # not themselves a scan root count as reference datasets.
                if not entry.is_dir() or entry.name in REF_SKIP:
                    continue
            n, size, trunc = _dir_stats(entry)
            out.append(Dataset(key=entry.name, root=root_name, kind=kind, path=str(entry),
                               exists=True, n_entries=n, size_bytes=size, size_truncated=trunc))
    return out


def _iter_text_files(roots: list[Path], suffixes: set[str]):
    for r in roots:
        if not r.exists():
            continue
        for p in r.rglob("*"):
            if p.suffix not in suffixes or not p.is_file():
                continue
            if "__pycache__" in p.parts:
                continue
            try:
                if p.stat().st_size > MAX_TEXT_BYTES:
                    continue
                yield p, p.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue


def build_reverse_index(keys: set[str], roots: list[Path], suffixes: set[str]) -> dict[str, list[str]]:
    """key -> sorted relative paths of files whose tokens contain that key.

    Tokenised rather than substring-matched: a substring scan makes `cache` match `embeddings_cache`
    and every dataset then looks universally used, which would make the field worthless.
    """
    idx: dict[str, set[str]] = {k: set() for k in keys}
    for path, text in _iter_text_files(roots, suffixes):
        toks = set(_TOKEN.findall(text))
        if not toks:
            continue
        try:
            rel = str(path.relative_to(ROOT)).replace("\\", "/")
        except ValueError:
            rel = str(path).replace("\\", "/")
        if is_self_reference(rel):
            continue  # the inventory's own files are not evidence that anything uses a dataset
        for k in keys & toks:
            idx[k].add(rel)
    return {k: sorted(v) for k, v in idx.items()}


def attach_uses(datasets: list[Dataset]) -> list[Dataset]:
    keys = {d.key for d in datasets if d.exists}
    code = build_reverse_index(keys, CODE_ROOTS, CODE_SUFFIXES)
    arts = build_reverse_index(keys, ARTIFACT_ROOTS, ARTIFACT_SUFFIXES)
    for d in datasets:
        d.consumed_by = code.get(d.key, [])
        d.reported_in = arts.get(d.key, [])
    return datasets


def attach_notes(datasets: list[Dataset]) -> tuple[list[Dataset], list[str]]:
    """Returns (datasets, stale_note_keys). A note whose dataset is absent is the staleness signal."""
    try:
        from dna_decode.data.inventory_notes import NOTES
    except Exception:
        return datasets, []
    present = {d.key for d in datasets if d.exists}
    for d in datasets:
        if d.key in NOTES:
            d.note = dict(NOTES[d.key])
    stale = sorted(k for k in NOTES if k not in present)
    return datasets, stale


def _human(n: int) -> str:
    for unit in ("B", "K", "M", "G"):
        if n < 1024 or unit == "G":
            return f"{n:.0f}{unit}" if unit == "B" else f"{n/1:.1f}{unit}" if False else f"{n}{unit}"
        n //= 1024
    return f"{n}G"


def size_str(d: Dataset) -> str:
    n = d.size_bytes
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            s = f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
            return f">={s}" if d.size_truncated else s
        n /= 1024
    return "?"


def build_report(datasets: list[Dataset], stale: list[str]) -> dict:
    present = [d for d in datasets if d.exists]
    documented = [d for d in present if d.documented]
    orphans = [d for d in present if d.orphan]
    missing_roots = [d for d in datasets if not d.exists]
    return {
        "schema": "data-inventory-v1",
        "date": str(_date.today()),
        "purpose": "Answer 'do we already have this, and has anything already used it?' without a "
                   "context window that holds the whole tree.",
        "derived_fields": ["consumed_by", "reported_in", "n_entries", "size_bytes"],
        "curated_fields": ["note"],
        "counts": {
            "datasets_present": len(present),
            "documented": len(documented),
            "undocumented": len(present) - len(documented),
            "orphans_no_reference_anywhere": len(orphans),
            "roots_missing": len(missing_roots),
        },
        "coverage_note": ("`documented` counts datasets with a curated semantic note. The remainder are "
                          "listed too -- an inventory that hid its own gaps would be the failure it "
                          "exists to prevent."),
        "known_limits": [
            "Matching is by DIRECTORY NAME as a whole token, so a dataset whose name is also an "
            "ordinary English or technical word picks up prose mentions as false 'uses' (several such "
            "names exist here; naming one in this sentence would itself create the false use). "
            "Over-counting a use is the SAFE direction here -- it makes a dataset look more used than "
            "it is, never less -- but treat a low count as the reliable signal and a high count on a "
            "common-word name as soft.",
            "Counts are keyed on the basename, so two datasets sharing a name under different roots "
            "(cache vs processed) share one count and cannot be told apart.",
            "A path built dynamically (f-string) rather than written literally is invisible to the "
            "index. Spot-checks found the hardcoded-literal form dominant, but this is a real floor: "
            "`orphan` means 'no file names this directory', not 'provably unused'.",
        ],
        "stale_notes": stale,
        "roots_missing": [{"root": d.root, "path": d.path} for d in missing_roots],
        "datasets": [
            {"key": d.key, "root": d.root, "kind": d.kind, "path": d.path,
             "n_entries": d.n_entries, "size_bytes": d.size_bytes, "size_truncated": d.size_truncated,
             "n_consumed_by": len(d.consumed_by), "n_reported_in": len(d.reported_in),
             "consumed_by": d.consumed_by[:12], "reported_in": d.reported_in[:12],
             "orphan": d.orphan, "note": d.note}
            for d in sorted(present, key=lambda x: (x.root, x.key))
        ],
    }


def render_md(rep: dict) -> str:
    c = rep["counts"]
    L: list[str] = []
    L.append(f"# Data inventory — what we have, and what already used it ({rep['date']})\n")
    L.append("**Generated. Do not hand-edit** — run `uv run python scripts/data_inventory.py`.\n")
    L.append("Read this before claiming a dataset is missing, or that a piece of work would be "
             "expensive because the data is not on disk. `consumed_by` / `reported_in` are DERIVED "
             "from the real tree, so they answer *has anything already used this?* — the question "
             "whose wrong answer causes rediscovery.\n")
    L.append(f"- datasets present: **{c['datasets_present']}** "
             f"({c['documented']} with a curated note, {c['undocumented']} undocumented)")
    L.append(f"- referenced nowhere in code or artifacts: **{c['orphans_no_reference_anywhere']}**")
    if rep["roots_missing"]:
        L.append(f"- **roots not reachable right now: {c['roots_missing']}** — "
                 + ", ".join(f"`{m['path']}`" for m in rep["roots_missing"])
                 + " (an unplugged drive is not an absent dataset)")
    if rep["stale_notes"]:
        L.append(f"- ⚠ curated notes naming datasets that no longer exist: "
                 + ", ".join(f"`{k}`" for k in rep["stale_notes"]))
    L.append("")

    by_root: dict[str, list[dict]] = {}
    for d in rep["datasets"]:
        by_root.setdefault(d["root"], []).append(d)

    titles = {"raw": "Cohorts / raw inputs (`data/raw`)",
              "processed": "Processed artifacts (`data/processed`)",
              "ref": "Reference databases (`data/`)",
              "cache": "External caches (`D:/dna_decode_cache`) — heavy, gitignored, regenerable"}
    for root in ("raw", "processed", "ref", "cache"):
        rows = by_root.get(root)
        if not rows:
            continue
        L.append(f"## {titles.get(root, root)}\n")
        L.append("| dataset | size | used by | reported in | what it is |")
        L.append("|---|---|---|---|---|")
        for d in rows:
            note = d.get("note") or {}
            what = note.get("what", "—")
            if not note:
                what = "_(undocumented)_"
            used = str(d["n_consumed_by"]) if d["n_consumed_by"] else "**0**"
            rep_n = str(d["n_reported_in"]) if d["n_reported_in"] else "**0**"
            size = f"{d['size_bytes']/1048576:.1f}MB" + (">" if d["size_truncated"] else "")
            L.append(f"| `{d['key']}` | {size} | {used} | {rep_n} | {what} |")
        L.append("")

    L.append("## Datasets nothing references\n")
    orph = [d for d in rep["datasets"] if d["orphan"]]
    if not orph:
        L.append("None — every present dataset is referenced by at least one script or artifact.\n")
    else:
        L.append("Either dead weight or a forgotten asset. This is the bucket rediscovery comes from.\n")
        for d in orph:
            L.append(f"- `{d['key']}` ({d['root']}, {d['size_bytes']/1048576:.1f}MB) — `{d['path']}`")
        L.append("")
    return "\n".join(L)


def self_check(rep: dict) -> tuple[bool, list[str]]:
    """The gate. Fails on a stale curated note or an empty scan."""
    problems = []
    if rep["stale_notes"]:
        problems.append(f"curated notes reference datasets that do not exist: {rep['stale_notes']}")
    if rep["counts"]["datasets_present"] == 0:
        problems.append("scan found zero datasets; the inventory would be vacuous")
    return (not problems), problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-check", action="store_true",
                    help="gate: non-zero exit on a stale curated note or a vacuous scan")
    ap.add_argument("--out-md", type=Path, default=ROOT / "wiki" / "data_inventory.md")
    ap.add_argument("--out-json", type=Path, default=ROOT / "wiki" / "data_inventory.json")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)

    datasets = attach_uses(scan_datasets())
    datasets, stale = attach_notes(datasets)
    rep = build_report(datasets, stale)

    if not a.no_write:
        a.out_json.write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
        a.out_md.write_text(render_md(rep), encoding="utf-8")

    c = rep["counts"]
    print(f"datasets present: {c['datasets_present']}  documented: {c['documented']}  "
          f"undocumented: {c['undocumented']}  orphans: {c['orphans_no_reference_anywhere']}")
    if rep["roots_missing"]:
        print(f"roots not reachable: {[m['path'] for m in rep['roots_missing']]}")
    if not a.no_write:
        print(f"-> {a.out_md}\n-> {a.out_json}")

    if a.self_check:
        ok, problems = self_check(rep)
        for p in problems:
            print(f"SELF-CHECK FAIL: {p}")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
