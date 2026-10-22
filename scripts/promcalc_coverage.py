"""Print how many tile_ids stage B has actually completed, across the main file and every shard.

Separate from the watchdog that calls it so the watchdog needs no embedded Python. (The first version
inlined this as a PowerShell here-string and the closing `'@` sat indented inside a function, which is
a parse error -- the watchdog silently never started and wrote no log.)

Run it only while stage B is IDLE. It decompresses the whole output, and doing that during a live run
is the recorded trap that turned a ~1 h job into a 5.0 h projection: observing the sweep slowed the
sweep, and the slowed rate was then quoted as the sweep's rate.

Reuses `promcalc_decompose.completed_tiles` rather than re-counting, so the number means exactly what
`--resume` will skip -- including its rule that a truncated shard's LAST tile is dropped, because a
tile killed mid-write has only some of its ~73 rows and "has at least one row" would read as complete.

    python scripts/promcalc_coverage.py <in_dir>

Prints a bare integer and exits 0, or prints nothing and exits 2. Stdlib only, so it runs under the
isolated stage-B interpreter with no project dependencies.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def completed(in_dir: Path) -> int:
    spec = importlib.util.spec_from_file_location("_sb", HERE / "promcalc_decompose.py")
    if not (spec and spec.loader):
        raise RuntimeError("cannot load promcalc_decompose.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return len(mod.completed_tiles(in_dir))


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    try:
        print(completed(Path(args[0])))
    except Exception as exc:  # noqa: BLE001 - the caller must not read a failure as a count
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
