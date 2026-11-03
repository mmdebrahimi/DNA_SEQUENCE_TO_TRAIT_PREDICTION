"""Mirror the load-bearing part of the PromoterAtlas clone off the volatile temp directory.

Three declared arms depend on a clone that lives at
`C:/Users/.../AppData/Local/Temp/pa_probe/PromoterAtlas`. A reboot or a temp sweep deletes it, and
two unattended windows -- one overnight, one ten days -- cannot depend on that. This copies the part
the runs actually need and VERIFIES the copy by hash, because a destination with the right filenames
is not evidence of the right contents.

**What is load-bearing, measured rather than guessed.** `load_atlas` needs exactly two things: the
base checkpoint and `src/promoter_atlas/models/dna_transformer.py`, which it loads standalone through
`importlib.util.spec_from_file_location`. That module imports ONLY torch -- no package siblings -- so
it genuinely stands alone. `atlas_heads.load_head` needs no source at all; it reads state dicts. The
whole `src/promoter_atlas/` tree is copied anyway because it is small and a narrower copy would break
silently if anything else is ever reached for.

**THE TARGET IS C:, WHICH DELIBERATELY INVERTS THIS PROJECT'S STANDING D: RULE, AND THE ARGUMENT IS
MEASURED RATHER THAN TIDY.** `run_tile_conv_detached.bat` records six deaths of long unattended jobs
on this host: four were `STATUS_IN_PAGE_ERROR (0xc0000006)` with torch memory-mapping 4.33 GB of
weights off the USB volume, and two produced no event on any of 118 logs and remain unexplained. The
single long job that ever reached `rc=0` here, at roughly 12 hours, kept its data on C:.

The D: rule exists to keep a chronically-full C: clear of HEAVY artifacts. This mirror is ~88 MB of
weights plus a small source tree -- about 0.3% of C:'s free space -- so the rule's purpose is not
violated, and the bulk caches stay on D:. The two facts are written into the artifact so the
exception is explicit rather than silent, and so a future reader can disagree with it knowingly.
This is favourable-but-not-proven, in the same register the existing launcher uses; it is not a
guarantee, which is why the runs that depend on it still checkpoint.

**No code change is needed to USE the mirror.** `atlas_embed.DEFAULT_PA_REPO` already reads
`os.environ.get("PA_REPO", <temp path>)`, so the launchers export the value this script records.

Read-only with respect to the source clone. Idempotent. Exits non-zero on any verification failure.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import time
from pathlib import Path

from dna_decode.glm.atlas_embed import DEFAULT_PA_REPO, AtlasEmbedError, audit_checkpoint
from scripts.expression_heads_probe import EXPECTED_BYTES, EXPECTED_SHA256
from scripts.positional_arm_probe import sha256_of

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/durable_clone_2026-11-03.json"

#: Default mirror location. On C: deliberately -- see the module docstring.
DEFAULT_DEST = REPO / "data/pa_clone_durable"

#: Subtrees copied, in the order they are verified. `trained_weights` first because it is the part
#: whose integrity the arms' numbers depend on.
SUBTREES = ("trained_weights", "src/promoter_atlas")

#: The one source file `load_atlas` loads standalone. Checked for presence explicitly, because the
#: whole mirror is useless without it and a missing file here would surface as a confusing
#: AtlasEmbedError hours into an unattended run.
REQUIRED_SOURCE = "src/promoter_atlas/models/dna_transformer.py"

BASE_CHECKPOINT = "trained_weights/base_model/promoteratlas-base.pt"


class DurableCloneError(RuntimeError):
    """The mirror cannot be made, or cannot be shown to match its source."""


# ---------------------------------------------------------------------------------------------------
# pure helpers
# ---------------------------------------------------------------------------------------------------
def copy_manifest(src_root: Path | str, subtrees=SUBTREES) -> list[dict]:
    """Every file the mirror should contain, as `{relpath, bytes}`.

    **Ordering: DETERMINISTIC, but not globally lexicographic.** Files are sorted within each
    subtree and the subtrees come in `SUBTREES` order, which puts `trained_weights` first
    deliberately -- it is the part whose integrity the arms' numbers depend on. Stated precisely
    because the looser word "sorted" is what made me assert global sortedness in a test and get it
    wrong.

    Pure apart from reading the filesystem: no copying, no hashing. Refuses an absent subtree rather
    than silently producing a short manifest, because a short manifest copies a partial clone and
    reports success.
    """
    root = Path(src_root)
    rows: list[dict] = []
    for sub in subtrees:
        d = root / sub
        if not d.is_dir():
            raise DurableCloneError("source subtree absent: %s" % d)
        for p in sorted(d.rglob("*")):
            if p.is_file() and p.suffix != ".pyc" and "__pycache__" not in p.parts:
                rows.append({"relpath": p.relative_to(root).as_posix(),
                             "bytes": p.stat().st_size})
    if not rows:
        raise DurableCloneError("manifest is empty; nothing would be mirrored")
    return rows


def verification_summary(rows: list[dict]) -> dict:
    """Roll a per-file verification list into a verdict. Pure.

    `ok` requires EVERY row to match. A summary that reported a percentage would let a partial
    mirror read as a qualified success, which for a checkpoint is the same as a failure.
    """
    total = len(rows)
    matched = [r for r in rows if r.get("digest_matches") is True]
    mismatched = [r["relpath"] for r in rows if r.get("digest_matches") is False]
    unchecked = [r["relpath"] for r in rows if r.get("digest_matches") is None]
    return {"n_files": total, "n_matched": len(matched),
            "mismatched": sorted(mismatched), "unchecked": sorted(unchecked),
            "ok": bool(total) and len(matched) == total,
            "total_bytes": sum(r.get("bytes", 0) for r in rows),
            "rule": ("ok requires EVERY file to match; a percentage would let a partial mirror read "
                     "as a qualified success, and for a checkpoint that is just a failure")}


def needs_copy(src: Path, dest: Path) -> bool:
    """True when the destination is absent or differs in size. Cheap pre-filter before hashing."""
    if not dest.is_file():
        return True
    return src.stat().st_size != dest.stat().st_size


# ---------------------------------------------------------------------------------------------------
# the mirror
# ---------------------------------------------------------------------------------------------------
def mirror(*, src_root: Path | str | None = None, dest_root: Path | str | None = None,
           verbose: bool = True) -> dict:
    src = Path(src_root if src_root is not None else DEFAULT_PA_REPO)
    dest = Path(dest_root if dest_root is not None else DEFAULT_DEST)
    if not src.is_dir():
        raise DurableCloneError(
            "source clone absent: %s. The volatile temp directory has gone, so the mirror cannot be "
            "made and the arms that depend on it are blocked until the clone is re-fetched." % src)
    if (src / REQUIRED_SOURCE).is_file() is False:
        raise DurableCloneError("source is missing %s, without which load_atlas cannot work"
                                % REQUIRED_SOURCE)

    t0 = time.time()
    rows = copy_manifest(src)
    n_copied, n_skipped = 0, 0

    for r in rows:
        s, d = src / r["relpath"], dest / r["relpath"]
        if needs_copy(s, d):
            d.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(s, d)
            n_copied += 1
            r["action"] = "copied"
        else:
            n_skipped += 1
            r["action"] = "already_present"

    # verify AFTER copying, by hashing both sides. A size match is not a content match.
    for r in rows:
        s, d = src / r["relpath"], dest / r["relpath"]
        if not d.is_file():
            r["digest_matches"] = False
            r["why"] = "destination absent after the copy pass"
            continue
        sd, dd = sha256_of(s), sha256_of(d)
        r["source_sha256"], r["dest_sha256"] = sd, dd
        r["digest_matches"] = sd == dd
        if verbose and not r["digest_matches"]:
            print("  MISMATCH %s" % r["relpath"], flush=True)

    summary = verification_summary(rows)

    # the base checkpoint additionally matches the PINNED identity every arm inherits
    base_dest = dest / BASE_CHECKPOINT
    pinned: dict = {"path": str(base_dest), "expected_sha256": EXPECTED_SHA256,
                    "expected_bytes": EXPECTED_BYTES}
    if base_dest.is_file():
        pinned.update({"sha256": sha256_of(base_dest), "bytes": base_dest.stat().st_size})
        pinned["matches_pin"] = (pinned["sha256"] == EXPECTED_SHA256
                                 and pinned["bytes"] == EXPECTED_BYTES)
    else:
        pinned.update({"sha256": None, "bytes": None, "matches_pin": False})

    # every checkpoint in the MIRROR re-passes the static pickle audit, so the mirror is not merely
    # byte-equal to the source but independently loadable under the same guarantee
    audits: dict = {}
    for p in sorted((dest / "trained_weights").rglob("*.pt")) if (dest / "trained_weights").is_dir() else []:
        key = p.relative_to(dest).as_posix()
        try:
            audits[key] = {"passed": True, "globals": audit_checkpoint(p)}
        except AtlasEmbedError as e:
            audits[key] = {"passed": False, "error": str(e), "kind": "off_allowlist_global"}
        except Exception as e:                                 # noqa: BLE001
            # A truncated or non-archive .pt raises `zipfile.BadZipFile`, NOT AtlasEmbedError, and a
            # narrow catch here let that crash the mirror with a confusing traceback instead of
            # recording a failed audit. This project has already been bitten by files of the right
            # NAME and the wrong CONTENTS (an assembly cohort where every cached genome was a
            # truncated gzip or an HTTP error page), so an unreadable checkpoint is a FAILED AUDIT,
            # never an exception that escapes.
            audits[key] = {"passed": False, "error": "%s: %s" % (type(e).__name__, e),
                           "kind": "unreadable_or_not_an_archive"}

    all_audited = bool(audits) and all(v["passed"] for v in audits.values())
    required_present = (dest / REQUIRED_SOURCE).is_file()

    rec = {
        "schema": "durable-clone-v1",
        "date": "2026-11-03",
        "source_root": str(src),
        "dest_root": str(dest),
        "PA_REPO": str(dest),
        "how_to_use": ("export PA_REPO=<dest_root>; atlas_embed.DEFAULT_PA_REPO already reads it, so "
                       "no code change is needed"),
        "host": {"platform": platform.platform(), "cpu_count": os.cpu_count()},
        "subtrees": list(SUBTREES),
        "required_source_present": required_present,
        "required_source": REQUIRED_SOURCE,
        "n_copied": n_copied, "n_already_present": n_skipped,
        "idempotent": n_copied == 0,
        "verification": summary,
        "pinned_base_checkpoint": pinned,
        "checkpoint_audits": audits,
        "all_checkpoints_audited": all_audited,
        "target_is_C_drive_deliberately": {
            "inverts": "the standing project rule that heavy artifacts live on D:, never on C:",
            "why": ("run_tile_conv_detached.bat records SIX deaths of long unattended jobs on this "
                    "host: four STATUS_IN_PAGE_ERROR (0xc0000006) with torch memory-mapping 4.33 GB "
                    "of weights off the USB volume, and two with no event on any of 118 logs. The "
                    "single long job that ever reached rc=0 kept its data on C:."),
            "why_the_rule_is_not_violated": ("the D: rule exists to keep a chronically-full C: clear "
                                             "of HEAVY artifacts; this mirror is %.1f MB, and the "
                                             "bulk caches stay on D:"
                                             % (summary["total_bytes"] / 1e6)),
            "confidence": ("favourable-but-NOT-proven, in the register the existing launcher uses. "
                           "Not a guarantee -- which is why the runs depending on this mirror still "
                           "checkpoint per unit."),
        },
        "what_is_load_bearing": {
            "load_atlas_needs": [BASE_CHECKPOINT, REQUIRED_SOURCE],
            "measured": ("dna_transformer.py imports ONLY torch -- no package siblings -- so it "
                         "genuinely stands alone under importlib spec loading. atlas_heads.load_head "
                         "needs no source at all; it reads state dicts."),
            "why_copy_more": ("the whole src/promoter_atlas tree is copied because it is small and a "
                              "narrower copy would break silently if anything else is ever reached "
                              "for"),
        },
        "wall_seconds": round(time.time() - t0, 1),
    }
    rec["ok"] = bool(summary["ok"] and pinned["matches_pin"] and all_audited and required_present)
    rec["files"] = rows
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", default=None, help="source clone root (default: PA_REPO or the temp dir)")
    ap.add_argument("--dest", default=None, help="mirror root (default: data/pa_clone_durable on C:)")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)

    try:
        rec = mirror(src_root=args.src, dest_root=args.dest)
    except (DurableCloneError, AtlasEmbedError) as e:
        print("REFUSED: %s" % e, flush=True)
        return 2
    except Exception as e:                                     # noqa: BLE001
        print("MIRROR FAILED (%s): %s" % (type(e).__name__, e), flush=True)
        return 1

    if not rec["ok"]:
        v = rec["verification"]
        print("VERIFICATION FAILED: %d/%d files matched; mismatched=%r; pin=%s; audits=%s"
              % (v["n_matched"], v["n_files"], v["mismatched"],
                 rec["pinned_base_checkpoint"]["matches_pin"], rec["all_checkpoints_audited"]),
              flush=True)
        return 3

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    v = rec["verification"]
    print("MIRROR   %s -> %s" % (rec["source_root"], rec["dest_root"]))
    print("FILES    %d verified (%d copied, %d already present), %.1f MB"
          % (v["n_matched"], rec["n_copied"], rec["n_already_present"], v["total_bytes"] / 1e6))
    print("PIN      base sha256 %s... %d bytes -> matches"
          % (rec["pinned_base_checkpoint"]["sha256"][:16], rec["pinned_base_checkpoint"]["bytes"]))
    print("AUDIT    %d checkpoints re-pass the static pickle audit" % len(rec["checkpoint_audits"]))
    print("PA_REPO  %s" % rec["PA_REPO"])
    print("-> %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
