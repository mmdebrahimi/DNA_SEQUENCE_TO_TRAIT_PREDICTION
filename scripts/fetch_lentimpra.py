#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fetch the Agarwal 2025 human lentiMPRA deposit, verifying every byte against the publisher's MD5.

    uv run python scripts/fetch_lentimpra.py --list
    uv run python scripts/fetch_lentimpra.py --select tables
    uv run python scripts/fetch_lentimpra.py --select all --yes

WHY VERIFICATION IS THE POINT OF THIS SCRIPT, not a nicety.

A first hand-rolled `curl` of two of these files wrote **a 20,899-byte Zenodo
"application-unavailable" page and a 150-byte nginx "502 Bad Gateway" page, under genuine `.xlsx`
filenames**. Zenodo was returning 502/503/504 intermittently. Fed to a parser those give an error
that is indistinguishable from "this table has no usable rows" -- corrupt input wearing the costume
of a biological result. This repo has already paid for that lesson once: the cached GPS assembly
cohort was 100% unusable (HTTP 403 pages and truncated gzips under genome filenames) and it took a
separate investigation to notice.

**A HEAD request or an opening-bytes sniff is the WRONG guard**, and that is measured here too: a
truncated gzip decompresses its opening blocks perfectly and only raises at the END, so reading the
first 2 KB passes exactly the files that later fail. The only sufficient check is a hash over the
WHOLE file, and Zenodo publishes an MD5 per file, so this script uses it. Anything that fails is
DELETED, never left on disk to be found later and trusted.

SELECTION MATTERS -- the deposit is 7.68 GB and most of it is not data. Measured by reading both
zip central directories over HTTP range requests (no download): `final_dump.zip` is 303 entries of
which **270 are model checkpoints totalling 3.83 GB**, leaving 10 MB of TSVs that are the authors'
own model PREDICTIONS; `final_joint_dump.zip` is the same shape at 4.33 GB. So `--select tables` is
the default: the supplementary tables carry the library design, and the 7.5 GB of LegNet weights are
opt-in rather than the thing you get by asking for the data.

Lands under `D:/dna_decode_cache/mpra_human/` -- the D: rule, and namespace-separate from the
E. coli `mpra/` directory so a human table can never be mistaken for a GSE144621 one.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

RECORD = "10558183"
API = "https://zenodo.org/api/records/%s" % RECORD
DEFAULT_DEST = Path("D:/dna_decode_cache/mpra_human")

#: Expected licence. Checked, not assumed -- a licence change is a reason to stop and ask.
EXPECTED_LICENSE = "cc-by-4.0"

#: Zenodo returns 502/503/504 under load. Transient, so retry with backoff rather than fail.
RETRIES = 6
BACKOFF_BASE = 4

SELECTIONS = {
    "tables": lambda k: k.lower().endswith(".xlsx"),
    "code": lambda k: k in ("sequence_cnn_models-master.zip", "human_legnet-main.zip"),
    "dumps": lambda k: k in ("final_dump.zip", "final_joint_dump.zip"),
    "all": lambda k: True,
}

#: Bodies that are error pages rather than payload. Cheap pre-check so the FAILURE is legible;
#: it is NOT the integrity guard -- the MD5 is.
HTML_SIGNATURES = (b"<html", b"<!DOCTYPE html", b"<?xml version")


class FetchError(RuntimeError):
    pass


def _curl_json(url: str, out: Path) -> dict:
    last = None
    for attempt in range(1, RETRIES + 1):
        r = subprocess.run(
            ["curl", "-sL", "--max-time", "90", url, "-o", str(out),
             "-w", "%{http_code}"],
            capture_output=True, text=True,
        )
        code = (r.stdout or "").strip()
        if code == "200" and out.is_file() and out.stat().st_size > 0:
            try:
                return json.loads(out.read_text(encoding="utf-8"))
            except json.JSONDecodeError as exc:
                last = "manifest was not JSON (%s)" % exc
        else:
            last = "http=%s" % code
        wait = BACKOFF_BASE * attempt
        print("  manifest attempt %d/%d failed (%s); retrying in %ds"
              % (attempt, RETRIES, last, wait))
        time.sleep(wait)
    raise FetchError("could not fetch the manifest: %s" % last)


def manifest(cache: Path) -> dict:
    """The publisher's own file list: name, size, md5. The source of truth for verification."""
    cache.mkdir(parents=True, exist_ok=True)
    raw = _curl_json(API, cache / "_zenodo_record.json")
    md = raw.get("metadata", {})
    lic = md.get("license")
    lic_id = lic.get("id") if isinstance(lic, dict) else lic
    files = []
    for f in raw.get("files", []):
        chk = f.get("checksum", "")
        algo, _, digest = chk.partition(":")
        files.append({
            "key": f["key"],
            "size": int(f["size"]),
            "algo": algo or "md5",
            "digest": digest,
            "url": f["links"]["self"],
        })
    return {
        "record": RECORD,
        "title": md.get("title"),
        "publication_date": md.get("publication_date"),
        "license": lic_id,
        "doi": raw.get("doi"),
        "files": files,
    }


def file_digest(path: Path, algo: str = "md5", chunk: int = 1 << 22) -> str:
    h = hashlib.new(algo)
    with path.open("rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def looks_like_error_page(path: Path) -> bool:
    try:
        head = path.open("rb").read(256).lstrip()
    except OSError:
        return False
    return any(head.startswith(s) for s in HTML_SIGNATURES)


def verify(path: Path, entry: dict) -> tuple[bool, str]:
    """Size AND whole-file digest. Both, because either alone has a real failure mode.

    A truncated transfer can coincidentally match neither; an HTML error page matches neither; but a
    size check alone would pass a same-length corruption and a digest alone would waste minutes
    hashing an obviously-short file. Size is the cheap gate, the digest is the verdict.
    """
    if not path.is_file():
        return False, "absent"
    got = path.stat().st_size
    if got != entry["size"]:
        extra = " (an HTML error page, not payload)" if looks_like_error_page(path) else ""
        return False, "size %d != expected %d%s" % (got, entry["size"], extra)
    if not entry["digest"]:
        return False, "publisher supplied no digest; refusing to call this verified"
    actual = file_digest(path, entry["algo"])
    if actual != entry["digest"]:
        return False, "%s %s != expected %s" % (entry["algo"], actual[:16], entry["digest"][:16])
    return True, "ok"


def fetch_one(entry: dict, dest: Path) -> tuple[bool, str]:
    path = dest / entry["key"]
    ok, why = verify(path, entry)
    if ok:
        return True, "already present and verified"

    for attempt in range(1, RETRIES + 1):
        # -C - resumes a partial transfer; a resumed body is still verified whole-file below.
        r = subprocess.run(
            ["curl", "-sL", "-C", "-", "--max-time", "3600", entry["url"],
             "-o", str(path), "-w", "%{http_code}"],
            capture_output=True, text=True,
        )
        code = (r.stdout or "").strip()
        ok, why = verify(path, entry)
        if ok:
            return True, "verified (%s, attempt %d)" % (entry["algo"], attempt)
        # A body that is an error page must not be left to be resumed from.
        if looks_like_error_page(path) or code not in ("200", "206"):
            path.unlink(missing_ok=True)
        if attempt < RETRIES:
            wait = BACKOFF_BASE * attempt
            print("    attempt %d/%d: http=%s %s; retrying in %ds"
                  % (attempt, RETRIES, code, why, wait))
            time.sleep(wait)

    path.unlink(missing_ok=True)
    return False, "FAILED after %d attempts: %s (nothing kept)" % (RETRIES, why)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dest", default=str(DEFAULT_DEST))
    ap.add_argument("--select", default="tables", choices=sorted(SELECTIONS))
    ap.add_argument("--list", action="store_true", help="show the manifest and exit")
    ap.add_argument("--yes", action="store_true",
                    help="required for a selection over 1 GB, so 7.5 GB is never an accident")
    a = ap.parse_args(argv)

    dest = Path(a.dest)
    print("record %s -> %s" % (RECORD, dest))
    try:
        man = manifest(dest)
    except FetchError as exc:
        print("REFUSED: %s" % exc, file=sys.stderr)
        return 2

    print("  title   : %s" % man["title"])
    print("  license : %s" % man["license"])
    if man["license"] != EXPECTED_LICENSE:
        print("REFUSED: licence is %r, expected %r. A licence change is a user decision, not a "
              "retry." % (man["license"], EXPECTED_LICENSE), file=sys.stderr)
        return 4

    picked = [f for f in man["files"] if SELECTIONS[a.select](f["key"])]
    total = sum(f["size"] for f in picked)
    print("  select  : %s -> %d of %d files, %.2f MB"
          % (a.select, len(picked), len(man["files"]), total / 1e6))

    if a.list:
        for f in man["files"]:
            mark = "*" if SELECTIONS[a.select](f["key"]) else " "
            print("   %s %-52s %10.2f MB  %s:%s"
                  % (mark, f["key"][:52], f["size"] / 1e6, f["algo"], f["digest"][:16]))
        return 0

    if total > 1e9 and not a.yes:
        print("REFUSED: selection is %.2f GB. Re-run with --yes if that is intended.\n"
              "  Note: the two dumps are ~99.7%% model checkpoints (measured by reading their zip\n"
              "  central directories over HTTP range requests); the DATA in them is ~10 MB of the\n"
              "  authors' own model predictions." % (total / 1e9), file=sys.stderr)
        return 5

    dest.mkdir(parents=True, exist_ok=True)
    results = []
    for f in picked:
        print("\n%s  (%.2f MB)" % (f["key"], f["size"] / 1e6))
        ok, why = fetch_one(f, dest)
        print("  %s %s" % ("OK  " if ok else "FAIL", why))
        results.append((f["key"], ok, why))

    kept = [k for k, ok, _ in results if ok]
    (dest / "_fetch_manifest.json").write_text(
        json.dumps({"record": RECORD, "license": man["license"], "doi": man["doi"],
                    "title": man["title"], "selection": a.select,
                    "verified": [dict(f, verified=True) for f in picked
                                 if f["key"] in kept]}, indent=2) + "\n",
        encoding="utf-8")

    print("\n%d/%d verified; manifest at %s" % (len(kept), len(picked),
                                                dest / "_fetch_manifest.json"))
    return 0 if len(kept) == len(picked) else 1


if __name__ == "__main__":
    raise SystemExit(main())
