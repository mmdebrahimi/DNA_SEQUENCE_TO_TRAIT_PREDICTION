"""Build the Mash reference sketch for the closed-set organism-identification router.

Read-mostly: enumerates locally cached genomes, labels each from its OWN assembly metadata, and emits
ONE `mash sketch` to D: plus a committed manifest. Writes nothing into the repo except the manifest.

    uv run python scripts/build_identify_reference.py [--limit N] [--dry-run]

TWO DEVIATIONS FROM THE PLAN, both measured, both recorded here rather than silently applied.

1. LABEL SOURCE IS THE ASSEMBLY'S OWN `ORGANISM` FIELD, NOT THE COHORT DIRECTORY NAME.
   The plan specified "cohort directory -> organism". Each cached genome ships an `annotations.gbk`
   whose GenBank header carries NCBI's own `  ORGANISM` line, and 1,111 of 1,129 genomes in the shared
   cache have one. That is strictly better: it is per-genome (so a mislabelled cohort cannot propagate
   to every genome under it) and it is the assembly's own metadata rather than an inference from a
   folder name. It is still NOT wet-lab identification -- submitter-provided organism names can be
   wrong -- so the evidence tier is unaffected.

2. NEAR-DUPLICATE REMOVAL IS NOT DONE HERE; IT MOVES TO VALIDATION.
   The plan de-duplicated at construction time via greedy-representative clustering. But duplicates do
   not hurt a nearest-neighbour router in DEPLOYMENT (a redundant neighbour changes no call) -- they
   only inflate LEAVE-ONE-OUT, because a held-out genome's near-identical twin stands in for it. That
   is a validation concern, and shrinking the deployed reference to fix a validation artifact would
   trade real coverage for a cleaner number. So the reference stays complete and Step 6 owns the
   exclusion (held-out genome AND anything within a lineage threshold of it).

MEASURED CORPUS as of 2026-10-04 (re-derive with --dry-run; do not trust this comment):
    in-set     2,026 genomes / 10 organisms, E. coli 1,379 (68%), K. oxytoca n=1
    out-of-set    48 genomes /  4 species, dominated by K. aerogenes (43) -- a CONGENER of an in-set
                  species, so the out-of-set control is the hard near-neighbour case, not a trivial one
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import sys
from datetime import date
from pathlib import Path

from dna_decode.data.organism_vocab import UnknownOrganism, canonical

REPO = Path(__file__).resolve().parent.parent
SHARED_CACHE = Path("D:/dna_decode_cache/refseq")
SKETCH_DIR = Path("D:/dna_decode_cache/identify")
SKETCH_PATH = SKETCH_DIR / "reference.msh"
MASH_IMAGE = "quay.io/biocontainers/mash:2.3--hb105d93_10"

#: `mash sketch` parameters, stamped into the manifest. k=21/s=1000 are Mash's own defaults for
#: genome-scale comparison and are what `eval/phylogeny.py` relies on implicitly; recording them makes
#: a later mismatch detectable (a sketch built at a different k cannot be compared against this one).
SKETCH_K = 21
SKETCH_S = 1000

#: A reference covering fewer than this many organisms cannot discriminate between organisms at all.
MIN_ORGANISMS = 2


class ReferenceBuildError(RuntimeError):
    """Raised instead of writing a reference that cannot do its job."""


def organism_from_genbank(gbk: Path, max_header_lines: int = 60) -> str | None:
    """The assembly's own `  ORGANISM` value, or None.

    Reads only the header. Returns None rather than raising on a missing/short/unreadable file, so one
    bad genome cannot abort a 2,000-genome enumeration -- but the caller MUST count the Nones, because
    a silently-shrinking corpus is the failure mode here.
    """
    try:
        with gbk.open(encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh):
                if line.startswith("  ORGANISM"):
                    got = line.split("ORGANISM", 1)[1].strip()
                    return got or None
                if i > max_header_lines:
                    return None
    except OSError:
        return None
    return None


def genus_species(organism_name: str) -> str:
    """First two whitespace tokens of a GenBank organism name, underscore-joined.

    `Escherichia coli 'BL21-Gold(DE3)pLysS AG'` -> `Escherichia_coli`. Strain suffixes are dropped
    because the vocabulary resolves at genus/species granularity, matching AMRFinder's `-O` list.
    """
    toks = [t for t in re.split(r"\s+", organism_name.strip()) if t]
    if not toks:
        return ""
    return "_".join(toks[:2]) if len(toks) >= 2 else toks[0]


def refseq_roots(repo: Path = REPO, shared: Path = SHARED_CACHE) -> list[Path]:
    """Every directory that holds `<ACCESSION>/genome.fna` -- the shared cache plus per-cohort caches.

    The shared cache alone is NOT enough and that was measured: it is 1,108 E. coli and 2 Klebsiella,
    so a shared-cache-only reference would cover 2 organisms and trip MIN_ORGANISMS.
    """
    roots = [shared] if shared.is_dir() else []
    # TWO on-disk layouts exist and the second was MISSED until 2026-10-04: a cohort may cache its
    # genomes under `refseq/` OR under `genomes/`. Measured repo-wide at the time of the fix:
    # 1,124 genomes under `data/raw/*/refseq/*/` and 8 under `data/raw/*/genomes/*/` -- the latter
    # being the ONLY local Candida auris genomes, so globbing one layout silently cost a whole
    # organism. The glob is enumerated from the filesystem rather than assumed, so a third layout
    # would show up as a count mismatch rather than vanishing.
    for layout in ("refseq", "genomes"):
        roots += sorted(p for p in (repo / "data" / "raw").glob(f"*/{layout}") if p.is_dir())
    return roots


def enumerate_labelled_genomes(roots: list[Path]) -> tuple[dict[str, list[tuple[str, Path]]], dict]:
    """Map canonical organism -> [(accession, fasta)], plus a stats dict.

    De-duplicates by ACCESSION across roots (the same genome is cached under several cohorts), keeping
    the first occurrence, so one genome never enters the reference twice.
    """
    per: dict[str, list[tuple[str, Path]]] = collections.defaultdict(list)
    out_of_set: collections.Counter = collections.Counter()
    seen: set[str] = set()
    stats = {"n_dirs_scanned": 0, "n_no_fasta": 0, "n_no_organism": 0, "n_duplicate_accession": 0}

    for root in roots:
        for p in sorted(root.iterdir()):
            if not p.is_dir():
                continue
            stats["n_dirs_scanned"] += 1
            fasta = p / "genome.fna"
            if not fasta.exists():
                stats["n_no_fasta"] += 1
                continue
            if p.name in seen:
                stats["n_duplicate_accession"] += 1
                continue
            seen.add(p.name)
            name = organism_from_genbank(p / "annotations.gbk")
            if not name:
                stats["n_no_organism"] += 1
                continue
            try:
                per[canonical(genus_species(name))].append((p.name, fasta))
            except UnknownOrganism:
                out_of_set[genus_species(name)] += 1

    stats["out_of_set_counts"] = dict(out_of_set.most_common())
    stats["n_out_of_set"] = sum(out_of_set.values())
    return dict(per), stats


def build_manifest(per: dict[str, list[tuple[str, Path]]], stats: dict) -> dict:
    """The committed sidecar. Records composition imbalance EXPLICITLY -- it bounds every downstream
    number, so it is reported, not hidden."""
    counts = {o: len(v) for o, v in sorted(per.items(), key=lambda kv: -len(kv[1]))}
    total = sum(counts.values())
    largest = max(counts.values()) if counts else 0
    return {
        "schema": "identify-reference-manifest-v1",
        "date": date.today().isoformat(),
        "sketch_path": str(SKETCH_PATH),
        "sketch_committed": False,
        "sketch_k": SKETCH_K,
        "sketch_s": SKETCH_S,
        "mash_image": MASH_IMAGE,
        "label_source": "assembly's own GenBank `  ORGANISM` header field (NCBI metadata), "
                        "NOT the cohort directory name",
        "label_is_wet_lab": False,
        "n_genomes": total,
        "n_organisms": len(counts),
        "per_organism": counts,
        "largest_organism_share": round(largest / total, 4) if total else None,
        "unscorable_under_leave_one_out": [o for o, n in counts.items() if n < 2],
        "out_of_set": {
            "n_genomes": stats.get("n_out_of_set", 0),
            "counts": stats.get("out_of_set_counts", {}),
            "note": "the Step 6 control. Dominated by congeners of in-set species "
                    "(e.g. Klebsiella aerogenes vs K. pneumoniae), so it is the HARD "
                    "near-neighbour case rather than a trivially distant one.",
        },
        "accessions_by_organism": {o: sorted(a for a, _ in v) for o, v in sorted(per.items())},
        "enumeration_stats": {k: v for k, v in stats.items()
                              if k not in ("out_of_set_counts", "n_out_of_set")},
    }


#: Host drive mounted into the container, and the mount point inside it. Every reference genome
#: resolves under D: -- the per-cohort ones via the `data/` junction
#: (`data/raw/...` -> `D:\dna_decode_data\raw\...`) and the shared cache directly
#: (``D:\dna_decode_cache\refseq\...`). VERIFIED by resolving real paths, not
#: assumed: the genomes sit in 40 separate root directories, and Docker needs ONE mount point, so a
#: single `D:/` mount is what makes a single `mash sketch` call possible at all.
HOST_DRIVE = Path("D:/")
CONTAINER_MOUNT = "/d"


def to_container_path(host: Path, drive: Path | None = None, mount: str | None = None) -> str:
    """Rewrite a resolved host path under `drive` to its in-container path.

    RAISES on a path outside the mounted drive. That is fail-closed on purpose: a genome on C: would
    simply not exist inside the container, and `mash` would skip it, so the sketch would quietly cover
    fewer genomes than the manifest claims -- a silently-shrinking reference.

    `drive`/`mount` default to None and are resolved from the module globals AT CALL TIME, not bound
    as default arguments. A default argument would freeze the import-time value into the signature,
    which makes the mount un-overridable and the guard untestable.
    """
    drive = drive if drive is not None else HOST_DRIVE
    mount = mount if mount is not None else CONTAINER_MOUNT
    r = host.resolve()
    try:
        rel = r.relative_to(drive.resolve())
    except ValueError as exc:
        raise ReferenceBuildError(
            f"{r} is not under the mounted drive {drive}; it would be invisible inside the container "
            f"and silently dropped from the sketch") from exc
    return mount + "/" + rel.as_posix()


def _mash_sketch(fastas: list[Path], out_prefix: Path) -> str:
    """One `mash sketch` call over every reference genome, via the no-shell Docker boundary.

    A file-of-filenames (`-l`) is used rather than thousands of argv entries: at 2,000+ genomes the
    command line exceeds the OS limit, and that failure looks like a Mash error rather than an argv
    problem. `-l` is a VERIFIED flag on the pinned image (`mash sketch -h`: "List input").

    Each input FILE becomes one sketch entry, so the `reference-ID` that `mash dist` later emits is
    the genome's path -- whose parent directory name is the accession, which the manifest maps to an
    organism.
    """
    from tools.docker_runner import run as docker_run

    SKETCH_DIR.mkdir(parents=True, exist_ok=True)
    paths = [to_container_path(f) for f in fastas]          # raises before any container starts
    fof = SKETCH_DIR / "reference_inputs.txt"
    # newline="\n" IS LOAD-BEARING, measured the hard way. On Windows `write_text` translates "\n" to
    # "\r\n", so every line of the list ends with a stray CR and mash tries to open "genome.fna\r" ->
    # `ERROR: could not open ... for reading` on a file that is demonstrably present and readable
    # inside the container. The CR also returns the terminal cursor, which mangles mash's own error
    # message and sends you looking at the mount instead of the list. Pinned by
    # tests/test_build_identify_reference.py::test_the_input_list_is_written_with_LF_endings.
    fof.write_text("\n".join(paths) + "\n", encoding="utf-8", newline="\n")

    res = docker_run(
        MASH_IMAGE,
        ["mash", "sketch", "-l", to_container_path(fof),
         "-o", to_container_path(out_prefix.with_suffix("")),
         "-k", str(SKETCH_K), "-s", str(SKETCH_S)],
        mounts={str(HOST_DRIVE): CONTAINER_MOUNT},
        capture_output=True, check=True, timeout=7200,
    )
    if not out_prefix.exists():
        raise ReferenceBuildError(
            f"mash reported success but {out_prefix} does not exist -- treat exit 0 as unverified and "
            f"check the mount. stdout/stderr: {getattr(res, 'stdout', '')!r} "
            f"{getattr(res, 'stderr', '')!r}")
    return str(out_prefix)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="enumerate + write the manifest, do NOT invoke mash")
    ap.add_argument("--limit", type=int, default=None,
                    help="cap genomes per organism (smoke runs only; stamped into the manifest)")
    ap.add_argument("--manifest-out", type=Path, default=None)
    a = ap.parse_args(argv)

    roots = refseq_roots()
    if not roots:
        print("NO refseq roots found -- nothing to enumerate", file=sys.stderr)
        return 2
    per, stats = enumerate_labelled_genomes(roots)

    if a.limit:
        per = {o: v[: a.limit] for o, v in per.items()}

    if len(per) < MIN_ORGANISMS:
        print(f"REFUSING: reference would cover {len(per)} organism(s); a reference with fewer than "
              f"{MIN_ORGANISMS} cannot discriminate between organisms at all.", file=sys.stderr)
        return 3

    manifest = build_manifest(per, stats)
    manifest["limit_per_organism"] = a.limit
    out = a.manifest_out or (REPO / "wiki" / f"identify_reference_manifest_{date.today().isoformat()}.json")
    out.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"roots={len(roots)}  in-set genomes={manifest['n_genomes']}  "
          f"organisms={manifest['n_organisms']}  "
          f"largest share={manifest['largest_organism_share']}")
    for o, n in manifest["per_organism"].items():
        print(f"   {o:<28} {n}")
    if manifest["unscorable_under_leave_one_out"]:
        print(f"   NOTE unscorable under leave-one-out (n<2): "
              f"{manifest['unscorable_under_leave_one_out']}")
    print(f"out-of-set control: {manifest['out_of_set']['n_genomes']} genomes, "
          f"{len(manifest['out_of_set']['counts'])} species")
    print(f"manifest: {out}")

    if a.dry_run:
        print("--dry-run: sketch NOT built")
        return 0

    try:
        _mash_sketch([f for v in per.values() for _, f in v], SKETCH_PATH)
    except ReferenceBuildError as exc:
        print(f"SKETCH NOT BUILT: {exc}", file=sys.stderr)
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
