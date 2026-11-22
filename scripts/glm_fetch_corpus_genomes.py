"""Fetch the genomes a corpus spec names, and VERIFY each one is the organism it was selected as.

    uv run python scripts/glm_fetch_corpus_genomes.py <spec.json> [--cache-root DIR] [--dry-run]

**A thin driver over `dna_decode/data/refseq.py`, which already does the hard part.** That module is an
NCBI Datasets v2 downloader with retry, atomic writes, idempotent caching and a completeness sentinel,
and it produced the 459 assemblies this project has been measuring. Nothing here re-implements any of
it; this adds the loop, the skip, and one gate it does not have.

**The gate, and why it is not optional.** `download_genome` takes an ACCESSION and nothing else. It
cannot know which organism a caller believed that accession to be, so a resolution pass that returned
the wrong assembly -- a synonym, a reclassified name, a different strain sharing a label, a transposed
digit -- produces a perfectly successful download of the wrong genome. It then enters the corpus under
the species label the spec asserted, and nothing downstream contradicts it: the sequence is real, the
annotation parses, the windows extract, the model trains. The error would surface, if ever, as a corpus
that does not behave like its own description.

So every entry must carry the organism name as NCBI reports it, and this script compares that against
the organism in the downloaded FASTA header. The comparison is GRADED rather than binary, because exact
string equality is the wrong bar -- headers vary in format and taxonomy genuinely gets revised:

  * **genus mismatch -> REFUSE.** Different genus is not a naming variation.
  * **genus matches, species epithet differs -> FLAG for review**, keep the file, exclude from the
    verified set. This is where a reclassification or a synonym lands, and a human should decide.
  * **no organism in the header -> UNVERIFIABLE**, which is reported as its own state and never counted
    as a pass. Some assemblies carry only an accession and a contig name.

**`accession: null` is refused, loudly, with the reason.** A corpus spec may legitimately be written
with species and strain and no accession -- but `refseq.py` has no taxon lookup, so there is nothing in
this repo that can resolve one. Failing with that stated beats failing with a confusing HTTP error.

**Disk is checked before fetching, not after.** This host runs at 97% full with roughly 19 GB free, and
this project carries a standing rule that heavy files live on an external drive which is currently
absent. Sixteen bacterial assemblies at about 5 MB each is trivial; a eukaryotic batch is not, and the
difference should stop a run rather than fill a disk.

Network: yes, to NCBI only, read-only, via the existing downloader. No GPU, no money.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CACHE = ROOT / "data" / "refseq_cache"
OUT_DIR = ROOT / "wiki"

MIN_FREE_GB = 2.0          # refuse rather than fill the last of a 97%-full disk
APPROX_MB_PER_BACTERIAL_GENOME = 6.0

_BINOMIAL = re.compile(r"\b([A-Z][a-z]{2,})\s+([a-z]{2,})\b")


def parse_organism_from_header(header: str) -> tuple[str, str] | None:
    """Pull `(genus, epithet)` out of a FASTA header, or None if no binomial is present.

    Headers in this cache look like
    `>DBQBFQ010000001.1 TPA_asm: Escherichia coli isolate BCW12650 SAMN4294...`,
    so the organism is embedded in free text rather than in a fixed field. The first capitalised word
    followed by a lowercase word is the binomial in every real example checked; `TPA_asm` does not match
    because the regex requires a space-separated lowercase second token of at least two letters.

    Returning None rather than guessing is the point: an unparsed header must become UNVERIFIABLE, not
    a silent pass.
    """
    m = _BINOMIAL.search(header)
    if not m:
        return None
    return m.group(1), m.group(2)


def verify_organism(expected: str, header: str) -> dict:
    """Graded comparison of a spec's expected organism against a downloaded header.

    Returns `{"state": VERIFIED | FLAG_SPECIES_MISMATCH | REFUSE_GENUS_MISMATCH | UNVERIFIABLE, ...}`.
    Graded because exact equality is the wrong bar -- taxonomy gets revised and header formats vary --
    while a genus mismatch is never a naming variation.
    """
    got = parse_organism_from_header(header)
    want = parse_organism_from_header(expected if expected[:1].isupper() else expected.capitalize())
    if want is None:
        return {"state": "UNVERIFIABLE", "reason": f"no binomial in expected organism {expected!r}",
                "expected": expected, "header_organism": None}
    if got is None:
        return {"state": "UNVERIFIABLE",
                "reason": "no binomial in the downloaded FASTA header; this assembly cannot be "
                          "verified by name and must not be counted as verified",
                "expected": expected, "header_organism": None}
    if got[0].lower() != want[0].lower():
        return {"state": "REFUSE_GENUS_MISMATCH",
                "reason": f"expected genus {want[0]}, header says {got[0]} -- a different genus is not "
                          "a naming variation, so this is the wrong assembly",
                "expected": expected, "header_organism": " ".join(got)}
    if got[1].lower() != want[1].lower():
        return {"state": "FLAG_SPECIES_MISMATCH",
                "reason": f"genus {got[0]} matches but epithet differs ({want[1]} vs {got[1]}); this is "
                          "where a reclassification or synonym lands, so a human should decide",
                "expected": expected, "header_organism": " ".join(got)}
    return {"state": "VERIFIED", "reason": "genus and epithet match",
            "expected": expected, "header_organism": " ".join(got)}


def _expected_organism(entry: dict) -> str:
    """The organism string to verify a download against, preferring NCBI AUTHORITATIVE over the spec.

    Field order matters and is not cosmetic. A spec author names a species when SELECTING, then resolves
    an accession and records what NCBI reports for it -- and those two can differ, because taxonomy gets
    revised at genus level (Clostridium to Clostridioides, Borrelia to Borreliella are both real) and
    because a resolver may return the current reference strain rather than the classic one. Verifying
    against the selector name instead of the resolved name would false-reject exactly the entries where
    a reclassification happened, which are the ones most worth keeping.

    An earlier version read `ncbi_organism` while the handed-over spec writes `ncbi_organism_name`, so it
    silently fell through to `species`. On that particular spec it was harmless -- the author had already
    used the modern genus in both fields, so 0 of 16 would have mismatched -- but it ignored the
    authoritative string the whole gate was built to compare against, and the next spec need not be so
    forgiving. Both spellings are accepted now so a field rename cannot reintroduce the silence.
    """
    return (entry.get("ncbi_organism_name")
            or entry.get("ncbi_organism")
            or entry.get("species")
            or "")


def check_disk(cache_root: Path, n_genomes: int) -> dict:
    """Refuse before fetching if the disk cannot take it. Checked up front, not discovered."""
    usage = shutil.disk_usage(cache_root if cache_root.exists() else cache_root.parent)
    free_gb = usage.free / 1e9
    need_gb = n_genomes * APPROX_MB_PER_BACTERIAL_GENOME / 1000
    return {
        "free_gb": round(free_gb, 2),
        "estimated_need_gb": round(need_gb, 3),
        "min_free_gb_required": MIN_FREE_GB,
        "ok": bool(free_gb - need_gb >= MIN_FREE_GB),
        "note": ("estimate assumes bacterial-scale assemblies; a eukaryotic batch is orders of "
                 "magnitude larger and this check is what should stop it"),
    }


def load_spec(path: Path) -> tuple[dict, list[dict]]:
    rec = json.loads(path.read_text(encoding="utf-8"))
    genomes = rec.get("genomes") or []
    if not genomes:
        raise SystemExit(f"REFUSE: no `genomes` array in {path}")
    return rec, genomes


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("spec", help="corpus species-selection spec JSON")
    p.add_argument("--cache-root", default=str(DEFAULT_CACHE))
    p.add_argument("--dry-run", action="store_true",
                   help="validate the spec and the disk, fetch nothing")
    p.add_argument("--allow-unverifiable", action="store_true",
                   help="keep assemblies whose organism cannot be checked by name")
    a = p.parse_args(argv)

    spec_path = Path(a.spec)
    if not spec_path.exists():
        print(f"REFUSE: spec not found at {spec_path}")
        return 2
    rec, genomes = load_spec(spec_path)
    cache_root = Path(a.cache_root)

    print(f"spec: {rec.get('spec', spec_path.name)}  ({len(genomes)} genomes)")
    print(f"cache root: {cache_root}")

    # --- the accession gate, before anything else ---
    missing = [g for g in genomes if not g.get("accession")]
    if missing:
        names = ", ".join(f"{g.get('species', '?')} {g.get('strain', '')}".strip()
                          for g in missing[:5])
        print(f"REFUSE: {len(missing)} of {len(genomes)} entries have no accession "
              f"(e.g. {names}{'...' if len(missing) > 5 else ''}).")
        print("  `dna_decode/data/refseq.py` has NO taxon lookup -- every entry point takes an "
              "accession, and `_build_dataset_url` hits /genome/accession/<acc>/download. There is "
              "nothing in this repo that resolves a species name to an assembly, so these cannot be "
              "fetched as written. Ask the spec author for resolved GCF/GCA accessions plus the "
              "organism name NCBI reports for each, which is what the verification gate compares "
              "against.")
        return 3

    missing_org = [g for g in genomes if not _expected_organism(g)]
    if missing_org:
        print(f"REFUSE: {len(missing_org)} entries carry neither `ncbi_organism` nor `species`, so "
              "the downloaded assembly could not be checked against what it was selected as.")
        return 3

    disk = check_disk(cache_root, len(genomes))
    print(f"disk: {disk['free_gb']} GB free, need about {disk['estimated_need_gb']} GB "
          f"-> {'OK' if disk['ok'] else 'INSUFFICIENT'}")
    if not disk["ok"]:
        print(f"REFUSE: fewer than {MIN_FREE_GB} GB would remain. This host runs near-full and the "
              "project's heavy-file drive is absent; send this batch to a host with room instead.")
        return 3

    if a.dry_run:
        print("\ndry run: spec and disk validated, nothing fetched.")
        return 0

    from dna_decode.data.refseq import download_genome, fasta_path, is_cache_complete, cache_dir_for

    t0 = time.time()
    rows, counts = [], {}
    for g in genomes:
        acc = g["accession"]
        expected = _expected_organism(g)
        row = {"accession": acc, "species": g.get("species"), "strain": g.get("strain"),
               "domain": g.get("domain"), "lineage_key": g.get("lineage_key"),
               "expected_organism": expected}
        cached = is_cache_complete(cache_dir_for(acc, cache_root))
        try:
            if cached:
                row["fetch"] = "SKIPPED_ALREADY_CACHED"
            else:
                download_genome(acc, cache_root)
                row["fetch"] = "DOWNLOADED"
        except Exception as exc:                       # noqa: BLE001 - report, do not abort the batch
            row["fetch"] = "FAILED"
            row["error"] = f"{type(exc).__name__}: {exc}"
            row["verification"] = {"state": "NOT_ATTEMPTED"}
            counts["FAILED"] = counts.get("FAILED", 0) + 1
            rows.append(row)
            print(f"  {acc:20s} FAILED  {row['error'][:70]}")
            continue

        try:
            header = fasta_path(acc, cache_root).open(encoding="utf-8",
                                                      errors="replace").readline().strip()
        except OSError as exc:
            header = ""
            row["header_read_error"] = str(exc)
        row["fasta_header"] = header[:200]
        v = verify_organism(expected, header)
        row["verification"] = v
        counts[v["state"]] = counts.get(v["state"], 0) + 1
        rows.append(row)
        print(f"  {acc:20s} {row['fetch']:24s} {v['state']:24s} {v.get('header_organism') or '-'}")

    verified = [r for r in rows if r["verification"]["state"] == "VERIFIED"]
    usable = list(verified)
    if a.allow_unverifiable:
        usable += [r for r in rows if r["verification"]["state"] == "UNVERIFIABLE"]

    out = OUT_DIR / f"glm_corpus_fetch_{date.today().isoformat()}.json"
    out.write_text(json.dumps({
        "record": "glm-corpus-fetch-v1",
        "date": date.today().isoformat(),
        "spec": str(spec_path),
        "spec_name": rec.get("spec"),
        "cache_root": str(cache_root),
        "disk": disk,
        "counts_by_verification_state": counts,
        "n_verified": len(verified),
        "n_usable": len(usable),
        "usable_accessions": [r["accession"] for r in usable],
        "rows": rows,
        "verification_policy": {
            "genus_mismatch": "REFUSE -- a different genus is not a naming variation",
            "species_epithet_mismatch": "FLAG, kept on disk but excluded from the usable set; this is "
                                        "where reclassifications and synonyms land",
            "no_binomial_in_header": "UNVERIFIABLE -- its own state, never counted as a pass; "
                                     "--allow-unverifiable opts it into the usable set explicitly",
            "why": ("download_genome takes an accession and cannot know which organism the caller "
                    "believed it to be, so a mis-resolved accession is a fully successful download of "
                    "the wrong genome that enters the corpus under the spec's label with nothing "
                    "downstream to contradict it"),
        },
        "honest_limits": [
            "Name verification checks the organism string in the FASTA header. It catches a "
            "mis-resolved accession; it does NOT verify the assembly is the right STRAIN, nor that the "
            "sequence is what the accession claims.",
            "The binomial parser takes the first capitalised-then-lowercase word pair in the header. "
            "That matched every real example checked in this cache, but header formats vary across "
            "submitters and a miss becomes UNVERIFIABLE rather than a wrong pass.",
            "The disk estimate assumes bacterial-scale assemblies at about 6 MB each. It is not a "
            "safe guide for a eukaryotic batch, which is orders of magnitude larger.",
            "A FLAG is not a rejection and a REFUSE is not a deletion: files stay on disk for a human "
            "to adjudicate, and only the usable set is offered to the corpus builder.",
        ],
    }, indent=2), encoding="utf-8")

    print(f"\nverification: {counts}")
    print(f"usable for the corpus: {len(usable)}/{len(rows)}")
    print(f"{round(time.time() - t0, 1)}s -> {out.relative_to(ROOT).as_posix()}")
    if counts.get("REFUSE_GENUS_MISMATCH"):
        print("REFUSE: at least one assembly is the wrong genus; fix the spec before building a corpus.")
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
