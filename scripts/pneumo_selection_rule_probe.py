"""Does the coverage-vs-identity selection rule change ANY pneumococcal call? Measure before changing.

THE QUESTION, and why it is open. Validating the E. coli serotype cell found a live defect: allele
selection by COVERAGE ONLY let cross-hybridizing fliC alleles win, and switching to identity-primary
lifted H accuracy 0.770 -> 0.926. A sweep then found the SAME coverage-first pattern in
`pneumoserotype` (`key = (coverage, identity)`) and `plasmid`. Those were deliberately NOT changed,
because the biology differs -- pneumococcal typing matches a WHOLE cps LOCUS against reference
sequences, not per-antigen alleles -- so identity-primary is not automatically correct there.

WHY A CHEAP PROBE FIRST. Settling it properly means re-running the 260-isolate Quellung cohort under
both rules, which costs ~235 ENA assembly fetches. Before spending that, this asks the far cheaper
question on the assemblies already cached on disk: **does the rule flip any call at all?**

  - flips 0 of N  -> weak evidence the rule is inert here; the expensive run is poorly justified
  - flips some    -> the rule is live for this cell and the full cohort run IS justified

A zero here is NOT proof the rule is harmless -- it is a small sample, and it is reported as such.
The probe also refuses to report if it cannot demonstrate it actually exercised both rules.

Offline: uses cached assemblies + the committed cps DB + native blastn.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.typing.blast_caller import call_alleles  # noqa: E402
from dna_decode.pneumoserotype.runner import serotype_of  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"


def best_under(per_allele: dict, identity_primary: bool):
    """The winning reference under one ordering. Returns (ref_id, serotype, pid, cov) or None."""
    best = None
    for ref_id, hit in per_allele.items():
        if not hit["called"]:
            continue
        key = ((hit["percent_identity"], hit["percent_coverage"]) if identity_primary
               else (hit["percent_coverage"], hit["percent_identity"]))
        if best is None or key > best[1]:
            best = (ref_id, key, hit)
    if best is None:
        return None
    ref_id, _, hit = best
    return (ref_id, serotype_of(ref_id), hit["percent_identity"], hit["percent_coverage"])


def fetch_assemblies(cohort: Path, dest: Path, n: int) -> list[Path]:
    """Re-fetch N assemblies from ENA, VALIDATING each before keeping it.

    The existing cache exists because the previous fetcher wrote whatever the server returned: 199-byte
    HTTP 403 pages, zero-byte files, and gzip streams truncated mid-download, all under genome
    filenames. So every download here is decompressed IN FULL and checked to be FASTA before it is
    written to disk -- a partial write is discarded rather than cached as a genome.
    """
    import gzip as _gz
    import urllib.request
    ua = {"User-Agent": "dna_decode/0.13 (research; genotype-phenotype validation)"}
    dest.mkdir(parents=True, exist_ok=True)
    kept: list[Path] = []
    rows = [ln.split("	") for ln in cohort.read_text(encoding="utf-8").splitlines() if ln.strip()]
    for r in rows:
        if len(kept) >= n:
            break
        ers = r[1].strip() if len(r) > 1 else ""
        if not ers:
            continue
        out = dest / f"{ers}.fetched.fa"
        if out.exists() and out.stat().st_size > 1024:
            kept.append(out)
            continue
        try:
            q = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ers}"
                 "&result=analysis&fields=generated_ftp&format=tsv")
            tsv = urllib.request.urlopen(urllib.request.Request(q, headers=ua),
                                         timeout=180).read().decode("utf8", "replace")
            ftp = [x for x in tsv.splitlines()[1:] if x.strip()]
            if not ftp:
                continue
            url = "https://" + ftp[0].split("	")[-1].strip()
            blob = urllib.request.urlopen(urllib.request.Request(url, headers=ua), timeout=600).read()
            body = _gz.decompress(blob)          # raises on a truncated stream -> discarded
            if not body.lstrip().startswith(b">"):
                continue
            out.write_bytes(body)
            kept.append(out)
            print(f"    fetched {ers} ({len(body)} bp uncompressed)", flush=True)
        except Exception as e:                   # noqa: BLE001
            print(f"    skip {ers}: {type(e).__name__}", flush=True)
    return kept


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--asm-dir", type=Path, default=Path("D:/dna_decode_cache/pneumo_gps/asm"))
    ap.add_argument("--cohort", type=Path,
                    default=Path("D:/dna_decode_cache/pneumo_gps/poland_quellung_cohort.tsv"))
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data" / "pneumoserotype_db" / "cps_references.fasta")
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--fetch", type=int, default=0,
                    help="re-fetch N assemblies from ENA (the cached ones are corrupt)")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"pneumo_selection_rule_probe_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    if not a.db.exists():
        print(f"cps DB absent at {a.db}", file=sys.stderr)
        return 2
    files = sorted(a.asm_dir.glob("*.fa.gz")) + sorted(a.asm_dir.glob("*.fa")) + \
        sorted(a.asm_dir.glob("*.fasta"))
    if not files:
        print(f"no cached assemblies in {a.asm_dir}", file=sys.stderr)
        return 2
    corrupt: list[dict] = []

    if a.fetch:
        files = list(fetch_assemblies(a.cohort, a.asm_dir, a.fetch)) or files

    # measured labels, keyed by the ERS accession embedded in the filename
    labels: dict[str, str] = {}
    if a.cohort.exists():
        for line in a.cohort.read_text(encoding="utf-8").splitlines():
            p = line.rstrip("\n").split("\t")
            # The cohort TSV is ERR<tab>ERS<tab>serotype. Assemblies are named by the SAMPLE (ERS)
            # accession, so keying on column 0 (the RUN accession) yields measured=None for EVERY
            # row -- which reads as "no label available" rather than "looked up the wrong key", and
            # would quietly make every comparison against the measured serotype impossible.
            if len(p) >= 3:
                labels[p[1].strip()] = p[2].strip()
            elif len(p) >= 2 and p[0].upper().startswith("ER"):
                labels[p[0].strip()] = p[1].strip()

    # INPUT VALIDATION FIRST. The cached .fa.gz files are NOT all assemblies: the ENA fetcher wrote
    # whatever the server returned, so 199-byte HTTP 403 pages and zero-byte files sit in the cache
    # under genome filenames. Feeding those to blastn yields a generic "invocation failed" that is
    # indistinguishable from a real no-match -- a corrupt input masquerading as a biological result.
    valid, corrupt = [], []
    for fp in files:
        try:
            if fp.stat().st_size < 1024:
                corrupt.append({"file": fp.name, "why": f"too small ({fp.stat().st_size} bytes)"})
                continue
            if fp.suffix == ".gz":
                import gzip as _gz
                # Decompress the WHOLE stream, not just a head. A truncated gzip decompresses its
                # first blocks perfectly and only raises EOFError at the end -- so a head-only check
                # passes exactly the files that will later fail, which is the worst kind of guard.
                with _gz.open(fp, "rb") as fh:
                    body = fh.read()
                if not body.lstrip().startswith(b">"):
                    corrupt.append({"file": fp.name, "why": "decompressed content is not FASTA"})
                    continue
            valid.append(fp)
        except Exception as e:                       # noqa: BLE001
            corrupt.append({"file": fp.name, "why": f"{type(e).__name__}: {str(e)[:80]}"})
    print(f"cached files: {len(files)} | usable assemblies: {len(valid)} | CORRUPT: {len(corrupt)}")
    for c in corrupt[:6]:
        print(f"    {c['file']}: {c['why']}")
    files = valid

    rows, flips, tie_at_max = [], 0, 0
    for fp in files:
        ers = fp.name.split(".")[0]
        try:
            # blastn cannot read a gzipped FASTA -- it exits non-zero and the caller reports a generic
            # "invocation failed", which is indistinguishable from a real no-match. Decompress first.
            target = fp
            tmp = None
            if fp.suffix == ".gz":
                import gzip
                import tempfile
                tmp = Path(tempfile.mkstemp(suffix=".fa")[1])
                with gzip.open(fp, "rb") as src, open(tmp, "wb") as dst:
                    dst.write(src.read())
                target = tmp
            try:
                res = call_alleles(target, a.db, identity_threshold=85.0, coverage_threshold=60.0,
                                   blastn_bin=a.blastn, timeout=600)
            finally:
                if tmp is not None:
                    tmp.unlink(missing_ok=True)
        except Exception as e:                       # noqa: BLE001
            rows.append({"ers": ers, "status": f"error:{type(e).__name__}", "error": str(e)[:150]})
            continue
        if res.get("status") != "ok":
            rows.append({"ers": ers, "status": res.get("status"), "reason": res.get("reason")})
            continue
        cov_first = best_under(res["per_allele"], identity_primary=False)
        id_first = best_under(res["per_allele"], identity_primary=True)
        rec = {"ers": ers, "status": "ok", "measured": labels.get(ers),
               "coverage_primary": (cov_first[1] if cov_first else None),
               "identity_primary": (id_first[1] if id_first else None),
               "cov_primary_metrics": (cov_first[2:] if cov_first else None),
               "id_primary_metrics": (id_first[2:] if id_first else None)}
        rec["flipped"] = rec["coverage_primary"] != rec["identity_primary"]
        if rec["flipped"]:
            flips += 1
        if cov_first and cov_first[2] >= 100.0 and cov_first[3] >= 100.0:
            tie_at_max += 1
        rows.append(rec)

    scored = [r for r in rows if r.get("status") == "ok"]
    # NON-VACUITY: if no assembly produced a call at all, a zero-flip result says nothing about the
    # rule -- it says the probe never exercised it. Refuse rather than report a hollow zero.
    if not scored:
        print("REFUSING: no assembly produced a cps call, so neither rule was exercised.",
              file=sys.stderr)
        return 3

    n = len(scored)
    print(f"cached assemblies scored: {n}")
    print(f"  winner already at 100 identity AND 100 coverage (rule cannot matter): {tie_at_max}")
    print(f"  calls that FLIP between the two rules: {flips}/{n}")
    for r in scored:
        if r["flipped"]:
            print(f"    {r['ers']} measured={r['measured']} "
                  f"coverage-primary={r['coverage_primary']} -> identity-primary={r['identity_primary']}")

    # A flip is only INTERESTING if it changes a wrong answer into a right one. Counting flips alone
    # would call a rule "live" when it merely swaps one miss for another, which is what happened here.
    def _sg(s):
        if not s:
            return None
        t = "".join(ch for ch in str(s) if ch.isdigit() or ch.isalpha())
        i = 0
        while i < len(t) and t[i].isdigit():
            i += 1
        return t[:i] or t
    improved = worsened = neutral = 0
    for r in scored:
        if not r["flipped"] or not r.get("measured"):
            continue
        m = _sg(r["measured"])
        cov_ok, id_ok = _sg(r["coverage_primary"]) == m, _sg(r["identity_primary"]) == m
        if id_ok and not cov_ok:
            improved += 1
        elif cov_ok and not id_ok:
            worsened += 1
        else:
            neutral += 1
    print(f"  of those flips (serogroup level vs measured): improved={improved} "
          f"worsened={worsened} neither-right={neutral}")

    if flips and improved == 0 and worsened == 0:
        verdict = "RULE_FLIPS_BUT_NEVER_IMPROVES_ON_SAMPLE"
        why = (f"{flips} of {n} calls change with the ordering, but NONE of the flips turns a wrong "
               f"serogroup into a right one ({neutral} swap one miss for another). On this sample "
               "there is no evidence identity-primary helps the pneumococcal cell, which is the "
               "opposite of the E. coli result and consistent with its misses having a DIFFERENT "
               "cause (a single-best-reference v0 ceiling, not cross-hybridization).")
    elif flips == 0:
        verdict = "RULE_INERT_ON_CACHED_SAMPLE"
        why = (f"neither call changed on {n} cached assemblies, so on this sample the coverage-vs-"
               "identity ordering is inert for the pneumococcal cell. That is WEAK evidence from a "
               "small sample, NOT proof the rule is harmless -- but it does not justify the ~235 ENA "
               "fetches a full both-rules cohort run would cost.")
    else:
        verdict = "RULE_IS_LIVE_FULL_RUN_JUSTIFIED"
        why = (f"{flips} of {n} calls change with the ordering, so the rule is live for this cell and "
               "the full Quellung-cohort run under both rules IS justified -- only measured labels can "
               "say which ordering is right.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "pneumo-selection-rule-probe-v1", "date": _date.today().isoformat(),
           "question": "does coverage-primary vs identity-primary change any pneumococcal serotype call?",
           "n_cached_scored": n, "n_flipped": flips, "n_winner_at_max_both_axes": tie_at_max,
           "flip_outcomes_vs_measured_serogroup": {"improved": improved, "worsened": worsened,
                                                   "neither_right": neutral},
           "cache_integrity": {"n_files": len(files) + len(corrupt), "n_usable": len(files),
                               "n_corrupt": len(corrupt), "corrupt": corrupt,
                               "note": "the ENA fetcher wrote server responses verbatim, so HTTP 403 "
                                       "pages and zero-byte files are cached under genome filenames"},
           "rows": rows, "verdict": verdict, "why": why,
           "honest_limits": [
               "SMALL SAMPLE: only the assemblies still cached on disk, not the 260-isolate cohort. A "
               "zero-flip result is weak evidence of inertness, never proof.",
               "This measures whether the CALL changes, not whether it changes for the BETTER. Only "
               "the measured Quellung labels can decide that, and that needs the full run.",
               "The pneumococcal cell's documented misses are WITHIN-serogroup (9A/9V, 6B/6E, 15B/15C) "
               "and are attributed to a single-best-reference v0 ceiling -- a different cause from the "
               "cross-hybridization defect that the E. coli fix addressed.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
