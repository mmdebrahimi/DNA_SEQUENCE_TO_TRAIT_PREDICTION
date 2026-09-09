"""Score the serovar caller against the REFERENCE TOOL, pinned — not against an undocumented field.

THE LIMIT THIS CLOSES. The 2026-09-04 salmserovar validation measured our caller at 0.783 against a
comparator of 0.925 and called that a -0.142 delta "vs the in-silico incumbent". That comparator was
NCBI Pathogen Detection's published `computed_types` -- a PRODUCTION call of undocumented tool, version
and configuration. A delta against a production field cannot separate two very different readings:

    (a) our caller is worse than the reference METHOD, or
    (b) our caller merely diverges from one undocumented NCBI pipeline.

The artifact recorded that limit and said resolving it needed SeqSero2 installed and pinned locally.
This is that run: SeqSero2 1.3.2 from the pinned biocontainer, on the SAME 200 isolates, against the
SAME wet-lab agglutination labels.

FAIRNESS IS STRUCTURAL. All three callers -- ours, SeqSero2, and NCBI's field -- are scored through the
SAME `equivalence.equivalent`, which grants only notation normalisation plus the committed
Kauffmann-White formula table and refuses fuzzy near-misses. No caller can be granted leniency another
is denied. Abstentions (`no_call`) are counted SEPARATELY from misses for every caller, because
abstention is not error and pooling them would flatter whichever caller abstains most.

WHY THE WET-LAB LABEL STILL ANCHORS THIS. Salmonella's gold standard is slide agglutination -- an
antisera reaction, not a computation -- so all three in-silico callers are being scored against
something none of them produced.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.salmserovar.equivalence import equivalent, load_formula_index  # noqa: E402
from dna_decode.salmserovar.runner import call_serovar  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
SEROVAR_DB = ROOT / "data" / "salmserovar_db"
IMAGE = "quay.io/biocontainers/seqsero2:1.3.2--pyhdfd78af_0"
# A Windows-form path is REQUIRED for the bind mount: a Git-Bash `/tmp/...` source silently produces an
# anonymous volume instead, and the container then sees an empty /data while the command still exits 0.
STAGE = Path("D:/dna_decode_cache/ss2_run")

NON_SPECIFIC = {"", "-", "na", "n/a", "none", "unknown", "undetermined", "pending",
                "not applicable", "not determined", "untypeable", "untypable"}


def run_seqsero2(fasta: Path, stage: Path, image: str, timeout: int = 900) -> dict:
    """One assembly -> SeqSero2's per-axis + serotype call. Never raises.

    SeqSero2 BASENAMES its input and chdirs, so the file must be staged into the mount and the
    container's working directory set to it; passing an absolute /data/x.fna fails with a
    FileNotFoundError naming the bare basename.
    """
    stage.mkdir(parents=True, exist_ok=True)
    # Clear ONLY what this function owns. An earlier version wiped the whole staging dir, which also
    # deleted the append-open checkpoint living beside it -- destroying the resume record, and on
    # Windows raising PermissionError because the handle was still open.
    shutil.rmtree(stage / "out", ignore_errors=True)
    (stage / "in.fna").unlink(missing_ok=True)
    staged = stage / "in.fna"
    shutil.copyfile(fasta, staged)
    env = dict(os.environ, MSYS_NO_PATHCONV="1")
    cmd = ["docker", "run", "--rm", "-v", f"{stage.as_posix()}:/data", "-w", "/data", image,
           "SeqSero2_package.py", "-m", "k", "-t", "4", "-i", "in.fna", "-d", "out"]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=timeout, env=env)
    except Exception as e:                                   # noqa: BLE001
        return {"error": type(e).__name__}
    tsv = stage / "out" / "SeqSero_result.tsv"
    if not tsv.exists():
        return {"error": f"no_result(exit={p.returncode})"}
    rows = list(csv.DictReader(tsv.read_text(encoding="utf-8", errors="replace").splitlines(),
                               delimiter="\t"))
    if not rows:
        return {"error": "empty_result"}
    r = rows[0]
    return {"serotype": (r.get("Predicted serotype") or "").strip(),
            "profile": (r.get("Predicted antigenic profile") or "").strip(),
            "O": (r.get("O antigen prediction") or "").strip(),
            "H1": (r.get("H1 antigen prediction(fliC)") or "").strip(),
            "H2": (r.get("H2 antigen prediction(fljB)") or "").strip()}


def score(pred: str | None, truth: str, idx) -> str:
    """hit / miss / no_call — the SAME judgement for every caller."""
    if not pred or pred.strip().lower() in NON_SPECIFIC:
        return "no_call"
    ok, _why = equivalent(pred, truth, idx)
    return "hit" if ok else "miss"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", type=Path,
                    default=ROOT / "wiki" / "salmserovar_cohort_2026-09-04.json")
    ap.add_argument("--asm-root", type=Path, default=Path("D:/dna_decode_cache/salm_asm"))
    ap.add_argument("--db-dir", type=Path, default=SEROVAR_DB)
    ap.add_argument("--image", default=IMAGE)
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/ss2_checkpoint/rows.jsonl"))
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"salmserovar_seqsero2_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    cohort = json.loads(a.cohort.read_text(encoding="utf-8"))
    isolates = cohort["isolates"]
    if a.limit:
        isolates = isolates[:a.limit]
    idx = load_formula_index(a.db_dir / "serovar_table.tsv")

    a.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for ln in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                r = json.loads(ln)
                done[r["asm_acc"]] = r
        print(f"resuming: {len(done)} already scored")
    fh = open(a.checkpoint, "a", encoding="utf-8")

    print(f"{len(isolates)} isolates; SeqSero2 image {a.image}\n")
    for n, iso in enumerate(isolates, 1):
        acc = iso["asm_acc"]
        if acc in done:
            continue
        fa = next((p for p in (a.asm_root / acc).glob("*.fna")), None) if (a.asm_root / acc).is_dir() \
            else None
        rec = {"asm_acc": acc, "truth": iso["serovar_label"], "ncbi": iso.get("computed_serotype")}
        if fa is None:
            rec["status"] = "assembly_missing"
        else:
            rec["status"] = "ok"
            try:
                ours = call_serovar(fa, a.db_dir, blastn_bin=a.blastn)
                rec["ours"] = ours.get("serovar")
                rec["ours_formula"] = ours.get("antigenic_formula")
            except Exception as e:                            # noqa: BLE001
                rec["ours"] = None
                rec["ours_error"] = type(e).__name__
            rec["seqsero2"] = run_seqsero2(fa, STAGE, a.image)
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        done[acc] = rec
        if n % 20 == 0:
            print(f"  [{n}/{len(isolates)}] {acc} ours={rec.get('ours')} "
                  f"ss2={(rec.get('seqsero2') or {}).get('serotype')}", flush=True)
    fh.close()

    rows = [done[i["asm_acc"]] for i in isolates if i["asm_acc"] in done]
    scored = [r for r in rows if r.get("status") == "ok"]
    missing = [r for r in rows if r.get("status") == "assembly_missing"]
    ss2_err = [r for r in scored if (r.get("seqsero2") or {}).get("error")]
    # A caller that errored on a genome cannot be scored on it; comparing the others there would be
    # scoring three callers on different denominators.
    comparable = [r for r in scored if not (r.get("seqsero2") or {}).get("error")]

    if not comparable:
        print(f"\nREFUSING: {len(rows)} rows, {len(scored)} with an assembly, {len(comparable)} "
              "comparable. Nothing was scored.", file=sys.stderr)
        return 3

    def tally(getter) -> dict:
        c = {"hit": 0, "miss": 0, "no_call": 0}
        for r in comparable:
            c[score(getter(r), r["truth"], idx)] += 1
        n = len(comparable)
        return {**c, "n": n, "accuracy": c["hit"] / n,
                "accuracy_on_called": c["hit"] / (c["hit"] + c["miss"]) if (c["hit"] + c["miss"]) else None,
                "abstention_rate": c["no_call"] / n}

    ours = tally(lambda r: r.get("ours"))
    ss2 = tally(lambda r: (r.get("seqsero2") or {}).get("serotype"))
    ncbi = tally(lambda r: r.get("ncbi"))

    print(f"\n=== {len(comparable)} comparable isolates "
          f"({len(missing)} assembly-missing, {len(ss2_err)} SeqSero2 errors, excluded) ===")
    for name, t in (("ours", ours), ("SeqSero2 1.3.2", ss2), ("NCBI computed_types", ncbi)):
        print(f"  {name:22s} acc {t['accuracy']:.4f}  "
              f"(hit {t['hit']} / miss {t['miss']} / no_call {t['no_call']}; "
              f"abstain {t['abstention_rate']:.1%})")

    d_tool = ours["accuracy"] - ss2["accuracy"]
    d_ncbi = ours["accuracy"] - ncbi["accuracy"]
    print(f"\n  DELTA vs the REFERENCE TOOL (SeqSero2): {d_tool:+.4f}   <- the number this run adds")
    print(f"  delta vs NCBI computed_types          : {d_ncbi:+.4f}   (the 2026-09-04 comparator)")

    verdict = ("OURS_TRAILS_THE_REFERENCE_TOOL" if d_tool < -0.02 else
               "OURS_MATCHES_THE_REFERENCE_TOOL" if abs(d_tool) <= 0.02 else
               "OURS_LEADS_THE_REFERENCE_TOOL")
    why = (f"against SeqSero2 1.3.2 -- a pinned, versioned reference implementation rather than an "
           f"undocumented production field -- our caller scores {ours['accuracy']:.4f} vs "
           f"{ss2['accuracy']:.4f} ({d_tool:+.4f}) on {len(comparable)} isolates with wet-lab "
           "agglutination labels. This resolves the ambiguity the previous artifact recorded: the "
           "earlier delta was measured against NCBI's `computed_types` and could not distinguish "
           "'worse than the reference method' from 'diverges from an undocumented pipeline'.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "salmserovar-seqsero2-v1", "date": _date.today().isoformat(),
           "closes": ("the 2026-09-04 limit that the comparator was NCBI-PD `computed_types`, a "
                      "production call of undocumented tool/version/config"),
           "reference_tool": {"name": "SeqSero2", "version": "1.3.2", "image": a.image,
                              "mode": "-m k -t 4 (k-mer, assembly input)"},
           "label": ("wet-lab agglutination serovar from the committed 2026-09-04 cohort; the gold "
                     "standard is an antisera reaction, so no in-silico caller produced it"),
           "fairness": ("all three callers scored through the SAME equivalence.equivalent (notation "
                        "normalisation + the committed Kauffmann-White formula table, no fuzzy "
                        "near-misses); abstentions counted separately from misses for each"),
           "n_cohort": len(rows), "n_assembly_missing": len(missing),
           "n_seqsero2_errors": len(ss2_err), "n_comparable": len(comparable),
           "ours": ours, "seqsero2": ss2, "ncbi_computed_types": ncbi,
           "delta_vs_reference_tool": d_tool, "delta_vs_ncbi_computed_types": d_ncbi,
           "verdict": verdict, "why": why,
           "honest_limits": [
               "SeqSero2 is a TOOL, not the wet-lab assay. It is scored against the same agglutination "
               "labels as our caller, so this is a like-for-like comparison of two in-silico callers "
               "against a shared external truth -- not a claim that SeqSero2 is correct.",
               "Per-serovar cap 12 in the cohort flattens prevalence, so these are per-isolate "
               "accuracies on a deliberately diverse mix, NOT population-weighted rates.",
               "Isolates whose assembly is absent, or where SeqSero2 errored, are EXCLUDED so all "
               "three callers share one denominator; the excluded counts are reported.",
               "Residual label circularity is bounded, not eliminated: per-isolate agglutination "
               "provenance is unprovable from metadata. Contamination would inflate every caller, so "
               "the DELTA is more trustworthy than the absolute levels.",
               "SeqSero2 was run in k-mer mode on assemblies (-m k -t 4). Its read-based or "
               "microassembly modes could score differently.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
