"""Is the O-antigen the systematic weak axis of the typing suite? Test it on a SECOND cell.

THE HYPOTHESIS, and where it came from. The Salmonella serovar cell was measured against a wet-lab
label on 2026-09-04 and scored 0.702 against an in-silico incumbent's 0.925. Its failure was NOT
uniform: 33 of 59 no-calls had an undetected phase-2 flagellin, and where a formula was produced the
H antigen was frequently right while the O antigen was unresolved (`O?`) or mis-grouped (Typhi called
1,3,19 rather than 9,12). That suggests something narrower and more useful than "the caller is weak":
**the O-antigen axis may be the systematic weakness, shared across serotyping cells.**

E. coli O:H serotype is the natural second test. It is an INDEPENDENT cell with its own allele DB, its
gold standard is likewise slide agglutination, and -- crucially -- its label decomposes into exactly the
two axes under test, so O and H accuracy can be scored SEPARATELY on the same isolates. A hypothesis
that reproduces on a second cell with a different database is a property of the approach; one that does
not is a property of the first caller.

WHAT IS DIFFERENT HERE, AND IT COSTS US SOMETHING. NCBI-PD carries `computed_types` for Salmonella but
NOT for E. coli -- the field is NULL. So there is **no in-silico incumbent** to score beside us, and the
absolute accuracy here is therefore LESS interpretable than the Salmonella number was. That is stated
rather than papered over. What survives without a comparator is the WITHIN-CELL contrast: O accuracy vs
H accuracy is measured on the same isolates, by the same caller, against the same labels, so the
comparison between axes is internally controlled even though its level is not externally anchored.

LABEL HYGIENE IS MOST OF THE WORK. The E. coli `serovar` field is a mixture: real serotypes
('O157:H7'), species names ('Escherichia coli'), partials ('O157'), Shigella serotypes on a different
scheme ('2A', '3a'), and placeholders ('E. coli O26:Pending'). Only a strict O:H shape is admitted.

Network + blastn. Writes wiki/serotype_oh_validation_<date>.json.
"""
from __future__ import annotations

import argparse
import collections
import json
import random
import re
import sys
import traceback
import urllib.request
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.data.refseq import download_genome, fasta_path  # noqa: E402
from dna_decode.serotype.runner import call_serotype  # noqa: E402
from source_diverse_validate import MAX_SOURCE_SHARE  # noqa: E402
from gentamicin_rmt_specificity_hunt import latest_metadata_url  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"

# A real E. coli serotype label: O-group then H-type. H may be '-' / 'NM' (non-motile) / 'H-' -- those
# are MEANINGFUL (flagellar antigen absent), not missing, and are kept.
_OH = re.compile(r"^\s*O(\d+)\s*:\s*(H\d+|H-|HNM|NM|H\?)\s*$", re.IGNORECASE)
# O-only labels ('O157') are admitted for the O axis alone -- discarding them would throw away the
# largest block of genuine agglutination labels, but they can never score the H axis.
_O_ONLY = re.compile(r"^\s*O(\d+)\s*$", re.IGNORECASE)

PLACEHOLDER = ("pending", "unknown", "not typed", "untypeable", "rough", "ont", "onut")


def parse_label(raw: str) -> tuple[str | None, str | None]:
    """-> (O, H) as normalised tokens; either may be None. Junk yields (None, None)."""
    t = (raw or "").strip().strip('"')
    if not t or any(p in t.lower() for p in PLACEHOLDER):
        return None, None
    m = _OH.match(t)
    if m:
        h = m.group(2).upper()
        h = "H-" if h in ("H-", "HNM", "NM") else h
        return f"O{int(m.group(1))}", (None if h == "H?" else h)
    m = _O_ONLY.match(t)
    if m:
        return f"O{int(m.group(1))}", None
    return None, None


def split_serotype(s: str | None) -> tuple[str | None, str | None]:
    """'O157:H7' -> ('O157','H7'); 'O?:H7' -> (None,'H7'). The caller emits ONE combined string."""
    if not s or ":" not in str(s):
        return None, None
    o, _, h = str(s).partition(":")
    o, h = o.strip(), h.strip()
    return (None if o in ("", "O?", "?") else o), (None if h in ("", "H?", "?") else h)


def norm_call(v: str | None, axis: str) -> str | None:
    """Normalise a caller's axis output to the label vocabulary; unresolved -> None."""
    if not v:
        return None
    t = str(v).strip().upper()
    if t in ("", "-", "?", "O?", "H?", "NONE", "NA", "UNKNOWN"):
        return None
    if axis == "O":
        m = re.search(r"O(\d+)", t)
        return f"O{int(m.group(1))}" if m else None
    if t in ("H-", "HNM", "NM"):
        return "H-"
    m = re.search(r"H(\d+)", t)
    return f"H{int(m.group(1))}" if m else None


def build_cohort(max_rows: int, target: int, seed: int) -> tuple[list[dict], dict]:
    url = latest_metadata_url("Escherichia_coli_Shigella")
    print(f"streaming {url}")
    r = urllib.request.urlopen(url, timeout=900)
    cols = r.readline().decode("utf8", "replace").rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(cols)}
    funnel = collections.Counter()
    pool: list[dict] = []
    for line in r:
        f = line.decode("utf8", "replace").rstrip("\n").split("\t")
        if len(f) < len(cols):
            continue
        funnel["rows"] += 1
        if funnel["rows"] > max_rows:
            break
        o, h = parse_label(f[ix["serovar"]])
        if not o:
            continue
        funnel["has_O_label"] += 1
        if h:
            funnel["has_OH_label"] += 1
        asm = f[ix["asm_acc"]].strip()
        if not asm or asm.upper() in ("NULL", "NA"):
            continue
        funnel["has_assembly"] += 1
        pool.append({"asm_acc": asm, "O": o, "H": h,
                     "raw": f[ix["serovar"]].strip(),
                     "bioproject": f[ix["bioproject_acc"]].strip() or "NO_BP"})

    rng = random.Random(seed)
    rng.shuffle(pool)
    # Cap per O-group AND prefer full O:H labels, so the H axis is actually scoreable and the cohort is
    # not 90% O157 (the most-sequenced serotype on earth).
    by_o: collections.Counter = collections.Counter()
    by_bp: collections.Counter = collections.Counter()
    cap_o, cap_bp = max(2, target // 10), max(2, target // 6)
    cohort: list[dict] = []
    for prefer_h in (True, False):
        for e in pool:
            if len(cohort) >= target:
                break
            if prefer_h and not e["H"]:
                continue
            if by_o[e["O"]] >= cap_o or by_bp[e["bioproject"]] >= cap_bp:
                continue
            if any(c["asm_acc"] == e["asm_acc"] for c in cohort):
                continue
            cohort.append(e); by_o[e["O"]] += 1; by_bp[e["bioproject"]] += 1
    src = collections.Counter(c["bioproject"] for c in cohort)
    top, ntop = src.most_common(1)[0] if src else ("NO_BP", 0)
    meta = {"source": url, "funnel": dict(funnel), "n": len(cohort),
            "n_distinct_O": len(by_o), "n_with_H_label": sum(1 for c in cohort if c["H"]),
            "n_bioprojects": len(src), "largest_source_share": (ntop / len(cohort)) if cohort else None,
            "passes_source_diversity_bar": bool(cohort) and (ntop / len(cohort)) <= MAX_SOURCE_SHARE}
    return cohort, meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-rows", type=int, default=250000)
    ap.add_argument("--target", type=int, default=120)
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--db", type=Path, default=ROOT / "data" / "serotypefinder_db" / "serotypefinder.fsa")
    ap.add_argument("--asm-dir", type=Path, default=Path("D:/dna_decode_cache/ecoli_sero_asm"))
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/ecoli_sero_asm/results.jsonl"))
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"serotype_oh_validation_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    cohort, meta = build_cohort(a.max_rows, a.target, a.seed)
    print(f"\nfunnel: {meta['funnel']}")
    print(f"cohort: {meta['n']} | {meta['n_distinct_O']} O-groups | {meta['n_with_H_label']} with an H "
          f"label | {meta['n_bioprojects']} BioProjects | largest {meta['largest_source_share']:.3f} "
          f"{'PASSES' if meta['passes_source_diversity_bar'] else 'FAILS'}")
    if not cohort:
        print("empty cohort", file=sys.stderr)
        return 2

    a.asm_dir.mkdir(parents=True, exist_ok=True)
    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for line in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                rec = json.loads(line)
                done[rec["asm_acc"]] = rec

    fh = open(a.checkpoint, "a", encoding="utf-8")
    for n, iso in enumerate(cohort, 1):
        acc = iso["asm_acc"]
        if acc in done:
            continue
        rec = {"asm_acc": acc, "O_label": iso["O"], "H_label": iso["H"], "raw": iso["raw"]}
        try:
            download_genome(acc, a.asm_dir)
            fa = fasta_path(acc, a.asm_dir)
            if not Path(fa).exists():
                rec["status"] = "assembly_missing"
            else:
                call = call_serotype(fa, a.db, blastn_bin=a.blastn)
                rec["status"] = "ok"
                rec["serotype"] = call.get("serotype")
                # The caller returns ONE combined 'O?:H7' string, NOT separate O_type/H_type keys.
                # Reading keys that do not exist yielded None for every axis and made a broken run look
                # like a clean 0.0-vs-0.0 finding -- so split the string the caller actually emits.
                o_c, h_c = split_serotype(rec["serotype"])
                rec["O_call"] = call.get("O_type") or call.get("O") or o_c
                rec["H_call"] = call.get("H_type") or call.get("H") or h_c
        except Exception as e:                       # noqa: BLE001 - recorded, never hidden
            rec["status"] = f"error:{type(e).__name__}"
            rec["error"] = str(e)[:200]
            rec["trace"] = traceback.format_exc()[-300:]
        fh.write(json.dumps(rec) + "\n"); fh.flush()
        done[acc] = rec
        if n % 10 == 0 or rec["status"] != "ok":
            print(f"  [{n}/{len(cohort)}] {acc} {rec['status']} "
                  f"O:{rec.get('O_call')}/{rec['O_label']} H:{rec.get('H_call')}/{rec['H_label']}",
                  flush=True)
    fh.close()

    rows = [done[c["asm_acc"]] for c in cohort if c["asm_acc"] in done]
    ok = [r for r in rows if r.get("status") == "ok"]
    axis: dict[str, collections.Counter] = {"O": collections.Counter(), "H": collections.Counter()}
    misses = []
    for r in ok:
        for ax in ("O", "H"):
            lab = r.get(f"{ax}_label")
            if not lab:
                continue                                    # unscoreable on this axis, not a miss
            call = norm_call(r.get(f"{ax}_call"), ax)
            if call is None:
                axis[ax]["no_call"] += 1
            elif call == lab:
                axis[ax]["hit"] += 1
            else:
                axis[ax]["miss"] += 1
                if len(misses) < 25:
                    misses.append({"asm": r["asm_acc"], "axis": ax, "label": lab, "call": call})

    def acc_of(c):
        s = c["hit"] + c["miss"]
        return (c["hit"] / s) if s else None

    def resolved_of(c):
        t = c["hit"] + c["miss"] + c["no_call"]
        return ((t - c["no_call"]) / t) if t else None

    print(f"\n=== {len(ok)} isolates called ===")
    for ax in ("O", "H"):
        c = axis[ax]
        print(f"  {ax}-antigen: hit={c['hit']:<4} miss={c['miss']:<4} no_call={c['no_call']:<4} "
              f"accuracy={acc_of(c) if acc_of(c) is None else round(acc_of(c),3)}  "
              f"resolved={resolved_of(c) if resolved_of(c) is None else round(resolved_of(c),3)}")

    o_res, h_res = resolved_of(axis["O"]), resolved_of(axis["H"])
    o_acc, h_acc = acc_of(axis["O"]), acc_of(axis["H"])
    gap_res = None if (o_res is None or h_res is None) else h_res - o_res
    gap_acc = None if (o_acc is None or h_acc is None) else h_acc - o_acc

    # PLUMBING GUARD. Zero resolution on BOTH axes is the signature of a broken extraction, not of a
    # caller that resolves nothing -- and the first run of this script printed a confident
    # "hypothesis not supported" from exactly that state. A verdict computed from nothing is worse
    # than no verdict, so refuse rather than report.
    if (o_res == 0.0 and h_res == 0.0) or (o_res is None and h_res is None):
        print("\nREFUSING to emit a verdict: BOTH axes resolved 0 calls. That is a plumbing "
              "signature (extraction reading keys the caller does not emit), not a result.",
              file=sys.stderr)
        return 3

    if gap_res is not None and gap_res > 0.10:
        verdict = "O_AXIS_IS_THE_WEAK_ONE_REPRODUCES"
        why = (f"the O antigen resolves on {o_res:.1%} of scoreable isolates against the H antigen's "
               f"{h_res:.1%} (gap {gap_res:+.3f}) -- the same asymmetry seen in the Salmonella cell, on "
               "an independent cell with a different allele database. The weakness travels with the "
               "O-antigen axis, not with one caller.")
    elif gap_res is not None and gap_res < -0.10:
        verdict = "H_AXIS_IS_THE_WEAK_ONE_HYPOTHESIS_INVERTED"
        why = (f"H resolves WORSE than O here ({h_res:.1%} vs {o_res:.1%}) -- the Salmonella asymmetry "
               "does NOT generalise, and is a property of that caller rather than of the O axis.")
    else:
        verdict = "NO_AXIS_ASYMMETRY_HYPOTHESIS_NOT_SUPPORTED"
        why = (f"O and H resolve comparably ({o_res} vs {h_res}) -- the Salmonella O-weakness does not "
               "reproduce here, so it is specific to that cell rather than systematic.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {
        "schema": "serotype-oh-validation-v1",
        "date": _date.today().isoformat(),
        "cell": "typing:Escherichia_coli:serotype",
        "hypothesis": "the O-antigen axis is the systematic weak point of the serotyping cells, as "
                      "suggested by the Salmonella serovar result (H often right where O failed)",
        "cohort": meta,
        "n_called": len(ok),
        "statuses": dict(collections.Counter(r.get("status", "?") for r in rows)),
        "O_axis": dict(axis["O"]), "H_axis": dict(axis["H"]),
        "O_accuracy": o_acc, "H_accuracy": h_acc,
        "O_resolved_rate": o_res, "H_resolved_rate": h_res,
        "resolved_gap_H_minus_O": gap_res, "accuracy_gap_H_minus_O": gap_acc,
        "verdict": verdict, "why": why,
        "sample_misses": misses,
        "honest_limits": [
            "NO IN-SILICO INCUMBENT: NCBI-PD does not populate `computed_types` for E. coli, so unlike "
            "the Salmonella cell there is no comparator and the ABSOLUTE accuracy here is much less "
            "interpretable. The internally-controlled O-vs-H contrast is what this run supports.",
            "The label is PD's submitter `serovar` parsed to a strict O:H shape. E. coli O:H serotyping "
            "is traditionally agglutination-based, but per-isolate provenance is unprovable here and "
            "-- unlike the Salmonella run -- there is no incumbent score to use as a circularity probe.",
            "O-only labels are scored on the O axis alone; they cannot score H. The two axes therefore "
            "have DIFFERENT denominators, which is why resolved-rate and accuracy are reported per axis "
            "rather than pooled.",
            "'H-' (non-motile) is a MEANINGFUL label, not a missing value, and is scored as such.",
            "Per-O-group and per-BioProject caps deliberately flatten prevalence away from O157.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
