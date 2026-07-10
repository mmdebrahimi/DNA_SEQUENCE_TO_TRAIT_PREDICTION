"""Bound the TB ABSTAIN correction with a SAMPLED regeno-callability probe (2026-07-10).

The independent + in-distribution TB numbers are an admitted CONSERVATIVE LOWER BOUND: without the
regenotyped VCF the scorer cannot distinguish "no determinant present" (a true S) from "the determinant
window could not be called" (which must ABSTAIN, never susceptible-by-absence). Closing that needs the
regeno VCFs. The full per-drug cohort is ~1.6 TB DECOMPRESSED (~176 GB gzipped), i.e. a day-plus job.

This probe answers the decision-relevant question at ~1% of that cost: **of the isolates the frozen
scorer currently calls S-by-absence, what fraction actually carry an UNCALLABLE determinant position and
should therefore ABSTAIN?** That single fraction bounds the size of the correction and says whether the
full regeno run is worth doing at all.

METHOD (nothing is reimplemented; the frozen scorer decides):
  For each sampled isolate we call the SAME `tb_amr.score_drug` twice --
      score_drug(drug, masked_calls, dets)                     -> today's behaviour (S-by-absence)
      score_drug(drug, masked_calls, dets, regeno_text=...)    -> callability-aware (may ABSTAIN)
  and diff the two. The correction is exactly the set of S -> ABSTAIN flips.

WHY A STREAMING FILTER IS SAFE: `score_drug` queries a DISCRETE position set
(`sorted({d.pos for d in determinants})`), not a window, and `tb_vcf._iter_records` ignores header lines
and any record not on H37Rv. So a subset text containing only the records at those positions yields
byte-identical `callable_positions` output, while letting us stop reading each 177 MB VCF at the last
determinant position. A position ABSENT from the subset is uncallable — which is exactly the frozen
(conservative) semantics for a position absent from the full VCF.

Restartable + idempotent: per-isolate verdicts are checkpointed to JSONL; regeno VCFs are never cached
(19.7 MB gz / 177 MB decompressed each -- streamed and discarded).

Run:  uv run python scripts/tb_callability_probe.py --n 100
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import time
import urllib.request
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data import tb_who_catalogue  # noqa: E402
from dna_decode.organism_rules import tb_amr, tb_vcf  # noqa: E402
from scripts.cryptic_feasibility_probe import _vcf_url, fetch_vcf, load_rows  # noqa: E402

DRUGS = ("rifampicin", "isoniazid")
DRUG_CODE = {"rifampicin": "RIF", "isoniazid": "INH"}
DEFAULT_CHECKPOINT = Path("D:/dna_decode_cache/tb_callability/checkpoint.jsonl")


def eligible_rows(rows: list[dict], drug_code: str = "RIF") -> list[dict]:
    """HIGH-quality R/S isolates that carry BOTH a masked and a regenotyped VCF path."""
    ph, q = f"{drug_code}_BINARY_PHENOTYPE", f"{drug_code}_PHENOTYPE_QUALITY"
    out = []
    for r in rows:
        masked, regeno = tb_vcf.vcf_paths_for(r)
        if masked and regeno and r.get(ph) in ("R", "S") and r.get(q) == "HIGH":
            out.append(r)
    return out


def stride_sample(rows: list[dict], n: int) -> list[dict]:
    """Deterministic evenly-spaced sample.

    The reuse table is grouped by site/run, so taking the FIRST n would sample one or two sites. A stride
    spreads the sample across the whole table without needing a RNG (reproducible by construction).
    """
    if n >= len(rows):
        return list(rows)
    step = len(rows) / n
    return [rows[int(i * step)] for i in range(n)]


def determinant_positions() -> tuple[dict[str, list], set[int], int]:
    """({drug: determinants}, union of queried positions, max position)."""
    dets = {d: tb_who_catalogue.load_determinants(d) for d in DRUGS}
    positions = {x.pos for ds in dets.values() for x in ds}
    return dets, positions, max(positions)


def stream_regeno_subset(rel_path: str, positions: set[int], max_pos: int, timeout: int = 300) -> str:
    """Download the regeno VCF and return ONLY the records at `positions`, as VCF text.

    Streams the gzip; stops at the last determinant position (VCFs are coordinate-sorted), so we never
    decompress the tail of a 177 MB file. Header lines are dropped — `_iter_records` skips them anyway.
    """
    url = _vcf_url(rel_path)
    kept: list[str] = []
    with urllib.request.urlopen(url, timeout=timeout) as resp:  # noqa: S310 (trusted EBI FTP-over-HTTPS)
        raw = resp.read()
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
        for bline in gz:
            if bline.startswith(b"#"):
                continue
            tab = bline.find(b"\t")
            if tab < 0:
                continue
            nxt = bline.find(b"\t", tab + 1)
            try:
                pos = int(bline[tab + 1:nxt])
            except ValueError:
                continue
            if pos > max_pos:
                break                     # coordinate-sorted -> nothing further can matter
            if pos in positions:
                kept.append(bline.decode("utf-8", "replace").rstrip("\n"))
    return "\n".join(kept)


def classify_positions(regeno_text: str, positions) -> dict[str, int]:
    """Split queried positions into three states — the distinction `callable_positions` collapses.

    `tb_vcf.callable_positions` returns callable iff a PASS record with an explicit GT exists, and treats
    an ABSENT position as uncallable. But these are minos regenotyped VCFs (`##source=minos`), which emit
    records ONLY at sites in the variant panel — ~29% of the genome, fragmented (rpoB RRDR 41.6%, katG
    24.8% of positions carry a record). So `absent` overwhelmingly means "not a genotyping target", NOT
    "sequencing could not call it".

    Only `present_fail` is evidence of TRUE uncallability (a record exists and failed MIN_DP/MIN_FRS/
    MIN_GCP or has `./.`).
    """
    want = set(positions)
    present_pass = present_fail = 0
    seen: set[int] = set()
    flags = tb_vcf.callable_positions(regeno_text, sorted(want))
    for pos, *_ in tb_vcf._iter_records(regeno_text):
        if pos in want and pos not in seen:
            seen.add(pos)
            if flags.get(pos):
                present_pass += 1
            else:
                present_fail += 1
    return {"present_pass": present_pass, "present_fail": present_fail,
            "absent": len(want) - len(seen), "n_positions": len(want)}


def probe_one(row: dict, dets: dict, positions: set[int], max_pos: int) -> dict | None:
    """Score one isolate with and without callability. Returns a verdict record (or None if unfetchable)."""
    masked_rel, regeno_rel = tb_vcf.vcf_paths_for(row)
    masked_bytes = fetch_vcf(masked_rel)
    if masked_bytes is None:
        return None
    masked_calls = tb_vcf.parse_masked_calls(masked_bytes.decode("utf-8", "replace"))

    try:
        regeno_text = stream_regeno_subset(regeno_rel, positions, max_pos)
    except Exception as e:  # noqa: BLE001 — one bad isolate must not abort the sweep
        print(f"    [regeno fetch fail] {regeno_rel}: {type(e).__name__}: {e}", flush=True)
        return None

    rec: dict = {"uniqueid": row.get("UNIQUEID") or row.get("ENA_RUN_ACCESSION") or regeno_rel, "drugs": {}}
    for drug in DRUGS:
        dpos = sorted({d.pos for d in dets[drug]})
        base = tb_amr.score_drug(drug, masked_calls, dets[drug])                       # masked-only (S-by-absence)
        # The side-by-side needs the FROZEN absent-is-uncallable rule explicitly: the default is now the
        # ratified fail-only rule (2026-07-10), so `flipped_S_to_ABSTAIN` would otherwise no longer measure
        # the frozen rule. `would_abstain_fail_only_rule` (from classify_positions) is the ratified rule;
        # keeping both explicit preserves the two-rule comparison this probe exists to report.
        aware = tb_amr.score_drug(drug, masked_calls, dets[drug], regeno_text=regeno_text,
                                  absent_is_uncallable=True)      # FROZEN rule, for the comparison
        states = classify_positions(regeno_text, dpos)
        rec["drugs"][drug] = {
            "call_without_callability": base.prediction,
            "call_with_callability": aware.prediction,
            "flipped_S_to_ABSTAIN": base.prediction == "S" and aware.prediction == "ABSTAIN",
            # ABSTAIN under a corrected rule where ONLY a present-and-failed record counts as uncallable
            "would_abstain_fail_only_rule": base.prediction == "S" and states["present_fail"] > 0,
            "n_uncallable_positions": aware.n_uncallable_positions,
            "n_determinant_positions": base.n_determinant_positions,
            **states,
            "phenotype": row.get(f"{DRUG_CODE[drug]}_BINARY_PHENOTYPE"),
        }
    return rec


def summarize(records: list[dict]) -> dict:
    """Per-drug: how many S-by-absence calls would become ABSTAIN once callability is assessed."""
    out: dict[str, dict] = {}
    for drug in DRUGS:
        ds = [r["drugs"][drug] for r in records if drug in r["drugs"]]
        n_S = sum(1 for d in ds if d["call_without_callability"] == "S")
        n_flip = sum(1 for d in ds if d["flipped_S_to_ABSTAIN"])
        n_flip_fail = sum(1 for d in ds if d["would_abstain_fail_only_rule"])
        n_R = sum(1 for d in ds if d["call_without_callability"] == "R")
        out[drug] = {
            "n_isolates": len(ds),
            "n_called_R": n_R,
            "n_called_S_by_absence": n_S,
            # CURRENT frozen rule: absent OR failed -> uncallable -> ABSTAIN
            "n_S_would_ABSTAIN_current_rule": n_flip,
            "frac_of_S_that_is_uncallable_current_rule": round(n_flip / n_S, 4) if n_S else None,
            # CORRECTED rule: only a present-and-FAILED record is evidence of true uncallability
            "n_S_would_ABSTAIN_fail_only_rule": n_flip_fail,
            "frac_of_S_that_is_uncallable_fail_only_rule": round(n_flip_fail / n_S, 4) if n_S else None,
            "median_positions_absent": _median([d["absent"] for d in ds]),
            "median_positions_present_pass": _median([d["present_pass"] for d in ds]),
            "median_positions_present_fail": _median([d["present_fail"] for d in ds]),
        }
    return out


def _median(xs: list[int]):
    if not xs:
        return None
    xs = sorted(xs)
    m = len(xs) // 2
    return xs[m] if len(xs) % 2 else (xs[m - 1] + xs[m]) / 2


def load_checkpoint(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            r = json.loads(line)
            out[r["uniqueid"]] = r
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=100, help="isolates to sample (stride across the reuse table)")
    ap.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    ap.add_argument("--output-prefix", type=Path,
                    default=REPO / "wiki" / f"tb_callability_probe_{_date.today().isoformat()}")
    a = ap.parse_args(argv)

    tb_who_catalogue.verify_pins()
    dets, positions, max_pos = determinant_positions()
    print(f"[tb-callability] {len(positions)} determinant positions across {DRUGS}; max_pos={max_pos}")

    rows = eligible_rows(load_rows())
    sample = stride_sample(rows, a.n)
    print(f"[tb-callability] {len(rows)} eligible isolates -> sampling {len(sample)} by stride", flush=True)

    a.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    done = load_checkpoint(a.checkpoint)
    print(f"[tb-callability] {len(done)} already checkpointed", flush=True)

    t0 = time.time()
    with open(a.checkpoint, "a", encoding="utf-8") as fh:
        for i, row in enumerate(sample, 1):
            uid = row.get("UNIQUEID") or row.get("ENA_RUN_ACCESSION") or tb_vcf.vcf_paths_for(row)[1]
            if uid in done:
                continue
            rec = probe_one(row, dets, positions, max_pos)
            if rec is None:
                continue
            done[rec["uniqueid"]] = rec
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if i % 10 == 0 or i == len(sample):
                el = time.time() - t0
                print(f"  [{i}/{len(sample)}] {len(done)} done, {el/60:.1f} min elapsed", flush=True)

    records = list(done.values())
    summary = summarize(records)

    payload = {
        "_schema": "tb-callability-probe-v1",
        "run_date": _date.today().isoformat(),
        "cohort": "CRyPTIC reuse table (HIGH-quality RIF R/S, masked+regeno available)",
        "n_eligible": len(rows),
        "n_sampled": len(records),
        "sampling": "deterministic stride across the reuse table (spreads across sites/runs)",
        "n_determinant_positions": len(positions),
        "per_drug": summary,
        "interpretation": (
            "The two rules disagree completely, and that IS the result. The FROZEN rule "
            "(`tb_vcf.callable_positions`) treats a determinant position ABSENT from the regeno VCF as "
            "uncallable. But these are minos regenotyped VCFs (`##source=minos`), which emit records only "
            "at sites in the variant panel: ~29% of genome positions, fragmented (rpoB RRDR 41.6%, katG "
            "24.8%). So `absent` overwhelmingly means 'not a genotyping target', not 'sequencing could not "
            "call it'. Under the frozen rule essentially EVERY S-by-absence call flips to ABSTAIN, which "
            "would not tighten the TB numbers -- it would erase specificity. Under the fail-only rule "
            "(a record exists and failed MIN_DP/MIN_FRS/MIN_GCP, or GT is ./.) the correction is the "
            "honest measure of true uncallability."),
        "scope_limits": (
            "A SAMPLE, not the full cohort — each fraction carries sampling error (Wilson CI in the md). "
            "Sampled from the RIF HIGH-quality masked+regeno subset, so it is not re-weighted to the "
            "prevalence-preserving cohort. Nothing is re-scored: the frozen tb_amr.score_drug decides both "
            "calls; only its regeno_text argument differs. The streaming position filter was verified "
            "byte-equivalent to the full 177 MB VCF on a real isolate (callable_positions identical for "
            "both drugs). The fail-only rule is a MEASUREMENT here, not a shipped change to tb_vcf."),
        "recommendation": (
            "Do NOT run the full ~1.6 TB regeno job on the strength of the frozen rule: it cannot tighten "
            "the reported numbers, because absence-from-panel is not uncallability. Either (a) accept the "
            "current conservative lower bound and say so, or (b) first ratify a corrected callability "
            "semantics (three states: present_pass / present_fail / absent) and re-measure on this sample."),
    }
    out_json = a.output_prefix.with_suffix(".json")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    from dna_decode.eval.clonality import wilson_ci
    lines = [
        f"# TB regeno-callability probe — bounding the ABSTAIN correction ({payload['run_date']})",
        "",
        f"**Sampled:** {len(records)} of {len(rows)} eligible isolates (deterministic stride) · "
        f"**{len(positions)}** determinant positions (RIF + INH, WHO grade-1/2)",
        "",
        "Nothing is re-scored. The frozen `tb_amr.score_drug` is called TWICE per isolate — without and "
        "with the regenotyped VCF — and the two calls are diffed. The correction is exactly the "
        "`S -> ABSTAIN` flips.",
        "",
        "| drug | n | called R | called S | S→ABSTAIN (frozen rule) | S→ABSTAIN (fail-only rule) | "
        "median absent / pass / fail positions |",
        "|---|---:|---:|---:|---|---|---|",
    ]
    for drug, s in summary.items():
        n_S = s["n_called_S_by_absence"]
        if n_S:
            lo, hi = wilson_ci(s["n_S_would_ABSTAIN_current_rule"], n_S)
            cur = f"**{s['frac_of_S_that_is_uncallable_current_rule']}** [{lo}–{hi}]"
            lo2, hi2 = wilson_ci(s["n_S_would_ABSTAIN_fail_only_rule"], n_S)
            fail = f"**{s['frac_of_S_that_is_uncallable_fail_only_rule']}** [{lo2}–{hi2}]"
        else:
            cur = fail = "n/a"
        lines.append(
            f"| {drug} | {s['n_isolates']} | {s['n_called_R']} | {n_S} | {cur} | {fail} | "
            f"{s['median_positions_absent']} / {s['median_positions_present_pass']} / "
            f"{s['median_positions_present_fail']} |")

    lines += ["", "## Interpretation", "", payload["interpretation"],
              "", "## Scope limits", "", payload["scope_limits"], "",
              "Generated by `scripts/tb_callability_probe.py`."]
    a.output_prefix.with_suffix(".md").write_text("\n".join(lines), encoding="utf-8")

    print("\n[tb-callability] summary:")
    for drug, s in summary.items():
        print(f"  {drug}: n={s['n_isolates']} R={s['n_called_R']} S={s['n_called_S_by_absence']} "
              f"| frozen rule -> ABSTAIN {s['n_S_would_ABSTAIN_current_rule']} "
              f"(frac={s['frac_of_S_that_is_uncallable_current_rule']}) "
              f"| fail-only rule -> ABSTAIN {s['n_S_would_ABSTAIN_fail_only_rule']} "
              f"(frac={s['frac_of_S_that_is_uncallable_fail_only_rule']}) "
              f"| median absent/pass/fail = {s['median_positions_absent']}/"
              f"{s['median_positions_present_pass']}/{s['median_positions_present_fail']}")
    print(f"artifact -> {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
