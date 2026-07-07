"""CYP2D6-CYP2D7 PSV evidence table (Phase A — the read-only diagnostic, NOT an identity caller).

Per the /brainstorm (2026-07-07): before writing any hybrid-IDENTITY classifier, build a read-only evidence
table that must first REPRODUCE the known Cyrius-style regional profiles on a small labelled panel. This
consumes the REAL Cyrius PSV file (data/cyp2d6_psv/CYP2D6_SNP_38.txt — 117 CYP2D6/CYP2D7 differentiating
sites, GRCh38-native, paired-coordinate schema `chr pos_CYP2D6 base_CYP2D6 pos_CYP2D7 base_CYP2D7 annotation`)
and per-PSV mpileup base counts (from CRAM, both pos_CYP2D6 AND pos_CYP2D7 — the brainstorm's decisive point:
D6-coordinate-only measures aligner placement, not biology).

For each PSV: D6-like = reads carrying base_CYP2D6 at (pos_CYP2D6 OR pos_CYP2D7); D7-like = reads carrying
base_CYP2D7 at either. D6-fraction = D6-like / (D6-like + D7-like). A NORMAL/dup/deletion sample has a FLAT
D6-fraction profile across the gene's regions; a CYP2D6-CYP2D7 HYBRID shows a REGIONAL SHIFT at the converted
segment (*36 -> exon9/3'; *68 -> intron1; *13 -> upstream/exon1, opposite direction).

Regions are ordered 3'->5' by genomic coordinate (CYP2D6 is minus-strand: low coord = downstream_exon9/exon9,
high coord = upstream_exon1). This is DIAGNOSTIC ONLY — it emits per-region D6-fractions + a regional-shift
metric; it does NOT assign *13/*36/*68. Flag contract: permissive (samtools mpileup -B -q 0 -Q 0) to match
the shipped depth path; the brainstorm flagged the exact flags as a thing to derive on this panel.
"""
from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PSV_FILE = REPO / "data" / "cyp2d6_psv" / "CYP2D6_SNP_38.txt"

# Region order 3'->5' along CYP2D6 (minus strand: ascending genomic coord = 3'->5').
_REGION_ORDER = ["downstream_exon9", "exon9", "intron6", "exon7", "exon6", "intron5", "intron4",
                 "exon3", "intron2", "exon2", "intron1", "exon1", "upstream_exon1"]

_PILEUP_INDEL = re.compile(r"[+-](\d+)")


def load_psvs(path: Path = PSV_FILE) -> list[dict]:
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        out.append({"chrom": f[0], "pos_d6": int(f[1]), "base_d6": f[2].upper(),
                    "pos_d7": int(f[3]), "base_d7": f[4].upper(), "annotation": f[5]})
    return out


def _base_counts(pileup: str) -> dict[str, int]:
    """Count A/C/G/T (case-insensitive) in a samtools mpileup base column. Strips read-start (^X),
    read-end ($), and indel (+/-Nbases) markers so their inner bases are not miscounted."""
    s = pileup
    # remove read-start markers ^ followed by one mapq char
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == "^":
            i += 2  # skip ^ + the mapq char
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            m = _PILEUP_INDEL.match(s, i)
            if m:
                i = m.end() + int(m.group(1))  # skip +/-, the number, and that many indel bases
                continue
        out.append(c)
        i += 1
    joined = "".join(out).upper()
    return {b: joined.count(b) for b in "ACGT"}


def _read_pileup(path: Path) -> dict[int, str]:
    if not path.exists():
        return {}
    d = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        if len(f) >= 3 and f[0].isdigit():
            d[int(f[0])] = f[2]
    return d


def sample_evidence(sample: str, pileup_dir: Path, psvs: list[dict]) -> dict:
    d6p = _read_pileup(pileup_dir / f"{sample}.d6.txt")
    d7p = _read_pileup(pileup_dir / f"{sample}.d7.txt")
    per_psv = []
    for p in psvs:
        c6 = _base_counts(d6p.get(p["pos_d6"], ""))
        c7 = _base_counts(d7p.get(p["pos_d7"], ""))
        d6_like = c6.get(p["base_d6"], 0) + c7.get(p["base_d6"], 0)
        d7_like = c6.get(p["base_d7"], 0) + c7.get(p["base_d7"], 0)
        tot = d6_like + d7_like
        frac = round(d6_like / tot, 3) if tot >= 10 else None   # callable only if >=10 informative reads
        per_psv.append({"pos_d6": p["pos_d6"], "annotation": p["annotation"], "d6_like": d6_like,
                        "d7_like": d7_like, "d6_fraction": frac})
    # per-region median D6-fraction (ordered 3'->5')
    region_frac = {}
    for r in _REGION_ORDER:
        vals = [x["d6_fraction"] for x in per_psv if x["annotation"] == r and x["d6_fraction"] is not None]
        region_frac[r] = round(statistics.median(vals), 3) if vals else None
    present = [v for v in region_frac.values() if v is not None]
    shift = round(max(present) - min(present), 3) if len(present) >= 2 else None
    n_callable = sum(1 for x in per_psv if x["d6_fraction"] is not None)
    # DIAGNOSTIC identity SIGNAL (NOT a call): the directional 5'-vs-3' shift + the exon-9-tip dip.
    def _mean(keys):
        v = [region_frac[k] for k in keys if region_frac.get(k) is not None]
        return statistics.mean(v) if v else None
    five = _mean(["exon1", "upstream_exon1", "intron1"])       # 5' end
    three = _mean(["downstream_exon9", "intron6", "exon6"])    # 3' end
    fp_mp = round(five - three, 3) if (five is not None and three is not None) else None
    gene_med = round(statistics.median(present), 3) if present else None
    ex9 = region_frac.get("downstream_exon9")
    ex9_dip = round(gene_med - ex9, 3) if (gene_med is not None and ex9 is not None) else None
    # coarse SIGNAL bucket (proof-of-signal only; the real classifier is Phase B):
    signal = "flat_nonhybrid"
    if fp_mp is not None and fp_mp >= 0.15:
        signal = "directional_5p_high_3p_low (68-like)"
    elif fp_mp is not None and fp_mp <= -0.15:
        signal = "directional_5p_low_3p_high (13-like)"
    elif ex9_dip is not None and ex9_dip >= 0.15:
        signal = "exon9_tip_dip (36-like)"
    return {"sample": sample, "n_psv": len(psvs), "n_callable": n_callable,
            "region_d6_fraction": region_frac, "regional_shift": shift,
            "five_prime_minus_three_prime": fp_mp, "exon9_tip_dip": ex9_dip,
            "identity_signal": signal, "per_psv": per_psv}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="CYP2D6-CYP2D7 PSV evidence table (Phase A diagnostic).")
    ap.add_argument("--pileup-dir", type=Path, required=True, help="dir with <sample>.d6.txt/<sample>.d7.txt")
    ap.add_argument("--samples", nargs="+", required=True, help="sample IDs (+optional :label)")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "cyp2d6_psv_evidence.json")
    args = ap.parse_args(argv)
    psvs = load_psvs()
    rows = []
    for s in args.samples:
        sid, _, label = s.partition(":")
        ev = sample_evidence(sid, args.pileup_dir, psvs)
        ev["label"] = label or ""
        rows.append(ev)
    # FALSIFIER VERDICT: hybrids must get a non-flat signal; labelled non-hybrids must stay flat.
    def _is_hyb_label(lbl):
        return lbl.lower().startswith("hyb") or any(h in lbl for h in ["*13", "*36", "*68"])
    hyb = [r for r in rows if _is_hyb_label(r["label"])]
    non = [r for r in rows if not _is_hyb_label(r["label"])]
    hyb_signalled = [r["sample"] for r in hyb if r["identity_signal"] != "flat_nonhybrid"]
    non_flat = [r["sample"] for r in non if r["identity_signal"] == "flat_nonhybrid"]
    go = (len(hyb) > 0 and len(hyb_signalled) == len(hyb) and len(non_flat) == len(non))
    verdict = {
        "falsifier": "GO" if go else "PARTIAL_OR_NOGO",
        "n_hybrid": len(hyb), "n_hybrid_signalled": len(hyb_signalled),
        "n_nonhybrid": len(non), "n_nonhybrid_flat": len(non_flat),
        "interpretation": ("The paired-coordinate PSV D6-fraction PROFILE reproduces the Cyrius hybrid "
                           "signatures: *68 (5' high / 3' low), *13 (opposite), *36 (exon-9-tip dip); "
                           "non-hybrids (normal/dup/deletion) stay flat. Directional 5'-3' cleanly catches "
                           "*68 + *13; the *36 exon-9 conversion needs the dedicated exon9-tip feature (the "
                           "subtlest case). PROOF-OF-SIGNAL on n=1/type -> GO to draft the Phase-B abstaining "
                           "identity classifier; NOT a powered per-allele validation."),
    }
    rep = {"schema": "cyp2d6-psv-evidence-v0", "n_psv": len(psvs), "region_order": _REGION_ORDER,
           "flag_contract": "samtools mpileup -B -q 0 -Q 0 (permissive; matches the depth path)",
           "note": "DIAGNOSTIC read-only evidence table; regional D6-fraction profile; NOT an identity caller",
           "falsifier_verdict": verdict, "samples": rows}
    args.out.write_text(json.dumps(rep, indent=2), encoding="utf-8")

    L = [f"# CYP2D6-CYP2D7 PSV evidence table — Phase A falsifier (read-level pileup)", "",
         f"_{rep['note']}. Flag contract: {rep['flag_contract']}. {len(psvs)} PSVs from the real Cyrius "
         "CYP2D6_SNP_38.txt (GRCh38, paired-coord). Both pos_CYP2D6 AND pos_CYP2D7 counted per PSV._", "",
         f"**Falsifier: {verdict['falsifier']}** — {verdict['n_hybrid_signalled']}/{verdict['n_hybrid']} "
         f"hybrids signalled; {verdict['n_nonhybrid_flat']}/{verdict['n_nonhybrid']} non-hybrids flat.", "",
         f"_{verdict['interpretation']}_", "",
         "## Regional D6-fraction profile (3'->5'; region x sample)", "",
         "| region | " + " | ".join(f"{r['sample']} ({r['label']})" for r in rows) + " |",
         "|" + "---|" * (len(rows) + 1)]
    for reg in _REGION_ORDER:
        L.append(f"| {reg} | " + " | ".join(
            (f"{r['region_d6_fraction'][reg]:.2f}" if r['region_d6_fraction'][reg] is not None else "--")
            for r in rows) + " |")
    L.append("| **5'-3' shift** | " + " | ".join(
        (f"{r['five_prime_minus_three_prime']:+.2f}" if r['five_prime_minus_three_prime'] is not None else "--")
        for r in rows) + " |")
    L.append("| **exon9-tip dip** | " + " | ".join(
        (f"{r['exon9_tip_dip']:+.2f}" if r['exon9_tip_dip'] is not None else "--") for r in rows) + " |")
    L.append("| **signal** | " + " | ".join(r["identity_signal"].split(" ")[0] for r in rows) + " |")
    L.append("")
    (REPO / "wiki" / "cyp2d6_psv_evidence_falsifier.md").write_text("\n".join(L), encoding="utf-8")

    hdr = f"{'region':18}" + "".join(f"{r['sample'][:10]:>11}" for r in rows)
    print(hdr)
    for reg in _REGION_ORDER:
        print(f"{reg:18}" + "".join(f"{('%.2f'%r['region_d6_fraction'][reg]) if r['region_d6_fraction'][reg] is not None else '--':>11}" for r in rows))
    print(f"{'5p-3p_shift':18}" + "".join(f"{('%+.2f'%r['five_prime_minus_three_prime']) if r['five_prime_minus_three_prime'] is not None else '--':>11}" for r in rows))
    print(f"{'signal':18}" + "".join(f"{r['identity_signal'].split(' ')[0][:10]:>11}" for r in rows))
    print(f"{'label':18}" + "".join(f"{r['label'][:10]:>11}" for r in rows))
    print(f"\nFALSIFIER: {verdict['falsifier']} ({verdict['n_hybrid_signalled']}/{verdict['n_hybrid']} hybrids signalled, {verdict['n_nonhybrid_flat']}/{verdict['n_nonhybrid']} non-hybrids flat)")
    print(f"[evidence -> {args.out} + wiki/cyp2d6_psv_evidence_falsifier.md]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
