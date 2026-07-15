"""Pre-resistance "one-nt-from-resistance" base-rate census (the /idea-validation-council go/no-go test).

The council verdict on the "escape / pre-resistance forecaster" idea was **Pursue test**: four independent
lenses converged on one crux — catalogued point-mutations are single-nt-accessible from wildtype BY
CONSTRUCTION, so "self-validation by recovering catalogued codon-neighbours" is circular, and the tool only
adds information beyond the catalogue if REAL genomes carry a NON-WT, NON-resistant INTERMEDIATE codon that
sits STRICTLY CLOSER to a resistant codon than wildtype does. This script measures both:

  PART A  accessibility census: fraction of catalogued DRMs whose wildtype reference codon is ONE nt from the
          nearest resistant codon. Expected >=90% -> confirms the self-validation is vacuous (single-nt by
          construction). (No genome needed — reference CDS + genetic code + catalogue only.)

  PART B  intermediate-carrier scan: the tool's value lives ONLY at multi-step DRMs (wildtype >=2 nt from
          resistance). For those, an aa I (I!=WT, I!=R) is a "closer intermediate" iff some codon of I is
          STRICTLY closer to a resistant codon than the WT codon is. Scan the real HIVDB isolate panels for
          genomes carrying such a residue at a catalogued position. PASS iff >=5 distinct intermediate-
          carriers exist; ~0 => NO-GO (the forecaster is a codon-table restatement of the catalogue).

Substrate: HIV (the cleanest — mutant-level WT->R catalogues NNRTI-RT + CAI-capsid, committed HXB2 ref CDS,
and position-columnar per-isolate aa panels in data/raw/hiv/*_DataSet.txt). HONEST SCOPE: the panels are
amino-acid-level, so Part B is a BEST-CASE UPPER BOUND on true intermediate-carriers (a residue is counted if
ANY of its codons is strictly closer to R than WT — the isolate's exact codon is unknown at aa resolution).
An upper bound <5 is therefore a decisive NO-GO; an upper bound >=5 would flag codon-level confirmation
(TB/cipro nucleotide data) as the follow-up. Deterministic, offline, no GPU. Frozen decoder surface untouched.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

HIV_RAW = REPO / "data" / "raw" / "hiv"
HIV_REF = REPO / "data" / "hiv_ref"

# Standard genetic code (DNA codons -> 1-letter aa; '*' = stop).
CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
AAS = sorted(set(CODON.values()) - {"*"})
CODONS_FOR = {aa: [c for c, a in CODON.items() if a == aa] for aa in set(CODON.values())}


def nt_hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def min_codon_dist(aa_from_codons: list[str], aa_to: str) -> int:
    """Minimum nt distance from any codon in `aa_from_codons` to any codon of `aa_to`."""
    tos = CODONS_FOR[aa_to]
    return min(nt_hamming(cf, ct) for cf in aa_from_codons for ct in tos)


def read_cds(path: Path) -> str:
    seq = "".join(l.strip() for l in path.read_text().splitlines() if not l.startswith(">"))
    return seq.upper().replace("U", "T")


def ref_codon(cds: str, aa_pos: int) -> str:
    i = (aa_pos - 1) * 3
    return cds[i:i + 3]


# ---- catalogues (mutant-level WT->R pairs only; position-based classes have no single R aa) -------------
def load_drms() -> list[dict]:
    """Return catalogued DRMs as {gene, pos, wt, res, cds_file, dataset, source}. Mutant-level only."""
    from dna_decode.data.hiv_amr import NNRTI_RT_MAJOR_DRMS, _RT_WT
    from dna_decode.data.hiv_amr import CAI_CLASS  # capsid mutant-level

    drms = []
    for m in sorted(NNRTI_RT_MAJOR_DRMS):
        wt, pos, res = m[0], int(m[1:-1]), m[-1]
        assert _RT_WT[pos] == wt, f"catalogue WT mismatch {m} vs _RT_WT"
        drms.append({"gene": "RT", "pos": pos, "wt": wt, "res": res, "drm": m,
                     "cds_file": "HIV1_RT_HXB2_cds.fna", "dataset": "NNRTI_DataSet.txt", "cls": "NNRTI"})
    for m in sorted(CAI_CLASS.major_drms):
        wt, pos, res = m[0], int(m[1:-1]), m[-1]
        assert CAI_CLASS.wt[pos] == wt, f"CAI WT mismatch {m}"
        drms.append({"gene": "CA", "pos": pos, "wt": wt, "res": res, "drm": m,
                     "cds_file": "HIV1_CA_HXB2_cds.fna", "dataset": "CAI_DataSet.txt", "cls": "CAI"})
    return drms


def load_panel(dataset: str) -> tuple[list[str], dict[int, int]]:
    """Return (rows, pos->col_index) for a HIVDB position-columnar panel. Header P<n> => residue col."""
    path = HIV_RAW / dataset
    with open(path, encoding="utf-8", newline="") as f:
        rdr = csv.reader(f, delimiter="\t")
        header = next(rdr)
        pos_col = {int(h[1:]): i for i, h in enumerate(header) if h.startswith("P") and h[1:].isdigit()}
        rows = [row for row in rdr]
    return rows, pos_col, header[0]


def cell_residues(cell: str) -> set[str]:
    """Residues at a panel cell. '-'/'.' = consensus(WT); a letter = substitution; mixtures = each letter."""
    cell = (cell or "").strip().upper()
    if cell in ("-", ".", "", "*", "X", "NA"):
        return set()
    return {c for c in cell if c.isalpha() and c in AAS}


def main() -> int:
    drms = load_drms()
    cds_cache = {f: read_cds(HIV_REF / f) for f in {d["cds_file"] for d in drms}}

    # All catalogued resistant residues at each position (a residue that IS a DRM aa is 'resistant', not an
    # intermediate — exclude them so Part B counts only SUSCEPTIBLE-but-primed genomes, not already-resistant).
    resistant_aas_by_pos = {}
    for d in drms:
        resistant_aas_by_pos.setdefault((d["gene"], d["pos"]), set()).add(d["res"])

    # ---------- PART A: accessibility + closer-intermediate residue sets ----------
    partA = []
    closer_by_key = {}      # (gene,pos) -> {residues strictly closer to R than WT AND not themselves a DRM aa}
    for d in drms:
        cds = cds_cache[d["cds_file"]]
        cW = ref_codon(cds, d["pos"])
        tr = CODON.get(cW, "?")
        ref_ok = (tr == d["wt"])
        d_wt = min_codon_dist([cW], d["res"]) if ref_ok else None
        r_at_pos = resistant_aas_by_pos[(d["gene"], d["pos"])]
        closer = set()
        if ref_ok:
            for I in AAS:
                if I == d["wt"] or I in r_at_pos:          # skip WT and ANY catalogued resistant residue
                    continue
                if min_codon_dist(CODONS_FOR[I], d["res"]) < d_wt:
                    closer.add(I)
        partA.append({**{k: d[k] for k in ("gene", "pos", "wt", "res", "drm", "cls")},
                      "ref_codon": cW, "ref_translates": tr, "ref_ok": ref_ok,
                      "wt_dist_to_resistance": d_wt, "single_nt_accessible": (d_wt == 1) if ref_ok else None,
                      "n_closer_intermediate_aas": len(closer)})
        closer_by_key.setdefault((d["gene"], d["pos"]), set()).update(closer)

    ok = [r for r in partA if r["ref_ok"]]
    n_acc = sum(1 for r in ok if r["single_nt_accessible"])
    acc_frac = n_acc / len(ok) if ok else float("nan")
    multi_step = [r for r in ok if r["wt_dist_to_resistance"] >= 2]

    # ---------- PART B: intermediate-carrier scan over real isolate panels ----------
    # Only DRMs with >=1 closer-intermediate aa can yield a meaningful carrier (all at multi-step positions).
    scan_keys = {(d["gene"], d["pos"], d["dataset"]) for d in drms if closer_by_key[(d["gene"], d["pos"])]}
    broad_keys = {(d["gene"], d["pos"], d["dataset"], d["wt"]) for d in drms}

    carriers_meaningful = set()   # (dataset,seqid) carrying a strictly-closer intermediate at a scan position
    carriers_broad = set()        # (dataset,seqid) carrying ANY non-WT non-R residue at a catalogued position
    per_position = {}
    panels = {}
    for ds in {k[2] for k in scan_keys} | {k[2] for k in broad_keys}:
        panels[ds] = load_panel(ds)

    for (gene, pos, ds) in scan_keys:
        rows, pos_col, _ = panels[ds]
        col = pos_col.get(pos)
        closer = closer_by_key[(gene, pos)]
        if col is None:
            per_position[f"{gene}{pos}"] = {"status": "POSITION_COLUMN_ABSENT", "dataset": ds}
            continue
        hits = 0
        for row in rows:
            if col >= len(row):
                continue
            res = cell_residues(row[col])
            if res & closer:
                carriers_meaningful.add((ds, row[0]))
                hits += 1
        per_position[f"{gene}{pos}"] = {"status": "scanned", "dataset": ds,
                                        "closer_intermediate_aas": sorted(closer),
                                        "n_isolates_carrying_closer_intermediate": hits}

    for (gene, pos, ds, wt) in broad_keys:
        rows, pos_col, _ = panels[ds]
        col = pos_col.get(pos)
        if col is None:
            continue
        r_at_pos = resistant_aas_by_pos[(gene, pos)]
        for row in rows:
            if col >= len(row):
                continue
            res = cell_residues(row[col])
            if res - {wt} - r_at_pos:        # any non-WT, non-(catalogued-resistant) residue present
                carriers_broad.add((ds, row[0]))

    n_meaningful = len(carriers_meaningful)
    verdict = "GO" if n_meaningful >= 5 else "NO_GO"

    result = {
        "cell": "pre_resistance_base_rate_census",
        "council_decision_under_test": "Pursue test (escape/pre-resistance forecaster)",
        "substrate": "HIV mutant-level catalogues (NNRTI-RT + CAI-capsid) x HIVDB isolate panels",
        "honest_scope": ("Part B is an aa-level UPPER BOUND on intermediate-carriers (isolate exact codon "
                         "unknown at aa resolution); upper bound <5 => decisive NO-GO"),
        "part_A_accessibility": {
            "n_drms_total": len(partA), "n_drms_ref_ok": len(ok),
            "n_single_nt_accessible": n_acc, "accessibility_fraction": round(acc_frac, 4),
            "expectation": ">=0.90 confirms self-validation is vacuous (single-nt by construction)",
            "confirmed_vacuous": acc_frac >= 0.90,
            "n_multi_step_drms": len(multi_step),
            "multi_step_drms": [r["drm"] for r in multi_step],
        },
        "part_B_intermediate_carriers": {
            "n_scan_positions": len({(k[0], k[1]) for k in scan_keys}),
            "n_intermediate_carriers_meaningful": n_meaningful,
            "n_carriers_broad_any_nonwt_nonR": len(carriers_broad),
            "pass_threshold": 5,
            "per_position": per_position,
        },
        "verdict": verdict,
        "interpretation": (
            f"Accessibility {acc_frac:.1%} of catalogued DRMs are single-nt from resistance by construction; "
            f"only {len(multi_step)} multi-step DRM(s) leave room for a 'closer intermediate'. "
            f"Meaningful intermediate-carriers found: {n_meaningful} (broad non-WT-at-DRM-position carriers: "
            f"{len(carriers_broad)}). Verdict {verdict} at the >=5 bar."),
    }
    out = REPO / "wiki" / f"pre_resistance_base_rate_census_{_date.today().isoformat()}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    a = result["part_A_accessibility"]; b = result["part_B_intermediate_carriers"]
    print(f"[census] PART A  accessibility: {a['n_single_nt_accessible']}/{a['n_drms_ref_ok']} "
          f"= {a['accessibility_fraction']:.1%} single-nt  (vacuous>={a['confirmed_vacuous']})")
    print(f"[census]         multi-step DRMs (WT>=2nt): {a['n_multi_step_drms']}  {a['multi_step_drms']}")
    print(f"[census] PART B  meaningful intermediate-carriers: {b['n_intermediate_carriers_meaningful']} "
          f"(broad non-WT@DRM: {b['n_carriers_broad_any_nonwt_nonR']})  bar>=5")
    for k, v in b["per_position"].items():
        if v.get("status") == "scanned":
            print(f"           {k}: closer={v['closer_intermediate_aas']} "
                  f"carriers={v['n_isolates_carrying_closer_intermediate']}")
        else:
            print(f"           {k}: {v['status']}")
    print(f"[census] VERDICT: {verdict}")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
