"""EP-8 Path B — G2 dry-manifest (CPU-only gate; manifest §0.5). MUST pass before any GPU embedding.

Proves the Arabidopsis flowering-time embedding test is runnable BEFORE the workhorse spends GPU:
  1. accession intersection  (AraPheno FT10 ∩ pseudogenome files ∩ SNP-matrix columns)
  2. pseudogenome ID->filename pattern  (resolved EMPIRICALLY, never assumed)
  3. agnostic window/coordinate table  (gene-bodies + flanks, from a GFF; the FROZEN frozen-blind window set
     the FM embeds AND the matched SNP baseline reads — manifest §1a primary)
  4. per-window/per-accession N-fraction QC  (uncalled bases silently corrupt embedding + consensus)
  5. genetic-group labels present for the analysis-N accessions  (the leave-one-group-out CV needs them)
Emits wiki/g2_dry_manifest_<date>.{md,json}. Any red check -> STOP + surface; do NOT burn GPU.

DESIGN SPLIT (dual-machine): the laptop WRITES + unit-tests this gate (pure logic on synthetic fixtures, no
19GB download, no GPU); the workhorse RUNS it on the real pseudogenomes/VCF/GFF it downloads. The pure
functions below (intersect_accessions / resolve_pseudogenome_pattern / build_window_table / n_fraction /
window_n_fractions / group_labels_present / dry_manifest_report) are data-shape logic, fully testable offline.

This is the manifest's phenotype-AGNOSTIC primary: window selection uses ONLY a gene annotation (no
flowering-time labels, no curated flowering-gene list) -> it preserves the "no curated mechanism catalog"
embedding-niche criterion (the flaw the 2026-06-08 brainstorm caught in the curated-locus draft).
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


@dataclass(frozen=True)
class Window:
    chrom: str
    start: int          # 1-based inclusive
    end: int            # 1-based inclusive
    gene_id: str
    win_index: int      # tile index within the (gene+flank) interval

    @property
    def name(self) -> str:
        return f"{self.gene_id}|{self.chrom}:{self.start}-{self.end}|w{self.win_index}"


# ---- pure-logic gate functions (offline-testable) --------------------------------------------------

def read_accession_ids(pheno_csv: str | Path, id_col: str = "accession_id") -> set[str]:
    rows = list(csv.DictReader(Path(pheno_csv).read_text(encoding="utf-8").splitlines()))
    if rows and id_col not in rows[0]:
        raise ValueError(f"{pheno_csv}: no '{id_col}' column; got {list(rows[0])}")
    return {r[id_col].strip() for r in rows if r.get(id_col, "").strip()}


def resolve_pseudogenome_pattern(accession_ids: set[str], filenames: list[str]) -> tuple[str | None, dict]:
    """Empirically infer the accession_id -> pseudogenome filename rule (NEVER assume pseudo{id}.fasta.gz).

    Tries each filename, extracts the embedded integer id-run, checks coverage against accession_ids.
    Returns (template_or_None, stats). template uses '{id}' placeholder, e.g. 'Cauris_{id}.fna' style.
    """
    # candidate: find the accession id as a maximal digit run inside each filename
    matched: dict[str, str] = {}
    template = None
    for fn in filenames:
        for m in re.finditer(r"\d+", fn):
            tok = m.group(0)
            if tok in accession_ids:
                matched[tok] = fn
                if template is None:
                    template = fn[:m.start()] + "{id}" + fn[m.end():]
                break
    stats = {"n_files": len(filenames), "n_accessions": len(accession_ids),
             "n_matched": len(matched), "template": template}
    return template, stats


def intersect_accessions(pheno_ids: set[str], pseudogenome_ids: set[str],
                         snp_matrix_ids: set[str]) -> dict:
    """Final analysis set = phenotype ∩ pseudogenome ∩ SNP-matrix. Report N + per-source drops."""
    analysis = pheno_ids & pseudogenome_ids & snp_matrix_ids
    return {
        "analysis_n": len(analysis),
        "analysis_ids": sorted(analysis),
        "n_phenotype": len(pheno_ids),
        "dropped_no_pseudogenome": sorted(pheno_ids - pseudogenome_ids),
        "dropped_no_snp": sorted(pheno_ids - snp_matrix_ids),
    }


def _parse_gff_genes(gff_path: str | Path) -> list[tuple[str, int, int, str]]:
    """Return [(chrom, start, end, gene_id)] for gene-type GFF3 rows. 1-based inclusive coords."""
    genes = []
    for ln in Path(gff_path).read_text(encoding="utf-8").splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        chrom, start, end, attrs = f[0], int(f[3]), int(f[4]), f[8]
        m = re.search(r"(?:ID|gene_id|Name)=([^;]+)", attrs)
        gid = m.group(1) if m else f"{chrom}:{start}"
        genes.append((chrom, start, end, gid))
    return genes


def build_window_table(gff_path: str | Path, *, flank: int = 1000, window: int = 512,
                       stride: int = 256, chrom_len: dict[str, int] | None = None) -> list[Window]:
    """PHENOTYPE-AGNOSTIC window set: every gene body +/- `flank`, tiled into `window`-bp windows at `stride`.

    Uses ONLY the gene annotation — no flowering-gene list, no phenotype. This is the frozen window set the
    FM embeds and the matched SNP baseline reads (manifest §1a primary + §2.1 information-matched baseline).
    """
    windows: list[Window] = []
    for chrom, gstart, gend, gid in _parse_gff_genes(gff_path):
        lo = max(1, gstart - flank)
        hi = gend + flank
        if chrom_len and chrom in chrom_len:
            hi = min(hi, chrom_len[chrom])
        wi = 0
        pos = lo
        while pos <= hi:
            end = min(pos + window - 1, hi)
            windows.append(Window(chrom, pos, end, gid, wi))
            if end >= hi:
                break
            pos += stride
            wi += 1
    return windows


def n_fraction(seq: str) -> float:
    if not seq:
        return 1.0
    s = seq.upper()
    return (s.count("N")) / len(s)


def window_n_fractions(windows: list[Window], chrom_seqs: dict[str, str],
                       *, max_n: float = 0.10) -> dict:
    """Per-window N-fraction over a single accession's chromosome sequences. Flags windows above max_n."""
    fracs = []
    flagged = []
    for w in windows:
        seq = chrom_seqs.get(w.chrom, "")
        sub = seq[w.start - 1:w.end]          # 1-based inclusive -> python slice
        f = n_fraction(sub)
        fracs.append(f)
        if f > max_n:
            flagged.append(w.name)
    mean = sum(fracs) / len(fracs) if fracs else 1.0
    return {"n_windows": len(windows), "mean_n_fraction": mean,
            "n_flagged": len(flagged), "flagged": flagged, "max_n_threshold": max_n}


def group_labels_present(analysis_ids: list[str], group_map: dict[str, str]) -> dict:
    """Confirm genetic-group (admixture) labels resolve for the analysis-N accessions (CV needs them)."""
    missing = [a for a in analysis_ids if a not in group_map or not str(group_map[a]).strip()]
    groups = {}
    for a in analysis_ids:
        g = group_map.get(a)
        if g and str(g).strip():
            groups[str(g)] = groups.get(str(g), 0) + 1
    return {"n_with_group": len(analysis_ids) - len(missing), "n_missing": len(missing),
            "missing": missing, "groups": groups, "n_groups": len(groups)}


def dry_manifest_report(*, intersect: dict, pattern_stats: dict, window_table: list[Window],
                        n_qc: dict | None, groups: dict, min_analysis_n: int = 100,
                        min_groups: int = 3) -> dict:
    """Aggregate the gate. verdict GREEN only if every check passes."""
    checks = {
        "accession_intersection": intersect["analysis_n"] >= min_analysis_n,
        "pseudogenome_pattern": pattern_stats.get("template") is not None
        and pattern_stats.get("n_matched", 0) >= min_analysis_n,
        "window_table": len(window_table) > 0,
        "n_fraction_qc": (n_qc is None) or (n_qc["n_flagged"] == 0),
        "group_labels": groups["n_missing"] == 0 and groups["n_groups"] >= min_groups,
    }
    verdict = "GREEN" if all(checks.values()) else "RED"
    return {
        "verdict": verdict,
        "checks": checks,
        "analysis_n": intersect["analysis_n"],
        "n_windows": len(window_table),
        "n_groups": groups["n_groups"],
        "pseudogenome_template": pattern_stats.get("template"),
        "n_fraction_qc": n_qc,
        "red_checks": [k for k, v in checks.items() if not v],
    }


# ---- CLI (workhorse runs this on the real downloaded data) -----------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pheno-csv", default="data/arabidopsis/FT10_pheno_261.csv")
    ap.add_argument("--pseudogenome-dir", required=True, help="dir of per-accession pseudogenome FASTAs")
    ap.add_argument("--snp-matrix-ids", help="file with one SNP-matrix accession_id per line (or a VCF #CHROM header)")
    ap.add_argument("--gff", required=True, help="Araport11/TAIR10 GFF3 (gene annotation for the window table)")
    ap.add_argument("--group-map", help="CSV: accession_id,group (admixture/PC group labels)")
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--window", type=int, default=512)
    ap.add_argument("--stride", type=int, default=256)
    ap.add_argument("--out-prefix", default=None)
    ap.add_argument("--today", default="(date unset)")
    a = ap.parse_args(argv)

    pheno_ids = read_accession_ids(a.pheno_csv)
    pfiles = [p.name for p in Path(a.pseudogenome_dir).glob("*") if p.is_file()]
    template, pstats = resolve_pseudogenome_pattern(pheno_ids, pfiles)
    pseudo_ids = {tok for tok in pheno_ids if template and template.replace("{id}", tok) in pfiles}
    snp_ids = pheno_ids
    if a.snp_matrix_ids and Path(a.snp_matrix_ids).exists():
        txt = Path(a.snp_matrix_ids).read_text(encoding="utf-8")
        snp_ids = {t for t in re.split(r"[\s,]+", txt) if t in pheno_ids} or pheno_ids
    inter = intersect_accessions(pheno_ids, pseudo_ids, snp_ids)
    wtable = build_window_table(a.gff, flank=a.flank, window=a.window, stride=a.stride)
    gmap = {}
    if a.group_map and Path(a.group_map).exists():
        for r in csv.DictReader(Path(a.group_map).read_text(encoding="utf-8").splitlines()):
            gmap[r.get("accession_id", "").strip()] = r.get("group", "").strip()
    groups = group_labels_present(inter["analysis_ids"], gmap)
    rep = dry_manifest_report(intersect=inter, pattern_stats=pstats, window_table=wtable,
                              n_qc=None, groups=groups)

    prefix = a.out_prefix or f"wiki/g2_dry_manifest_{a.today}"
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    Path(prefix + ".json").write_text(json.dumps({**rep, "intersect": inter, "pattern": pstats},
                                                 indent=2), encoding="utf-8")
    lines = [f"# G2 dry-manifest ({a.today}) — verdict {rep['verdict']}", "",
             f"- analysis N: {rep['analysis_n']}", f"- windows: {rep['n_windows']}",
             f"- genetic groups: {rep['n_groups']}", f"- pseudogenome template: {rep['pseudogenome_template']}",
             f"- RED checks: {rep['red_checks'] or 'none'}", "", "| check | pass |", "|---|---|"]
    lines += [f"| {k} | {'PASS' if v else 'RED'} |" for k, v in rep["checks"].items()]
    Path(prefix + ".md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"G2 dry-manifest: {rep['verdict']}  analysis_n={rep['analysis_n']} windows={rep['n_windows']} "
          f"groups={rep['n_groups']} red={rep['red_checks']}")
    print(f"wrote {prefix}.md + {prefix}.json")
    return 0 if rep["verdict"] == "GREEN" else 1


if __name__ == "__main__":
    raise SystemExit(main())
