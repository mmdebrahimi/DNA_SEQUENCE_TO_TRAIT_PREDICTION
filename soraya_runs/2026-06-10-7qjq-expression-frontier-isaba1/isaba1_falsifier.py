"""EXPRESSION-frontier falsifier: does ISAba1-upstream-of-blaOXA-51-family recover the carbapenem-R
strains that gene-PRESENCE misses (the documented Acinetobacter FN ceiling)?

Hypothesis (from /hypothesise #3 + the wider-AMR EXPRESSION boundary): in A. baumannii, an ISAba1 insertion
immediately upstream of the intrinsic blaOXA-51-family gene provides a hybrid promoter that OVEREXPRESSES it
-> meropenem-R. This is invisible to gene-presence (every isolate HAS OXA-51), so it is the dominant cause
of the deterministic rule's false negatives. If detecting the ISAba1->OXA-51 junction in the assembly
SEPARATES the intrinsic-only-R strains from S, it is a real expression-context signal that could cross the
abstain floor.

Method (per cached assembly): blastn ISAba1 ref + OXA-51-family ref vs the genome; for each OXA-51 hit,
check whether an ISAba1 hit sits within UPSTREAM_BP (5' of the gene, strand-aware) on the same contig.
Cross-tabulate junction-positive vs R/S label AND vs strong-acquired-carbapenemase presence (the FN set =
R strains with NO strong acquired carbapenemase).

KILL if the junction does NOT enrich in R / does NOT mark the strong-acquired-negative R strains.
Pure: blastn over cached assemblies; no money, no Docker.
"""
import subprocess, sys, tempfile
from collections import Counter
from pathlib import Path
sys.path.insert(0, ".")
from scripts.fungal_erg11_caller import _find

UPSTREAM_BP = 400
ISABA1 = "data/isaba1_ref/ISAba1_ref.fna"
OXA = "data/isaba1_ref/OXA51fam_ref.fna"
BASE = Path("data/raw/acinetobacter_meropenem")
STRONG = ("OXA-23", "OXA-24", "OXA-40", "OXA-72", "OXA-25", "OXA-26", "OXA-143", "OXA-235",
          "NDM", "IMP", "VIM", "KPC")


def _blast(query, genome, blastn, makedb):
    with tempfile.TemporaryDirectory() as td:
        db = str(Path(td) / "g")
        subprocess.run([makedb, "-in", genome, "-dbtype", "nucl", "-out", db], check=True, capture_output=True)
        out = subprocess.run([blastn, "-query", query, "-db", db, "-outfmt",
                              "6 sseqid sstart send pident length bitscore", "-max_target_seqs", "20"],
                             check=True, capture_output=True, text=True).stdout
    hits = []
    for ln in out.strip().splitlines():
        f = ln.split("\t")
        if float(f[3]) >= 85 and int(f[4]) >= 120:     # decent identity + length
            ss, se = int(f[1]), int(f[2])
            hits.append({"contig": f[0], "lo": min(ss, se), "hi": max(ss, se),
                         "strand": "+" if ss < se else "-", "bits": float(f[5])})
    return hits


def junction_positive(genome, blastn, makedb):
    oxa = _blast(OXA, genome, blastn, makedb)
    isa = _blast(ISABA1, genome, blastn, makedb)
    for o in oxa:
        # 5' (upstream) boundary of the OXA gene, strand-aware
        if o["strand"] == "+":
            win_lo, win_hi = o["lo"] - UPSTREAM_BP, o["lo"]
        else:
            win_lo, win_hi = o["hi"], o["hi"] + UPSTREAM_BP
        for i in isa:
            if i["contig"] == o["contig"] and i["hi"] >= win_lo and i["lo"] <= win_hi:
                return True, len(oxa), len(isa)
    return False, len(oxa), len(isa)


def strong_acquired(acc):
    mt = BASE / "amrfinder_runs" / acc / "main.tsv"
    if not mt.exists():
        return None
    txt = mt.read_text(encoding="utf-8", errors="replace")
    return any(tok in txt for tok in STRONG)


def main():
    blastn, makedb = _find("blastn"), _find("makeblastdb")
    if not blastn or not makedb:
        print("NO BLAST"); return
    labels = {}
    for ln in (BASE / "selected.tsv").read_text().splitlines():
        if "\t" in ln:
            a, rs = ln.split("\t"); labels[a] = rs.strip()
    rows = []
    for acc, lab in labels.items():
        g = BASE / "refseq" / acc / "genome.fna"
        if not g.exists():
            continue
        jpos, n_oxa, n_isa = junction_positive(str(g), blastn, makedb)
        sa = strong_acquired(acc)
        rows.append({"acc": acc, "label": lab, "junction": jpos, "strong_acquired": sa,
                     "n_oxa": n_oxa, "n_isaba1": n_isa})
        print(f"  {acc} {lab} junction={jpos} strong_acq={sa} (OXA={n_oxa} ISAba1={n_isa})")
    # cross-tabs
    R = [r for r in rows if r["label"] == "R"]; S = [r for r in rows if r["label"] == "S"]
    jR = sum(r["junction"] for r in R); jS = sum(r["junction"] for r in S)
    print(f"\nJunction-positive: R {jR}/{len(R)}  S {jS}/{len(S)}")
    # the FN set: R strains with NO strong acquired carbapenemase (what gene-presence misses)
    fn = [r for r in R if r["strong_acquired"] is False]
    fn_j = sum(r["junction"] for r in fn)
    print(f"Intrinsic-only R (no strong acquired = the FN ceiling): {len(fn)}; of these junction-positive: {fn_j}")
    # would junction RESCUE any: R strains that are junction+ AND strong-acquired-negative
    rescued = [r["acc"] for r in fn if r["junction"]]
    # false-rescue risk: S strains that are junction+ (would become wrongly-R)
    false_rescue = [r["acc"] for r in S if r["junction"]]
    print(f"\nRESCUE candidates (intrinsic-only R, junction+): {rescued}")
    print(f"FALSE-rescue risk (S, junction+): {false_rescue}")
    import json
    Path("soraya_runs/2026-06-10-7qjq-expression-frontier-isaba1/isaba1_falsifier_result.json").write_text(
        json.dumps({"rows": rows, "junction_R": [jR, len(R)], "junction_S": [jS, len(S)],
                    "intrinsic_only_R": len(fn), "intrinsic_only_R_junction_pos": fn_j,
                    "rescued": rescued, "false_rescue": false_rescue}, indent=2))
    verdict = ("SURVIVES" if (fn_j >= 1 and len(false_rescue) <= len(S) * 0.25 and jR > jS)
               else "KILLED")
    print(f"\nVERDICT: {verdict}  (junction enriches in R: {jR}>{jS}? ; recovers intrinsic-only-R: {fn_j}/{len(fn)}; "
          f"false-rescue S: {len(false_rescue)}/{len(S)})")


if __name__ == "__main__":
    main()
