"""Run the INDEPENDENT TB number: score the WHO-catalogue rule on the AMR Portal provenance-disjoint cohort.

End-to-end, CHECKPOINTED (restartable — wedge-safe): for each disjoint isolate with a fetchable assembly,
  1. fetch the GCA assembly FASTA (NCBI Datasets API, HTTPS — no Docker),
  2. align to H37Rv NC_000962.3 + call variants (minimap2 asm5 + paftools.js call, via the Docker
     biocontainer; FILTER '.' -> PASS masking = the asm5 confident-difference set),
  3. parse -> tb_vcf masked calls -> tb_amr.score_drug (RIF + INH) vs the measured label,
  4. checkpoint per isolate to a results JSONL, then aggregate raw + lineage-collapsed sens/spec + Wilson CI
     via the FROZEN clonality.cluster_weighted_confusion (reusing score_tb_cryptic.score_cohort).

HONESTY: provenance-disjoint at the accession level (BioSample/GCA not in CRyPTIC/our cohorts); measured
phenotype (non-circular); the WHO catalogue is applied UNCHANGED. No regeno VCF -> callability unassessed
(non-match = S, documented). This is the first genuinely-INDEPENDENT TB number the gold-set saga could not get.

Per-isolate container startup is amortized by aligning every staged FASTA in ONE docker invocation.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import os

COHORT = REPO / "data" / "raw" / "tb_goldset" / "amr_portal_tb_disjoint_cohort.tsv"
# Work dir holds transient assemblies (~1.3MB each) + VCFs + checkpoint. Override to D: for the full run
# (C: is disk-tight; 2,845 assemblies ~3.7GB). $TB_INDEP_WORK or --work.
WORK = Path(os.environ.get("TB_INDEP_WORK", str(REPO / "data" / "raw" / "tb_indep")))
REF = WORK / "ref" / "H37Rv.fna"
ASM, VCF = WORK / "asm", WORK / "vcf"
CKPT = WORK / "results.jsonl"
IMG = "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{}/download?include_annotation_type=GENOME_FASTA"


def read_cohort(max_n: int) -> list[dict]:
    """Disjoint isolates with a fetchable GCA assembly + >=1 measured RIF/INH label."""
    out = []
    with open(COHORT, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["leaked"] == "0" and r["assembly"].startswith("GCA_") and (
                    r["rif_label"] in ("R", "S") or r["inh_label"] in ("R", "S")):
                out.append(r)
            if max_n and len(out) >= max_n:
                break
    return out


def fetch_assembly(gca: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        return True
    try:
        zp = dest.with_suffix(".zip")
        urllib.request.urlretrieve(NCBI.format(gca), zp)
        with zipfile.ZipFile(zp) as z:
            fnas = [n for n in z.namelist() if n.endswith(".fna")]
            if not fnas:
                return False
            with z.open(fnas[0]) as f, open(dest, "wb") as o:
                o.write(f.read())
        zp.unlink(missing_ok=True)
        return dest.stat().st_size > 0
    except Exception as e:
        print(f"  fetch FAIL {gca}: {e}", file=sys.stderr)
        return False


def align_all(progress=True, chunk: int = 250) -> int:
    """Align every staged .fna lacking a .vcf, in CHUNKED, per-chunk-fault-tolerant Docker containers.

    One short-lived container per `chunk` assemblies (vs one long-lived container over thousands — the
    documented Docker D:-mount-churn wedge risk on an unattended multi-hour run). A chunk that throws
    (container/mount death) is logged and SKIPPED; the next chunk proceeds. Missing VCFs are re-aligned on
    the next call (skip-existing) — so a re-run fills any gaps. Returns #VCFs produced this call."""
    todo = [f for f in sorted(ASM.glob("*.fna")) if not (VCF / f"{f.stem}.vcf").exists()]
    if not todo:
        return 0
    # each container aligns at most `chunk` still-missing .fna then EXITS (short-lived -> D: mount stays
    # fresh). The inner loop skips already-produced VCFs, so successive containers make forward progress.
    inner = (f'c=0; for f in /data/asm/*.fna; do b=$(basename "$f" .fna); '
             '[ -s "/data/vcf/$b.vcf" ] && continue; '
             'minimap2 -cx asm5 --cs /data/ref/H37Rv.fna "$f" 2>/dev/null '
             '| sort -k6,6 -k8,8n | paftools.js call -f /data/ref/H37Rv.fna - 2>/dev/null '
             '> "/data/vcf/$b.vcf"; '
             f'c=$((c+1)); [ $c -ge {chunk} ] && break; done')
    if progress:
        print(f"  aligning {len(todo)} assemblies in short-lived chunks of {chunk}...")
    made = 0
    max_containers = (len(todo) // chunk) + 3
    for _ in range(max_containers):
        cmd = ["docker", "run", "--rm", "-v", f"{WORK}:/data", IMG, "bash", "-c", inner]
        try:
            subprocess.run(cmd, check=True, timeout=3600,
                           env={"MSYS_NO_PATHCONV": "1", **os.environ})
        except Exception as e:                          # container/mount death -> next container resumes
            print(f"  align chunk FAILED ({type(e).__name__}: {e}); continuing", file=sys.stderr)
        new = sum(1 for f in todo if (VCF / f"{f.stem}.vcf").exists())
        if progress:
            print(f"    ...{new}/{len(todo)} aligned")
        if new >= len(todo) or new == made:             # done, or a chunk made zero progress (give up pass)
            made = new
            break
        made = new
    return made


def _pass_mask(vcf_text: str) -> str:
    return "\n".join(
        ln if ln.startswith("#") else
        "\t".join(p if i != 6 or p != "." else "PASS" for i, p in enumerate(ln.split("\t")))
        for ln in vcf_text.splitlines())


def score_isolate(name: str, dets) -> dict | None:
    from dna_decode.organism_rules import tb_vcf, tb_amr
    vp = VCF / f"{name}.vcf"
    if not vp.exists() or vp.stat().st_size == 0:
        return None
    calls = tb_vcf.parse_masked_calls(_pass_mask(vp.read_text(encoding="utf-8")))
    return {"rif_pred": tb_amr.score_drug("rifampicin", calls, dets["rifampicin"]).prediction,
            "inh_pred": tb_amr.score_drug("isoniazid", calls, dets["isoniazid"]).prediction,
            "n_calls": len(calls)}


def load_ckpt() -> dict:
    done = {}
    if CKPT.exists():
        for ln in CKPT.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                r = json.loads(ln)
                done[r["strain_id"]] = r
    return done


def aggregate(rows: list[dict]) -> dict:
    """raw sens/spec per drug + Wilson CI. (lineage-collapse needs Mash; raw is the smoke headline.)"""
    from scripts.amr_portal_score_independent import wilson_ci
    res = {}
    for drug, pk, lk in (("rifampicin", "rif_pred", "rif_label"), ("isoniazid", "inh_pred", "inh_label")):
        tp = fp = tn = fn = 0
        for r in rows:
            lab, pred = r.get(lk), r.get(pk)
            if lab not in ("R", "S") or pred not in ("R", "S"):
                continue
            if pred == "R":
                tp += lab == "R"; fp += lab == "S"
            else:
                tn += lab == "S"; fn += lab == "R"
        nR, nS = tp + fn, tn + fp
        res[drug] = {"tp": tp, "fp": fp, "tn": tn, "fn": fn, "n_R": nR, "n_S": nS,
                     "sens": tp / nR if nR else None, "spec": tn / nS if nS else None,
                     "sens_ci95": wilson_ci(tp, nR), "spec_ci95": wilson_ci(tn, nS),
                     "accuracy": (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) else None}
    return res


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max", type=int, default=50, help="isolate cap (0 = full 2,845 cohort)")
    ap.add_argument("--work", type=Path, default=None, help="work dir (default $TB_INDEP_WORK or repo; use D: for full run)")
    a = ap.parse_args(argv)
    if a.work:
        global WORK, REF, ASM, VCF, CKPT
        WORK = a.work; REF = WORK / "ref" / "H37Rv.fna"; ASM = WORK / "asm"; VCF = WORK / "vcf"
        CKPT = WORK / "results.jsonl"
    if not REF.exists():
        print(f"ERROR: H37Rv reference missing at {REF}", file=sys.stderr); return 2
    for d in (ASM, VCF):
        d.mkdir(parents=True, exist_ok=True)
    from dna_decode.data import tb_who_catalogue as cat
    cat.verify_pins()
    dets = {"rifampicin": cat.load_determinants("rifampicin"), "isoniazid": cat.load_determinants("isoniazid")}

    cohort = read_cohort(a.max)
    done = load_ckpt()
    print(f"[tb-indep] cohort={len(cohort)} (max={a.max or 'full'}) | already scored={len(done)}")
    # 1. fetch
    label_by = {r["strain_id"]: r for r in cohort}
    n_fetch = 0
    for r in cohort:
        if r["strain_id"] in done:
            continue
        if fetch_assembly(r["assembly"], ASM / f"{r['strain_id']}.fna"):
            n_fetch += 1
    print(f"[tb-indep] assemblies present (fetched this run: {n_fetch})")
    # 2. align (one container)
    align_all()
    # 3. score + checkpoint
    n_new = 0
    with open(CKPT, "a", encoding="utf-8") as ck:
        for r in cohort:
            sid = r["strain_id"]
            if sid in done:
                continue
            sc = score_isolate(sid, dets)
            if sc is None:
                continue
            row = {"strain_id": sid, "rif_label": r["rif_label"], "inh_label": r["inh_label"], **sc}
            ck.write(json.dumps(row) + "\n"); ck.flush()
            done[sid] = row; n_new += 1
    print(f"[tb-indep] scored this run: {n_new} | total: {len(done)}")
    # 4. aggregate (only the cohort subset we're running)
    rows = [done[r["strain_id"]] for r in cohort if r["strain_id"] in done]
    agg = aggregate(rows)
    out = REPO / "wiki" / "tb_independent_amr_portal_scores.json"
    out.write_text(json.dumps({"n_isolates": len(rows), "drugs": agg}, indent=2, default=str), encoding="utf-8")
    for drug, m in agg.items():
        print(f"  {drug}: n_R={m['n_R']} n_S={m['n_S']} sens={m['sens']} spec={m['spec']} acc={m['accuracy']}")
        for k in ("sens", "spec"):
            if m[f"{k}_ci95"]:
                print(f"    {k} 95% CI: [{m[f'{k}_ci95'][0]:.3f}, {m[f'{k}_ci95'][1]:.3f}]")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
