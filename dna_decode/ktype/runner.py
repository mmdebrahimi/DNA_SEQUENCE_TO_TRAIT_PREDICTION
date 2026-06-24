"""Klebsiella K-antigen (capsule) typing via wzi -- the serotype sibling.

`wzi` is a conserved capsule-locus (cps) gene whose allele predicts the Klebsiella capsule (K) type. This
is the SAME shape as the E. coli serotype caller: blastn a curated allele DB vs an assembly -> best wzi
allele -> the K-locus (KL) type via the wzi->KL map. DB = the BIGSdb Pasteur wzi scheme as bundled by
Kleborate (`wzi.fasta` alleles + `wzi.txt` ST/wzi/KL map).

HONESTY (load-bearing): the wzi->K-type relationship is NOT one-to-one (~94% predictive; isolates with
distinct K-types can share a wzi allele -- Brisse 2013 JCM). So the call is a PREDICTED KL type with that
ceiling, FAITHFUL-TO-TOOL (the Kleborate/BIGSdb wzi method), NOT an independent baseline. Full-locus typing
(Kaptive) is more accurate; this single-gene v0 is the "smallest credible slice" (serotype pattern).
Offline-safe via the shared blast engine (missing blastn/DB -> status 'unavailable').
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.typing.blast_caller import call_alleles

WZI_IDENTITY_THRESHOLD = 90.0    # wzi typing is near-exact allele matching
WZI_COVERAGE_THRESHOLD = 80.0


def load_wzi_kl_map(wzi_txt: str | Path) -> dict[str, str]:
    """Parse wzi.txt (ST<tab>wzi<tab>KL) -> {wzi_number: primary_KL}. 'KL2 (KL30)' -> 'KL2' (alternatives in
    parens are dropped for the primary call)."""
    m: dict[str, str] = {}
    for line in Path(wzi_txt).read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 3 or parts[1].strip().lower() == "wzi":   # skip header
            continue
        wzi, kl = parts[1].strip(), parts[2].strip()
        if wzi and kl:
            m[wzi] = kl.split()[0]
    return m


def _wzi_num(allele_id: str) -> str | None:
    """'wzi_1' -> '1'."""
    return allele_id.split("_", 1)[1] if (allele_id.lower().startswith("wzi_") and "_" in allele_id) else None


def call_ktype(fasta: str | Path, db_dir: str | Path, *,
               identity_threshold: float = WZI_IDENTITY_THRESHOLD,
               coverage_threshold: float = WZI_COVERAGE_THRESHOLD,
               blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn the wzi allele DB vs `fasta`; best wzi allele -> predicted Klebsiella K-locus (KL) type."""
    db_dir = Path(db_dir)
    wzi_fasta, wzi_txt = db_dir / "wzi.fasta", db_dir / "wzi.txt"
    if not wzi_fasta.exists() or not wzi_txt.exists():
        return {"status": "unavailable", "wzi_allele": None, "kl_type": None,
                "reason": f"wzi DB not found in {db_dir} (need wzi.fasta + wzi.txt)"}
    res = call_alleles(fasta, wzi_fasta, identity_threshold=identity_threshold,
                       coverage_threshold=coverage_threshold, blastn_bin=blastn_bin, timeout=timeout)
    if res["status"] != "ok":
        return {"status": "unavailable", "tool": res.get("tool"), "wzi_allele": None, "kl_type": None,
                "reason": res.get("reason")}
    kl_map = load_wzi_kl_map(wzi_txt)
    best = None  # (allele_id, (identity, coverage), hit) -- best CALLED wzi allele
    for aid, hit in res["per_allele"].items():
        if not hit["called"]:
            continue
        key = (hit["percent_identity"], hit["percent_coverage"])
        if best is None or key > best[1]:
            best = (aid, key, hit)
    base = {"status": "ok", "tool": "blastn", "method": "wzi_blastn_v0",
            "parameters": {"identity_threshold": identity_threshold, "coverage_threshold": coverage_threshold}}
    if best is None:
        return {**base, "wzi_allele": None, "kl_type": None, "predicted_k": None,
                "note": "no wzi allele called at threshold (novel/partial/absent wzi locus)"}
    aid, _, hit = best
    kl = kl_map.get(_wzi_num(aid) or "")
    return {**base, "wzi_allele": aid, "kl_type": kl, "predicted_k": kl,
            "percent_identity": hit["percent_identity"], "percent_coverage": hit["percent_coverage"]}
