"""gnomAD frequency-benign source — supply the missing BENIGN class for constrained clinical genes.

The clinical-gene census (2026-07-22) found 18/28 genes SINGLE_CLASS because curated ClinVar missense is
pathogenic-dominated for most disease genes (benign missense are rarely submitted). Population-genetics gives
an independent benign proxy: a missense observed in gnomAD ABOVE the maximum credible frequency for a
penetrant pathogenic allele is benign (the ACMG BS1/BA1 frequency principle). This fetches those.

**CRITICAL CIRCULARITY (load-bearing, must be surfaced):** AlphaMissense (Cheng 2023) was TRAINED using
population-common variants as its benign weak-label. So scoring AM's benign-calling on a gnomAD-frequency-benign
set is CIRCULAR — AM will trivially call them benign because that is what it was trained on. gnomAD-benign is
therefore a FAIR test ONLY for decoders that never saw population frequency: the DMS-itself ceiling (wet-lab),
ESM2 (self-supervised MLM, no benign/pathogenic labels), and BLOSUM62. It is NOT a fair test for AlphaMissense.
Consumers MUST flag AM as circular on this benign source and headline the DMS/ESM2/BLOSUM numbers.

Benign threshold: AF >= 1e-4 (0.01%) is the primary proxy — above the max credible frequency for most
penetrant pathogenic missense (ACMG-adjacent). Highly-constrained oncogenes (KRAS) may still lack enough
common variants even here — reported honestly, not forced.

  uv run python scripts/gnomad_benign.py --gene LDLR          # probe one gene
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_CACHE = Path("D:/dna_decode_cache/gnomad")
DEFAULT_AF_MIN = 1e-4  # ACMG-adjacent frequency-benign floor

_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
        "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
        "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}

_QUERY = ("query($g:String!){gene(gene_symbol:$g,reference_genome:GRCh38)"
          "{variants(dataset:gnomad_r4){consequence hgvsp genome{af} exome{af}}}}")


def parse_hgvsp(hgvsp: str) -> tuple[str, int, str] | None:
    """`p.Gly324Ser` -> ('G', 324, 'S'). None for synonymous/nonsense/non-single-missense. PURE."""
    if not hgvsp or not hgvsp.startswith("p."):
        return None
    b = hgvsp[2:].strip()
    if b[:3] not in _AA3:
        return None
    i = 3
    while i < len(b) and b[i].isdigit():
        i += 1
    if i == 3 or i >= len(b):
        return None
    pos = int(b[3:i])
    alt3 = b[i:]
    if alt3 in ("Ter", "=") or alt3 not in _AA3:
        return None
    return _AA3[b[:3]], pos, _AA3[alt3]


def _variant_af(v: dict) -> float:
    a = [(v.get("exome") or {}).get("af"), (v.get("genome") or {}).get("af")]
    a = [x for x in a if x]
    return max(a) if a else 0.0


def fetch_gnomad_benign(gene: str, af_min: float = DEFAULT_AF_MIN, use_cache: bool = True) -> dict[tuple[str, int, str], float]:
    """{(wt,pos,alt): af} for gnomAD missense with max(exome,genome) AF >= af_min. Cached per gene to D:.
    The RAW per-gene variant list is cached (all missense + AF) so a threshold change needs no re-fetch."""
    cache = GNOMAD_CACHE / f"{gene}.json"
    raw: list[dict] = []
    if use_cache and cache.exists():
        raw = json.loads(cache.read_text(encoding="utf-8"))
    else:
        body = json.dumps({"query": _QUERY, "variables": {"g": gene}}).encode()
        for attempt in range(5):
            try:
                req = urllib.request.Request(GNOMAD_API, data=body, headers={"Content-Type": "application/json"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = json.loads(r.read().decode())
                vs = (((data.get("data") or {}).get("gene") or {}).get("variants")) or []
                raw = [{"hgvsp": v.get("hgvsp"), "af": _variant_af(v)}
                       for v in vs if v.get("consequence") == "missense_variant" and v.get("hgvsp")]
                break
            except urllib.error.HTTPError as e:
                if e.code in (429, 500, 502, 503) and attempt < 4:
                    time.sleep(2 * (2 ** attempt))
                    continue
                raise
            except Exception:  # noqa: BLE001
                if attempt < 4:
                    time.sleep(2)
                    continue
                raise
        if use_cache:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(raw), encoding="utf-8")
    out: dict[tuple[str, int, str], float] = {}
    for v in raw:
        if v["af"] < af_min:
            continue
        p = parse_hgvsp(v["hgvsp"])
        if p:
            out[p] = max(out.get(p, 0.0), v["af"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--af-min", type=float, default=DEFAULT_AF_MIN)
    a = ap.parse_args()
    b = fetch_gnomad_benign(a.gene, a.af_min)
    print(f"{a.gene}: {len(b)} gnomAD frequency-benign missense (AF >= {a.af_min})")
    for k, v in sorted(b.items(), key=lambda kv: -kv[1])[:10]:
        print(f"  {k[0]}{k[1]}{k[2]}  AF={v:.5f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
