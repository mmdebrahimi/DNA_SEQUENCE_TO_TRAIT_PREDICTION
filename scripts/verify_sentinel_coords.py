"""verify_sentinel_coords -- the EXECUTABLE anti-fabrication rail for PGx SENTINELS population.

A wrong GRCh38 coordinate on a SentinelVariant is invisible: it produces a silent no-op sentinel (never
fires, the non-core leak persists) or a false-withhold (fires on the wrong site). This tool makes the
"source coords verbatim, never fabricate" discipline MACHINE-CHECKED: for each proposed sentinel it
cross-checks (rsid -> GRCh38 chrom:pos, ref) against Ensembl REST and FAILS CLOSED on any mismatch.

Policy (grounded in the Ensembl variation endpoint contract):
- Query GET /variation/human/<rsid> (JSON). Consider ONLY `mappings` with:
    assembly_name == "GRCh38"  AND  coord_system == "chromosome"  AND  seq_region_name == expected chrom.
- MULTIPLE qualifying chromosome mappings  -> FAIL (multi-mapping / ambiguous; do not guess).
- ZERO qualifying mappings                 -> FAIL (no canonical GRCh38 chromosome mapping).
- start != expected pos                     -> FAIL (coordinate mismatch).
- expected ref not consistent with the mapping allele_string -> FAIL (ref-base mismatch).
- expected alt given but not among the mapping's ALT alleles -> WARN (Ensembl allele_string may be
    strand-oriented / not list every ALT); reported, not a hard fail (the caller's "*" wildcard tolerates it).
- rsID merges/synonyms: Ensembl returns the CURRENT `name`; if it differs from the queried rsid, the queried
    id is a MERGED/synonym id -> reported as `merged_into` (not a fail, but surfaced so the catalog can be
    updated to the current id).
- Network unreachable / non-200 / parse error -> status UNVERIFIED (LOUD), never a silent pass.

Offline-safe + testable: the network call is injected (`fetch=`), so unit tests run with a mock and never
touch the network. CLI does a real fetch. NON-frozen tooling; reads catalogs read-only, edits nothing.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ENSEMBL = "https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"


def _norm_chrom(c: str) -> str:
    return c[3:] if str(c).lower().startswith("chr") else str(c)


def _http_fetch(rsid: str, timeout: float = 10.0) -> dict:
    """Real Ensembl REST fetch. Raises urllib errors on network/HTTP failure (caught by verify_sentinel)."""
    req = urllib.request.Request(ENSEMBL.format(rsid=rsid), headers={"User-Agent": "dna_decode-pgx-verify"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


@dataclass
class VerifyResult:
    rsid: str
    status: str                 # "OK" | "MISMATCH" | "UNVERIFIED"
    detail: str
    merged_into: str | None = None
    alt_warning: bool = False

    @property
    def ok(self) -> bool:
        return self.status == "OK"


def verify_sentinel(rsid: str, chrom: str, pos: int, ref: str, alt: str | None = None, *,
                    fetch=_http_fetch) -> VerifyResult:
    """Verify one (rsid -> GRCh38 chrom:pos, ref[, alt]) against Ensembl. fetch() is injected for tests."""
    want_chrom = _norm_chrom(chrom)
    try:
        data = fetch(rsid)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError, ValueError) as e:
        return VerifyResult(rsid, "UNVERIFIED", f"Ensembl fetch failed ({type(e).__name__}): {e}")

    if not isinstance(data, dict) or "mappings" not in data:
        return VerifyResult(rsid, "UNVERIFIED", "no 'mappings' in Ensembl response")

    merged = None
    current = data.get("name")
    if current and current != rsid:
        merged = current   # queried id is a merged/synonym -> surface the current id

    quals = [m for m in data.get("mappings", [])
             if m.get("assembly_name") == "GRCh38" and m.get("coord_system") == "chromosome"
             and _norm_chrom(m.get("seq_region_name", "")) == want_chrom]
    if not quals:
        return VerifyResult(rsid, "MISMATCH",
                            f"no GRCh38 chromosome mapping on chr{want_chrom}", merged_into=merged)
    if len(quals) > 1:
        return VerifyResult(rsid, "MISMATCH",
                            f"{len(quals)} GRCh38 chr{want_chrom} mappings (ambiguous multi-map)",
                            merged_into=merged)
    m = quals[0]
    if int(m.get("start", -1)) != int(pos):
        return VerifyResult(rsid, "MISMATCH",
                            f"pos {m.get('start')} != expected {pos}", merged_into=merged)

    allele_string = m.get("allele_string", "") or ""
    alleles = allele_string.split("/")
    if ref and alleles and alleles[0] != ref:
        return VerifyResult(rsid, "MISMATCH",
                            f"ref {alleles[0]!r} != expected {ref!r} (allele_string {allele_string!r})",
                            merged_into=merged)
    alt_warning = bool(alt) and alt != "*" and alt not in alleles[1:]
    return VerifyResult(rsid, "OK",
                        f"GRCh38 chr{want_chrom}:{pos} {allele_string}", merged_into=merged,
                        alt_warning=alt_warning)


def verify_catalog(sentinels, *, fetch=_http_fetch) -> list[VerifyResult]:
    """Verify a whole SENTINELS list; returns one VerifyResult per entry."""
    return [verify_sentinel(s.rsid, s.chrom, s.pos, s.ref, getattr(s, "alt", None), fetch=fetch)
            for s in sentinels]


def _load_gene_sentinels(gene: str):
    import importlib
    mod = importlib.import_module(f"dna_decode.pgx.{gene}_catalog")
    return getattr(mod, "SENTINELS", [])


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(prog="verify_sentinel_coords",
                                 description="Machine-check PGx SENTINELS coords against Ensembl GRCh38 "
                                             "(anti-fabrication rail; fail-closed).")
    ap.add_argument("--rsid"); ap.add_argument("--chrom"); ap.add_argument("--pos", type=int)
    ap.add_argument("--ref"); ap.add_argument("--alt")
    ap.add_argument("--catalog", help="verify a gene's whole SENTINELS list (e.g. cyp2c9, tpmt)")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    results: list[VerifyResult] = []
    if args.catalog:
        results = verify_catalog(_load_gene_sentinels(args.catalog))
    elif args.rsid and args.chrom and args.pos is not None and args.ref:
        results = [verify_sentinel(args.rsid, args.chrom, args.pos, args.ref, args.alt)]
    else:
        print("error: give --catalog <gene> OR --rsid --chrom --pos --ref [--alt]", file=sys.stderr)
        return 2

    out = [{"rsid": r.rsid, "status": r.status, "detail": r.detail,
            "merged_into": r.merged_into, "alt_warning": r.alt_warning} for r in results]
    if args.json:
        print(json.dumps(out, indent=2))
    else:
        for r in results:
            tag = {"OK": "OK  ", "MISMATCH": "FAIL", "UNVERIFIED": "????"}[r.status]
            extra = (f"  [merged_into {r.merged_into}]" if r.merged_into else "") + \
                    ("  [ALT not in allele_string -- confirm strand/rep]" if r.alt_warning else "")
            print(f"[{tag}] {r.rsid}: {r.detail}{extra}")
    # exit 1 on any MISMATCH; 3 if any UNVERIFIED (and none mismatched); 0 all OK.
    if any(r.status == "MISMATCH" for r in results):
        return 1
    if any(r.status == "UNVERIFIED" for r in results):
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
