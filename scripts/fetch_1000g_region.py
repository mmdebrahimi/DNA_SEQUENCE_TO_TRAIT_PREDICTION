#!/usr/bin/env python
"""Docker-free 1000G-region VCF fetcher — pure-Python tabix(.tbi)-over-HTTP-range + BGZF block decode.

The prior PGx pipeline sliced 1000G regions with Docker bcftools (htslib is unavailable natively on this
Windows host). Docker is flaky here (WSL mount corruption on churn). This does the same slice with ZERO
external tools: parse the remote `.tbi` tabix index, HTTP-range-fetch only the compressed BGZF blocks that
cover the target region + the header, decode them in pure Python, and write a plain single-region VCF the
existing `scripts/pgx_getrm_concordance.py` harness consumes. Reusable for ALL future PGx genes ("~13 more
PGx genes for free" — wiki/data_source_research_human_animal_2026-06-25).

Usage:
    uv run python scripts/fetch_1000g_region.py --chrom chr10 --start 95036000 --end 95069000 \
        --out data/pgx_1000g/cyp2c8_1000g.vcf
"""
from __future__ import annotations

import argparse
import gzip
import os
import struct
import subprocess
import sys
import tempfile
import zlib
from pathlib import Path

# Canonical 1000G 30x phased panel (NYGC/EBI). Per-chrom bgzipped + tabixed.
BASE = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
        "working/20220422_3202_phased_SNV_INDEL_SV")
VCF_TMPL = BASE + "/1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def _http_range(url: str, start: int, end: int | None = None, timeout: int = 90) -> bytes:
    """HTTP Range GET [start, end] inclusive (end=None -> open-ended). Via curl.

    curl is used instead of urllib: on this (Windows) host urllib intermittently HANGS on the large 1000G
    chr VCFs while curl fetches the same range in <1s (empirically verified 2026-07-05). curl is a hard
    dependency of the fetch path (already required for the Docker-free posture)."""
    rng = f"{start}-" + ("" if end is None else str(end))
    last = None
    for _ in range(4):
        fd, tmp = tempfile.mkstemp(suffix=".bin")
        os.close(fd)
        try:
            r = subprocess.run(["curl", "-sSL", "--max-time", str(timeout), "-r", rng, url, "-o", tmp],
                               capture_output=True, text=True)
            data = Path(tmp).read_bytes()
            if r.returncode == 0 and data:
                return data
            last = RuntimeError(f"curl rc={r.returncode} bytes={len(data)}: {r.stderr[:160]}")
        except Exception as e:  # noqa: BLE001 — retry transient failures
            last = e
        finally:
            try:
                os.unlink(tmp)
            except OSError:
                pass
    raise last or RuntimeError("curl range fetch failed")


def _http_get(url: str, timeout: int = 90) -> bytes:
    return _http_range(url, 0, None, timeout)


# ---- tabix (.tbi) index parse (little-endian; spec: samtools tabix / SAMv1 §5.2) ----
def _reg2bins(beg: int, end: int) -> list[int]:
    """UCSC binning: bins that may overlap [beg, end) (0-based)."""
    end -= 1
    out = [0]
    for start, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        out += list(range(start + (beg >> shift), start + (end >> shift) + 1))
    return out


def _parse_tbi(tbi_gz: bytes):
    """Return (names->ref_id, per-ref {bin->[(cbeg,cend)]}, per-ref linear-index list)."""
    buf = gzip.decompress(tbi_gz)  # .tbi is a BGZF (gzip-compatible) stream
    if buf[:4] != b"TBI\x01":
        raise ValueError("not a TBI index")
    off = 4
    (n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm) = struct.unpack_from("<8i", buf, off)
    off += 32
    names_blob = buf[off:off + l_nm]; off += l_nm
    names = [n.decode() for n in names_blob.split(b"\x00") if n]
    name2id = {n: i for i, n in enumerate(names)}
    ref_bins = []
    ref_linear = []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", buf, off); off += 4
        bins: dict[int, list[tuple[int, int]]] = {}
        for _b in range(n_bin):
            (bin_id, n_chunk) = struct.unpack_from("<Ii", buf, off); off += 8
            chunks = []
            for _c in range(n_chunk):
                (cbeg, cend) = struct.unpack_from("<QQ", buf, off); off += 16
                chunks.append((cbeg, cend))
            bins[bin_id] = chunks
        (n_intv,) = struct.unpack_from("<i", buf, off); off += 4
        linear = list(struct.unpack_from("<%dQ" % n_intv, buf, off)) if n_intv else []
        off += 8 * n_intv
        ref_bins.append(bins)
        ref_linear.append(linear)
    return name2id, ref_bins, ref_linear


def _bgzf_blocks(data: bytes):
    """Yield decompressed payloads of every COMPLETE BGZF block at the front of `data`."""
    i = 0
    n = len(data)
    while i + 18 <= n:
        if data[i:i + 4] != b"\x1f\x8b\x08\x04":
            break
        xlen = struct.unpack_from("<H", data, i + 10)[0]
        # find the BC subfield (SI1='B'(66) SI2='C'(67)) -> BSIZE
        extra = data[i + 12:i + 12 + xlen]
        bsize = None
        j = 0
        while j + 4 <= len(extra):
            si1, si2, slen = extra[j], extra[j + 1], struct.unpack_from("<H", extra, j + 2)[0]
            if si1 == 66 and si2 == 67:
                bsize = struct.unpack_from("<H", extra, j + 4)[0]
                break
            j += 4 + slen
        if bsize is None:
            break
        block_len = bsize + 1
        if i + block_len > n:
            break  # incomplete trailing block
        block = data[i:i + block_len]
        try:
            yield gzip.decompress(block)
        except (OSError, EOFError, zlib.error):
            break
        i += block_len


def _voff_coffset(voff: int) -> int:
    return voff >> 16


def fetch_region(chrom: str, start: int, end: int, out: Path, *, verbose: bool = True) -> dict:
    url = VCF_TMPL.format(chrom=chrom)
    if verbose:
        print(f"[fetch] {chrom}:{start}-{end}  <- {url}", flush=True)
    # 1) header: first chunk of the file holds ## meta + the #CHROM sample line. 400KB is ample even at
    # 3202 samples (~26KB #CHROM line, within the first ~2 bgzf blocks); a smaller GET is far less
    # stall-prone than a 3MB one under EBI throttling.
    head_raw = _http_range(url, 0, 400_000, timeout=60)
    header_lines: list[str] = []
    done_header = False
    carry = ""
    for payload in _bgzf_blocks(head_raw):
        carry += payload.decode("utf-8", "replace")
        while "\n" in carry:
            line, carry = carry.split("\n", 1)
            if line.startswith("#"):
                header_lines.append(line)
                if line.startswith("#CHROM"):
                    done_header = True
                    break
            else:
                done_header = True
                break
        if done_header:
            break
    if not header_lines or not header_lines[-1].startswith("#CHROM"):
        raise RuntimeError("failed to recover #CHROM header from the first 3MB")

    # 2) region chunks via the tabix index
    tbi = _http_get(url + ".tbi")
    name2id, ref_bins, ref_linear = _parse_tbi(tbi)
    if chrom not in name2id:
        # some indexes drop the 'chr' prefix
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt not in name2id:
            raise RuntimeError(f"{chrom} not in tabix names: {list(name2id)[:5]}...")
        chrom_key = alt
    else:
        chrom_key = chrom
    rid = name2id[chrom_key]
    bins = ref_bins[rid]
    linear = ref_linear[rid]
    lin_min = linear[start >> 14] if (start >> 14) < len(linear) else (linear[-1] if linear else 0)
    chunks = []
    for b in _reg2bins(start, end + 1):
        for cbeg, cend in bins.get(b, []):
            if cend >= lin_min:
                chunks.append((cbeg, cend))
    if not chunks:
        raise RuntimeError("no tabix chunks overlap the region (wrong coords/assembly?)")
    cmin = min(_voff_coffset(c[0]) for c in chunks)
    cmax = max(_voff_coffset(c[1]) for c in chunks)
    # +65536 to include the final block that cmax points into
    region_raw = _http_range(url, cmin, cmax + 65536, timeout=90)
    if verbose:
        print(f"[fetch] header {len(header_lines)} lines; region bytes {cmin}-{cmax} "
              f"({len(region_raw)//1024} KB compressed)", flush=True)

    # 3) decode region blocks, keep records in [start, end]
    data_lines: list[str] = []
    carry = ""
    for payload in _bgzf_blocks(region_raw):
        carry += payload.decode("utf-8", "replace")
        while "\n" in carry:
            line, carry = carry.split("\n", 1)
            if not line or line.startswith("#"):
                continue
            tab = line.find("\t")
            tab2 = line.find("\t", tab + 1)
            try:
                pos = int(line[tab + 1:tab2])
            except ValueError:
                continue
            if start <= pos <= end:
                data_lines.append(line)

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(header_lines) + "\n" + "\n".join(data_lines) + "\n", encoding="utf-8")
    return {"records": len(data_lines), "samples": len(header_lines[-1].split("\t")) - 9,
            "out": str(out), "region": f"{chrom}:{start}-{end}"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Docker-free 1000G region VCF slice (tabix-over-HTTP).")
    ap.add_argument("--chrom", required=True, help="e.g. chr10")
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    try:
        rep = fetch_region(args.chrom, args.start, args.end, args.out)
    except Exception as e:  # noqa: BLE001
        print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
        return 1
    print(f"[done] {rep['records']} records x {rep['samples']} samples -> {rep['out']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
