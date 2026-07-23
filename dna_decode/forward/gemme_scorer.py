"""Evolution modality for the forward hybrid — GEMME (a deterministic evolutionary-conservation model).

GEMME (Laine et al. 2019) scores a variant by how well it fits the evolutionary history in an MSA — a global
epistatic model with ZERO learned parameters. On ProteinGym it is a top-tier evolution predictor (~0.455)
and, combined with sequence (ESM2) + structure (ProSST), forms the sweep-top 3-way
(`wiki/forward_modality_hybrid_2026-07-17.md`: `ESM2+GEMME+ProSST` 0.547).

This module is the GEMME component of the hybrid, mirroring `prosst_scorer` / `structure_scorer`:
- `gemme_table_from_column(rows, col="GEMME")` — the DEPLOYABLE-NOW path: adapt a precomputed GEMME column
  (e.g. ProteinGym's `pg_zeroshot`) into a `{mutation: score}` table, oriented higher=preserved, ready for
  `variant_effect.rank_average_hybrid`. Because GEMME is DETERMINISTIC (0 params), a precomputed column IS
  canonical GEMME output — using it is the same move as using ProteinGym's pre-quantized ProSST structures.
- `run_gemme(msa_path, query_seq)` — the NOVEL-protein path: run GEMME locally. GEMME needs JET2 (Java) + R +
  an MSA; it is Windows-hostile (the same class as ESM-IF/ProSST's `torch_geometric` quantizer). So it
  LAZY-checks the toolchain and raises `GemmeUnavailable` when absent; run it on a Linux host (feed the MSA
  from `msa_fetch.fetch_msa`). Deferred — the seam is complete, the real run needs the toolchain.

Sign: GEMME's native score is a NEGATIVE conservation delta (more negative = more deleterious). ProteinGym's
column is already oriented higher=fitter, so the adapter passes it through; a raw GEMME run is oriented in
`run_gemme`.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

# GEMME ships a self-contained Docker image (JET2 + R + MUSCLE + python2.7 all bundled) -- the deployable
# runner path (www.lcqb.upmc.fr/GEMME/download.html). This host runs it via Docker Desktop (the same route as
# Mash/AMRFinder/Bakta in tools/docker_runner). `$GEMME_IMAGE` overrides the pinned tag.
GEMME_IMAGE = os.environ.get("GEMME_IMAGE", "elodielaine/gemme:gemme")

_AA = "ACDEFGHIKLMNPQRSTVWY"


class GemmeUnavailable(RuntimeError):
    """Raised when the GEMME toolchain (Docker image, or JET2/Java + R) is absent (this host)."""


# GEMME score tiers on the ProteinGym-oriented scale (higher = preserved); coarse, like the other modalities.
_GEMME_PRESERVED = -0.5
_GEMME_DAMAGING = -4.0


def gemme_table_from_column(rows, col: str = "GEMME", mutant_key: str = "mutant") -> dict[str, float]:
    """Adapt a precomputed GEMME column into {mutation: score} (higher=preserved). `rows` is an iterable of
    dict rows (e.g. csv.DictReader over ProteinGym's pg_zeroshot). Skips rows with a missing/NA GEMME cell.

    GEMME is deterministic (0 params) -> a precomputed column is canonical GEMME output, so this is the
    deployable evolution component for any protein ProteinGym (or a GEMME run) has already scored."""
    out: dict[str, float] = {}
    for r in rows:
        m = r.get(mutant_key)
        v = r.get(col)
        if not m or v in (None, "", "NA"):
            continue
        try:
            out[m] = float(v)
        except (TypeError, ValueError):
            continue
    if not out:
        raise ValueError(f"no usable {col!r} values in the supplied rows")
    return out


def find_jet2() -> str | None:
    """Resolve the JET2 launcher (GEMME's structure/MSA engine) — env override -> PATH."""
    import os
    env = os.environ.get("JET2_BIN")
    if env and Path(env).exists():
        return env
    return shutil.which("JET2") or shutil.which("jet2")


def a3m_to_aligned_fasta(a3m_path: str | Path, out_fasta: str | Path) -> tuple[str, int]:
    """Convert a ColabFold a3m MSA to a GEMME-ready aligned FASTA. GEMME needs every row the SAME length
    (query length); a3m encodes homolog INSERTS as lowercase relative to the query, so keeping only
    uppercase + '-' (the match columns) yields a query-length alignment. Query row goes FIRST (GEMME's
    focus). Returns (query_seq, n_sequences). PURE (file IO only)."""
    recs: list[tuple[str, str]] = []
    hdr, seq = None, []
    for ln in Path(a3m_path).read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if hdr is not None:
                recs.append((hdr, "".join(seq)))
            hdr, seq = ln[1:].strip() or f"seq{len(recs)}", []
        elif ln.strip():
            seq.append(ln.strip())
    if hdr is not None:
        recs.append((hdr, "".join(seq)))
    if not recs:
        raise ValueError(f"no sequences in {a3m_path}")

    def match_only(s: str) -> str:
        return "".join(c for c in s if c.isupper() or c == "-")

    query_match = match_only(recs[0][1])
    L = len(query_match)
    out_lines = []
    for i, (h, s) in enumerate(recs):
        m = match_only(s)
        if len(m) != L:               # ragged row (rare) -> pad/trim to query length, keeps GEMME happy
            m = (m + "-" * L)[:L]
        # GEMME/JET2 want clean headers; sanitize to alnum+underscore
        safe = "".join(ch if ch.isalnum() else "_" for ch in h)[:40] or f"seq{i}"
        out_lines.append(f">{safe}")
        out_lines.append(m)
    Path(out_fasta).write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return query_match.replace("-", ""), len(recs)


def _parse_gemme_matrix(pred_path: Path, query_seq: str) -> dict[str, float]:
    """Parse a GEMME prediction matrix (`*_normPred_evolCombi.txt`) -> {wt{pos}alt: score, higher=preserved}.

    GEMME's matrix is rows=the 20 amino acids, columns=query positions; a cell is the predicted effect of
    mutating that position TO that amino acid. GEMME's native scale is NEGATIVE = more deleterious, so it is
    ALREADY oriented higher=preserved (the ProteinGym-column convention). A '-inf'/'NA' cell (the WT itself,
    or an unscored position) is skipped."""
    lines = [ln.rstrip("\n") for ln in pred_path.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"empty GEMME output {pred_path}")
    # header row: quoted position labels (or blank first cell); data rows start with an aa label like "A"
    rows = [ln.split() for ln in lines]
    header = rows[0]
    # column j (0-based within the header, minus the leading label cell) -> query position (1-based)
    data = {}
    for r in rows[1:]:
        aa = r[0].strip('"').upper()
        if aa not in _AA:
            continue
        for j, cell in enumerate(r[1:]):
            pos = j + 1
            c = cell.strip('"')
            if c in ("NA", "-inf", "inf", "nan", ""):
                continue
            try:
                data[(aa, pos)] = float(c)
            except ValueError:
                continue
    out: dict[str, float] = {}
    for (alt, pos), val in data.items():
        if 1 <= pos <= len(query_seq):
            wt = query_seq[pos - 1].upper()
            if wt != alt and wt in _AA:
                out[f"{wt}{pos}{alt}"] = val
    if not out:
        raise ValueError(f"no variants parsed from GEMME matrix {pred_path} (header={header[:5]})")
    return out


def run_gemme(msa_a3m_or_fasta: str | Path, query_seq: str, *, work_dir: str | Path | None = None,
              image: str = GEMME_IMAGE, timeout: int = 3600) -> dict[str, float]:
    """Run GEMME on a real MSA -> {mutation: score} for every single substitution (higher=preserved).

    Deployable NOVEL-protein path via the self-contained GEMME Docker image (JET2 + R + MUSCLE + python2.7).
    Feed an a3m (from `msa_fetch.fetch_msa`, auto-converted) or an already-aligned FASTA. Runs
    `python2.7 $GEMME_PATH/gemme.py ali.fasta -r input -f ali.fasta` in the image with the work dir bind-
    mounted at /project, then parses the normPred matrix. Raises GemmeUnavailable if Docker/the image is
    absent (use `gemme_table_from_column` for the precomputed/validation path)."""
    if shutil.which("docker") is None:
        raise GemmeUnavailable("GEMME needs Docker (the elodielaine/gemme image) or a native JET2+R toolchain.")
    src = Path(msa_a3m_or_fasta)
    wd = Path(work_dir) if work_dir else src.parent
    wd.mkdir(parents=True, exist_ok=True)
    ali = wd / "gemme_ali.fasta"
    if src.suffix.lower() in (".a3m", ".a2m"):
        a3m_to_aligned_fasta(src, ali)
    else:
        ali.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    # MSYS_NO_PATHCONV so Git Bash does not mangle the container-side /project path.
    env = {**os.environ, "MSYS_NO_PATHCONV": "1"}
    cmd = ["docker", "run", "--rm", "--mount",
           f"type=bind,source={wd.resolve()},target=/project", image, "bash", "-lc",
           "cd /project && python2.7 $GEMME_PATH/gemme.py gemme_ali.fasta -r input -f gemme_ali.fasta"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    except FileNotFoundError as e:
        raise GemmeUnavailable(f"docker not runnable: {e}") from e
    if proc.returncode != 0:
        raise GemmeUnavailable(f"GEMME run failed (rc={proc.returncode}): {proc.stderr[-800:]}")
    preds = sorted(wd.glob("*normPred_evolCombi.txt")) or sorted(wd.glob("*normPred*.txt")) or sorted(wd.glob("*pred*.txt"))
    if not preds:
        raise GemmeUnavailable(f"GEMME produced no prediction matrix in {wd} (stdout tail: {proc.stdout[-400:]})")
    return _parse_gemme_matrix(preds[0], query_seq)


def gemme_tier(score: float) -> str:
    """GEMME (ProteinGym-oriented, higher=preserved) -> forward-cell tier."""
    if score >= _GEMME_PRESERVED:
        return "preserved"
    if score <= _GEMME_DAMAGING:
        return "damaging"
    return "uncertain"
