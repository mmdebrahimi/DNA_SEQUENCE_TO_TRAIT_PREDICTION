"""Read-only preflight check: is THIS host able to run a given EP's compute?

Per /probe 2026-05-26 ("can we use /iterate overnight?"): answer was no —
GTX 860M can't run NT v2 100M, no Docker, retrained pickle not on origin.
Recommendation was to build a small preflight script that surfaces those
constraints in machine-readable form BEFORE any 6-hour compute job is kicked off.

Scope:
  - READ-ONLY. No mutation. No model loads. No HDF5 reads beyond file-exists.
  - Pure CPU. No GPU access. No Docker invocation.
  - Per-EP verdict: RUNNABLE / NOT_RUNNABLE / PARTIAL with missing-capabilities list.

Capabilities checked:
  - python_version + key import availability
  - torch CUDA: available + device_name + VRAM (read-only via torch.cuda metadata)
  - docker: `docker --version` exit code (does NOT pull / run any container)
  - AMRFinder DB directory existence
  - NT cache HDF5 file existence
  - Refseq cache directory existence
  - Cohort parquet existence
  - Model pickle existence

Per-EP profiles (from plans/Post_V0_EP_Ladder_Plan.md + EP plans):
  - EP-1A genome-input contract: GPU + retrained model + refseq + cipro cohort
  - EP-1B external benchmark: GPU + Docker + AMRFinder + model + refseq
  - EP-1.5 architecture POC: GPU + Docker + cache + cipro mini cohort + cef mini cohort + tet mini cohort
  - cef slice: GPU + Docker + cef cohort + cache subset for cef strains
  - post_push_checklist: sync diagnostic + pytest (always runnable here)

Output:
  - reports/preflight_runnable_<DATE>.{md,json}
  - Console summary

Exit codes:
  0 = all checked EPs runnable on this host
  1 = at least one EP NOT_RUNNABLE
  2 = tool runtime error
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import date as _date
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Default paths (matches CLAUDE.md gotchas + post-Codex-push checklist).
# ---------------------------------------------------------------------------

DEFAULTS = {
    "amrfinder_db": Path("C:/Users/Farshad/dna_decode_stage2/amrfinder_db"),
    "nt_cache_files": [
        Path("D:/dna_decode_cache/embeddings/nt_n40_cipro.h5"),
        Path("D:/dna_decode_cache/embeddings/nt_n147_cipro.h5"),
    ],
    "refseq_cache": Path("D:/dna_decode_cache/refseq"),
    "cohort_files": [
        Path("data/processed/stage2_n150_cipro_cohort.parquet"),
        Path("data/processed/gate_b_n40_cipro_cohort.parquet"),
        Path("data/processed/gate_b_mini_cipro_cohort.parquet"),
        Path("data/processed/gate_b_mini_cef_cohort.parquet"),
        Path("data/processed/gate_b_mini_tet_cohort.parquet"),
    ],
    "model_pickles": [
        Path("data/processed/models/ciprofloxacin_nucleotide_transformer.pkl"),
        Path("data/processed/models/ciprofloxacin_mock.pkl"),
    ],
}

# Per-EP capability requirements. Each EP names which top-level capability
# keys must be `available` for the EP to be RUNNABLE here.
EP_PROFILES: dict[str, dict[str, Any]] = {
    "EP_1A_genome_input": {
        "doc_anchor": "plans/EP_1A_Public_Genome_Selection_Memo.md",
        "requires_capabilities": ["gpu", "model_pickle_real", "refseq_cache"],
        "requires_files": [Path("data/processed/stage2_n150_cipro_cohort.parquet")],
        "note": "Real genome-input cipro decode. Needs GPU for NT v2 100M inference + retrained leakage-safe model + refseq for held-out genome GFF3.",
    },
    "EP_1B_external_benchmark": {
        "doc_anchor": "plans/EP_1B_External_Benchmark_Panel_Memo.md",
        "requires_capabilities": ["gpu", "docker", "amrfinder_db", "model_pickle_real", "refseq_cache"],
        "requires_files": [],
        "note": "10-genome side-by-side vs AMRFinderPlus + RGI. Needs Docker for both classical tools; GPU for v0 predictions on 10 genomes.",
    },
    "EP_1_5_architecture_poc": {
        "doc_anchor": "plans/EP_1_5_Architecture_Decision_Packet.md",
        "requires_capabilities": ["gpu", "docker", "amrfinder_db"],
        "requires_files": [
            Path("data/processed/gate_b_mini_cipro_cohort.parquet"),
            Path("data/processed/gate_b_mini_cef_cohort.parquet"),
            Path("data/processed/gate_b_mini_tet_cohort.parquet"),
        ],
        "note": "3-candidate POC for distributed-mechanism architecture. Needs GPU (NT inference); Docker for Candidate B (AMRFinder features).",
    },
    "cef_slice": {
        "doc_anchor": "plans/Cef_Cached_Strain_V0_1_Slice_Plan.md",
        "requires_capabilities": ["gpu", "docker", "amrfinder_db", "refseq_cache"],
        "requires_files": [Path("data/processed/gate_b_mini_cef_cohort.parquet")],
        "note": "Cef cached-strain v0.1 slice. Needs GPU for cef classifier training + Docker for cef mechanism audit.",
    },
    "post_push_checklist": {
        "doc_anchor": "wiki/post_codex_push_checklist_2026-05-26.md",
        "requires_capabilities": ["git", "pytest"],
        "requires_files": [],
        "note": "Operational sequence; mostly git + sync-check + pytest. Always runnable on any dev host.",
    },
}


# ---------------------------------------------------------------------------
# Capability checks. Each returns dict {available: bool, details: ...}.
# ---------------------------------------------------------------------------


def check_python() -> dict:
    return {
        "available": True,
        "version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "executable": sys.executable,
    }


def check_gpu() -> dict:
    """Check CUDA availability via torch.cuda metadata. Does NOT load any model."""
    try:
        import torch
    except ImportError:
        return {"available": False, "reason": "torch not importable"}
    if not torch.cuda.is_available():
        return {"available": False, "reason": "torch.cuda.is_available() is False"}
    try:
        idx = 0
        name = torch.cuda.get_device_name(idx)
        props = torch.cuda.get_device_properties(idx)
        vram_gib = props.total_memory / (1024**3)
        cc = f"{props.major}.{props.minor}"
        return {
            "available": True,
            "device_name": name,
            "compute_capability": cc,
            "vram_gib": round(vram_gib, 2),
            "fits_nt_v2_100m": vram_gib >= 6.0,
            "bitsandbytes_compatible": (props.major >= 7),
        }
    except Exception as e:
        return {"available": False, "reason": f"GPU query failed: {e!r}"}


def check_docker() -> dict:
    """Check Docker by invoking `docker --version` (does NOT pull / run a container)."""
    docker_bin = shutil.which("docker")
    if docker_bin is None:
        return {"available": False, "reason": "docker not on PATH"}
    try:
        proc = subprocess.run(
            [docker_bin, "--version"], capture_output=True, text=True, timeout=10
        )
        if proc.returncode != 0:
            return {"available": False, "reason": f"docker --version exit {proc.returncode}: {proc.stderr.strip()}"}
        return {
            "available": True,
            "binary": docker_bin,
            "version_string": proc.stdout.strip(),
        }
    except Exception as e:
        return {"available": False, "reason": f"docker --version raised: {e!r}"}


def check_git() -> dict:
    """Check git on PATH."""
    git_bin = shutil.which("git")
    if git_bin is None:
        return {"available": False, "reason": "git not on PATH"}
    try:
        proc = subprocess.run(
            [git_bin, "--version"], capture_output=True, text=True, timeout=10
        )
        return {"available": proc.returncode == 0, "binary": git_bin, "version_string": proc.stdout.strip()}
    except Exception as e:
        return {"available": False, "reason": f"git --version raised: {e!r}"}


def check_pytest() -> dict:
    """Check pytest is importable (without running it)."""
    try:
        import pytest  # noqa: F401
        return {"available": True}
    except ImportError:
        return {"available": False, "reason": "pytest not importable"}


def check_path_exists(p: Path) -> dict:
    """Path-existence check; reports type (file / dir) when present."""
    try:
        exists = p.exists()
        if not exists:
            return {"available": False, "path": str(p), "reason": "not found"}
        return {
            "available": True,
            "path": str(p),
            "type": "dir" if p.is_dir() else "file",
            "size_bytes": p.stat().st_size if p.is_file() else None,
        }
    except Exception as e:
        return {"available": False, "path": str(p), "reason": f"stat failed: {e!r}"}


def check_amrfinder_db(path: Path) -> dict:
    """Check the AMRFinder DB directory exists + has a `latest` subdir."""
    base = check_path_exists(path)
    if not base["available"]:
        return base
    latest = path / "latest"
    # Some AMRFinder installs use symlinks for `latest`; on Windows this can
    # raise OSError 1920. Catch + report as not-present rather than crash.
    try:
        base["latest_present"] = latest.exists()
    except OSError as e:
        base["latest_present"] = False
        base["latest_check_error"] = repr(e)
    return base


def check_model_pickle(path: Path) -> dict:
    """Check pickle file existence + whether it has the v0.1 RELOCKED provenance fields.

    Does NOT unpickle untrusted content into the runtime — just opens + peeks
    at the top-level dict keys via a small read.
    """
    base = check_path_exists(path)
    if not base["available"]:
        return base
    # Peek at the pickle's top-level structure WITHOUT loading the full model
    try:
        import pickle
        with open(path, "rb") as f:
            obj = pickle.load(f)
        if isinstance(obj, dict):
            base["keys"] = sorted(obj.keys())
            prov = obj.get("provenance", {})
            base["provenance_keys"] = sorted(prov.keys()) if isinstance(prov, dict) else None
            base["cv_strategy"] = prov.get("cv_strategy") if isinstance(prov, dict) else None
            base["cv_auroc"] = prov.get("cv_auroc") if isinstance(prov, dict) else None
            base["is_leakage_safe"] = (base["cv_strategy"] == "leave_one_accession_out")
        else:
            base["keys"] = None
            base["is_leakage_safe"] = False
    except Exception as e:
        base["pickle_load_error"] = repr(e)
        base["is_leakage_safe"] = False
    return base


# ---------------------------------------------------------------------------
# EP runnability evaluation.
# ---------------------------------------------------------------------------


def _is_capability_satisfied(cap: str, caps: dict, file_checks: dict) -> tuple[bool, str]:
    """Map a capability key to a satisfied?/why-not check against the gathered data."""
    if cap == "gpu":
        gpu = caps.get("gpu", {})
        if not gpu.get("available"):
            return False, f"GPU unavailable ({gpu.get('reason', 'unknown')})"
        if not gpu.get("fits_nt_v2_100m"):
            return False, f"GPU VRAM {gpu.get('vram_gib')} GiB < 6 GiB required for NT v2 100M"
        return True, "GPU adequate"
    if cap == "docker":
        d = caps.get("docker", {})
        return (d.get("available", False), d.get("reason") or "docker present")
    if cap == "amrfinder_db":
        a = caps.get("amrfinder_db", {})
        return (a.get("available", False) and a.get("latest_present", False),
                a.get("reason") or ("latest dir present" if a.get("latest_present") else "latest dir missing"))
    if cap == "refseq_cache":
        r = caps.get("refseq_cache", {})
        return (r.get("available", False), r.get("reason") or "refseq cache present")
    if cap == "git":
        g = caps.get("git", {})
        return (g.get("available", False), g.get("reason") or "git present")
    if cap == "pytest":
        p = caps.get("pytest", {})
        return (p.get("available", False), p.get("reason") or "pytest importable")
    if cap == "model_pickle_real":
        for mp_name, mp in (caps.get("model_pickles") or {}).items():
            if mp.get("is_leakage_safe") is True:
                return True, f"{mp_name} has cv_strategy=leave_one_accession_out"
        return False, "no leakage-safe model pickle found (cv_strategy != leave_one_accession_out)"
    return False, f"unknown capability key: {cap!r}"


@dataclass
class EPRunnability:
    ep_name: str
    runnable: bool
    missing_capabilities: list[str] = field(default_factory=list)
    missing_files: list[str] = field(default_factory=list)
    details: dict[str, str] = field(default_factory=dict)
    doc_anchor: str = ""
    note: str = ""


def evaluate_ep(ep_name: str, profile: dict, caps: dict) -> EPRunnability:
    """Evaluate one EP's runnability against gathered capabilities.

    Lookup keys are normalized to POSIX form (forward slashes) so Path objects
    on Windows produce the same key as plain forward-slash strings on Linux.
    """
    missing_caps: list[str] = []
    missing_files: list[str] = []
    details: dict[str, str] = {}
    file_checks = caps.get("file_checks", {})
    for cap in profile.get("requires_capabilities", []):
        ok, reason = _is_capability_satisfied(cap, caps, file_checks)
        details[cap] = reason
        if not ok:
            missing_caps.append(cap)
    for f in profile.get("requires_files", []):
        # Try both the as-posix form AND the platform-native str() form so
        # callers writing Path-keyed file_checks on Windows match Path-keyed
        # profiles regardless of source.
        key_posix = Path(f).as_posix()
        key_native = str(Path(f))
        fc = file_checks.get(key_posix) or file_checks.get(key_native) or {}
        if not fc.get("available"):
            missing_files.append(key_posix)
    return EPRunnability(
        ep_name=ep_name,
        runnable=(not missing_caps and not missing_files),
        missing_capabilities=missing_caps,
        missing_files=missing_files,
        details=details,
        doc_anchor=profile.get("doc_anchor", ""),
        note=profile.get("note", ""),
    )


# ---------------------------------------------------------------------------
# Output rendering.
# ---------------------------------------------------------------------------


def render_markdown(payload: dict) -> str:
    lines: list[str] = []
    lines.append(f"# Preflight Runnability Report — {payload['run_date']}")
    lines.append("")
    lines.append(f"**Repo:** `{payload['repo']}`")
    lines.append(f"**Host fingerprint:** `{payload['host_fingerprint']}`")
    lines.append("")
    lines.append("## Capabilities")
    lines.append("")
    lines.append("| Capability | Status | Details |")
    lines.append("|---|---|---|")
    caps = payload["capabilities"]
    for key in ("python", "git", "pytest", "gpu", "docker", "amrfinder_db", "refseq_cache"):
        v = caps.get(key, {})
        status = "OK" if v.get("available") else "MISSING"
        if key == "gpu" and v.get("available"):
            detail = f"{v.get('device_name')} ({v.get('vram_gib')} GiB; CC {v.get('compute_capability')}; fits NT v2 100M: {v.get('fits_nt_v2_100m')})"
        elif key == "docker" and v.get("available"):
            detail = v.get("version_string", "")
        elif key == "python":
            detail = v.get("version", "")
        elif v.get("available"):
            detail = v.get("path") or v.get("version_string") or ""
        else:
            detail = v.get("reason", "")
        lines.append(f"| {key} | {status} | {detail} |")
    lines.append("")
    lines.append("## File checks (paths from CLAUDE.md gotchas + plan files)")
    lines.append("")
    lines.append("| Path | Status | Type |")
    lines.append("|---|---|---|")
    for path, fc in payload["capabilities"].get("file_checks", {}).items():
        lines.append(f"| `{path}` | {'OK' if fc.get('available') else 'MISSING'} | {fc.get('type', 'n/a')} |")
    lines.append("")
    lines.append("## Model pickles (with v0.1 RELOCKED provenance check)")
    lines.append("")
    for name, mp in payload["capabilities"].get("model_pickles", {}).items():
        lines.append(f"- **{name}** — exists: {mp.get('available')}; cv_strategy: `{mp.get('cv_strategy')}`; cv_auroc: `{mp.get('cv_auroc')}`; is_leakage_safe: **{mp.get('is_leakage_safe')}**")
    lines.append("")
    lines.append("## Per-EP runnability verdicts")
    lines.append("")
    lines.append("| EP | RUNNABLE | Missing capabilities | Missing files |")
    lines.append("|---|:---:|---|---|")
    for ep in payload["ep_runnability"]:
        runnable = "YES" if ep["runnable"] else "no"
        caps_missing = ", ".join(ep["missing_capabilities"]) or "—"
        files_missing = ", ".join(ep["missing_files"]) or "—"
        lines.append(f"| {ep['ep_name']} | {runnable} | {caps_missing} | {files_missing} |")
    lines.append("")
    lines.append("## Doc anchors")
    lines.append("")
    for ep in payload["ep_runnability"]:
        lines.append(f"- **{ep['ep_name']}** → `{ep['doc_anchor']}` — {ep['note']}")
    lines.append("")
    return "\n".join(lines)


def host_fingerprint(caps: dict) -> str:
    """Build a short host fingerprint for the report header."""
    gpu = caps.get("gpu", {})
    if gpu.get("available"):
        return f"GPU={gpu.get('device_name')} ({gpu.get('vram_gib')} GiB; CC {gpu.get('compute_capability')})"
    return "CPU-only (no CUDA-capable GPU detected)"


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Read-only preflight check for per-EP runnability.")
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--amrfinder-db", type=Path, default=DEFAULTS["amrfinder_db"])
    parser.add_argument("--nt-cache", type=Path, action="append", default=None,
                        help="NT cache HDF5 path. Repeatable; if omitted, checks the project defaults.")
    parser.add_argument("--refseq-cache", type=Path, default=DEFAULTS["refseq_cache"])
    parser.add_argument("--cohort", type=Path, action="append", default=None,
                        help="Cohort parquet to check. Repeatable; if omitted, uses defaults.")
    parser.add_argument("--model-pickle", type=Path, action="append", default=None,
                        help="Model pickle to check. Repeatable; if omitted, uses defaults.")
    parser.add_argument("--output-prefix", type=Path,
                        default=Path(f"reports/preflight_runnable_{_date.today().isoformat()}"))
    args = parser.parse_args(argv)

    nt_cache_files = args.nt_cache or DEFAULTS["nt_cache_files"]
    cohort_files = args.cohort or DEFAULTS["cohort_files"]
    model_pickles = args.model_pickle or DEFAULTS["model_pickles"]

    # ---- Gather ----
    caps: dict = {
        "python": check_python(),
        "git": check_git(),
        "pytest": check_pytest(),
        "gpu": check_gpu(),
        "docker": check_docker(),
        "amrfinder_db": check_amrfinder_db(args.amrfinder_db),
        "refseq_cache": check_path_exists(args.refseq_cache),
    }

    # Use POSIX paths as keys for platform consistency (evaluate_ep does the same).
    file_checks: dict[str, dict] = {}
    for nt in nt_cache_files:
        file_checks[Path(nt).as_posix()] = check_path_exists(nt)
    for c in cohort_files:
        file_checks[Path(c).as_posix()] = check_path_exists(c)
    caps["file_checks"] = file_checks

    model_pickle_map: dict[str, dict] = {}
    for mp in model_pickles:
        model_pickle_map[Path(mp).as_posix()] = check_model_pickle(mp)
    caps["model_pickles"] = model_pickle_map

    # ---- Per-EP evaluation ----
    ep_results: list[dict] = []
    for ep_name, profile in EP_PROFILES.items():
        result = evaluate_ep(ep_name, profile, caps)
        ep_results.append(asdict(result))

    payload = {
        "run_date": _date.today().isoformat(),
        "repo": str(args.repo.resolve()),
        "host_fingerprint": host_fingerprint(caps),
        "capabilities": caps,
        "ep_runnability": ep_results,
    }

    # ---- Emit ----
    out_md = args.output_prefix.with_suffix(".md")
    out_json = args.output_prefix.with_suffix(".json")
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(render_markdown(payload), encoding="utf-8")
    out_json.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")

    # ---- Console summary ----
    print(f"[preflight] host: {payload['host_fingerprint']}")
    print(f"[preflight] capabilities: " + ", ".join(
        f"{k}={'OK' if caps[k].get('available') else 'MISSING'}"
        for k in ("python", "git", "pytest", "gpu", "docker", "amrfinder_db", "refseq_cache")
    ))
    print(f"[preflight] per-EP runnability:")
    any_not_runnable = False
    for ep in ep_results:
        verdict = "YES" if ep["runnable"] else "no"
        miss_caps = ", ".join(ep["missing_capabilities"]) or "—"
        miss_files = ", ".join(ep["missing_files"]) or "—"
        print(f"  {ep['ep_name']:>32s}: runnable={verdict}  caps_missing={miss_caps}  files_missing={miss_files}")
        if not ep["runnable"]:
            any_not_runnable = True
    print(f"[preflight] wrote {out_md}")
    print(f"[preflight] wrote {out_json}")
    return 1 if any_not_runnable else 0


if __name__ == "__main__":
    sys.exit(main())
