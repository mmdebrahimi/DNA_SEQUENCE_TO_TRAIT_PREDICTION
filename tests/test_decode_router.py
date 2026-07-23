"""Offline tests for the input-aware decode router (dna_decode/decode_router.py).

Pure file-sniff + grounded table + registry read. No network, no decoder execution.
"""
from __future__ import annotations

from dna_decode.decode_router import (
    detect_input_kind, applicable_decoders, render_decode_plan, run_decode_plan, DECODERS,
)


def test_detect_vcf(tmp_path):
    v = tmp_path / "s.vcf"
    v.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n10\t1\t.\tC\tT\n", encoding="utf-8")
    assert detect_input_kind(v) == "vcf_human"


def test_detect_protein_fasta(tmp_path):
    f = tmp_path / "p.fasta"
    f.write_text(">TP53\nMEEPQSDPSVEPPLSQETFSDLWKLLPENN\n", encoding="utf-8")  # has E/F/L/P -> protein
    assert detect_input_kind(f) == "protein_fasta"


def test_detect_nucleotide_fasta(tmp_path):
    f = tmp_path / "g.fna"
    f.write_text(">contig\nATGAGCATTCAACATTTCCGTGTCGCCCTTATT\n", encoding="utf-8")  # only ACGT -> nucleotide
    assert detect_input_kind(f) == "nucleotide_fasta"


def test_detect_nucleotide_with_iupac_ambiguity(tmp_path):
    f = tmp_path / "g.fna"
    f.write_text(">c\nACGTNRYSWKM\n", encoding="utf-8")  # IUPAC nt ambiguity codes -> still nucleotide
    assert detect_input_kind(f) == "nucleotide_fasta"


def test_detect_unknown(tmp_path):
    f = tmp_path / "x.txt"
    f.write_text("this is not a sequence file\n", encoding="utf-8")
    assert detect_input_kind(f) == "unknown"


def test_applicable_protein_is_forward_and_inverse():
    routes = {d.route for d in applicable_decoders("protein_fasta")}
    assert routes == {"dna-forward", "dna-inverse"}


def test_applicable_vcf_covers_the_orphaned_human_cells():
    # the 3 human cells that were orphaned from the unified dna-decode entry must be routed here
    routes = {d.route for d in applicable_decoders("vcf_human")}
    assert {"dna-pgx", "dna-clinvar", "dna-hla"} <= routes


def test_applicable_nucleotide_covers_amr_typing_and_finder():
    routes = {d.route for d in applicable_decoders("nucleotide_fasta")}
    assert {"dna-amr", "dna-mlst", "dna-serotype", "dna-plasmid", "dna-coloc", "dna-forward"} <= routes


def test_every_decoder_example_names_its_route():
    # each example command must invoke its own route (grounded, not a guessed runbook)
    for kind, decs in DECODERS.items():
        for d in decs:
            assert d.example.startswith(d.route), f"{kind}:{d.route} example does not start with the route"


def test_render_plan_lists_commands_and_tier(tmp_path):
    f = tmp_path / "p.fasta"
    f.write_text(">P\nMEEPQSDPSVEPP\n", encoding="utf-8")
    out = render_decode_plan(f)
    assert "protein FASTA" in out
    assert "dna-forward" in out and "run:" in out


def test_render_plan_unknown_is_honest(tmp_path):
    f = tmp_path / "x.txt"
    f.write_text("nonsense\n", encoding="utf-8")
    out = render_decode_plan(f)
    assert "No decoder recognizes" in out


class _RecordRunner:
    """Injectable runner: records (trait, argv) calls instead of executing a real decoder."""
    def __init__(self, rc=0):
        self.calls = []
        self.rc = rc

    def __call__(self, trait, argv):
        self.calls.append((trait, argv))
        return self.rc


def test_run_protein_invokes_inverse(tmp_path):
    f = tmp_path / "p.fasta"
    f.write_text(">P\nMEEPQSDPSVEPP\n", encoding="utf-8")
    r = _RecordRunner()
    rc = run_decode_plan(f, runner=r)
    assert rc == 0
    traits = [t for t, _ in r.calls]
    assert "inverse" in traits          # inverse auto-runs offline
    assert "forward" not in traits      # forward needs a --mutation -> report-only, not run
    inv_argv = next(a for t, a in r.calls if t == "inverse")
    assert "--protein-fasta" in inv_argv and "--target-percentile" in inv_argv


def test_run_nucleotide_delegates_to_profile(tmp_path):
    f = tmp_path / "g.fna"
    f.write_text(">c\nATGAGCATTCAACAT\n", encoding="utf-8")
    r = _RecordRunner()
    run_decode_plan(f, runner=r)
    traits = [t for t, _ in r.calls]
    assert traits == ["profile"]                       # reuse profile, no per-decoder re-run
    assert "--sample-id" in r.calls[0][1]


def test_run_vcf_is_report_only(tmp_path):
    f = tmp_path / "s.vcf"
    f.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n10\t1\t.\tC\tT\n", encoding="utf-8")
    r = _RecordRunner()
    rc = run_decode_plan(f, runner=r)
    assert rc == 0
    assert r.calls == []                               # never auto-guess a --gene


def test_run_unknown_is_honest_and_runs_nothing(tmp_path):
    f = tmp_path / "x.txt"
    f.write_text("nonsense\n", encoding="utf-8")
    r = _RecordRunner()
    rc = run_decode_plan(f, runner=r)
    assert rc == 0 and r.calls == []


def test_run_propagates_nonzero_exit(tmp_path):
    f = tmp_path / "p.fasta"
    f.write_text(">P\nMEEPQSDPSVEPP\n", encoding="utf-8")
    r = _RecordRunner(rc=3)
    assert run_decode_plan(f, runner=r) == 3           # an auto-run decoder's failure surfaces


def test_cli_decode_run_flag_calls_runner(tmp_path, monkeypatch):
    # the unified CLI's `decode <file> --run` must route to run_decode_plan
    import dna_decode.cli as uni
    import dna_decode.decode_router as router
    seen = {}

    def _stub(p, **k):
        seen["path"] = str(p)
        return 0

    monkeypatch.setattr(router, "run_decode_plan", _stub)
    f = tmp_path / "p.fasta"
    f.write_text(">P\nMEEPQSDPSVEPP\n", encoding="utf-8")
    rc = uni.main(["decode", str(f), "--run"])
    assert rc == 0 and seen["path"] == str(f)


def test_protein_plan_shows_forward_capability_hint(tmp_path):
    # slice 2 seam: when the forward route applies, the router surfaces which learned methods run HERE
    f = tmp_path / "p.fasta"
    f.write_text(">P\nMEEPQSDPSVEPP\n", encoding="utf-8")
    out = render_decode_plan(f)
    assert "forward:" in out
    assert "--capabilities" in out
    # blosum62 is always runnable, so it must appear in the runnable list
    assert "blosum62" in out
