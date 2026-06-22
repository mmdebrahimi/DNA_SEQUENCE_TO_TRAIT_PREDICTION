# result — tb-goldset  (verdict: gold-set-help mvp-reached; independent number = external wall; indel-norm scoped-deferred)
THREAD 1 (gold-set help) DONE + verified:
- Feasibility found: BV-BRC has 0 measured TB DST (taxon 1773) -> the golden baseline is NOT an API fetch;
  it needs a specific non-CRyPTIC TB DST+genome source.
- Built: tb_goldset.cryptic_accessions + assert_independent (the provably-not-in-CRyPTIC check, verified on
  the REAL CRyPTIC table: ERR4810489 correctly excluded); scripts/build_tb_goldset_manifest.py (source-TSV
  -> leakage-checked per-drug manifest the existing scorer consumes); wiki/tb_goldset_howto_2026-06-22.md
  (the playbook for a SWE). 4 tests + existing arm. The scored INDEPENDENT number is gated on the user
  obtaining a source (external/data wall, precisely characterized + a one-command path once data is in hand).
THREAD 2 (indel-norm): INVESTIGATED -> NOT shipped. CRyPTIC indels are unanchored (`2155741_del_c`) vs WHO
  VCF-anchored (`2152935 CACTCG>C`); a quick matcher would produce WRONG matches (worse than the honest
  SNV+MNV lower bound). Correct paths (scoped, deferred): genomic left-alignment normalization, OR a
  CRyPTIC-EFFECTS protein-form match with a circularity check. Per the don't-ship-passing-looking-but-wrong
  rail, deferred rather than risk wrong numbers.
Full suite 1597 passed; frozen surface + TB leak-guard byte-unchanged.
