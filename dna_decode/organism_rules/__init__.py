"""Non-frozen, organism-routed decoder extensions.

This package holds decoder cells for organisms OUTSIDE the frozen E. coli/Klebsiella/S. aureus/
C. auris AMR surface (reproducibility freeze 2026-06-13). The first member is the M. tuberculosis
AMR cell (RIF + INH) scored on CRyPTIC against the WHO mutation catalogue — a labelled in-distribution
KNOWLEDGE_BASELINE, deterministic, no learned model. Mirrors the non-frozen overlay pattern of
`dna_decode/data/experimental_drug_rules.py` (TMP-SMX): scorer-local, branded, frozen helpers reused
but never edited.
"""
