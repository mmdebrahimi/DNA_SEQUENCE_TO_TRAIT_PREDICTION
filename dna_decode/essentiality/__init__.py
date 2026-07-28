"""Essentiality decoder (NEW cell, NON-frozen). Single-gene KO -> essential/non-essential.

v0 = the R1 deterministic conserved-core decoder: essential genes are dominated by a universal
functional core (translation, replication, transcription, cell-envelope/division). Label-independent
by construction (it reads gene FUNCTION, not an essentiality label), so it is buildable + testable
without the (externally-walled) gold-standard essential-gene labels. See
`plans/New_Phenotype_Direction_Ideation_And_Plan_2026-07-28.md`.
"""
