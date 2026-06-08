"""Cross-tool AMR concordance analysis — compare the two INDEPENDENT acquired-gene callers.

dna-amr (AMRFinder DB) vs dna-resfinder (ResFinder DB) on the same genome: which acquired AMR genes both
call, which only one calls. Surfaces the independent cross-check resfinder was built for. An ANALYSIS over
existing decoders (not a new trait/DB). Gene naming differs between the two tools, so the fair comparison is
at the gene-FAMILY level (allele-variant suffix stripped). Pure core (offline-testable) + a CLI.
"""
