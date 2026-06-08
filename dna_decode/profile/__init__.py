"""Unified genome profile — run every applicable assembly-FASTA decoder in one command.

`dna-decode profile <assembly.fna>` runs pathotype + plasmid + serotype + resfinder (all blastn / pure-stdlib,
no Docker) on one assembly and emits a single per-genome report. The "tell me everything" UX over the shipped
decoders. Each decoder degrades independently (missing DB / no blastn -> that section status 'unavailable',
the rest still run). amr genome-mode (AMRFinder/Docker) is folded in only when an --amrfinder-run is given.
"""
