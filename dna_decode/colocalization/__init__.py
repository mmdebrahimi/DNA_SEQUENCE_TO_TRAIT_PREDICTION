"""Resistance-gene × plasmid-replicon co-localization analysis.

Answers the concrete question the AMR + plasmid decoders only hinted at separately: is THIS acquired
resistance gene plasmid-borne? — by checking whether it sits on the SAME assembly contig as a called
plasmid replicon. Composes resfinder + plasmid via the engine's positions mode (subject contig per hit).
A same-contig co-location is suggestive (a circularized plasmid often assembles into one contig), NOT proof
of a single plasmid — stated as a caveat. Pure core (offline-testable) + a CLI.
"""
