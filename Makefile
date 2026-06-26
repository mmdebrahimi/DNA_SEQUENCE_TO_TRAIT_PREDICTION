# dna-decode — common commands (agent- and human-friendly). Recipes use TABs.
.PHONY: install test smoke example-amr example-observed list

install:
	pip install -e ".[dev]"

test:
	pytest tests/ -q

smoke list:
	dna-decode list

# Verified zero-dependency example (wheel-only; no Docker / no GPU / no DB)
example-observed:
	dna-decode amr --drug efavirenz --observed RT:K103N

# Verified example from committed AMRFinder run data (bacterial path; no Docker needed)
example-amr:
	dna-decode amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_000492655.1
