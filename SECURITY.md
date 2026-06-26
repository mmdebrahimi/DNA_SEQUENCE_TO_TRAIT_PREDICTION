# Security & responsible-use policy

## Not a clinical tool

`dna-decode` is **research-use software**. It is **not** a clinical diagnostic, not a medical device, and
its output must **not** be used to guide patient treatment, antimicrobial therapy, infection control, or
any clinical or regulatory decision. Every call carries a validation tier and named blind spots precisely
so that scope is never overstated; honor them.

## Scientific-integrity expectations

- Do not present calls outside their stated validation tier or supported organism/drug surface
  (`dna-decode list` is authoritative).
- Do not average distinct validation tiers into a single aggregate "score" — the project deliberately
  emits no aggregate headline.
- Resistance calls have known blind spots (e.g. an S call cannot rule out efflux / porin-loss / regulatory
  resistance the curated determinant databases do not capture). These are reported by the tool; preserve them.

## Reporting a vulnerability or a scientific defect

Two distinct channels, both via the GitHub repository
(<https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION>):

- **Software vulnerability** (code execution, dependency, supply-chain): open a **private security
  advisory** (GitHub → Security → Report a vulnerability) rather than a public issue.
- **Scientific defect** (a wrong call, a mis-stated validation claim, an over-broad supported-surface
  entry): open a public issue with the input, the command, and the observed vs expected output. Wrong
  *science* in a tool like this is a safety issue, not just a bug.

## Supply-chain notes

- The deterministic decoder imports no `torch` / `transformers` / GPU stack (the heavy foundation-model
  extra is optional and unused by the shipped CLI).
- External tool integrations (AMRFinder, BLAST+) run pinned container images / local binaries; see
  `docs/quickstart.md`. No telemetry is collected; the CLI runs fully offline for the committed-data paths.
