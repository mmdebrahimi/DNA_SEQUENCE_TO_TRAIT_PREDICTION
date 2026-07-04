# Soraya --until-mvp: frozen-spec local conservation (own the number)

Verdict: mvp-reached (4/4), target EXCEEDED. The project computes its OWN independent-sites conservation (pinned+hashed spec; ProteinGym MSAs+weights). Function median 0.476 (>= 0.43 target); weighted subset (7) reproduces ProteinGym Site-Independent to 0.0145 (implementation validated, number owned not cited); unweighted fallback (12) 0.480. Closes the purist path. verify-in-batch caught 2 bugs (focus-coordinate map: 0.047->0.451; run weight-gate). Commit HEAD. No-resume: bounded attempts per active session only.
