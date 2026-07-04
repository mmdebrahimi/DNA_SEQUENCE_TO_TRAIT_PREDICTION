# Soraya --until-mvp: DGRP learned-path test

Verdict: mvp-reached (4/4). The expected 5th de-confounded NEGATIVE, confirmed + informative. DGRP2 x StarvationRes, N=197: naive cv_r2 0.0154 / de-confounded within-clade r2 0.0349 (perm floor ~0.02) -> LEARNED_PATH_NEGATIVE_UNDER_DECONFOUNDING. Learned/embedding path now 0-for-5 under de-confounding. New lesson: DGRP lacks discrete structure (191/197 one clade) -> clade-de-confounding degenerate; used PCA-continuous structure. Perf trap fixed (skip-before-parse: 8min->47s). No-resume: bounded attempts per active session only.
