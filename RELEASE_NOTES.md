# Version 1.0.0

**Certified whole-domain bounds and equality barriers for the odd theta–Weil form**  
Leonard van Hemert

This release presents explicit bounds on the full odd form domain, with complete physical residuals, arithmetic tails and a transparent loss/retained-energy tradeoff. It includes a quantitative obstruction showing why the fixed positive comparison cannot reach zero loss by perfect trial enrichment alone.

The least scalar loss among the positive-retention certificates in the 231 tested candidates is `0.031324476155738682727829`, retaining `0.000118321042` of the positive energy. Balanced bounds retain `0.132584532412` and `0.212458845389`, with losses below `0.034762` and `0.039482`. These are certified tested values, not a proved global optimum.

Assets include the main paper, complete supplement, exact certificate inputs and bounded read-only verifier, and a pinned Lean project. Formal verification covers the stated abstract proof chain and rational consequences; the theta–arithmetic analytical and numerical specialization inputs are explicitly identified in the verification reports.

RH, zero-loss positivity, the even-sector argument and the full Li transfer remain unresolved.

## Publication

Published on 13 October 2026. Cite the versioned preprint at [10.5281/zenodo.23337494](https://doi.org/10.5281/zenodo.23337494). The [Zenodo archive](https://doi.org/10.5281/zenodo.23337494) contains both PDFs, the complete source ZIP and original checksums. The [GitHub release](https://github.com/LeonardSEO/odd-theta-weil-bounds/releases/tag/v1.0.0) also includes the alternative tar.gz archive.

The [publication CI run](https://github.com/LeonardSEO/odd-theta-weil-bounds/actions/runs/38374755999) passed the Lean and certificate jobs. Detailed scope and evidence are in [the verification report](verification/VERIFICATION_REPORT.md).

Subsequent author-attribution, citation and documentation updates are recorded on the default branch. The v1.0.0 tag, PDFs and release archives remain unchanged.
