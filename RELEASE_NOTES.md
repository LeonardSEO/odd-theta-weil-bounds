# Version 1.0.0

**Certified whole-domain bounds and equality barriers for the odd theta–Weil form**  
Leonard van Hemert

This release presents explicit bounds on the full odd form domain, with complete physical residuals, arithmetic tails and a transparent loss/retained-energy tradeoff. It includes a quantitative obstruction showing why the fixed positive comparison cannot reach zero loss by perfect trial enrichment alone.

The least scalar loss among the positive-retention certificates in the 231 tested candidates is `0.031324476155738682727829`, retaining `0.000118321042` of the positive energy. Balanced bounds retain `0.132584532412` and `0.212458845389`, with losses below `0.034762` and `0.039482`. These are certified tested values, not a proved global optimum.

Assets include the main paper, complete supplement, exact certificate inputs and bounded read-only verifier, and a pinned Lean project. Formal verification covers the stated abstract proof chain and rational consequences; the theta–arithmetic analytical and numerical specialization inputs are explicitly identified in the verification reports.

RH, zero-loss positivity, the even-sector argument and the full Li transfer remain unresolved.

## Suggested GitHub repository description

Certified odd theta–Weil bounds, equality barriers, reproducible certificates and a pinned Lean proof project.

## Publication handoff

1. Review the PDFs, `verification/` reports, `NOTICE`, metadata and checksums. Confirm the author and preparation date. The public repository is [LeonardSEO/odd-theta-weil-bounds](https://github.com/LeonardSEO/odd-theta-weil-bounds).
2. Publish this readable source tree, including `lean/`, its lock file and verification logs. Run the provided GitHub Actions workflow. Resolve any platform-specific failure before release.
3. Tag the reviewed, verified source commit `v1.0.0`. Create a matching GitHub release and attach `paper.pdf`, `supplement.pdf`, the versioned source archives and `SHA256SUMS`. Do not change assets after freezing the release.
4. Archive that exact release through Zenodo (connect the repository before publishing the release, or upload the matching reviewed assets manually). Confirm the resulting record corresponds to the tag and retain the actual version DOI and concept DOI separately.
5. In an explicit follow-up metadata commit or new version, record the issued version DOI in `CITATION.cff`, `.zenodo.json` and `README.md`. Update the Zenodo record as appropriate. Do not invent a DOI, relabel the earlier work's DOI, or rewrite the frozen source tag. If changed source/PDF bytes are to be released, use a new version with new archives and checksums.

The follow-up instruction authorizes GitHub publication. Zenodo archival and DOI minting remain a separate handoff.
