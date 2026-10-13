# Certified whole-domain bounds and equality barriers for the odd theta–Weil form

**Leonard van Hemert · Version 1.0.0 · 13 October 2026**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23337494.svg)](https://doi.org/10.5281/zenodo.23337494)

[Read the paper](https://github.com/LeonardSEO/odd-theta-weil-bounds/releases/download/v1.0.0/paper.pdf) · [Supplement](https://github.com/LeonardSEO/odd-theta-weil-bounds/releases/download/v1.0.0/supplement.pdf) · [Download v1.0.0](https://github.com/LeonardSEO/odd-theta-weil-bounds/releases/tag/v1.0.0) · [Verification](https://github.com/LeonardSEO/odd-theta-weil-bounds/actions/workflows/verify.yml)

This publication proves explicit lower bounds on the full complex odd form domain of a theta-weighted Riemann–Weil form. It combines an admissible residual-enriched trial, complete physical residual coupling and certified arithmetic tails. A quantitative equality-vector obstruction explains a limitation of the fixed positive comparison.

The **14-page paper** develops the main theorem and its significance. The **17-page supplement** supplies the domain arguments, complete arithmetic tails, matrix estimates and certificate derivations. This repository contains both manuscripts, the required certificate data and scripts, and the complete Lean proof project.

## Main result

Under the exact hypotheses of the main theorem, with `z = (I − P_K)u` on the original odd form domain and the original positive energy `E₈`, the certified family contains:

| Retained energy coefficient | Scalar loss α (exact decimal endpoint) |
|---:|---:|
| 0.000118321042 | 0.031324476155738682727829 |
| 0.132584532412 | 0.034761534500246499521894 |
| 0.212458845389 | 0.0394810814425857399223805 |

Each row means `a_W(u) ≥ coefficient · E₈(z) − α · ‖z‖²`. The first loss is the least among the **positive-retention certificates** in 231 tested candidates, not a proved global optimum. The bounds must be compared as whole energy expressions; the smaller scalar loss retains less energy.

The equality barrier proves that the fixed positive midpoint comparison overshoots an actual equality vector by more than `7/5000` in the original Hilbert norm. This obstructs zero loss within that comparison even with exact residual enrichment. It is not a negative direction of the signed form. The supplement also records exact signed cancellation, controlled signed action, five-dimensional positivity and nontrivial coupling outside that subspace.

RH, zero-loss positivity, the even-sector argument and the full Li transfer remain unresolved.

## Reproduce the verification

The [release verification report](verification/VERIFICATION_REPORT.md) records the executed checks. The frozen publication commit passed both the Lean and certificate jobs in [GitHub Actions](https://github.com/LeonardSEO/odd-theta-weil-bounds/actions/runs/38374755999).

For the exact published source, download a source archive from [v1.0.0](https://github.com/LeonardSEO/odd-theta-weil-bounds/releases/tag/v1.0.0), extract it and run the commands below from its root.

### Lean

The complete readable project is in [lean/](lean/), pinned to **Lean 4.32.0** and mathlib commit `81a5d257c8e410db227a6665ed08f64fea08e997`, with dependency revisions in `lake-manifest.json`. The checked proof chain includes the domain-explicit full-residual resolvent identity and upper bound, complete action-error inequalities, coefficient-metric transfer on the form domain, and exact rational certificate consequences.

This is **conditional formal verification of the abstract proof chain and its rational consequences**, not formal verification of the entire theta–Weil realization. Analytical construction of the operators/forms, theta identities, domain admissibility, verified numerical error enclosures and the physical specialization remain inputs. The precise theorem mapping and remaining analytical inputs are recorded in [verification/FORMAL_CORRESPONDENCE.md](verification/FORMAL_CORRESPONDENCE.md). Build and transitive-axiom evidence accompany the project; see [verification/VERIFICATION_REPORT.md](verification/VERIFICATION_REPORT.md).

With the pinned Lean toolchain available through Elan:

```sh
cd lean
lake exe cache get
lake build
lake env lean Axioms.lean > axioms-current.txt
python3 ../scripts/check_axioms.py axioms-current.txt
cd ..
```

The axiom check accepts only Lean's standard logical axioms `propext`, `Classical.choice` and `Quot.sound`. It rejects `sorryAx`, custom axioms and native-evaluation trust extensions. Explicit hypotheses in the theorem types are mathematical premises, not hidden axioms; they are listed in the correspondence report.

### Certificates

Use Python 3.14 and an isolated virtual environment:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-certificate.txt
.venv/bin/python scripts/check_certificates.py
```

Run without Python's `-O` flag: assertions perform the checks. The checker is read-only. It verifies recorded dependency hashes, exact rational endpoints, finite interval spectral consequences, the published residual prices from cached full matrices, and selected signed and barrier consequences. It performs no new parameter search or expensive integration/zero/DFT replay. [verification/CERTIFICATE_SCOPE.md](verification/CERTIFICATE_SCOPE.md) identifies precisely which upstream enclosures remain accepted premises.

Historical generators in `computations/` may overwrite certificate files and perform expensive work; do not use them as the routine verifier. Optional generation dependencies are recorded separately in `requirements-generation.txt`. Numerical regeneration may depend on BLAS locator choices; cross-platform byte-identical regeneration is not claimed.

### Manuscripts

With Tectonic 0.17.0 (used for this release), run from the repository root:

```sh
tectonic paper.tex
tectonic supplement.tex
```

Both sources are standalone and include their bibliography. TeX Live with `latexmk -pdf paper.tex` and `latexmk -pdf supplement.tex` is an alternative, not tested for this release. The supplement is linked by filename; keep both PDFs together.

## Downloads and integrity

The [Zenodo record](https://doi.org/10.5281/zenodo.23337494) archives the paper, supplement, complete source ZIP and original checksums. The [GitHub v1.0.0 release](https://github.com/LeonardSEO/odd-theta-weil-bounds/releases/tag/v1.0.0) also provides the alternative tar archive:

- `paper.pdf` and `supplement.pdf`.
- Complete source archives in `.tar.gz` and `.zip` formats, including readable Lean proofs and the dependency lock.
- `SHA256SUMS` for the four files above.

To check downloaded release assets, place both PDFs and `SHA256SUMS` in one directory and both archives in its `dist/` subdirectory, then run:

```sh
shasum -a 256 -c SHA256SUMS
```

Each source archive also contains `MANIFEST.sha256`. From the extracted source directory, verify its contents with:

```sh
shasum -a 256 -c MANIFEST.sha256
```

The tagged release and its downloadable assets are fixed. The default branch includes subsequent citation and documentation updates; its source manifest describes that checkout. Installed dependencies and build caches are excluded from the source archives.

## Citation

Leonard van Hemert. *Certified whole-domain bounds and equality barriers for the odd theta–Weil form*. Version 1.0.0, Zenodo, 2026. [10.5281/zenodo.23337494](https://doi.org/10.5281/zenodo.23337494).

Machine-readable citation metadata is available in [CITATION.cff](CITATION.cff).

## Documentation and rights

[Release notes](RELEASE_NOTES.md) describe version 1.0.0. [Provenance](PROVENANCE.md) records the mathematical sources and release lineage. The [formal correspondence](verification/FORMAL_CORRESPONDENCE.md) and [certificate scope](verification/CERTIFICATE_SCOPE.md) identify the verification boundary in detail.

Copyright © 2026 Leonard van Hemert. All rights reserved; see [NOTICE](NOTICE). Public availability does not grant an open-content or open-source licence. Third-party dependency attribution is recorded in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
