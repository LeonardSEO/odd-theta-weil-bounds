# Certified whole-domain bounds and equality barriers for the odd theta–Weil form

**Leonard van Hemert · Version 1.0.0 · 13 October 2026**

This publication proves explicit lower bounds on the full complex odd form domain of a theta-weighted Riemann–Weil form. It combines an admissible residual-enriched trial, complete physical residual coupling and certified arithmetic tails. A quantitative equality-vector obstruction explains a limitation of the fixed positive comparison.

Read [paper.pdf](paper.pdf) for the main argument and [supplement.pdf](supplement.pdf) for the detailed analytical and certificate construction. Their sources are `paper.tex` and `supplement.tex`. The existing rights notice is [NOTICE](NOTICE).

## Mathematical scope

Under the exact hypotheses of the main theorem, with `z = (I − P_K)u` on the original odd form domain and the original positive energy `E₈`, the certified family contains:

| Retained energy coefficient | Scalar loss α (exact decimal endpoint) |
|---:|---:|
| 0.000118321042 | 0.031324476155738682727829 |
| 0.132584532412 | 0.034761534500246499521894 |
| 0.212458845389 | 0.0394810814425857399223805 |

Each row means `a_W(u) ≥ coefficient · E₈(z) − α · ‖z‖²`. The first loss is the least among the **positive-retention certificates** in 231 tested candidates, not a proved global optimum. The bounds must be compared as whole energy expressions; the smaller scalar loss retains less energy.

The equality barrier proves that the fixed positive midpoint comparison overshoots an actual equality vector by more than `7/5000` in the original Hilbert norm. This obstructs zero loss within that comparison even with exact residual enrichment. It is not a negative direction of the signed form. The supplement also records exact signed cancellation, controlled signed action, five-dimensional positivity and nontrivial coupling outside that subspace.

RH, zero-loss positivity, the even-sector argument and the full Li transfer remain unresolved.

## Lean verification

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

## Certificate checks

Use Python 3.14 and an isolated virtual environment:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-certificate.txt
.venv/bin/python scripts/check_certificates.py
```

Run without Python's `-O` flag: assertions perform the checks. The checker is read-only. It verifies recorded dependency hashes, exact rational endpoints, finite interval spectral consequences, the published residual prices from cached full matrices, and selected signed and barrier consequences. It performs no new parameter search or expensive integration/zero/DFT replay. [verification/CERTIFICATE_SCOPE.md](verification/CERTIFICATE_SCOPE.md) identifies precisely which upstream enclosures remain accepted premises.

Historical generators in `computations/` may overwrite certificate files and perform expensive work; do not use them as the routine verifier. Optional generation dependencies are recorded separately in `requirements-generation.txt`. Numerical regeneration may depend on BLAS locator choices; cross-platform byte-identical regeneration is not claimed.

## Build the PDFs

With Tectonic 0.17.0 (used for this release), run from the repository root:

```sh
tectonic paper.tex
tectonic supplement.tex
```

Both sources are standalone and include their bibliography. TeX Live with `latexmk -pdf paper.tex` and `latexmk -pdf supplement.tex` is an alternative, not tested for this release. The supplement is linked by filename; keep both PDFs together.

## Release assets and checksums

`MANIFEST.sha256` fingerprints publication/source files, including the Lean project and PDFs. `SHA256SUMS` fingerprints the two PDFs and the versioned source archives in `dist/`. The archives include the source manifest and readable Lean proof sources and dependency lock, excluding installed dependencies, caches and `.git`. `SHA256SUMS` accompanies the archives rather than being embedded in them, avoiding a circular checksum.

```sh
python3 scripts/package_release.py
shasum -a 256 -c SHA256SUMS
shasum -a 256 -c MANIFEST.sha256
```

Packaging verifies archive membership and byte contents before reporting success. Rebuilding PDFs or modifying any source changes release bytes; regenerate checksums and use a new version after a release is frozen.

## GitHub and Zenodo handoff

[RELEASE_NOTES.md](RELEASE_NOTES.md) contains the repository description, release notes and concrete publication sequence. Publish the reviewed source and matching tagged GitHub release, archive that exact release through Zenodo, and then add the **actual issued version DOI** through an explicit metadata update. No DOI is invented in this package. Do not silently rewrite a frozen tag or reuse the earlier paper's DOI.

The prepared GitHub Actions workflow builds Lean and checks its axioms and certificates. It has not been executed on GitHub during local preparation. [PROVENANCE.md](PROVENANCE.md) records the earlier versioned publication and the source lineage; [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) records dependency attribution. [ANNOUNCEMENT.md](ANNOUNCEMENT.md) is a proposed announcement, not an external message.
