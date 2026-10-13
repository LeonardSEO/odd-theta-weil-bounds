# Verification report — version 1.0.0

**Leonard van Hemert · 13 October 2026**

Scope: a 14-page main article and 17-page supplement, with the readable pinned Lean project and portable finite certificate closure. The title is *Certified whole-domain bounds and equality barriers for the odd theta–Weil form*.

## Mathematical and formal correspondence

The publication concentrates the certified whole-domain odd bounds, the full-residual mechanism, the retained-energy tradeoff and the quantitative equality barrier. Secondary signed-action, finite positivity and outside-coupling results are in the supplement. Standard factorization/residual/Schur arguments are attributed; the claimed contribution is the explicit theta–arithmetic realization and its controlled obstruction. Prior publication v1.0.0 and its version DOI were confirmed against the public GitHub tag/release and Zenodo record; no DOI has been assigned to this new work.

The independent reviewer inspected the actual formal definitions, proof sources, fully explicit elaborated types, manuscript mapping and transitive axioms. In particular, `enriched_certified_retained_energy` composes the full action decomposition and separate errors with the finite RE5 price, derives the resolvent comparison, transfers using the original coefficient metric N to the full form domain, and derives the retained-energy/scalar-loss expression. The target does not assume that final expression or substitute the physical Gram for N. Complex real-scalar restriction preserves every complex vector. See `FORMAL_CORRESPONDENCE.md`.

## Checks actually executed

| Check | Result and limits |
|---|---|
| Pinned whole Lean build | **PASS**, Lean 4.32.0, mathlib `81a5d257c8e410db227a6665ed08f64fea08e997`; 3005 jobs. Actual log: `lean/build-output.txt`. One harmless `unnecessarySimpa` style warning remains. |
| Transitive axiom inspection | **PASS**, all 22 declarations in `lean/Axioms.lean`; actual output in `lean/axioms-output.txt`. Every closure contains only `propext`, `Classical.choice`, `Quot.sound`. No `sorryAx`, custom mathematical axiom or native-evaluation trust extension. |
| Strict axiom-log checker | **PASS**, all expected reports present. Its rejection self-check covers missing/duplicate reports, `sorryAx`, a native-evaluation extension and a custom axiom. |
| Independent intended-statement/type audit | **PASS**, separate Lean invocation checked the domain structure, composed target, metric/complex targets and public numeric types, and printed selected transitive axiom sets. |
| Exact endpoints | **PASS**, generated source consistency for all 231 rows; independent literal-to-JSON comparison also passed, with 173 positive-retention rows. Three displayed numeric consequences are kernel-checked conditional statements. No global parameter optimum is asserted. |
| Bounded numerical certificate replay | **PASS**, performed from a different working directory: 61 file fingerprints, 158 embedded digest links, original coupled spectral guards, all 231 rational endpoint rows, signed five-dimensional physical/signed spectral and outside-coupling guards. Four published residual prices were reconstructed from cached full matrices; no new search. |
| Fixed-price robust-capture check | **PASS**, fresh finite inverse/Rayleigh calculation gives `0.1478676174248644… > 0.14`, within the stated sufficient-certificate cone. |
| Provenance-only portability repair | **PASS**, original versus sanitized norm-cell payloads compare identically after excluding exactly the historical command/cwd fields. The hash-only closure check then passed; numerical replay was not repeated for unchanged mathematics. Historical digest alias and exact adaptation are explicit in the manifest. |
| PDF compilation | **PASS**, Tectonic 0.17.0. Final logs have no overfull boxes or undefined references. Main 14 pages, supplement 17 pages. |
| PDF inspection | Every page rendered and inspected for layout; corrected overflowing expressions and one visible spacing-command typo. PDF text contains no unresolved reference markers. PDF titles and author metadata match the sources. |
| Release metadata | CFF and workflow YAML and Zenodo JSON parsed; author, title, version, preparation date and preprint metadata checked. CFF uses a software-package record with the paper as preferred article citation. No new DOI present. |

The detailed numerical trust boundary is in `CERTIFICATE_SCOPE.md`. The authoritative read-only command is `python3 scripts/check_certificates.py`, with its pinned dependencies. The generator sources are archival and may be expensive or overwrite certificates; they are not the routine verifier.

## Analytical boundary

The Lean proof is a **verified conditional abstract core**, not a formalized theta–Weil theorem. Construction of the theta operators and closed forms, inverse existence, actual trial admissibility, Fourier/prime/zero tails, native interval matrix evidence and original-form/kernel specialization are written analytical/certificate inputs. Upstream integrations, complete zero verification, Fourier DFT and large-Gram/action construction were reused unchanged, not recomputed in this publication-preparation task.

RH, zero-loss positivity, the even-sector argument and the full Li transfer remain unresolved.

## Packaging and publication boundary

`scripts/package_release.py` fingerprints the publication source, creates both versioned source archives, verifies their exact membership and contents (including all readable Lean modules and the lockfile), and checks PDF/archive SHA-256 values. Its successful output and the final checksum verification are reported at handoff. `MANIFEST.sha256` records the source set; `SHA256SUMS` records the release assets. Archives contain the source manifest but exclude the asset checksum file to avoid circularity.

The provided GitHub Actions workflow has been prepared and parsed, but not executed on GitHub. A clean fresh dependency install, cross-platform numerical regeneration and an external independent Lean kernel were not tested. No upload, push, public release or DOI minting occurred. The active research archive and prior frozen publications were not changed.
