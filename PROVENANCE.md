# Provenance and rights

This is version 1.0.0 of **Certified whole-domain bounds and equality barriers for the odd theta–Weil form**, by Leonard van Hemert. It is a new publication, distinct from the earlier *Coupling estimates* publication. The publication date is the manuscript's preparation date. GitHub publication is recorded in [the repository history](https://github.com/LeonardSEO/odd-theta-weil-bounds); this work is archived at [10.5281/zenodo.23337494](https://doi.org/10.5281/zenodo.23337494).

## Mathematical source

The source is the established odd-form research and its publication candidate. The main paper concentrates the full-domain residual-enriched bounds and the equality barrier. The included supplement supplies the public proof details. The source research archive and frozen earlier publications are preserved separately; their directories are not runtime dependencies.

The earlier publication is cited at its version-specific DOI, [10.5281/zenodo.23298499](https://doi.org/10.5281/zenodo.23298499), and [GitHub v1.0.0](https://github.com/LeonardSEO/theta-coupling-estimates/tree/v1.0.0). Its recorded release commit is `a8cf3560d5fbec587b23a45c93550681eeb27cd7`. The recorded earlier source archive SHA-256 is `538fe25a41df103c3ca73228986fa601e88d2ed9907b901845b0e8c1fae4613e`. These identifiers refer to that earlier work; they are not a DOI or checksum for this release.

## Certificates and formal source

`CERTIFICATE_MANIFEST.json` records certificate/source byte hashes, embedded input links and the explicit portability adaptation. Executed-source snapshots and all but one certificate JSON retain their original bytes. In the norm-cell record, only the historical command and working-directory text were sanitized to remove machine paths; the manifest records both digests and canonical equality of every other field. Downstream historical digests remain unchanged and are resolved through an explicit provenance alias, not presented as hashes of the sanitized file. Historical route identifiers inside those immutable records are provenance strings, not hidden filesystem dependencies. The public analytical account is `supplement.tex`. See `verification/CERTIFICATE_SCOPE.md` for what is replayed and what is accepted from the cached verified calculations.

`lean/lean-toolchain` fixes Lean 4.32.0. `lean/lakefile.toml` pins mathlib to commit `81a5d257c8e410db227a6665ed08f64fea08e997`; `lean/lake-manifest.json` locks its dependencies. Readable project source is included; installed compilers, dependencies and caches are excluded. Formal hypotheses and their relationship to the manuscript are documented in `verification/FORMAL_CORRESPONDENCE.md` and the Lean scope documentation.

## Rights

The existing `NOTICE` applies to this publication and packaging. No new open-content license is granted by public availability. Third-party dependencies are fetched under their own licenses and are not vendored here; see `THIRD_PARTY_NOTICES.md`. Public literature is cited, not redistributed as third-party PDFs. AI-assisted drafting, proof development and code preparation were used; the release reports executable checks and their actual limits rather than treating assistance as verification.

## Release identity

`MANIFEST.sha256` identifies the tracked publication files in the corresponding checkout or extracted archive. `SHA256SUMS` identifies the PDF and versioned source archive assets. Neither checksum file asserts mathematical truth. A subsequent DOI metadata update must be an explicit new commit or version; never silently replace the bytes of a frozen tagged release.

## Post-publication metadata and author attribution

On 13 October 2026, the Zenodo record [10.5281/zenodo.23337494](https://doi.org/10.5281/zenodo.23337494) was published with the two PDFs, complete source ZIP and original `SHA256SUMS`. The public API confirmed the title, author, version, publication state and all four file sizes and MD5 checksums against the frozen local assets. The alternative tar.gz archive remains available on the matching GitHub release. The concept DOI is [10.5281/zenodo.23337493](https://doi.org/10.5281/zenodo.23337493).

The frozen `v1.0.0` tag remains at `aadbdb0affb9e97124f68d2c37e388fdd2304345`, whose [CI run](https://github.com/LeonardSEO/odd-theta-weil-bounds/actions/runs/38374755999) passed. At the author's request, the five commits on the default branch were reattributed to the LeonardSEO account, producing `5ab7c958bcca1100c408f49f7809a4001da67be9`; their trees were unchanged. The tag and release asset bytes were preserved. This follow-up documentation commit records the issued DOI and updates the public README and citation metadata. Its source manifest covers the updated checkout; the manifests inside the frozen archives remain unchanged.
