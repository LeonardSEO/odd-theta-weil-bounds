# Certificate verification and source closure

Run from any working directory:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-certificate.txt
.venv/bin/python scripts/check_certificates.py
```

The pinned verification runtime is NumPy 2.5.3 and python-flint 0.9.0 (the checked runtime uses FLINT 3.6.0). The checker sets Arb precision to 192 bits. Run Python without `-O`; assertions are the certificate checks. The compatibility entry point `computations/verify_candidate.py` runs the same checker. Neither entry point writes certificates.

## What the bounded check establishes

- SHA-256 and byte identity for every computational release file in the manifest; every embedded input, generator, executed-source, owner, and matrix digest resolves to included bytes or the explicitly declared provenance-only adaptation below.
- Exact-input Arb residual/polar spectral guards for the original 56-dimensional coupled matrices, with cached full-error allowances, and exact rational endpoint transfer.
- Exact rational endpoint arithmetic for the 231 existing enriched rows. It reconstructs the four published enriched rows (best scalar, best retained energy, and the two reported intermediate tradeoffs) from cached `GY`, `CY`, `YA`, `CA`, `AG`, and physical matrices and full saved errors, then proves the published rational thresholds directly: every interval LDL pivot of `(theta * coefficient_lower - complete_error) I - symmetric(comparison_center)` must be strictly positive. All original input and binary64 assembly error terms are retained. This strict finite Loewner check avoids a fresh eigensolver/polar surcharge for the four published rows; their very small endpoint margins had made that surcharge platform-sensitive. It does not rerun the parameter search or certify an optimum outside the recorded admissible finite candidates.
- Exact-input five-dimensional physical/signed spectral guards and outside coupling bounds from cached signed centers, and the small fixed-price scalar barrier arithmetic.
- A fresh finite ST2 inverse and Rayleigh witness check, including the declared Frobenius perturbation price. This is a barrier for a sufficient comparison method, not a negative direction of the actual signed form.

## Immutable accepted premises

The release includes the source and data dependency closure, but this bounded check does not establish every upstream numerical premise from scratch. The following remain explicitly accepted from the historical certificate records and the analytic arguments in `supplement.tex`:

1. The upstream whole-domain odd-theta identity, equality/score compression, zero counting and tail estimates cited by the manuscript. The upstream publication identifier and source digest are recorded in the manifest; that publication is a mathematical premise, not a silently downloaded runtime input.
2. Certified scalar quadratures and theta/Fourier tails, the physical frame and coefficient-metric error bounds, the 2469 zero enclosures and completeness claims, and the cached low-zero Taylor evaluations.
3. The large physical/zero Gram matrices and their full operator error allowances, continuous-to-midpoint errors, prime and zero tail allowances, and the complete enriched action enclosure, including the original native DFT and whole-action error calculations.
4. The saved signed-center error allowances and whole signed-action error, and the ST2 physical/mixed/Gram uncertainty reserves. The checker verifies their finite consequences and byte lineage; hashes alone do not prove these error bounds.

The JSON records name the precise input hashes. `certificate_hash_links` resolves every such digest, including abbreviated owner keys, to an included file. For the sanitized cell certificate, the explicit `provenance_adaptations` entry preserves its original digest as a historical alias; its release-byte digest is separately checked. Missing inputs are fatal. The tighter-resolvent JSON, generator and executed source were added because the pilot generator reads that JSON; its matrices were already included. No archived research tree is required at runtime.

## Source lineage and generation

Executed-source snapshots and all certificate JSON except the cell certificate are preserved byte for byte. In `complement_signed_symbol_norm_cells.json`, only `execution_provenance.command` and `execution_provenance.cwd` were replaced with portable descriptions. All other fields, including every numerical field, are unchanged. Downstream original input digests are preserved. The manifest records original and release digests plus a canonical payload digest excluding exactly those two fields. Release preparation compared the original and sanitized canonical payload bytes for equality; the checker recomputes that payload digest. The original private paths and original unredacted JSON are not shipped. Therefore this adaptation is an explicitly documented provenance equivalence, not a claim that the sanitized bytes satisfy the historical SHA-256 or a standalone cryptographic proof of the omitted original provenance. A certificate's historical `source_file`, `source`, `owning_note`, execution-log identifiers or reproduction text may refer to the original workspace or route notes. These are provenance text, not runtime imports or additional proof dependencies. The public analytical account is `supplement.tex`.

The one portable source adaptation replaces the ST2 check's hardcoded `li-positivity/certificates/` working-directory path with a path relative to its own release file. Its original bytes are retained as `certificates/complement_enriched_signed_tail_capture_executed_source.py`; `source_adaptations` explicitly records both hashes. The original certificate is not rewritten to pretend that the adapted script produced the old record.

`routes/complement-residual-enrichment.md` is retained only as an immutable analytical source snapshot whose hash is recorded by three certificates and whose path is read by their generators. It is not an additional public supplement. The barrier generator and its scalar JSON remain in `routes/` to preserve their source hash and relative path conventions. All other historical route reports are omitted after their mathematical content has been incorporated into the supplement.

The `computations/` generators are historical source, not the safe verification command. Original generator assertions referring to the cell certificate intentionally retain its historical digest; direct generator execution against the sanitized file requires an explicit provenance-aware adaptation. The authoritative checker uses the manifest's validated historical-digest resolver for finite assembly, without changing mathematical instructions or downstream certificate records. They can overwrite saved certificates and some perform expensive integration, zero, DFT or matrix work. Several generators also require SciPy (`scipy.linalg.eigh`, `scipy.special.digamma`, `scipy.signal.fftconvolve`); SciPy 1.18.1 is pinned separately in `requirements-generation.txt`; it is not needed by the read-only checker and is intentionally not installed by `requirements-certificate.txt`. These pins describe the existing local runtime, not a fresh-install test. The binary trial coefficients are the exact dyadic trial used for verification; a new floating-point generation run may choose different locator/trial bytes on a different BLAS platform. Full byte-identical regeneration across platforms has not been established.

For provenance or packaging changes alone, `python3 scripts/check_certificates.py --hash-only` verifies release bytes, all hash links, source adaptations, and the canonical payload commitment without rerunning numerical checks. The previously successful numerical pass remains valid because this adaptation changes no mathematical fields.
