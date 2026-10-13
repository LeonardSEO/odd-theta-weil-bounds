# Third-party software and literature

The manuscript and release packaging retain the rights in `NOTICE`. That notice does not replace any third-party license. Dependency source and installed binaries are not vendored in this release.

| Dependency | Role | License / upstream notice |
|---|---|---|
| Lean 4.32.0 | Formal proof elaboration and kernel | [Apache-2.0](https://github.com/leanprover/lean4/blob/v4.32.0/LICENSE) |
| mathlib, pinned commit `81a5d257c8e410db227a6665ed08f64fea08e997` | Standard mathematical library | [Apache-2.0](https://github.com/leanprover-community/mathlib4/blob/81a5d257c8e410db227a6665ed08f64fea08e997/LICENSE); transitive package revisions in `lean/lake-manifest.json` retain their upstream licenses |
| python-flint 0.9.0 | Arb interval arithmetic bindings | [MIT](https://github.com/flintlib/python-flint/blob/0.9.0/LICENSE) |
| FLINT / Arb | Rigorous real and matrix arithmetic through python-flint | [Upstream distribution and licensing](https://flintlib.org/); consult the installed distribution's bundled notices |
| NumPy 2.5.3 | Finite arrays and floating locator computations | BSD-3-Clause and bundled component licenses (0BSD, MIT, Zlib, CC0-1.0), as recorded by the installed package; [upstream](https://github.com/numpy/numpy) |
| SciPy 1.18.1 | Optional historical certificate generation | BSD-3-Clause and bundled notices; [upstream](https://github.com/scipy/scipy). Not required by the read-only verifier |
| Tectonic 0.17.0 | PDF compilation during preparation | [Upstream](https://github.com/tectonic-typesetting/tectonic); engine, fonts and TeX packages retain their own notices |

No third-party publication PDF is redistributed. The bibliography identifies the standard residual, factorization, Schur-complement and Weil-form literature. The theta–arithmetic contribution and its limits are stated in the main article.
