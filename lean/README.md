# Checked formal core and exact endpoint arithmetic

This directory is an executable Lean 4.32.0 project, with mathlib pinned to
`81a5d257c8e410db227a6665ed08f64fea08e997` and its dependency lockfile.
It formalizes the domain-sensitive residual/metric argument and exact
arithmetic of the 231 saved candidates. It does **not** formalize the
complete theta analysis or independently verify Arb's interval calculations.
All such mathematical inputs are visible theorem arguments, never axioms.

Run from this directory, using the pinned Lean installation:

```sh
lake exe cache get
lake build
lake env lean Axioms.lean
python3 generate_endpoints.py --check
```

The final Python command checks the generated module and mapping against the
saved certificate without rewriting files. See the generator's command-line
interface if generating after an intentionally revised source certificate.

## Statement correspondence

| Manuscript / route assertion | Declaration | Exact scope |
|---|---|---|
| Coercive inverse quadratic upper bound | `inverse_bound` | Derived from `s>0` and coercivity using a nonnegative ambient square; the inverse estimate is not a hypothesis. |
| Full residual completion of the square | `full_residual_identity` | `c` is any ambient vector; trial `y` belongs to the explicit operator domain. |
| Full residual upper bound | `full_residual_upper` | Combines the exact identity with the derived inverse upper bound. No residual projection. |
| Original coefficient-metric form transfer | `coefficient_metric_to_form` | Uses the original `N`, `Ninv`, coefficient comparison, and operator/form representation on the whole form domain. |
| Two Young bounds / action-error prices | `ActionError.vector_young`, `norm_young`, `nested_residual`, `inner_error`, `inner_cross_error` | Proven from inner-product nonnegativity and whole-vector norm bounds. Coefficients and physical vectors inhabit independent spaces. |
| RE5 enriched resolvent comparison | `enriched_resolvent_upper` | Derives the complete RE5 quadratic upper bound including both action errors and the full ambient residual. |
| Original complex Hilbert interpretation | `complex_quadratic_order_iff`, `complex_linear_restriction`, `complex_metric_restriction` | Real scalar restriction retains every complex vector and exactly the real parts of all Hermitian quadratics. |
| Composed full RE5 certificate to original-form bound | `enriched_certified_retained_energy` | Takes action decomposition, norm error bounds and finite RE5 price; calls RE5, metric transfer and retained-energy theorem. |
| Original-form retained-energy/loss transfer | `certified_retained_energy` | Derives the finite-rank bound from the coefficient certificate before substituting it in the original lower bound. |
| All saved rational endpoints | `Endpoints.all_endpoints_checked`, `endpoint_count` | Kernel-checked exact rational identities and pass predicates for all 231 original records, including failed rows. |
| Three displayed main bounds | `Endpoints.paper_first`, `paper_second`, `paper_third` | Exact real arithmetic with explicit original lower-bound and residual premises; first two losses rounded strictly upward. |

Names in the table are in namespace `OddThetaWeil`. `Axioms.lean` audits
each of the 22 listed declarations, including the aggregate whose dependency
closure contains all 231 per-record proofs.

## Domain and metric mapping

For the paper instantiate the physical space with the complex odd Hilbert
space (or its quotient after the verified score/kernel projection), endowed
with `InnerProductSpace.rclikeToReal ℂ H`. This does not restrict the vectors
to real coordinates. `ComplexBridge.lean` checks the exact quadratic
correspondence and the unchanged underlying functions of complex-linear maps.

`CoerciveResolvent.domain` is `Dom(E_Y)`, as a real submodule of the same
physical space. `op` is `T=E_Y+s`, and `inverse` is the resolvent into that
operator domain. The structure explicitly requires both inverse identities,
symmetry on that domain, and `s ||y||² ≤ <Ty,y>`. A selfadjoint nonnegative
`E_Y` and `s>0` supply these inputs analytically. The formal proof neither
asserts that an unbounded operator acts on all vectors nor assumes a bounded
operator on the whole Hilbert space.

`V` in the metric/transfer theorems is the full form domain of `E_Y`.
`hDV` is the inclusion of the operator domain. `form` is the real part of
the shifted closed form. `hrep` is the operator/form representation theorem
on an operator-domain first argument and arbitrary form-domain second
argument. The resulting inequality is therefore not restricted to
operator-domain vectors or trial vectors.

Use the coefficient space `EuclideanSpace ℂ (Fin 56)`, restricted to real
scalars, with `C=BZ`, `Cadj=C*`, **`N=Z*Z`**, and its actual inverse.
`N` is never replaced by `C*C` or the identity. `hNinv` is the right inverse
identity; `hadj` is the real-part adjoint identity. Positivity/invertibility
of the actual coefficient metric is an analytic/certificate input. The
algebraic implication needs only its stated right-inverse identity and
coefficient comparison, so it makes no additional silent metric assumption.

For `enriched_resolvent_upper`, at each arbitrary coefficient vector `a`,
substitute `c=C a`, `y=Y a`, `ahat=Ahat a`, `d=Delta_R(Y a)`,
`e=D_gamma a`, and `L=||Y||_upper`. The hypothesis `haction` is
`(E_Y+s)Y a = Ahat a + s Y a + Delta_R(Y a) + D_gamma a`.
The two norm hypotheses follow from the full RE2/RE4 map bounds.
Its output is exactly the real quadratic evaluation of RE5, including
`||C a-Ahat a-sY a||²`. The coefficient space `F` need not equal `H`.
The symmetric cross Gram is the same real quadratic as `<Ahat a,Y a>`.

The numerical finite Gram enclosure must establish that `11/10` times
this RE5 upper expression is at most `theta <a,N a>` for **every** complex
coefficient `a`. Applying `enriched_resolvent_upper` gives `hcoeff` for
`coefficient_metric_to_form` or `certified_retained_energy`. The composed theorem
`enriched_certified_retained_energy` performs this complete chain in Lean,
taking the finite RE5 upper-price inequality as input. This finite
matrix enclosure is not assumed to be the desired Weil bound: the subsequent
inverse/metric/form implication is proved in Lean. Nevertheless, its Arb
calculation is not itself replayed by the Lean kernel.

For `certified_retained_energy`, use `nu=11/10`,
`rho=979236349/2000000000000`, `s=mu-rho`, `energy=E8(z)`, and the actual
Weil form value, where `z=(I-P_K)u`. `henergy` identifies the shifted form,
and `hweil` is the separately established original Weil lower bound.
The closed-kernel projection and its invariance are analytic inputs; no
finite-dimensional kernel assumption is introduced.

## Exact data and trust boundary

`generate_endpoints.py` reads the saved shift JSON using exact `Fraction`
arithmetic and emits literal rational Lean statements. `endpoint-source.json`
records its SHA-256 and every original zero-based source row. The Lean proofs
use `norm_num`/linear and polynomial arithmetic; none uses `native_decide`.
The three manuscript records are rows 59, 210, and 221. Failed candidates
remain failed records; these proofs do not claim a minimum over all 231
losses or over untested shifts. The arithmetic theorem checks consistency of
saved rational endpoints, not the outward matrix certificate behind them.

Unformalized analytic inputs: theta and Fourier identities/tails, original
closed forms and selfadjoint operators, positivity/coercivity and resolvent
existence, admissibility of the actual trial columns, whole-action error
bounds from RE2/RE4, physical Gram interval enclosures, coefficient metric
positivity, and original-form/kernel transfer premises. They are stated in
the manuscript/routes and consumed explicitly as theorem arguments here.

`build-output.txt` records the successful pinned build. `axioms-output.txt`
records actual `#print axioms` output. Only Lean's standard logical axioms
`propext`, `Classical.choice`, and `Quot.sound` are permitted in the reported
dependency closures. This is a checked conditional formal core, not a formal
proof of RH or a claim that the entire manuscript has been formalized.
