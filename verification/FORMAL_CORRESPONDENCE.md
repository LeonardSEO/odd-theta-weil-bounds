# Formal correspondence and verification obligations

Verdict: the supplied abstract algebraic implication and exact endpoint arithmetic are formally verified in the pinned environment, under the explicit analytic/certificate interfaces described below. The complete theta–Weil theorem is not formalized. Independent intended-statement, source, elaborated-type and transitive-axiom review is complete.

## Independently recovered public statement

Sources: `paper.tex` (main theorem, matrix residual lemma, enrichment section), `routes/complement-residual-enrichment.md` (RE1–RE5), the action and shift audits, and their JSON certificates. These were read before inspecting any candidate Lean implementation.

The physical space is the complex odd subspace of L²(R,dx), equivalently the positive-half-line norm with measure 2dx. The public theorem quantifies over every original form-domain u and every closed verified K contained in ker A_W and containing the score q. With z=(I-P_K)u, it bounds a_W(u) below by c E8(z)−alpha ||z||². It imposes no additional finite orthogonality conditions. Formalizing a finite real matrix alone cannot prove this whole-domain statement.

The three advertised pairs derive from rho=979236349/2000000000000 and alpha=rho+theta(mu−rho):

| mu | beta | theta | retained energy | displayed loss |
|---|---|---|---|---|
| 1/25 | 3/2 | 867415467588/10^12 | 132584532412/10^12 | 0.034762 (upward) |
| 1/20 | 1277/512 | 787541154611/10^12 | 212458845389/10^12 | 0.039482 (upward) |
| 401/12800 | 1913/2048 | 499940839479/500000000000 | 59160521/500000000000 | 31324476155738682727829/10^24 (exact) |

The saved search has 231 candidate rows, including failed positive-energy candidates with theta≥1. “Least scalar loss” means among the passing candidates, not all rows or all possible shifts. In particular some failed rows have smaller alpha and negative retained energy.

## Required load-bearing correspondence

| Obligation | Exact intended content | Required formal boundary |
|---|---|---|
| Complex convention and domains | E is positive selfadjoint, potentially unbounded; s>0; C and Y map finite complex coefficient space into physical H, with Y columns in Dom(E). Form conclusions hold on Dom(E^(1/2)). | Record conjugate-linear/linear inner-product convention, real parts, operator versus form domains. A bounded or finite-dimensional specialization must be labeled. |
| Metric separation | N=Z*Z is positive on coefficient space; physical Gram H_C=C*C is distinct. C=BZ, nu=11/10. | Do not replace N by C*C, assume orthonormality of exact binary Z, or invert an enlarged trial Gram. |
| Shifted inverse | T=E+s has positive bounded inverse and <T^-1 r,r>≤||r||²/s. | Derive from positivity/coercivity, or explicitly classify an assumed inverse with its laws as an abstract interface not a constructed theta resolvent. |
| Completion of square | R=C−TY, K=C*Y+Y*C−Y*TY; C*T^-1 C=K+R*T^-1R. | Prove with actual maps/products and inverse laws. Assuming this identity as a premise only verifies downstream arithmetic. |
| Actual action bridge | EY=Ahat+Delta_R Y+D, Delta_R selfadjoint, ||Delta_R||≤epsilon_R, ||D||≤g, ||Y||≤y. | Preserve independent physical map bounds. Scalar/matrix Gram uncertainty is not a replacement for a physical map error. |
| Two Young estimates | ||Delta_R Yx+Dx||²≤(1+tau)epsilon_R²||Yx||²+(1+1/tau)g²||x||²; then split Rhat from that error with beta>0. | Prove both norm inequalities for complex vectors; tau,beta>0, g,epsilon_R,y≥0. |
| Enriched correction | K≤Khat+epsilon_R Y*Y+g y I, and RE5 combines the preceding bounds with factors 1/s. | Prove correction rather than assume the desired resolvent bound. The full residual Rhat=C−Ahat−sY is retained. |
| Shifted Gram expansion | Khat=CY+CY*−sG−(YA+YA*)/2; Rhat*Rhat=H_C−CA−CA*−s(CY+CY*)+AG+s(YA+YA*)+s²G. | AG=Ahat*Ahat is the whole physical action Gram; no projection onto a finite trial span. |
| Coefficient-to-form transfer | nu C*T^-1 C≤theta N implies nu CN^-1C*≤theta(E+s) as forms. | Prove via weighted duality/factorization; an assumed desired physical form inequality would omit this decisive step. |
| Numerical input lift | Complete matrix endpoint p<theta(1−10^-8), N≥(1−10^-8)I; every uncertainty multiplier is nonnegative. | Exact arithmetic may prove the implication from these analytic/finite-matrix premises. Actual matrices, all interval arithmetic, and their provenance remain external unless imported and checked as well. |
| Whole-domain floor and kernel transfer | a_W(z)≥E8(z)−<R56 z,z>−rho||z||² and a_W(u)=a_W(z). | An abstract scalar combination is useful but conditional; verify actual theta definitions, projection/domain preservation and kernel identity before claiming the public theorem itself is formalized. |

## Concrete realization that remains separate unless explicitly formalized

Theta positivity/smoothing and Fourier transform; the complete gamma/prime/contact/pole identity; q in Dom(G)∩ker A_W; all closed-domain/compression identities; exact B/V0/Z/Y realization; positive metric and physical rank 56; zero enclosures and infinite tails; strip and Poisson estimates; entire physical Fourier action/DFT Gram and full error ledger; native interval arithmetic and source/data hashes; the global rest price; and closed-kernel transfer. None is established merely by accepting a record of scalar bounds in Lean. The equality barrier, kernel-aware replacement, and signed-action cancellation are additional public results requiring their own targets to be called formalized.

## Actual declaration inspection

Inspected the complete sources of Resolvent, Metric, ActionError, Enriched, ComplexBridge and Transfer, plus the endpoint definitions and final public-number consequences.

| Checked source target | What its statement and proof actually establish |
|---|---|
| `CoerciveResolvent` | An explicit real-linear operator on a submodule domain, with a supplied two-sided inverse into that domain, symmetry and coercivity. It does not construct the resolvent of a theta operator or assume the desired inverse bound. |
| `inverse_bound` | Derives the inverse quadratic upper bound from coercivity and a nonnegative square. |
| `full_residual_identity`, `full_residual_upper` | Prove the exact square completion and upper bound for arbitrary ambient c and domain y. The residual is c−Ty in the whole physical space. |
| `ActionError.vector_young`, `norm_young`, `nested_residual` | Prove Young's inequality and its two-stage use from actual physical vector bounds. |
| `ActionError.inner_error`, `inner_cross_error` | Bound the real cross terms from the separate whole-action errors. |
| `enriched_resolvent_upper` | Composes the preceding results into RE5 at an arbitrary coefficient vector, with action identity Ty=Ahat+s y+d+e and physical norm bounds. It derives rather than assumes the inverse comparison. |
| `coefficient_metric_to_form` | Uses a supplied coefficient metric inverse Ninv, the actual adjoint identity, positive symmetric form and operator representation to derive the finite-rank bound on the form domain. It never identifies N with C*C. No finite-dimensional assumption is needed for this stronger abstract implication. |
| `complex_quadratic_order_iff`, `complex_linear_restriction`, `complex_metric_restriction` | Restriction of scalars retains all complex vectors and real parts of their Hermitian quadratic forms. The proofs do not restrict to real-coordinate vectors. A concrete complex theta instantiation is not supplied. |
| `certified_retained_energy` | Composes the coefficient comparison with the positive-form representation and original Weil floor into the energy/loss conclusion. The coefficient resolvent comparison is an explicit premise of this lower-level theorem. |
| `enriched_certified_retained_energy` | Composes the full chain: whole action split and norm error bounds + the explicit full RE5 quadratic upper-price certificate → enriched resolvent comparison → coefficient-metric-to-form duality → retained-energy loss. It does not accept the desired resolvent or physical finite-rank comparison as a premise. |
| `Endpoints.all_endpoints_checked` | Checks exact positivity of shifts/betas/thetas, retained/loss identities, and pass↔theta<1 for all 231 archived rows. This is not interval-matrix or theta verification. |
| `Endpoints.paper_first`, `paper_second`, `paper_third` | Prove the three exact public scalar consequences from the original floor and a supplied physical residual comparison; first two losses round upward. These are conditional scalar consequences, not declarations about a defined theta form. |

The explicit structure/argument premises were inspected in addition to the axiom checks. Inverse existence, form representation, the original floor and finite certificate premises are substantive external analytic inputs; none is silently proved by constructing a `CoerciveResolvent`. The initial modules do not define theta, Xi, E8, K, Z or their interval enclosures.

Independent standard-library check: parsed all 231 literal Lean records using exact rational arithmetic and compared each field and pass flag with its owning JSON row; all matched. The source SHA256 in `endpoint-source.json` also matched. There are 173 passing rows. This check did not regenerate Lean, rerun native interval calculations or use display floats.

## Execution and trust evidence

The initial umbrella build log records success under Lean 4.32.0 and mathlib commit `81a5d257c8e410db227a6665ed08f64fea08e997`, with one harmless `unnecessarySimpa` linter warning. The linked action-to-form composition was then added and independently inspected. The updated umbrella build also completed successfully. An independent check ran `lake env lean` on a temporary independent type-inspection file in the isolated build project, using the pinned 4.32.0 executable directory. It exited 0. The scratch file imported the release umbrella, printed the actual `CoerciveResolvent` structure, checked fully explicit types of the composed theorem, metric transfer, complex order equivalence and all three public-number consequences, then printed axioms of the composed theorem, aggregate endpoint check and three consequences. All five closures contain only `propext`, `Classical.choice`, and `Quot.sound`. The released `lean/axioms-output.txt` independently records the same allowed set for all 22 listed targets. No `sorryAx`, custom mathematical axiom or native-evaluation dependency occurs. Source inspection found no unsafe/kernel-bypass declaration. The only custom option raises elaboration recursion depth for the 231-row conjunction.

Required public wording: the formal result covers an abstract coercive-resolvent/action-error/coefficient-metric implication and exact endpoint arithmetic under explicit analytic/certificate premises. It does not formalize the entire public theta theorem, interval arithmetic, inverse existence, kernel construction, global spectral statement, or an equality barrier unless a separate checked barrier declaration is delivered.

Skipped: expensive numerical replay. Risk: full theta–Weil realization and native interval evidence remain outside the Lean kernel boundary even when all abstract declarations build without unsupported axioms.

## Exact remaining correspondence boundary

The composed theorem stops at an explicitly supplied *whole RE5 quadratic upper-price inequality*. The code does not define the cached matrix Grams and prove their polynomial expansion, replay polar spectral certificates, or prove that archived `price_ball` strings imply that inequality. The 231-row Lean check verifies rational loss identities and pass classifications; it does not prove the finite minimum/Pareto search assertion. Those are external certificate checks. The generic composed theorem and the three numeric scalar consequences are checked separately; no concrete theta object is instantiated inside Lean. This scope is mathematically meaningful and contains the decisive abstract action/resolvent/form chain, while leaving the actual analytical realization visible as an external obligation.

Final source correspondence check: `paper.tex`, “Formal verification scope,” now names the composed abstract chain, all 231 exact records and the 22 audited declarations, and explicitly excludes the complete theta specialization and equality barrier. The supplement's concluding replay paragraph distinguishes finite certificate replay from formalization. `lean/README.md` correctly explains domains, inverse existence, real restriction of the entire complex space, N versus C*C, whole physical residual, explicit finite RE5 input and the three numeric rows. These statements agree with the checked declarations. No equality-barrier Lean theorem was supplied or is now claimed.
