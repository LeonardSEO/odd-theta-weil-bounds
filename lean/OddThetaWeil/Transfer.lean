import OddThetaWeil.Metric
import OddThetaWeil.Enriched
import OddThetaWeil.Endpoints

namespace OddThetaWeil
open scoped RealInnerProductSpace
variable {H F : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]
  [NormedAddCommGroup F] [InnerProductSpace ℝ F]

/-- Complete form-domain transfer from the coefficient certificate. The original
Weil lower bound and identification of the shifted energy are explicit inputs;
the finite-rank bound itself is derived by coefficient_metric_to_form. -/
theorem certified_retained_energy {s ν θ ρ μ energy weil : ℝ}
    (T : CoerciveResolvent H s)
    (V : Submodule ℝ H) (hDV : T.domain ≤ V)
    (form : V →ₗ[ℝ] V →ₗ[ℝ] ℝ)
    (hsym : ∀ u v, form u v = form v u)
    (hpos : ∀ u, 0 ≤ form u u)
    (hrep : ∀ (y : T.domain) (u : V),
      form ⟨y, hDV y.property⟩ u = ⟪T.op y, (u : H)⟫)
    (C : F →ₗ[ℝ] H) (Cadj : H →ₗ[ℝ] F)
    (hadj : ∀ a u, ⟪C a, u⟫ = ⟪a, Cadj u⟫)
    (N Ninv : F →ₗ[ℝ] F) (hNinv : ∀ a, N (Ninv a) = a)
    (hν : 0 < ν) (hθ : 0 < θ)
    (hcoeff : ∀ a, ν * ⟪C a, (T.inverse (C a) : H)⟫ ≤ θ * ⟪a, N a⟫)
    (u : V) (henergy : form u u = energy + (μ - ρ) * ‖(u : H)‖ ^ 2)
    (hweil : energy - ρ * ‖(u : H)‖ ^ 2 -
      ν * ⟪C (Ninv (Cadj (u : H))), (u : H)⟫ ≤ weil) :
    (1 - θ) * energy - (ρ + θ * (μ - ρ)) * ‖(u : H)‖ ^ 2 ≤ weil := by
  have hr := coefficient_metric_to_form T V hDV form hsym hpos hrep
    C Cadj hadj N Ninv hNinv hν hθ hcoeff u
  rw [henergy] at hr
  exact Endpoints.retained_energy_loss hweil hr

/-- End-to-end conditional transfer: the finite input is the RE5 matrix upper
price, not the desired physical form bound or a resolvent estimate. -/
theorem enriched_certified_retained_energy {s ν θ ρ μ energy weil β τ ε g L : ℝ}
    (T : CoerciveResolvent H s)
    (V : Submodule ℝ H) (hDV : T.domain ≤ V)
    (form : V →ₗ[ℝ] V →ₗ[ℝ] ℝ)
    (hsym : ∀ u v, form u v = form v u)
    (hpos : ∀ u, 0 ≤ form u u)
    (hrep : ∀ (y : T.domain) (u : V),
      form ⟨y, hDV y.property⟩ u = ⟪T.op y, (u : H)⟫)
    (C : F →ₗ[ℝ] H) (Cadj : H →ₗ[ℝ] F)
    (hadj : ∀ a u, ⟪C a, u⟫ = ⟪a, Cadj u⟫)
    (N Ninv : F →ₗ[ℝ] F) (hNinv : ∀ a, N (Ninv a) = a)
    (Y : F → T.domain) (Ahat D Err : F → H)
    (hs : 0 < s) (hν : 0 < ν) (hθ : 0 < θ)
    (hβ : 0 < β) (hτ : 0 < τ) (hε : 0 ≤ ε) (hg : 0 ≤ g)
    (haction : ∀ a, T.op (Y a) = Ahat a + s • (Y a : H) + D a + Err a)
    (hd : ∀ a, ‖D a‖ ≤ ε * ‖(Y a : H)‖)
    (he : ∀ a, ‖Err a‖ ≤ g * ‖a‖)
    (hy : ∀ a, ‖(Y a : H)‖ ≤ L * ‖a‖)
    (hprice : ∀ a, ν * (
      (2 * ⟪C a, (Y a : H)⟫ - s * ‖(Y a : H)‖ ^ 2 - ⟪Ahat a, (Y a : H)⟫) +
      ε * ‖(Y a : H)‖ ^ 2 + g * L * ‖a‖ ^ 2 +
      ((1 + β) * ‖C a - Ahat a - s • (Y a : H)‖ ^ 2 +
       (1 + 1 / β) * ((1 + τ) * ε ^ 2 * ‖(Y a : H)‖ ^ 2 +
         (1 + 1 / τ) * g ^ 2 * ‖a‖ ^ 2)) / s) ≤ θ * ⟪a, N a⟫)
    (u : V) (henergy : form u u = energy + (μ - ρ) * ‖(u : H)‖ ^ 2)
    (hweil : energy - ρ * ‖(u : H)‖ ^ 2 -
      ν * ⟪C (Ninv (Cadj (u : H))), (u : H)⟫ ≤ weil) :
    (1 - θ) * energy - (ρ + θ * (μ - ρ)) * ‖(u : H)‖ ^ 2 ≤ weil := by
  have hcoeff : ∀ a, ν * ⟪C a, (T.inverse (C a) : H)⟫ ≤ θ * ⟪a, N a⟫ := by
    intro a
    have hu := enriched_resolvent_upper T (C a) (Y a) (Ahat a) (D a) (Err a) a
      hs hβ hτ hε hg (haction a) (hd a) (he a) (hy a)
    exact (mul_le_mul_of_nonneg_left hu (le_of_lt hν)).trans (hprice a)
  exact certified_retained_energy T V hDV form hsym hpos hrep C Cadj hadj N Ninv
    hNinv hν hθ hcoeff u henergy hweil

end OddThetaWeil
