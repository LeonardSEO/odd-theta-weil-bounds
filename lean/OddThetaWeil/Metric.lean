import OddThetaWeil.Resolvent

namespace OddThetaWeil
open scoped RealInnerProductSpace
variable {H F : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]
  [NormedAddCommGroup F] [InnerProductSpace ℝ F]

/-- Coefficient-metric comparison implies the full form-domain inequality.
`N` is the actual coefficient metric, not the physical Gram. The form
representation and existence of its coercive resolvent are analytic inputs. -/
theorem coefficient_metric_to_form {s ν θ : ℝ}
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
    (u : V) :
    ν * ⟪C (Ninv (Cadj (u : H))), (u : H)⟫ ≤ θ * form u u := by
  let a := Ninv (Cadj (u : H))
  let v : V := ⟨T.inverse (C a), hDV (T.inverse (C a)).property⟩
  have hvu : form v u = ⟪C a, (u : H)⟫ := by
    simpa [v, T.right_inverse] using hrep (T.inverse (C a)) u
  have hvv : form v v = ⟪C a, (T.inverse (C a) : H)⟫ := by
    simpa [v, T.right_inverse] using hrep (T.inverse (C a)) v
  have hmetric : ⟪a, N a⟫ = ⟪C a, (u : H)⟫ := by
    rw [hadj]
    simp only [a, hNinv]
  have hc := hcoeff a
  rw [hmetric, ← hvv, ← hvu] at hc
  have hp := hpos (ν • v - θ • u)
  simp only [map_sub, map_smul, LinearMap.sub_apply, LinearMap.smul_apply,
    smul_eq_mul] at hp
  rw [hsym u v] at hp
  have hb : ν * form v u ≤ θ * form u u := by
    nlinarith [mul_le_mul_of_nonneg_left hc (le_of_lt hν)]
  simpa only [hvu, a] using hb

end OddThetaWeil
