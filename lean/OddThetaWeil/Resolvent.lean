import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-! Domain-explicit coercive resolvent algebra. A complex Hilbert space is
regarded as a real inner product space: these are precisely real parts of
its Hermitian quadratic forms. No boundedness of the unbounded operator T
is asserted. Existence of its inverse is an explicit analytic input. -/
namespace OddThetaWeil
open scoped RealInnerProductSpace

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

structure CoerciveResolvent (H : Type*) [NormedAddCommGroup H]
    [InnerProductSpace ℝ H] (s : ℝ) where
  domain : Submodule ℝ H
  op : domain →ₗ[ℝ] H
  inverse : H →ₗ[ℝ] domain
  right_inverse : ∀ x, op (inverse x) = x
  left_inverse : ∀ y, inverse (op y) = y
  symmetric : ∀ (x y : domain), ⟪op x, (y : H)⟫ = ⟪(x : H), op y⟫
  coercive : ∀ (y : domain), s * ⟪(y : H), (y : H)⟫ ≤ ⟪op y, (y : H)⟫

variable {s : ℝ} (T : CoerciveResolvent H s)

/-- The inverse quadratic bound is derived from coercivity, not assumed. -/
theorem inverse_bound (hs : 0 < s) (r : H) :
    ⟪r, (T.inverse r : H)⟫ ≤ ⟪r,r⟫ / s := by
  have hc := T.coercive (T.inverse r)
  rw [T.right_inverse] at hc
  have hp := real_inner_self_nonneg (x := r - s • (T.inverse r : H))
  simp only [inner_sub_left, inner_sub_right, real_inner_smul_left,
    real_inner_smul_right] at hp
  rw [real_inner_comm r (T.inverse r : H)] at hp
  apply (le_div_iff₀ hs).2
  nlinarith [mul_le_mul_of_nonneg_left hc (le_of_lt hs)]

/-- Exact completion of the full ambient residual; y must be in the operator domain. -/
theorem full_residual_identity (c : H) (y : T.domain) :
    ⟪c, (T.inverse c : H)⟫ =
      2 * ⟪c, (y : H)⟫ - ⟪T.op y, (y : H)⟫ +
        ⟪c - T.op y, (T.inverse (c - T.op y) : H)⟫ := by
  have hy : T.inverse (c - T.op y) = T.inverse c - y := by
    rw [map_sub, T.left_inverse]
  rw [hy]
  simp only [Submodule.coe_sub, inner_sub_left, inner_sub_right]
  have hsym := T.symmetric y (T.inverse c)
  rw [T.right_inverse, real_inner_comm c (y : H)] at hsym
  rw [hsym]
  ring

/-- Rigorous upper bound retaining the residual outside every trial subspace. -/
theorem full_residual_upper (hs : 0 < s) (c : H) (y : T.domain) :
    ⟪c, (T.inverse c : H)⟫ ≤
      2 * ⟪c, (y : H)⟫ - ⟪T.op y, (y : H)⟫ +
        ⟪c - T.op y, c - T.op y⟫ / s := by
  rw [full_residual_identity T c y]
  exact add_le_add_right (inverse_bound T hs _) _

end OddThetaWeil
