import OddThetaWeil.Resolvent
import OddThetaWeil.ActionError

namespace OddThetaWeil
open scoped RealInnerProductSpace
variable {H F : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]
  [NormedAddCommGroup F]

/-- RE5 at an arbitrary coefficient vector a. Its whole residual and both
whole-action errors occur as vectors in H; no output projection is taken. -/
theorem enriched_resolvent_upper {s β τ ε g L : ℝ}
    (T : CoerciveResolvent H s) (c : H) (y : T.domain)
    (ahat d e : H) (a : F)
    (hs : 0 < s) (hβ : 0 < β) (hτ : 0 < τ)
    (hε : 0 ≤ ε) (hg : 0 ≤ g)
    (haction : T.op y = ahat + s • (y : H) + d + e)
    (hd : ‖d‖ ≤ ε * ‖(y : H)‖) (he : ‖e‖ ≤ g * ‖a‖)
    (hy : ‖(y : H)‖ ≤ L * ‖a‖) :
    ⟪c, (T.inverse c : H)⟫ ≤
      (2 * ⟪c, (y : H)⟫ - s * ‖(y : H)‖ ^ 2 - ⟪ahat, (y : H)⟫) +
      ε * ‖(y : H)‖ ^ 2 + g * L * ‖a‖ ^ 2 +
      ((1 + β) * ‖c - ahat - s • (y : H)‖ ^ 2 +
       (1 + 1 / β) * ((1 + τ) * ε ^ 2 * ‖(y : H)‖ ^ 2 +
         (1 + 1 / τ) * g ^ 2 * ‖a‖ ^ 2)) / s := by
  have hdinner := (abs_le.mp (ActionError.inner_error (y : H) d hd)).1
  have heinner := (abs_le.mp (ActionError.inner_cross_error (y : H) e a he)).1
  have heL : g * ‖(y : H)‖ * ‖a‖ ≤ g * L * ‖a‖ ^ 2 := by
    nlinarith [mul_le_mul_of_nonneg_left hy (mul_nonneg hg (norm_nonneg a))]
  have hk : 2 * ⟪c, (y : H)⟫ - ⟪T.op y, (y : H)⟫ ≤
      (2 * ⟪c, (y : H)⟫ - s * ‖(y : H)‖ ^ 2 - ⟪ahat, (y : H)⟫) +
      ε * ‖(y : H)‖ ^ 2 + g * L * ‖a‖ ^ 2 := by
    rw [haction]
    simp only [inner_add_left, real_inner_smul_left, real_inner_self_eq_norm_sq]
    rw [real_inner_comm (y : H) d, real_inner_comm (y : H) e]
    linarith
  have hr := ActionError.nested_residual (c - ahat - s • (y : H)) d e
    (y : H) a hβ hτ hε hg hd he
  have hid : c - T.op y = (c - ahat - s • (y : H)) - d - e := by
    rw [haction]; abel
  have hu := full_residual_upper T hs c y
  rw [real_inner_self_eq_norm_sq, hid] at hu
  exact hu.trans (add_le_add hk (div_le_div_of_nonneg_right hr (le_of_lt hs)))

end OddThetaWeil
