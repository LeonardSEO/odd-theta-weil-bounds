import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

namespace OddThetaWeil.ActionError

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]
variable {F : Type*} [NormedAddCommGroup F]

/-- Vector Young inequality, from the nonnegative square of `β • x - y`. -/
theorem vector_young (x y : H) {β : ℝ} (hβ : 0 < β) :
    inner ℝ (x + y) (x + y) ≤
      (1 + β) * inner ℝ x x + (1 + 1 / β) * inner ℝ y y := by
  have hs : 0 ≤ inner ℝ (β • x - y) (β • x - y) := by
    simpa using (inner_self_nonneg (𝕜 := ℝ) (x := β • x - y))
  rw [real_inner_sub_sub_self] at hs
  simp only [inner_smul_left, inner_smul_right, conj_trivial] at hs
  rw [real_inner_add_add_self]
  apply (mul_le_mul_iff_right₀ hβ).mp
  field_simp
  nlinarith

 theorem norm_young (x y : H) {β : ℝ} (hβ : 0 < β) :
    ‖x + y‖ ^ 2 ≤ (1 + β) * ‖x‖ ^ 2 + (1 + 1 / β) * ‖y‖ ^ 2 := by
  simpa only [real_inner_self_eq_norm_sq] using vector_young x y hβ

/-- Both error components are controlled before the full residual is squared. -/
theorem nested_residual (r d e y : H) (a : F) {β τ ε g : ℝ}
    (hβ : 0 < β) (hτ : 0 < τ) (hε : 0 ≤ ε) (hg : 0 ≤ g)
    (hd : ‖d‖ ≤ ε * ‖y‖) (he : ‖e‖ ≤ g * ‖a‖) :
    ‖r - d - e‖ ^ 2 ≤ (1 + β) * ‖r‖ ^ 2 +
      (1 + 1 / β) * ((1 + τ) * ε ^ 2 * ‖y‖ ^ 2 +
        (1 + 1 / τ) * g ^ 2 * ‖a‖ ^ 2) := by
  have hd2 : ‖d‖ ^ 2 ≤ ε ^ 2 * ‖y‖ ^ 2 := by
    simpa only [mul_pow] using (sq_le_sq₀ (norm_nonneg d) (mul_nonneg hε (norm_nonneg y))).2 hd
  have he2 : ‖e‖ ^ 2 ≤ g ^ 2 * ‖a‖ ^ 2 := by
    simpa only [mul_pow] using (sq_le_sq₀ (norm_nonneg e) (mul_nonneg hg (norm_nonneg a))).2 he
  have hde := norm_young d e hτ
  have hr := norm_young r (-(d + e)) hβ
  have h1 : 0 ≤ 1 + τ := by positivity
  have h2 : 0 ≤ 1 + 1 / τ := by positivity
  have h3 : 0 ≤ 1 + 1 / β := by positivity
  have hde' : ‖d + e‖ ^ 2 ≤
      (1 + τ) * ε ^ 2 * ‖y‖ ^ 2 + (1 + 1 / τ) * g ^ 2 * ‖a‖ ^ 2 := by
    nlinarith [mul_le_mul_of_nonneg_left hd2 h1, mul_le_mul_of_nonneg_left he2 h2]
  have herr := mul_le_mul_of_nonneg_left hde' h3
  have hid : r - d - e = r + -(d + e) := by abel
  rw [hid]
  simp only [norm_neg] at hr
  linarith

/-- Cross-term control uses the full Hilbert norm of the action error. -/
theorem inner_error (y d : H) {ε : ℝ} (hd : ‖d‖ ≤ ε * ‖y‖) :
    |inner ℝ y d| ≤ ε * ‖y‖ ^ 2 := by
  calc
    |inner ℝ y d| ≤ ‖y‖ * ‖d‖ := abs_real_inner_le_norm y d
    _ ≤ ‖y‖ * (ε * ‖y‖) := mul_le_mul_of_nonneg_left hd (norm_nonneg y)
    _ = ε * ‖y‖ ^ 2 := by ring

 theorem inner_cross_error (y e : H) (a : F) {g : ℝ} (he : ‖e‖ ≤ g * ‖a‖) :
    |inner ℝ y e| ≤ g * ‖y‖ * ‖a‖ := by
  calc
    |inner ℝ y e| ≤ ‖y‖ * ‖e‖ := abs_real_inner_le_norm y e
    _ ≤ ‖y‖ * (g * ‖a‖) := mul_le_mul_of_nonneg_left he (norm_nonneg y)
    _ = g * ‖y‖ * ‖a‖ := by ring

end OddThetaWeil.ActionError
