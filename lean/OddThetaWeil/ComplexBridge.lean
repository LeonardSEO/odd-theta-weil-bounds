import OddThetaWeil.Metric

namespace OddThetaWeil

noncomputable section ComplexBridge
variable {H F : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]
  [NormedAddCommGroup F] [InnerProductSpace ℂ F]
local instance : InnerProductSpace ℝ H := InnerProductSpace.rclikeToReal ℂ H
local instance : InnerProductSpace ℝ F := InnerProductSpace.rclikeToReal ℂ F

/-- Restricting scalars changes neither vectors nor the Hermitian quadratic order.
The quantifier still ranges over every complex vector, not a real-coordinate subset. -/
theorem complex_quadratic_order_iff (A B : H → H) :
    (∀ z : H, (inner ℂ z (A z)).re ≤ (inner ℂ z (B z)).re) ↔
    (∀ z : H, inner ℝ z (A z) ≤ inner ℝ z (B z)) := Iff.rfl

/-- Every complex coefficient is retained by the real-linear restriction. -/
theorem complex_linear_restriction (C : F →ₗ[ℂ] H) (a : F) :
    C.restrictScalars ℝ a = C a := rfl

/-- The restricted metric is exactly the real part of the original N quadratic. -/
theorem complex_metric_restriction (N : F →ₗ[ℂ] F) (a : F) :
    inner ℝ a (N.restrictScalars ℝ a) = (inner ℂ a (N a)).re := rfl
end ComplexBridge

end OddThetaWeil
