/-
Equations (3) and (4): in the limit of one slow mode the response is the rank-one matrix
`c₁ p_B q_A` and the covariance is the rank-one matrix `c₂ p_B p_A`.
-/
import Mathlib

set_option linter.unusedSectionVars false
set_option linter.unusedDecidableInType false
set_option linter.unusedFintypeInType false

noncomputable section

namespace Reciprocity

variable {n : Type*}

/-- Skew-symmetric part of the rank-one response. -/
theorem rankOne_response_skew (c : ℝ) (p q : n → ℝ) (a b : n) :
    (c * p b * q a - c * p a * q b) / 2 = c / 2 * (p b * q a - p a * q b) := by
  ring

/-- The rank-one covariance has no skew-symmetric part. -/
theorem rankOne_covariance_skew (c : ℝ) (p : n → ℝ) (a b : n) :
    (c * p b * p a - c * p a * p b) / 2 = 0 := by
  ring

/-- The skew-symmetric part of the rank-one response vanishes for every pair exactly when
`p` and `q` are linearly dependent. -/
theorem rankOne_skew_zero_iff (p q : n → ℝ) :
    (∀ a b, p b * q a - p a * q b = 0) ↔
      ∃ x y : ℝ, (x ≠ 0 ∨ y ≠ 0) ∧ ∀ i, x * p i + y * q i = 0 := by
  constructor
  · intro h
    by_cases hq : ∀ i, q i = 0
    · exact ⟨0, 1, Or.inr one_ne_zero, fun i => by rw [hq i]; ring⟩
    · obtain ⟨k, hk⟩ := not_forall.mp hq
      exact ⟨q k, -p k, Or.inl hk, fun i => by linear_combination h k i⟩
  · rintro ⟨x, y, hxy, h⟩ a b
    have ha := h a
    have hb := h b
    have h1 : x * (p b * q a - p a * q b) = 0 := by linear_combination q a * hb - q b * ha
    have h2 : y * (p b * q a - p a * q b) = 0 := by linear_combination p b * ha - p a * hb
    rcases hxy with hx | hy
    · exact (mul_eq_zero.mp h1).resolve_left hx
    · exact (mul_eq_zero.mp h2).resolve_left hy

/-- Equation (13): on the standardised scale the linear response ratio is the correlation. -/
theorem standardised_ratio (sBA sA sB : ℝ) (hA : sA ≠ 0) (hB : sB ≠ 0) :
    (sBA / (sA * sA)) * (sA / sB) = sBA / (sA * sB) := by
  field_simp

end Reciprocity
