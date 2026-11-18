/-
Equation (22): closed form of the leave-one-out prediction of ridge regression.

The regularised kernel `K + λ I` is written in block form for the left-out gene,
`[[a, bᵀ], [b, C]]`, with `C` the regularised kernel of the other training genes and `b` the
kernel values between the left-out gene and the others. The vector `(α₀, α')` solves the full
system with right-hand side `(y₀, y')`, so `α₀` is the entry of `(K + λ I)⁻¹ y` for the
left-out gene. The vector `(β₀, β')` solves the system with the unit vector as right-hand
side, so `β₀` is the diagonal entry of `(K + λ I)⁻¹`. The prediction of the model fitted
without the gene is `bᵀ C⁻¹ y'`.
-/
import Mathlib

open Matrix

set_option linter.unusedSectionVars false
set_option linter.unusedDecidableInType false
set_option linter.unusedFintypeInType false

noncomputable section

namespace Reciprocity

variable {m : Type*} [Fintype m] [DecidableEq m]

/-- Elimination of the block of the other genes (Schur complement). -/
theorem schur_first_row (C : Matrix m m ℝ) (hC : IsUnit C.det) (a x0 r0 : ℝ)
    (b x' r' : m → ℝ)
    (h1 : a * x0 + b ⬝ᵥ x' = r0) (h2 : x0 • b + C *ᵥ x' = r') :
    a * x0 + b ⬝ᵥ (C⁻¹ *ᵥ r') - x0 * (b ⬝ᵥ (C⁻¹ *ᵥ b)) = r0 := by
  have hx : x' = C⁻¹ *ᵥ (r' - x0 • b) := by
    have h : C *ᵥ x' = r' - x0 • b := by rw [← h2]; simp
    rw [← h, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul C hC, Matrix.one_mulVec]
  have hdot : b ⬝ᵥ x' = b ⬝ᵥ (C⁻¹ *ᵥ r') - x0 * (b ⬝ᵥ (C⁻¹ *ᵥ b)) := by
    rw [hx, Matrix.mulVec_sub, Matrix.mulVec_smul, dotProduct_sub, dotProduct_smul, smul_eq_mul]
  rw [← h1, hdot]; ring

/-- Equation (22): the prediction for a gene by the model fitted without it equals
`y₀ - α₀ / β₀`. -/
theorem loo_closed_form (C : Matrix m m ℝ) (hC : IsUnit C.det) (a y0 α0 β0 : ℝ)
    (b y' α' β' : m → ℝ)
    (h1 : a * α0 + b ⬝ᵥ α' = y0) (h2 : α0 • b + C *ᵥ α' = y')
    (h3 : a * β0 + b ⬝ᵥ β' = 1) (h4 : β0 • b + C *ᵥ β' = 0) :
    b ⬝ᵥ (C⁻¹ *ᵥ y') = y0 - α0 / β0 := by
  have hα := schur_first_row C hC a α0 y0 b α' y' h1 h2
  have hβ := schur_first_row C hC a β0 1 b β' 0 h3 h4
  rw [Matrix.mulVec_zero, dotProduct_zero, add_zero] at hβ
  have hβ0 : β0 ≠ 0 := by
    intro h; rw [h] at hβ; simp at hβ
  have hdiv : α0 / β0 = α0 * (a - b ⬝ᵥ (C⁻¹ *ᵥ b)) := by
    rw [div_eq_iff hβ0]; linear_combination (-α0) * hβ
  rw [hdiv]; linear_combination hα

end Reciprocity
