/-
Equations (1) and (5): Lyapunov equation, probability flux and static response.
-/
import Mathlib

open Matrix

set_option linter.unusedSectionVars false
set_option linter.unusedDecidableInType false
set_option linter.unusedFintypeInType false

noncomputable section

namespace Reciprocity

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- The flux matrix `Q = (J Σ - Σ Jᵀ) / 2` of equation (5). -/
def flux (J Sg : Matrix n n ℝ) : Matrix n n ℝ := (1 / 2 : ℝ) • (J * Sg - Sg * Jᵀ)

/-- (a) `Q` is skew-symmetric when `Σ` is symmetric. -/
theorem flux_transpose (J Sg : Matrix n n ℝ) (hS : Sgᵀ = Sg) :
    (flux J Sg)ᵀ = -flux J Sg := by
  unfold flux
  rw [transpose_smul, transpose_sub, transpose_mul, transpose_mul, transpose_transpose, hS,
    ← smul_neg, neg_sub]

/-- (b) The Lyapunov equation `J Σ + Σ Jᵀ = -2D` gives `-J Σ = D - Q`. -/
theorem neg_J_mul_Sigma (J Sg D : Matrix n n ℝ)
    (hL : J * Sg + Sg * Jᵀ = -((2 : ℝ) • D)) :
    -(J * Sg) = D - flux J Sg := by
  ext i j
  have h := congrFun (congrFun hL i) j
  unfold flux
  simp only [Matrix.add_apply, Matrix.sub_apply, Matrix.neg_apply, Matrix.smul_apply,
    smul_eq_mul] at h ⊢
  linarith

/-- (c) Inverse of the product: `(-J Σ)⁻¹ = -Σ⁻¹ J⁻¹`. -/
theorem inv_neg_J_mul_Sigma (J Sg : Matrix n n ℝ) (hJ : IsUnit J.det) (hSg : IsUnit Sg.det) :
    (-(J * Sg))⁻¹ = -(Sg⁻¹ * J⁻¹) := by
  apply Matrix.inv_eq_right_inv
  rw [neg_mul_neg, Matrix.mul_assoc, ← Matrix.mul_assoc Sg, Matrix.mul_nonsing_inv Sg hSg,
    Matrix.one_mul, Matrix.mul_nonsing_inv J hJ]

/-- (c) `-Σ⁻¹ J⁻¹ = (D - Q)⁻¹`. -/
theorem inv_D_sub_flux (J Sg D : Matrix n n ℝ)
    (hL : J * Sg + Sg * Jᵀ = -((2 : ℝ) • D)) (hJ : IsUnit J.det) (hSg : IsUnit Sg.det) :
    (D - flux J Sg)⁻¹ = -(Sg⁻¹ * J⁻¹) := by
  rw [← neg_J_mul_Sigma J Sg D hL]
  exact inv_neg_J_mul_Sigma J Sg hJ hSg

/-- (c) Equation (5): `-J⁻¹ = Σ (D - Q)⁻¹`. -/
theorem response_eq_Sigma_mul_inv (J Sg D : Matrix n n ℝ)
    (hL : J * Sg + Sg * Jᵀ = -((2 : ℝ) • D))
    (hJ : IsUnit J.det) (hP : IsUnit (D - flux J Sg).det) :
    -J⁻¹ = Sg * (D - flux J Sg)⁻¹ := by
  have hb := neg_J_mul_Sigma J Sg D hL
  have hSg : Sg = (-J⁻¹) * (D - flux J Sg) := by
    rw [← hb, neg_mul_neg, ← Matrix.mul_assoc, Matrix.nonsing_inv_mul J hJ, Matrix.one_mul]
  calc -J⁻¹ = (-J⁻¹) * ((D - flux J Sg) * (D - flux J Sg)⁻¹) := by
        rw [Matrix.mul_nonsing_inv _ hP, Matrix.mul_one]
    _ = ((-J⁻¹) * (D - flux J Sg)) * (D - flux J Sg)⁻¹ := by rw [Matrix.mul_assoc]
    _ = Sg * (D - flux J Sg)⁻¹ := by rw [← hSg]

/-- (d) Detailed balance, `Q = 0`: the response is `Σ D⁻¹`. -/
theorem response_detailed_balance (J Sg D : Matrix n n ℝ)
    (hL : J * Sg + Sg * Jᵀ = -((2 : ℝ) • D))
    (hJ : IsUnit J.det) (hD : IsUnit D.det) (hQ : flux J Sg = 0) :
    -J⁻¹ = Sg * D⁻¹ := by
  have h := response_eq_Sigma_mul_inv J Sg D hL hJ (by rw [hQ, sub_zero]; exact hD)
  rw [hQ, sub_zero] at h
  exact h

/-- (d) With `Q = 0` and isotropic noise `D = d • 1`, the response is `d⁻¹ • Σ` and symmetric. -/
theorem response_symmetric_of_detailed_balance (J Sg : Matrix n n ℝ) (d : ℝ) (hd : d ≠ 0)
    (hS : Sgᵀ = Sg)
    (hL : J * Sg + Sg * Jᵀ = -((2 : ℝ) • (d • (1 : Matrix n n ℝ))))
    (hJ : IsUnit J.det) (hQ : flux J Sg = 0) :
    -J⁻¹ = d⁻¹ • Sg ∧ (-J⁻¹)ᵀ = -J⁻¹ := by
  have hb := neg_J_mul_Sigma J Sg _ hL
  rw [hQ, sub_zero] at hb
  have h1 : -Sg = d • J⁻¹ := by
    have h : J⁻¹ * (-(J * Sg)) = J⁻¹ * (d • (1 : Matrix n n ℝ)) := by rw [hb]
    rw [Matrix.mul_neg, ← Matrix.mul_assoc, Matrix.nonsing_inv_mul J hJ, Matrix.one_mul,
      Matrix.mul_smul, Matrix.mul_one] at h
    exact h
  have h2 : -J⁻¹ = d⁻¹ • Sg := by
    have h3 : Sg = d • (-J⁻¹) := by rw [smul_neg, ← h1, neg_neg]
    conv_rhs => rw [h3]
    rw [smul_smul, inv_mul_cancel₀ hd, one_smul]
  refine ⟨h2, ?_⟩
  rw [h2, transpose_smul, hS]

end Reciprocity
