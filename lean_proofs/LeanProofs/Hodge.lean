/-
Equations (16) to (18): gene-level terms of the skew-symmetric and symmetric parts.
-/
import Mathlib
import LeanProofs.Frobenius

open Finset

set_option linter.unusedSectionVars false
set_option linter.unusedDecidableInType false
set_option linter.unusedFintypeInType false

noncomputable section

namespace Reciprocity

variable {n : Type*} [Fintype n]

/-- Gene-level term `p̂_i = (1/N) Σ_j A_ij` of equation (16). -/
def pot (A : n → n → ℝ) : n → ℝ := fun i => (∑ j, A i j) / (Fintype.card n : ℝ)

/-- The matrix `p_i - p_j` of a potential over genes. -/
def grad (p : n → ℝ) : n → n → ℝ := fun i j => p i - p j

/-- Pair-specific remainder `Ã_ij = A_ij - (p̂_i - p̂_j)`. -/
def resid (A : n → n → ℝ) : n → n → ℝ := fun i j => A i j - grad (pot A) i j

theorem sum_skew_zero (A : n → n → ℝ) (hA : ∀ i j, A j i = -A i j) : ∑ i, ∑ j, A i j = 0 := by
  have h1 : ∑ i, ∑ j, A i j = ∑ i, ∑ j, A j i := Finset.sum_comm
  have h2 : ∑ i, ∑ j, A j i = -∑ i, ∑ j, A i j := by
    rw [← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl fun j _ => hA i j
  linarith

/-- (a) The gene-level terms sum to zero. -/
theorem pot_sum_zero (A : n → n → ℝ) (hA : ∀ i j, A j i = -A i j) : ∑ i, pot A i = 0 := by
  unfold pot
  rw [← Finset.sum_div, sum_skew_zero A hA, zero_div]

theorem resid_skew (A : n → n → ℝ) (hA : ∀ i j, A j i = -A i j) (i j : n) :
    resid A j i = -resid A i j := by
  unfold resid grad
  rw [hA i j]; ring

theorem resid_row_sum [Nonempty n] (A : n → n → ℝ) (hA : ∀ i j, A j i = -A i j) (i : n) :
    ∑ j, resid A i j = 0 := by
  have hN : (Fintype.card n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  have hrow : ∑ j, A i j = (Fintype.card n : ℝ) * pot A i := by
    unfold pot; field_simp
  have h : ∑ j, resid A i j
      = ∑ j, A i j - (Fintype.card n : ℝ) * pot A i + ∑ j, pot A j := by
    unfold resid grad
    simp only [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    ring
  rw [h, pot_sum_zero A hA, hrow]; ring

theorem resid_col_sum [Nonempty n] (A : n → n → ℝ) (hA : ∀ i j, A j i = -A i j) (j : n) :
    ∑ i, resid A i j = 0 := by
  have h : ∑ i, resid A i j = -∑ i, resid A j i := by
    rw [← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl fun i _ => resid_skew A hA j i
  rw [h, resid_row_sum A hA j, neg_zero]

/-- A matrix with zero row sums and zero column sums is orthogonal to every potential term. -/
theorem ip_grad_zero (r : n → ℝ) (B : n → n → ℝ) (hr : ∀ i, ∑ j, B i j = 0)
    (hc : ∀ j, ∑ i, B i j = 0) : ip (grad r) B = 0 := by
  have h1 : ∑ i, ∑ j, r i * B i j = 0 := by
    apply Finset.sum_eq_zero; intro i _
    rw [← Finset.mul_sum, hr, mul_zero]
  have h2 : ∑ i, ∑ j, r j * B i j = 0 := by
    rw [Finset.sum_comm]
    apply Finset.sum_eq_zero; intro j _
    rw [← Finset.mul_sum, hc, mul_zero]
  have h3 : ip (grad r) B = ∑ i, ∑ j, r i * B i j - ∑ i, ∑ j, r j * B i j := by
    unfold ip grad
    simp only [← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    ring
  rw [h3, h1, h2, sub_zero]

/-- (b) The gene-level term and the pair-specific remainder are orthogonal. -/
theorem ip_grad_resid [Nonempty n] (A : n → n → ℝ) (hA : ∀ i j, A j i = -A i j) :
    ip (grad (pot A)) (resid A) = 0 :=
  ip_grad_zero _ _ (resid_row_sum A hA) (resid_col_sum A hA)

/-- (c) `‖A‖² = ‖G‖² + ‖Ã‖²`. -/
theorem nsq_eq_grad_add_resid [Nonempty n] (A : n → n → ℝ) (hA : ∀ i j, A j i = -A i j) :
    nsq A = nsq (grad (pot A)) + nsq (resid A) := by
  refine nsq_of_add A _ _ (fun i j => ?_) (ip_grad_resid A hA)
  unfold resid; ring

/-- (c) `p̂` minimises `Σ_ij (A_ij - (p_i - p_j))²` over all potentials. -/
theorem pot_minimises [Nonempty n] (A : n → n → ℝ) (hA : ∀ i j, A j i = -A i j) (p : n → ℝ) :
    nsq (resid A) ≤ nsq (fun i j => A i j - grad p i j) := by
  have h : nsq (fun i j => A i j - grad p i j)
      = nsq (grad (fun i => pot A i - p i)) + nsq (resid A) := by
    refine nsq_of_add _ _ _ (fun i j => ?_)
      (ip_grad_zero _ _ (resid_row_sum A hA) (resid_col_sum A hA))
    show A i j - grad p i j = grad (fun i => pot A i - p i) i j + resid A i j
    unfold resid grad; ring
  rw [h]
  linarith [nsq_nonneg (grad (fun i => pot A i - p i))]

/-! Equation (17): symmetric part, off-diagonal least squares. `S` has zero diagonal. -/

/-- `ĝ`, the mean off-diagonal entry. -/
def gHat (S : n → n → ℝ) : ℝ :=
  (∑ i, ∑ j, S i j) / ((Fintype.card n : ℝ) * ((Fintype.card n : ℝ) - 1))

/-- `û_i`. -/
def uHat (S : n → n → ℝ) : n → ℝ := fun i =>
  (∑ j, S i j - ((Fintype.card n : ℝ) - 1) * gHat S) / ((Fintype.card n : ℝ) - 2)

/-- The gene-level terms of the symmetric part sum to zero. -/
theorem uHat_sum_zero (S : n → n → ℝ) (hN : 2 < (Fintype.card n : ℝ)) : ∑ i, uHat S i = 0 := by
  have h0 : (Fintype.card n : ℝ) ≠ 0 := by intro h; linarith
  have h1 : (Fintype.card n : ℝ) - 1 ≠ 0 := by intro h; linarith
  have hg : (Fintype.card n : ℝ) * (((Fintype.card n : ℝ) - 1) * gHat S) = ∑ i, ∑ j, S i j := by
    unfold gHat; field_simp
  unfold uHat
  rw [← Finset.sum_div, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ,
    nsmul_eq_mul, hg, sub_self, zero_div]

/-- Normal equation of the off-diagonal least-squares fit for gene `i`: the remainder
`S_ij - ĝ - û_i - û_j` sums to zero over `j ≠ i`. The diagonal term, with `S_ii = 0`, is
subtracted from the full row sum. -/
theorem uHat_normal_equation (S : n → n → ℝ) (hN : 2 < (Fintype.card n : ℝ)) (i : n) :
    ∑ j, (S i j - gHat S - uHat S i - uHat S j) - (0 - gHat S - uHat S i - uHat S i) = 0 := by
  have h2 : (Fintype.card n : ℝ) - 2 ≠ 0 := by intro h; linarith
  have hu := uHat_sum_zero S hN
  have hrow : ∑ j, (S i j - gHat S - uHat S i - uHat S j)
      = ∑ j, S i j - (Fintype.card n : ℝ) * gHat S - (Fintype.card n : ℝ) * uHat S i
        - ∑ j, uHat S j := by
    simp only [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  have hui : ((Fintype.card n : ℝ) - 2) * uHat S i
      = ∑ j, S i j - ((Fintype.card n : ℝ) - 1) * gHat S := by
    unfold uHat; field_simp
  rw [hrow, hu]
  linarith

end Reciprocity
