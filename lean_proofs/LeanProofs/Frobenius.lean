/-
Equations (6) and (7): orthogonal decomposition into symmetric and skew-symmetric parts and
the bound on the variance explained by a symmetric prediction.
Matrices are functions `n → n → ℝ`. Sums run over all entries, which equals the sum over
`i ≠ j` for matrices with zero diagonal.
-/
import Mathlib

open Finset

set_option linter.unusedSectionVars false
set_option linter.unusedDecidableInType false
set_option linter.unusedFintypeInType false

noncomputable section

namespace Reciprocity

variable {n : Type*} [Fintype n]

/-- Frobenius inner product. -/
def ip (X Y : n → n → ℝ) : ℝ := ∑ i, ∑ j, X i j * Y i j

/-- Squared Frobenius norm. -/
def nsq (X : n → n → ℝ) : ℝ := ip X X

/-- Symmetric part. -/
def symPart (M : n → n → ℝ) : n → n → ℝ := fun i j => (M i j + M j i) / 2

/-- Skew-symmetric part. -/
def skewPart (M : n → n → ℝ) : n → n → ℝ := fun i j => (M i j - M j i) / 2

theorem nsq_nonneg (X : n → n → ℝ) : 0 ≤ nsq X :=
  Finset.sum_nonneg fun i _ => Finset.sum_nonneg fun j _ => mul_self_nonneg (X i j)

theorem sym_add_skew (M : n → n → ℝ) (i j : n) : M i j = symPart M i j + skewPart M i j := by
  unfold symPart skewPart; ring

theorem symPart_symm (M : n → n → ℝ) (i j : n) : symPart M j i = symPart M i j := by
  unfold symPart; ring

theorem skewPart_skew (M : n → n → ℝ) (i j : n) : skewPart M j i = -skewPart M i j := by
  unfold skewPart; ring

/-- A symmetric and a skew-symmetric matrix are orthogonal. -/
theorem ip_symm_skew (S A : n → n → ℝ) (hS : ∀ i j, S j i = S i j)
    (hA : ∀ i j, A j i = -A i j) : ip S A = 0 := by
  have h1 : ip S A = ∑ i, ∑ j, S j i * A j i := by
    unfold ip; exact Finset.sum_comm
  have h2 : ∑ i, ∑ j, S j i * A j i = -ip S A := by
    unfold ip
    rw [← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [hS i j, hA i j]; ring
  linarith

/-- Pythagoras for an orthogonal sum. -/
theorem nsq_of_add (Z X Y : n → n → ℝ) (h : ∀ i j, Z i j = X i j + Y i j)
    (hXY : ip X Y = 0) : nsq Z = nsq X + nsq Y := by
  have key : nsq Z = nsq X + 2 * ip X Y + nsq Y := by
    unfold nsq ip
    simp only [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [h i j]; ring
  rw [key, hXY]; ring

/-- (a) `⟨S, A⟩ = 0`. -/
theorem ip_symPart_skewPart (M : n → n → ℝ) : ip (symPart M) (skewPart M) = 0 :=
  ip_symm_skew _ _ (symPart_symm M) (skewPart_skew M)

/-- (b) `‖M‖² = ‖S‖² + ‖A‖²`. -/
theorem nsq_eq_sym_add_skew (M : n → n → ℝ) :
    nsq M = nsq (symPart M) + nsq (skewPart M) :=
  nsq_of_add M _ _ (sym_add_skew M) (ip_symPart_skewPart M)

/-- (c) For a symmetric prediction, `‖M - M̂‖² = ‖S - M̂‖² + ‖A‖²`. -/
theorem nsq_error_symmetric (M Mh : n → n → ℝ) (hMh : ∀ i j, Mh j i = Mh i j) :
    nsq (fun i j => M i j - Mh i j)
      = nsq (fun i j => symPart M i j - Mh i j) + nsq (skewPart M) := by
  refine nsq_of_add _ _ _ (fun i j => ?_) (ip_symm_skew _ _ (fun i j => ?_) (skewPart_skew M))
  · show M i j - Mh i j = symPart M i j - Mh i j + skewPart M i j
    unfold symPart skewPart; ring
  · show symPart M j i - Mh j i = symPart M i j - Mh i j
    rw [symPart_symm M i j, hMh i j]

/-- (c) The error of a symmetric prediction is at least `‖A‖²`. -/
theorem nsq_skew_le_error (M Mh : n → n → ℝ) (hMh : ∀ i j, Mh j i = Mh i j) :
    nsq (skewPart M) ≤ nsq (fun i j => M i j - Mh i j) := by
  rw [nsq_error_symmetric M Mh hMh]
  linarith [nsq_nonneg (fun i j => symPart M i j - Mh i j)]

/-- (d) Equation (7): the reciprocity bound. -/
theorem reciprocity_bound (M Mh : n → n → ℝ) (hMh : ∀ i j, Mh j i = Mh i j)
    (hpos : 0 < nsq M) :
    1 - nsq (fun i j => M i j - Mh i j) / nsq M
      ≤ nsq (symPart M) / (nsq (symPart M) + nsq (skewPart M)) := by
  have hM := nsq_eq_sym_add_skew M
  have hc := nsq_skew_le_error M Mh hMh
  rw [← hM, sub_le_iff_le_add, ← add_div, one_le_div hpos]
  linarith

/-- (d) The bound equals one minus the skew-symmetric share. -/
theorem sym_share_eq (s a : ℝ) (h : s + a ≠ 0) : s / (s + a) = 1 - a / (s + a) := by
  field_simp; ring

end Reciprocity
