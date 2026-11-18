/-
Equations (14) and (15): measurement error inflates the squared norms of both parts within
one screen by the same amount, and the cross product of two screens is unbiased.

The expectation is a linear functional `E` on random variables `Ω → ℝ` with `E 1 = 1`.
`part 1` is the symmetric part and `part (-1)` the skew-symmetric part. The sums run over all
entries and the error variance on the diagonal is zero, which is the sum over `i ≠ j`.
-/
import Mathlib

open Finset

set_option linter.unusedSectionVars false
set_option linter.unusedDecidableInType false
set_option linter.unusedFintypeInType false

noncomputable section

namespace Reciprocity

variable {n : Type*} [Fintype n] {Ω : Type*}

/-- `(X_ij + c X_ji) / 2`: symmetric part for `c = 1`, skew-symmetric part for `c = -1`. -/
def part (c : ℝ) (X : n → n → ℝ) : n → n → ℝ := fun i j => (X i j + c * X j i) / 2

theorem part_one (X : n → n → ℝ) (i j : n) : part 1 X i j = (X i j + X j i) / 2 := by
  unfold part; ring

theorem part_neg_one (X : n → n → ℝ) (i j : n) : part (-1) X i j = (X i j - X j i) / 2 := by
  unfold part; ring

/-- Expectation of a product of two noisy entries of one part. -/
theorem E_prod (E : (Ω → ℝ) →ₗ[ℝ] ℝ) (hE1 : E (fun _ => (1 : ℝ)) = 1)
    (c m1 m2 : ℝ) (e f e' f' : Ω → ℝ)
    (he : E e = 0) (hf : E f = 0) (he' : E e' = 0) (hf' : E f' = 0) :
    E (fun ω => (m1 + e ω + c * (m2 + f ω)) / 2 * ((m1 + e' ω + c * (m2 + f' ω)) / 2))
      = ((m1 + c * m2) / 2) ^ 2
        + (E (fun ω => e ω * e' ω) + c * E (fun ω => e ω * f' ω)
          + c * E (fun ω => f ω * e' ω) + c ^ 2 * E (fun ω => f ω * f' ω)) / 4 := by
  have h : (fun ω => (m1 + e ω + c * (m2 + f ω)) / 2 * ((m1 + e' ω + c * (m2 + f' ω)) / 2))
      = (((m1 + c * m2) / 2) ^ 2) • (fun _ => (1 : ℝ))
        + ((m1 + c * m2) / 4) • e + ((m1 + c * m2) * c / 4) • f
        + ((m1 + c * m2) / 4) • e' + ((m1 + c * m2) * c / 4) • f'
        + (1 / 4 : ℝ) • (fun ω => e ω * e' ω) + (c / 4) • (fun ω => e ω * f' ω)
        + (c / 4) • (fun ω => f ω * e' ω) + (c ^ 2 / 4) • (fun ω => f ω * f' ω) := by
    funext ω
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul]
    ring
  rw [h]
  simp only [map_add, map_smul, smul_eq_mul, hE1, he, hf, he', hf']
  ring

/-- The expectation commutes with a finite double sum. -/
theorem E_sum (E : (Ω → ℝ) →ₗ[ℝ] ℝ) (g : n → n → Ω → ℝ) :
    E (fun ω => ∑ i, ∑ j, g i j ω) = ∑ i, ∑ j, E (g i j) := by
  have h : (fun ω => ∑ i, ∑ j, g i j ω) = ∑ i, ∑ j, g i j := by
    funext ω
    simp only [Finset.sum_apply]
  rw [h, map_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [map_sum]

theorem sum_split (P s : n → n → ℝ) :
    ∑ i, ∑ j, (P i j + (s i j + s j i) / 4) = ∑ i, ∑ j, P i j + (1 / 2) * ∑ i, ∑ j, s i j := by
  have hswap : ∑ i, ∑ j, s j i = ∑ i, ∑ j, s i j := Finset.sum_comm
  have h : ∑ i, ∑ j, (P i j + (s i j + s j i) / 4)
      = ∑ i, ∑ j, P i j + (1 / 4) * ∑ i, ∑ j, s i j + (1 / 4) * ∑ i, ∑ j, s j i := by
    simp only [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    ring
  rw [h, hswap]; ring

/-- Equation (14). Within one screen the expected squared norm of each part is the true
squared norm plus half of the total error variance. -/
theorem single_screen_bias (E : (Ω → ℝ) →ₗ[ℝ] ℝ) (hE1 : E (fun _ => (1 : ℝ)) = 1)
    (c : ℝ) (hc : c ^ 2 = 1) (M : n → n → ℝ) (ε : n → n → Ω → ℝ) (σ2 : n → n → ℝ)
    (h0 : ∀ i j, E (ε i j) = 0)
    (hvar : ∀ i j, E (fun ω => ε i j ω * ε i j ω) = σ2 i j)
    (hcross : ∀ i j, i ≠ j → E (fun ω => ε i j ω * ε j i ω) = 0)
    (hdiag : ∀ i, σ2 i i = 0) :
    E (fun ω => ∑ i, ∑ j, part c (fun a b => M a b + ε a b ω) i j
        * part c (fun a b => M a b + ε a b ω) i j)
      = ∑ i, ∑ j, part c M i j * part c M i j + (1 / 2) * ∑ i, ∑ j, σ2 i j := by
  have key : ∀ i j, E (fun ω => part c (fun a b => M a b + ε a b ω) i j
        * part c (fun a b => M a b + ε a b ω) i j)
      = part c M i j * part c M i j + (σ2 i j + σ2 j i) / 4 := by
    intro i j
    have hx : E (fun ω => ε i j ω * ε j i ω) = 0 := by
      by_cases hij : i = j
      · subst hij; rw [hvar, hdiag]
      · exact hcross i j hij
    have hy : E (fun ω => ε j i ω * ε i j ω) = 0 := by
      by_cases hij : i = j
      · subst hij; rw [hvar, hdiag]
      · exact hcross j i (Ne.symm hij)
    have h := E_prod E hE1 c (M i j) (M j i) (ε i j) (ε j i) (ε i j) (ε j i)
      (h0 i j) (h0 j i) (h0 i j) (h0 j i)
    rw [hvar i j, hvar j i, hx, hy] at h
    simp only [part]
    rw [h]
    linear_combination (σ2 j i / 4) * hc
  rw [E_sum]
  simp only [key]
  exact sum_split _ _

/-- Equation (15). The cross product of the same part in two screens with independent errors
has the true squared norm as expectation. -/
theorem cross_screen_unbiased (E : (Ω → ℝ) →ₗ[ℝ] ℝ) (hE1 : E (fun _ => (1 : ℝ)) = 1)
    (c : ℝ) (M : n → n → ℝ) (ε₁ ε₂ : n → n → Ω → ℝ)
    (h1 : ∀ i j, E (ε₁ i j) = 0) (h2 : ∀ i j, E (ε₂ i j) = 0)
    (hind : ∀ i j k l, E (fun ω => ε₁ i j ω * ε₂ k l ω) = 0) :
    E (fun ω => ∑ i, ∑ j, part c (fun a b => M a b + ε₁ a b ω) i j
        * part c (fun a b => M a b + ε₂ a b ω) i j)
      = ∑ i, ∑ j, part c M i j * part c M i j := by
  have key : ∀ i j, E (fun ω => part c (fun a b => M a b + ε₁ a b ω) i j
        * part c (fun a b => M a b + ε₂ a b ω) i j)
      = part c M i j * part c M i j := by
    intro i j
    have h := E_prod E hE1 c (M i j) (M j i) (ε₁ i j) (ε₁ j i) (ε₂ i j) (ε₂ j i)
      (h1 i j) (h1 j i) (h2 i j) (h2 j i)
    rw [hind i j i j, hind i j j i, hind j i i j, hind j i j i] at h
    simp only [part]
    rw [h]
    ring
  rw [E_sum]
  simp only [key]

end Reciprocity
