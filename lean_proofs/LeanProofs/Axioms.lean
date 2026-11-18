/-
Axioms used by the main theorems. Only the three standard axioms of Lean and Mathlib
(`propext`, `Classical.choice`, `Quot.sound`) may appear.
-/
import LeanProofs.Lyapunov
import LeanProofs.SlowMode
import LeanProofs.Frobenius
import LeanProofs.Noise
import LeanProofs.Hodge
import LeanProofs.Ridge

open Reciprocity

#print axioms flux_transpose
#print axioms neg_J_mul_Sigma
#print axioms inv_D_sub_flux
#print axioms response_eq_Sigma_mul_inv
#print axioms response_detailed_balance
#print axioms response_symmetric_of_detailed_balance
#print axioms rankOne_skew_zero_iff
#print axioms rankOne_covariance_skew
#print axioms standardised_ratio
#print axioms ip_symPart_skewPart
#print axioms nsq_eq_sym_add_skew
#print axioms nsq_error_symmetric
#print axioms reciprocity_bound
#print axioms single_screen_bias
#print axioms cross_screen_unbiased
#print axioms pot_sum_zero
#print axioms ip_grad_resid
#print axioms nsq_eq_grad_add_resid
#print axioms pot_minimises
#print axioms uHat_sum_zero
#print axioms uHat_normal_equation
#print axioms loo_closed_form
