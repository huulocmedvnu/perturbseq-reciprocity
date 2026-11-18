# Machine-checked proofs of the identities in Methods

Lean 4 proofs, with Mathlib, of the algebraic identities behind the numbered equations of the
article. No proof uses `sorry`, and the theorems depend on the three standard axioms of Lean and
Mathlib only (`propext`, `Classical.choice`, `Quot.sound`), which `LeanProofs/Axioms.lean` prints
during the build.

| File | Equations | Statements |
|---|---|---|
| `Lyapunov.lean` | (1), (5) | `flux_transpose`: Q is skew-symmetric. `neg_J_mul_Sigma`: the Lyapunov equation gives -JΣ = D - Q. `inv_D_sub_flux`: (D - Q)⁻¹ = -Σ⁻¹J⁻¹. `response_eq_Sigma_mul_inv`: -J⁻¹ = Σ(D - Q)⁻¹. `response_detailed_balance`, `response_symmetric_of_detailed_balance`: with Q = 0 the response is ΣD⁻¹, and it is symmetric when D is a multiple of the identity. |
| `SlowMode.lean` | (4), (13) | `rankOne_skew_zero_iff`: the skew-symmetric part of the rank-one response p_B q_A vanishes for every pair exactly when p and q are linearly dependent. `rankOne_covariance_skew`: the rank-one covariance p_B p_A has no skew-symmetric part. `standardised_ratio`: the response ratio on the standardised scale is the correlation. |
| `Frobenius.lean` | (6), (7) | `ip_symPart_skewPart`: the two parts are orthogonal. `nsq_eq_sym_add_skew`: the squared norms add. `nsq_error_symmetric`: the error of a symmetric prediction is the error on the symmetric part plus the squared norm of the skew-symmetric part. `reciprocity_bound`: equation (7). |
| `Noise.lean` | (14), (15) | `single_screen_bias`: within one screen, each part gains half of the total error variance. `cross_screen_unbiased`: the cross product of two screens has the true squared norm as expectation. The expectation is a linear functional with E 1 = 1. |
| `Hodge.lean` | (16), (17) | `pot_sum_zero`, `ip_grad_resid`, `nsq_eq_grad_add_resid`, `pot_minimises`: the gene-level terms sum to zero, are orthogonal to the remainder and minimise the squared error. `uHat_sum_zero`, `uHat_normal_equation`: the terms of the symmetric part sum to zero and satisfy the normal equations of the off-diagonal fit. |
| `Ridge.lean` | (22) | `loo_closed_form`: the prediction for a gene by the ridge model fitted without it equals the closed form, with the regularised kernel in block form. |

Equations (2) and (3) are standard results that are stated and not formalised here. Equation (4)
is verified as the algebraic statement on the two rank-one matrices, not as a limit. Equations
(8) to (12), (18) to (21) and (23) are definitions.

## Build

Requires `elan`. The toolchain and the Mathlib revision are pinned in `lean-toolchain` and
`lake-manifest.json` (Lean 4.36.0, Mathlib v4.36.0).

```
cd lean_proofs
lake exe cache get     # precompiled Mathlib
lake build
```

The build ends with `Build completed successfully` and lists the axioms of each theorem.
