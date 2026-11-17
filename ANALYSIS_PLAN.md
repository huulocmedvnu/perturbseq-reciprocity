# Analysis plan

Definitions, thresholds and decision rules of the study, with the outcome of every rule.
Sections 1 to 4 were fixed before the two K562 screens were compared. Section 5 lists analyses
that were added after the primary results were seen. They are exploratory and are labelled as
such in the article.

## 1. Question and data
- Question: is the response of a gene to the knockdown of another gene reciprocal, as linear
  response theory predicts for a system near a steady state and as predictors built on
  co-expression assume?
- Screens (CRISPR interference, single-cell RNA sequencing): K562 genome-wide, K562
  essential-gene and RPE1 essential-gene (Replogle et al. 2022, scPerturb release), HepG2 and
  Jurkat essential-gene (Nadig et al. 2025, GEO GSE264667).

## 2. Definitions
- Expression of a cell: ln(1 + counts per 10,000). Every cell is standardised against the
  non-targeting control cells of its own batch (mean and standard deviation per gene, standard
  deviation not below 0.05, batches with fewer than 20 control cells not used).
- Response of gene B to the knockdown of gene A: mean standardised expression of B in the cells
  that carry the guides against A. t statistic: response divided by
  sqrt(variance across cells / number of cells + 1 / number of control cells).
- Square gene set of a pair of screens: genes that are knocked down and measured in both, with
  an on-target t of -3 or lower and a mean of at least 0.25 unique molecular identifiers (UMI)
  per control cell in both screens (sensitivity: 0.5).
- Response matrix M. Scaled form (primary): M_AB = response of B / |on-target response of A| x
  w_A, with w_A = t_AA^2 / (t_AA^2 + 1). Unscaled form (sensitivity): the response itself.
- Symmetric part S = (M + M') / 2 and skew-symmetric part A = (M - M') / 2, diagonal excluded.
- Replicable variance between two screens: V_sym = <S1, S2>, V_asym = <A1, A2>. Independent
  measurement error has expectation zero in both products.
- Skew-symmetric share = V_asym / (V_sym + V_asym).
- Active pair: |t| of 3 or more in at least one direction in one half of the cells of screen 1.
  Estimates for active pairs use the other half of screen 1 and all of screen 2.
- Within one screen, the two screens are disjoint quarters of the cells of each knockdown
  (quarters 1 and 2 against quarters 3 and 4, active pairs selected in quarter 1 and read in
  quarter 2).
- Uncertainty: 2,000 bootstrap samples of genes (95% percentile interval). Null distribution:
  10,000 permutations of the gene labels of screen 2.

## 3. Rules fixed before the comparison of the two K562 screens
1. The skew-symmetric part replicates: Pearson correlation of A between the screens of at
   least 0.20 among active pairs, above all 10,000 permutations.
2. The share is not negligible: skew-symmetric share of at least 10% among active pairs.
3. The direction has biological meaning: among active transcription factor and target pairs
   (DoRothEA confidence A and B, TRRUST), the larger effect is from factor to target in at
   least 60%, and subunits of one protein complex have a lower share than other pairs.
If rules 1 to 3 held in K562, the same analysis was to be repeated in RPE1, HepG2 and Jurkat,
and the prediction from co-expression was to be examined, with the expectation that its errors
concentrate in pairs with a large skew-symmetric part.

## 4. Outcome of the rules
1. Met. K562, active pairs: correlation 0.865 (95% interval 0.846 to 0.881), largest value
   among 10,000 permutations 0.037.
2. Met. Share 0.499 (0.496 to 0.502).
3. Met in part. Complex subunits (six complexes named by gene family, the definition in use
   when the rule was fixed) had a share of 0.32 (0.20 to 0.42). The direction rule reached 60%
   in K562 on 15 and 17 pairs only (10 of 15, 13 of 17), and it was met in HepG2 (39 of 49, 39
   of 55) and Jurkat (53 of 64, 45 of 53) but not in RPE1 (17 of 29, 16 of 31).
- Expectation on prediction errors: not supported. The squared error of the prediction had the
  same rank correlation with the size of the skew-symmetric part (0.17) as with the size of the
  symmetric part (0.18).

## 5. Analyses added after the primary results (exploratory)
- Protein complexes from CORUM in place of the six named complexes (the named complexes are
  kept as a second definition).
- CollecTRI as a second source of transcription factor and target pairs.
- Exclusion of gene pairs whose transcription start sites lie within 1 Mb on one chromosome.
- Calibration of the estimator on simulated screens with a known share.
- Gene-level terms: least-squares fit of A_ij = p_i - p_j and S_ij = g + u_i + u_j, and the
  share that remains after their removal.
- Counts of pairs with a replicated effect (|t| of 3 or more in both screens, same sign) in one
  direction or in both. Reverse direction absent: |t| below 2 in both screens.
- Effect sent and received by each gene.
- Prediction of held-out knockdowns. The knocked-down genes are divided into five groups. All
  responses to the knockdown of the genes of one group are withheld from the first screen, a
  model is fitted to the other knockdowns and its prediction is compared with the second screen
  for the whole response and its two parts. Comparators: co-expression with one scale, with the
  on-target effect as scale and with a scale fitted to the scored response, the mean of the
  training knockdowns, and the reverse direction of each pair.
  Outcome: the model reached 0.37 to 0.57 for the whole response and for the skew-symmetric
  part, against 0.10 to 0.19 for co-expression with one scale and 0.24 to 0.43 with a scale
  fitted to the scored response. With the response to the same knockdown in three other cell
  lines as a further input it reached 0.45 to 0.62 (ceiling 0.46 to 0.67).
- The same comparison on unnormalised counts.
  Outcome: co-expression with a scale fitted to the scored response reached 0.77 to 0.90, and
  the mean count of each gene with a scale fitted in the same way 0.81 to 0.91. The model
  trained on other knockdowns reached 0.41 to 0.54 and co-expression with one scale at most 0.20.
- Prediction of the response matrix of one cell line from the three others.
  Outcome: 0.49 to 0.54. Co-expression of the target line changed this by -0.04 to 0.02.
