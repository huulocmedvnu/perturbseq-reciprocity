# Non-reciprocal responses to gene knockdown in genome-scale Perturb-seq screens

Code, analysis plan and results for a test of reciprocity in the gene-by-gene response matrix of
five CRISPR interference Perturb-seq screens in four human cell lines.

Author: Thai Huu Loc Tran, Department of Obstetrics and Gynecology, Hong Hung Hospital, Tay Ninh City 80116,
Viet Nam (tranthaihuuloc2@gmail.com, ORCID 0009-0009-1810-010X).

## Question

Linear response theory predicts that, near a steady state and with every gene on the scale of its
variation in unperturbed cells, the response of gene B to the knockdown of gene A equals the response
of A to the knockdown of B. Predictors of perturbation response that are built on co-expression in
unperturbed cells inherit this symmetry. The study measures how much of the replicable response is
symmetric and how much is skew-symmetric.

## Method in brief

For the genes that are both knocked down and measured, the response matrix M is divided into a
symmetric part S = (M + M') / 2 and a skew-symmetric part A = (M - M') / 2. The variance of each part
that replicates is the inner product of the part between two independent screens, in which independent
measurement error has expectation zero. The skew-symmetric share is V_asym / (V_sym + V_asym).
Definitions, thresholds, the rules fixed in advance and their outcomes are in `ANALYSIS_PLAN.md`.

## Main results

| Comparison | Genes | Active pairs | Correlation, symmetric part | Correlation, skew-symmetric part | Skew-symmetric share (95% interval) |
|---|---|---|---|---|---|
| K562, two independent screens | 1,501 | 88,630 | 0.85 | 0.86 | 0.499 (0.496 to 0.502) |
| K562, cell quarters | 1,589 | 83,221 | 0.88 | 0.90 | 0.498 (0.495 to 0.501) |
| RPE1, cell quarters | 1,596 | 141,202 | 0.87 | 0.87 | 0.462 (0.453 to 0.471) |
| HepG2, cell quarters | 1,738 | 85,938 | 0.75 | 0.78 | 0.475 (0.467 to 0.482) |
| Jurkat, cell quarters | 1,643 | 56,149 | 0.85 | 0.88 | 0.496 (0.493 to 0.499) |

With all responses to a knockdown withheld, a model trained on the other knockdowns predicted them
with a correlation of 0.37 to 0.57 with the second screen (skew-symmetric part 0.37 to 0.57), against
0.10 to 0.19 for co-expression with one scale and 0.24 to 0.43 for co-expression with a scale fitted
to the scored response (`results/heldout_prediction.tsv`).

Simulated screens with a true share of 0, 0.10, 0.25 and 0.50 return 0.000, 0.100, 0.250 and 0.500.
All estimates are in `results/`, and the data behind every figure panel are in
`results/figure_source_data/`.

## Layout

- `ANALYSIS_PLAN.md` definitions, thresholds and decision rules, with the outcome of every rule.
- `src/lib.rs` readers and shared statistics. `src/bin/` one program per analysis:
  - `reduce.rs` single-cell counts to the square response summary of one screen (mean and variance of
    standardised expression per knockdown, cell halves and quarters, correlation in control cells).
  - `reciprocity.rs` symmetric and skew-symmetric parts between two screens, bootstrap intervals,
    permutation test, classes of gene pairs, direction of transcription factor and target pairs.
  - `calibration.rs` the estimator on simulated screens with a known share.
  - `gene_level.rs` gene-level terms, sign agreement of the two directions, effect sent and received.
  - `directed_pairs.rs` counts of pairs with a replicated effect in one direction or in both.
  - `linear_response.rs` prediction of the response from co-expression in control cells.
  - `heldout_prediction.rs` prediction of knockdowns that are held out of the training screen, on
    the standardised response and on unnormalised counts, against co-expression predictors.
  - `transfer_prediction.rs` prediction of the response matrix of one cell line from other lines.
- `scripts/download_annotations.sh` curated gene-pair annotations (DoRothEA, TRRUST, CollecTRI, CORUM,
  transcription start sites). `scripts/complex_families.sh` six complexes named by gene family.
- `scripts/build_reduced.sh` public single-cell files to `reduced/` and `reduced_raw/` (count shifts).
- `run_all.sh` every analysis of the article, from `reduced/` to `results/`.
- `scripts/plot_figures.py`, `scripts/make_flowchart.py` figures 1 to 7 and their source data.
- `results/` result tables. `figures/` figure files (PDF).

## Reproduce

Requirements: Rust 1.75 or later. For the figures, Python 3 with plotly 5, kaleido 0.2, pandas, numpy,
and `rsvg-convert`. For the step from raw data, the HDF5 library and Python with h5py.

```
bash scripts/download_annotations.sh     # annotations, a few megabytes
bash scripts/build_reduced.sh            # raw data to reduced screens (about 25 GB of downloads)
bash run_all.sh                          # all analyses
python3 scripts/make_flowchart.py && python3 scripts/plot_figures.py
```

The reduced screens (the output of `build_reduced.sh`) are attached to the release as
`reduced_screens.tar.zst` and `reduced_raw_screens.tar.zst`, which allows `run_all.sh` to be run without
the raw data. Unpack both archives in the repository root.

`run_all.sh` uses all processor threads (set `RAYON_NUM_THREADS` to change this). The bootstrap and
permutation steps draw from fixed random streams, and repeated runs give identical tables.

## Data sources

| Data | Source |
|---|---|
| K562 genome-wide, K562 essential-gene, RPE1 essential-gene screens | Replogle et al., Cell 2022, as distributed by scPerturb (Zenodo record 13350497) |
| HepG2 and Jurkat essential-gene screens | Nadig et al., Nature Genetics 2025 (GEO GSE264667) |
| Transcription factor and target pairs | DoRothEA confidence A and B, CollecTRI (through OmniPath), TRRUST v2 |
| Protein complexes | CORUM (through OmniPath) |
| Transcription start sites | NCBI RefSeq Select, hg38 (UCSC) |

Third-party annotation files are downloaded by the script and are not redistributed here.

## Licence

The code is released under the MIT licence (`LICENSE`). The result tables and figure source data in
`results/` are released under CC0 1.0.
